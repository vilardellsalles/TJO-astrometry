"""
Select and reduce images to be downloaded by observers
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os.path
import argparse
import warnings
from re import match
from glob import glob
from shutil import rmtree
from socket import gethostname
from traceback import format_exc
from collections import defaultdict

import astropy.coordinates as coord
from astropy import table
from astropy import log as astrolog
from astropy.time import Time
from astropy.units import day, second
import ccdproc
from numpy import cos, floor, ceil

import icat.pipelines as ppl
from icat import system
from icat import astrometry
from icat import photometry
from icat.logs import get_logger, LogOrProgressBar
from icat.config import config_option, config_section
from icat.notifications import send_mail
from icat.mysql_database import ImagesDB


imgqs_ppl = ppl.import_pipeline("gaia-imgqs")
catstore = ppl.import_pipeline("gaia-catstore")
astrolog.setLevel("ERROR")
logger = get_logger(__name__)

min_nbias = 5
min_ndarks = 5
min_nflats = 3

base_msg = "Dear TJO operator,\n\n{}\n\nSincerely,\n\nICAT reduction pipeline."
same_focus = config_option(__name__, "same_focus")
site = config_section("Location")
OAdM = coord.EarthLocation(lat=site["latitude"], lon=site["longitude"],
                           height=site["elevation"])


class MissingError(ValueError):
    """ Custom class to manage missing images that impede the reduction """
    pass


class TJODB(ImagesDB):
    """
    Class to read TJO specific information from an ICAT MySQL database
    """

    columns = {"id": 0, "original_ids": "", "filename": "",
               "naxis1": 0, "naxis2": 0, "object": "", "propcode": 0,
               "imgqs": "", "obstype": "", "instrume": "", "filter": "",
               "date-obs": "", "time-obs": "", "jd": 0,
               "elevatio": "", "ra": "", "dec": "", "domestat": "",
               "camtemp": 0, "exptime": 0, "defocus": 0, "focuspos": 0}

    @property
    def sql_columns(self):
        """
        Return columns suited for SQL queries
        """

        return ",".join("`{}`".format(col) for col in self.columns)

    def _one_or_many(self, field, elements):
        """
        Return an SQL comparison containting either 'IN' or '=',
        depending on the number of elements to be found
        """

        values = tuple(_value for _value in set(elements)
                       if _value is not None)

        num_values = len(values)
        if num_values > 1:
            return "{} IN {}".format(field, values)
        elif num_values == 1:
            return "{}={}".format(field, values).replace(",)", ")")
        else:
            return "True"

    def _masked(self, images):
        """
        Return a table with all None values masked with fill values
        specified in self.columns
        """

        masked_images = table.Table(images, masked=True)
        for col, fill in self.columns.iteritems():
            masked_images[col].mask = [images[col] == [None]]
            masked_images[col] = masked_images[col].filled(fill).tolist()
            masked_images[col].fill_value = fill
            masked_images[col].mask = [images[col] == [None]]

        null_column = [""] * len(masked_images)
        for col in ["bias", "darks", "flats", "invalid", "destination"]:
            masked_images[col] = null_column
            masked_images[col].mask = True

        return masked_images

    def proposals(self, propcode=None):
        """
        Return active proposals at the TJO
        """

        if propcode is None:
            query = ("SELECT prop.id, user.password "
                     "FROM ph1_proposals AS prop, ph0_users AS user "
                     "WHERE prop.user_id=user.id AND prop.admin_state_id>=4")
        else:
            try:
                valid = [int(value[1:]) if value.startswith("p")
                         else int(value) for value in propcode]
                valid_propcodes = self._one_or_many("prop.id", valid)
            except Exception:
                err_msg = "Invalid proposal code found in: {}"
                raise ppl.PipelineError(err_msg.format(propcode))

            query_msg = ("SELECT prop.id, user.password "
                         "FROM ph1_proposals AS prop, ph0_users AS user "
                         "WHERE prop.user_id=user.id AND {}")

            query = query_msg.format(valid_propcodes)

        return dict(self.execute(query))

    def get_images(self, image_list):
        """
        Look for images in the database from a list of images.
        An Astropy table is returnded with several values
        """

        log_msg = "Looking for {} images in database..."
        self.logger.debug(log_msg.format(len(image_list)))

        all_names = tuple(self._trim_prefix(name) for name in image_list)

        query_msg = "SELECT {} FROM {} WHERE {}"
        query = query_msg.format(self.sql_columns, self.images_table,
                                 self._one_or_many("filename", all_names))
        result = self.execute(query)

        images = None
        if result:
            images = table.Table(rows=result, names=self.columns)
            missing = [_image for _image in all_names
                       if _image not in images["filename"]]

        else:
            missing = all_names

        if missing:
            miss_len = len(missing)
            miss_plural = "s were" if miss_len > 1 else " was"
            err_msg = "{} image{} not found in the database"
            raise ppl.PipelineError(err_msg.format(miss_len, miss_plural))

        return self._masked(images)

    def constraint(self, column):
        """
        Return a string from an Astropy table.Column specifying the
        constraints to be used in a SQL query
        """

        return self._one_or_many("`{}`".format(column.name), column)

    def find_calibration(self, images):
        """
        Retrieve all possible calibration images from database
        """

        obstype = str("('Bias', 'Dark', 'Flat')")
        jd_min = int(min(images["jd"][~images["jd"].mask]-30))
        jd_max = int(max(images["jd"][~images["jd"].mask]+30))
        naxis1 = images["naxis1"][~images["naxis1"].mask]
        naxis2 = images["naxis2"][~images["naxis2"].mask]
        instrume = images["instrume"][~images["instrume"].mask].astype(str)

        criteria = ["`original_ids` IS NULL"]
        criteria += ["`obstype` IN {}".format(obstype)]
        criteria += ["`{0}`>{1} AND `{0}`<{2}".format("jd", jd_min, jd_max)]
        criteria += [self.constraint(naxis1)]
        criteria += [self.constraint(naxis2)]
        criteria += [self.constraint(instrume)]

        log_msg = "Looking for calibration images in database with {} "
        log_msg += "constraints"
        logger.debug(log_msg.format(len(criteria)))

        query_msg = "SELECT {} FROM {} WHERE {}"
        query = query_msg.format(self.sql_columns, self.images_table,
                                 " AND ".join(criteria))
        result = self.execute(query)

        if result:
            calibration_images = table.Table(rows=result, names=self.columns)
            log_msg = "{} calibration images found in database"
            self.logger.debug(log_msg.format(len(calibration_images)))
        else:
            log_msg = "No valid calibration images found in '{}' table"
            self.logger.warning(log_msg.format(self.images_table))
            return images[:0].copy()

        return self._masked(calibration_images)

    def cloud_coverage(self, images):
        """
        Return the sky minus ambient temperature from the cloud sensor
        in the given julian day range
        """

        try:
            science_mask = images["obstype"] == "Science"
            min_jd = min(images["jd"][science_mask].filled(1e20))
            max_jd = max(images["jd"][science_mask].filled(0))

            query = ("SELECT skyminusambienttemperature_cloud_sensor_jd, "
                     "skyminusambienttemperature_cloud_sensor "
                     "FROM oadm_data.meteo2 "
                     "WHERE skyminusambienttemperature_cloud_sensor_jd > {} "
                     "AND skyminusambienttemperature_cloud_sensor_jd < {}")

            results = self.execute(query.format(floor(min_jd), ceil(max_jd)))

            clouds = table.Table(rows=results, names=["jd", "temperature"],
                                 masked=True, dtype=(float, float))

            mask = [clouds["temperature"] < -900]
            clouds["jd"].mask = mask
            clouds["temperature"].mask = mask

            return clouds

        except Exception as err:
            err_msg = "An error occurred when retrieving cloud coverage: {}"
            self.logger.warning(err_msg.format(err))
            return


def limit_nights(images, night, min_num=-1):
    """
    If more than the minimum number of calibration images are found,
    limit the number of nights to the closest nights, if desired
    """

    num_images = len(images)
    selected_nights = [True] * num_images
    if num_images >= min_num:
        images["diff_jd"] = abs(images["jd"].filled().astype(int) - night)
        images.sort(["diff_jd"])
        max_diff = images[min_num-1]["diff_jd"]
        selected_nights = images["diff_jd"] <= max_diff

    image_ids = set(images["id"][selected_nights].astype(str))

    num_filtered = len(image_ids)
    if num_filtered >= min_num:
        log_msg = "%s '%s' filtered to those %s "
        log_msg += "having a maximum difference of %s nights"
        logger.debug(log_msg, num_images, images["obstype"][0], num_filtered,
                     max_diff)

        return "|{}|".format("|".join(image_ids))
    else:
        return "||"


def update_calibration_ids(image, column, value):
    """
    Add all the ids in a calibration field without truncating
    its content
    """

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            image[column] = value
    except Exception:
        string_len = len(value)
        image.table[column] = image.table[column].astype(("U", string_len))
        image[column] = value


def select_bias(image, calib):
    """
    Select suitable bias images for a given image
    """

    if image["bias"]:
        return

    if not calib:
        return update_calibration_ids(image, "bias", "||")

    camtemp = image["camtemp"] or calib["camtemp"].fill_value
    jd = int(image["jd"] or calib["jd"].fill_value)

    mask = calib["obstype"] == "Bias"
    mask *= calib["filter"] == "C"
    mask *= calib["exptime"].astype(int) == 0
    mask *= calib["instrume"] == image["instrume"]
    mask *= calib["naxis1"] == image["naxis1"]
    mask *= calib["naxis2"] == image["naxis2"]
    mask *= calib["original_ids"].mask
    mask *= abs(calib["camtemp"].filled() - camtemp) < 3
    mask *= abs(calib["jd"].filled().astype(int) - jd) < 30
    mask *= [match("^[02][0-9][0-9]$", imgqs) is not None
             for imgqs in calib["imgqs"].filled()]

    selected_bias = limit_nights(calib[mask], jd, min_nbias)
    update_calibration_ids(image, "bias", selected_bias)


def select_darks(image, bias, darks):
    """
    Select suitable dark images for a given image
    """

    if image["darks"]:
        return

    if not darks:
        return update_calibration_ids(image, "darks", "||")

    camtemp = image["camtemp"] or darks["camtemp"].fill_value
    jd = int(image["jd"] or darks["jd"].fill_value)

    mask = darks["obstype"] == "Dark"
    mask *= darks["filter"] == "C"
    mask *= darks["exptime"] >= image["exptime"]
    mask *= darks["instrume"] == image["instrume"]
    mask *= darks["naxis1"] == image["naxis1"]
    mask *= darks["naxis2"] == image["naxis2"]
    mask *= darks["original_ids"].mask
    mask *= abs(darks["camtemp"].filled() - camtemp) < 3
    mask *= abs(darks["jd"].filled().astype(int) - jd) < 30
    mask *= [match("^[02][0-9][0-9]$", imgqs) is not None
             for imgqs in darks["imgqs"].filled()]

    for dark, valid in zip(darks, mask):
        if not valid:
            continue

        select_bias(dark, bias)

    with_bias = (darks["bias"] != "||") & (darks["bias"] != "")
    selected_darks = limit_nights(darks[mask & with_bias], jd, min_ndarks)
    update_calibration_ids(image, "darks", selected_darks)


def select_flats(image, bias, darks, flats):
    """
    Select suitable flat images to reduce science images
    """

    if image["flats"]:
        return

    if not flats:
        return update_calibration_ids(image, "flats", "||")

    camtemp = image["camtemp"] or flats["camtemp"].fill_value
    jd = int(image["jd"] or flats["jd"].fill_value)
    focuspos = image["focuspos"] or flats["focuspos"].fill_value
    defocus = image["defocus"] or flats["defocus"].fill_value

    mask = flats["obstype"] == "Flat"
    mask *= flats["domestat"] == "Open"
    mask *= flats["filter"] == image["filter"]
    mask *= abs(flats["defocus"].filled() - defocus) < same_focus
    mask *= abs(flats["focuspos"].filled() - focuspos) < same_focus
    mask *= flats["instrume"] == image["instrume"]
    mask *= flats["naxis1"] == image["naxis1"]
    mask *= flats["naxis2"] == image["naxis2"]
    mask *= flats["original_ids"].mask
    mask *= abs(flats["camtemp"].filled() - camtemp) < 3
    mask *= abs(flats["jd"].filled().astype(int) - jd) < 30
    mask *= [match("^[02][0-9][0-9]$", imgqs) is not None
             for imgqs in flats["imgqs"].filled()]

    for flat, valid in zip(flats, mask):
        if not valid:
            continue

        select_bias(flat, bias)
        select_darks(flat, bias, darks)

    with_bias = (flats["bias"] != "||") & (flats["bias"] != "")
    with_darks = (flats["darks"] != "||") & (flats["darks"] != "")
    selected_flats = limit_nights(flats[mask & with_bias & with_darks],
                                  jd, min_nflats)
    update_calibration_ids(image, "flats", selected_flats)


def find_originals(selected, images):
    """
    Find all the original images required to reduce a list of images
    """

    all_ids = images["id"].astype(str).tolist()
    required = set(selected["id"].astype(str).tolist())
    for sel in selected:
        if sel["bias"]:
            required |= set(sel["bias"].strip("|").split("|"))

        if sel["darks"]:
            dark_ids = set(sel["darks"].strip("|").split("|"))
            sel_darks = images[[num in dark_ids for num in all_ids]]

            for dark in sel_darks:
                required |= set(dark["bias"].strip("|").split("|"))

            required |= dark_ids

    return "|{}|".format("|".join(sorted(required)))


def check_existence(reduced_images, ids, destination, imtype="image"):
    """
    Check whether a master frame with the same ids exists,
    either in the current proposal or in any other proposals
    """

    # Check whether it has been created for the current proposal

    same_path = [name for name, original in reduced_images.iteritems()
                 if original == ids and name.startswith(destination)]

    num_images = len(same_path)
    if num_images > 1:
        err_msg = "More than one {} possible in the same proposal: {}"
        raise ppl.PipelineError(err_msg.format(imtype, ", ".join(same_path)))

    elif num_images == 1:
        image_name = same_path[0]
        logger.debug("Using {}: {}".format(imtype, image_name))
        return image_name

    # Check whether it has been created for any other proposal

    all_paths = [name for name, original in reduced_images.iteritems()
                 if original == ids]
    different = set(os.path.basename(name) for name in all_paths)

    num_images = len(different)
    if num_images > 1:
        err_msg = "More than one {} possible in other proposals: {}"
        raise ppl.PipelineError(err_msg.format(imtype, ", ".join(all_paths)))

    elif num_images == 1:
        orig_image = all_paths[0]
        orig_ids = reduced_images[orig_image]

        system.new_copy(orig_image, destination)

        image_name = os.path.join(destination, different.pop())
        reduced_images[image_name] = orig_ids

        logger.debug("Using copied {}: {}".format(imtype, image_name))
        return image_name


def ccdred(*args, **kwargs):
    """
    Unfortunately, ccd_process in ccdproc (1.2.0) adds a lot of rubbish
    in header as HIERARCH keywords. This function removes these keywords
    """

    sub_data = ccdproc.ccd_process(*args, **kwargs)

    hierarch = [key for key in sub_data.header.iteritems()
                if len(key[0]) > 8 and key[1] in sub_data.header]

    for key, value in hierarch:
        del(sub_data.header[key])
        del(sub_data.header[value])

    return sub_data


def bias_combine(reduced_images, selected_bias, origin, destination):
    """
    Create master bias, if required
    """

    bias_ids = sorted(selected_bias["id"].astype(str).tolist())
    ids = "|{}|".format("|".join(bias_ids))

    bias_name = check_existence(reduced_images, ids, destination, "'Bias'")

    if bias_name:
        return bias_name

    names = [os.path.splitext(bias)[0] if bias.endswith(".gz") else bias
             for bias in selected_bias["filename"]]
    bias_list = [os.path.join(origin, os.path.basename(bias))
                 for bias in names]

    instrume = [name for name in bias_list[0].split("_") if len(name) == 1]
    bias_basename = "TJO{:.5f}_{}_cbs.fits".format(Time.now().jd,
                                                   instrume.pop())
    bias_name = os.path.join(destination, bias_basename)

    log_msg = "Using {} bias to create master bias '{}'"
    logger.debug(log_msg.format(len(bias_list), bias_name))

    ccdproc.combine(bias_list, output_file=bias_name, mem_limit=1e8,
                    method="median", unit="adu", dtype="float32")

    reduced_images[bias_name] = ids

    return bias_name


def dark_combine(reduced_images, selected_darks, images, origin, destination):
    """
    Create master dark, if required
    """

    ids = find_originals(selected_darks, images)

    dark_name = check_existence(reduced_images, ids, destination, "'Dark'")

    if dark_name:
        return dark_name

    names = [os.path.splitext(dark)[0] if dark.endswith(".gz") else dark
             for dark in selected_darks["filename"]]
    dark_list = [os.path.join(origin, os.path.basename(dark))
                 for dark in names]

    instrume = [name for name in dark_list[0].split("_") if len(name) == 1]
    dark_basename = "TJO{:.5f}_{}_cth.fits".format(Time.now().jd,
                                                   instrume.pop())
    dark_name = os.path.join(destination, dark_basename)

    log_msg = "Cleaning {} darks to create master dark '{}'"
    logger.debug(log_msg.format(len(dark_list), dark_name))

    scaling = []
    clean_images = []
    all_ids = images["id"].astype(str)
    for selected in selected_darks:
        requires = selected["bias"].strip("|").split("|")
        bias = images[[num in requires for num in all_ids]]

        master_bias = bias_combine(reduced_images, bias, origin, destination)

        raw_data = ccdproc.CCDData.read(selected["filename"], unit="adu")
        master_data = ccdproc.CCDData.read(master_bias, unit="adu")

        sub_data = ccdproc.subtract_bias(raw_data, master_data,
                                         add_keyword=False)

        clean_images += [sub_data]
        scaling += [1.0 / selected["exptime"]]

    master_dark = ccdproc.combine(clean_images, mem_limit=1e8, scale=scaling,
                                  method="median", unit="adu", dtype="float32")
    master_dark.header["EXPTIME"] = 1.0
    master_dark.write(dark_name)

    reduced_images[dark_name] = ids

    return dark_name


def flat_combine(reduced_images, selected_flats, images, origin, destination):
    """
    Create master flat, if required
    """

    ids = find_originals(selected_flats, images)

    flat_name = check_existence(reduced_images, ids, destination, "'Flat'")

    if flat_name:
        return flat_name

    names = [os.path.splitext(flat)[0] if flat.endswith(".gz") else flat
             for flat in selected_flats["filename"]]
    flat_list = [os.path.join(origin, os.path.basename(flat))
                 for flat in names]

    band = selected_flats["filter"].tolist().pop()
    instrume = [name for name in flat_list[0].split("_") if len(name) == 1]
    flat_basename = "TJO{:.5f}_{}_cf{}.fits".format(Time.now().jd,
                                                    instrume.pop(), band)
    flat_name = os.path.join(destination, flat_basename)

    log_msg = "Cleaning {} flats to create master flat '{}'"
    logger.debug(log_msg.format(len(flat_list), flat_name))

    scaling = []
    clean_images = []
    all_ids = images["id"].astype(str).tolist()
    for selected in selected_flats:
        requires = selected["bias"].strip("|").split("|")
        bias = images[[num in requires for num in all_ids]]

        master_bias = bias_combine(reduced_images, bias, origin, destination)

        requires = selected["darks"].strip("|").split("|")
        darks = images[[num in requires for num in all_ids]]

        master_dark = dark_combine(reduced_images, darks, images,
                                   origin, destination)

        raw_data = ccdproc.CCDData.read(selected["filename"], unit="adu")
        bias_data = ccdproc.CCDData.read(master_bias, unit="adu")
        dark_data = ccdproc.CCDData.read(master_dark, unit="adu")

        sub_data = ccdred(raw_data, master_bias=bias_data,
                          dark_frame=dark_data, add_keyword=False,
                          exposure_key="EXPTIME", exposure_unit=second)

        clean_images += [sub_data]

        scaling += [1.0 / sub_data.data.mean()]

    master_flat = ccdproc.combine(clean_images, mem_limit=1e8, scale=scaling,
                                  method="median", unit="adu", dtype="float32")
    master_flat.header["EXPTIME"] = 1.0
    master_flat.write(flat_name)

    reduced_images[flat_name] = ids

    return flat_name


def is_newer(one, other):
    """
    Check whether a file is newer than another file
    """

    if not os.path.isfile(other):
        return True

    # Apparently, modified time can be some milliseconds lower in the
    # copy than in the original file due to floating point issues

    one_modified = os.stat(one).st_mtime - 1E-5
    other_modified = os.stat(other).st_mtime

    return True if one_modified > other_modified else False


def log_sky_brightness(images):
    """
    Return the sky brightness at the time some images where taken
    """

    try:
        julian_days = Time(images["jd"].filled(0), format="jd")

        sun = coord.get_sun(julian_days)
        moon = coord.get_moon(julian_days)
        altaz = moon.transform_to(coord.AltAz(obstime=julian_days,
                                              location=OAdM))
        max_phase = max((1 - cos(moon.separation(sun).radian)) / 2)
        max_elevation = max(altaz.alt.degree)

        if max_phase < 0.3:
            logger.debug("Moon phase lower than 0.3: %s", max_phase)
            return "Sky brightness: Dark"
        elif max_elevation < 0:
            log_msg = "Dark: Moon phase %s and elevation %s"
            logger.debug(log_msg, max_phase, max_elevation)
            return "Sky brightness: Dark"
        else:
            log_msg = "Bright: Moon phase %s and elevation %s"
            logger.debug(log_msg, max_phase, max_elevation)
            return "Sky brightness: Bright"

    except Exception as err:
        if images:
            err_msg = "An error occurred when logging sky brightness: %s"
            logger.warning(err_msg, err, exc_info=True)
        else:
            logger.warning("Unknown sky brightness")

        return "Sky brightness: Bright | Dark"


def log_solar_elevation(images):
    """
    Return the solar elevation at the time some images where taken
    """

    try:
        julian_days = Time(images["jd"].filled(0), format="jd")

        sun = coord.get_sun(julian_days)
        altaz = sun.transform_to(coord.AltAz(obstime=julian_days,
                                             location=OAdM))
        max_elevation = max(altaz.alt.degree)

        if max_elevation > -0.8:
            return "Solar elevation: Day time"
        elif max_elevation > -6:
            return "Solar elevation: Sunset to sunrise"
        elif max_elevation > -12:
            return "Solar elevation: Civil twilight"
        elif max_elevation > -18:
            return "Solar elevation: Nautical twilight"
        else:
            return "Solar elevation: Night time"

    except Exception as err:
        if images:
            err_msg = "An error occurred when logging solar elevation: %s"
            logger.warning(err_msg, err, exc_info=True)
        else:
            logger.warning("Unknown solar elevation")

        return ("Solar elevation: Night time | Nautical twilight | "
                "Civil twilight | Sunset to sunrise | Day Time")


def match_or_error(reduction, existing, path):
    """
    Ensure that reduction images match existing images
    """

    prop_images = table.join(reduction, existing)

    num_selected = len(prop_images)
    num_existing = len(existing)
    if num_existing != num_selected:
        err_msg = "Number of selected images ({}) do not match "
        err_msg += "{} images in '{}'"
        raise ppl.PipelineError(err_msg.format(num_selected,
                                               num_existing, path))

    return prop_images


def create_rawdata_log(rawdata, seeing, cloud_coverage):
    """
    Create rawdata log file/s to be submitted to observers
    """

    raw_images = [os.path.basename(image) for image in rawdata["filename"]]
    rawdata["basename"] = [os.path.splitext(image)[0]
                           if image.endswith(".gz") else image
                           for image in raw_images]

    for path in set(rawdata["destination"][~rawdata["destination"].mask]):

        if not path:
            continue

        logger.debug("Creating rawdata logs for: %s", path)

        in_path = [os.path.basename(image)
                   for image in glob(os.path.join(path, "*"))]
        existing_images = table.Table([in_path], names=["basename"])

        prop_images = match_or_error(rawdata, existing_images, path)

        rawpath = os.path.basename(path)
        prop_images["filename"] = [os.path.join(rawpath, image)
                                   for image in prop_images["basename"]]
        prop_images.sort(["obstype", "filename"])

        cat_path = path.replace("_rawdata", "_catalog")
        try:
            os.makedirs(cat_path)
        except OSError:
            pass
        filename = os.path.join(cat_path, "{}.txt".format(rawpath))

        log_text = ["General remarks:", ""]

        # Only science images having invalid values should produce the
        # following text

        text = ("Some of the Science images do not have the minimum "
                "amount of required {0}. Additional {0} may be acquired "
                "in future nights, but observed time for the affected "
                "Sequences will not be computed.")

        science_mask = ~prop_images["obstype"].mask
        science_mask *= prop_images["obstype"] == "Science"
        invalid_mask = ~prop_images["invalid"].mask
        invalid = prop_images["invalid"][science_mask & invalid_mask]

        if any("bias" in image for image in invalid):
            log_text += [text.format("Bias"), ""]

        if any("darks" in image for image in invalid):
            log_text += [text.format("Darks"), ""]

        if any("flats" in image for image in invalid):
            log_text += [text.format("Sky flats"), ""]

        # Only science images using bias or darks from different nights
        # should produce the following text

        text = ("Some of the Science images use {} from different nights."
                "Those images could be reduced, but observed time for "
                "the affected Sequences will not be computed.")

        image_ids = prop_images["id"].astype(str)

        bias_mask = ~prop_images["bias"].mask
        bias_list = [num.strip("|").split("|")
                     for num in prop_images["bias"][science_mask & bias_mask]]
        bias_id = [num for calib in bias_list for num in calib]
        bias_jds = prop_images["jd"][[num in bias_id for num in image_ids]]

        if len(set(bias_jds.filled(0).astype(int))) > 1:
            log_text += [text.format("Bias"), ""]

        dark_mask = ~prop_images["darks"].mask
        dark_list = [num.strip("|").split("|")
                     for num in prop_images["darks"][science_mask & dark_mask]]
        dark_id = [num for calib in dark_list for num in calib]
        dark_jds = prop_images["jd"][[num in dark_id for num in image_ids]]

        if len(set(dark_jds.filled(0).astype(int))) > 1:
            log_text += [text.format("Darks"), ""]

        log_text += ["", "Observing conditions:", ""]

        night_images = prop_images[science_mask]
        all_nights = set(night_images["jd"].filled(0).astype(int))
        if len(all_nights) > 1:
            night_images = None
            log_msg = "More than one night found when computing observing "
            log_msg = "conditions for '%s'"
            logger.warning(log_msg, filename)

        log_text += [log_sky_brightness(night_images)]
        log_text += [log_solar_elevation(night_images)]

        log_text += ["", "", "Observed images:", ""]

        columns = ("filename", "obstype", "object", "filter",
                   "date-obs", "time-obs", "jd",  "elevatio", "ra", "dec",
                   "instrume", "naxis1", "naxis2",
                   "camtemp", "exptime", "defocus", "focuspos", "domestat")

        log_text += prop_images[columns].pformat(max_width=-1, max_lines=-1)

        if os.path.isfile(filename):
            log_msg = "Adding data to an existing file '{}'"
            logger.warning(log_msg.format(filename))

        with open(filename, "a") as logf:
            logf.write("\n".join(log_text))


def create_reddata_log(reddata):
    """
    Create reddata log file/s to be submitted to observers
    """

    for path in set(os.path.dirname(image) for image in reddata["reduced"]):

        existing_images = table.Table([glob(os.path.join(path, "*"))],
                                      names=["reduced"])

        prop_images = match_or_error(reddata, existing_images, path)

        red_path = os.path.basename(path)
        raw_path = path.replace("_reddata", "_rawdata")
        cat_path = path.replace("_reddata", "_catalog")

        science_images = prop_images[~prop_images["filename"].mask]

        rawnames = [os.path.basename(os.path.splitext(image)[0])
                    if image.endswith(".gz") else os.path.basename(image)
                    for image in science_images["filename"].filled("")]

        science_images["original"] = [os.path.join(raw_path, image)
                                      for image in rawnames]

        columns = ["original", "master_bias", "master_dark", "master_flat",
                   "reduced", "catalog"]
        for col in columns:
            basenames = [os.path.basename(image)
                         for image in science_images[col].filled("")]
            dirnames = [os.path.basename(os.path.dirname(image))
                        for image in science_images[col].filled("")]
            science_images[col] = [os.path.join(impath, image)
                                   for impath, image in
                                   zip(dirnames, basenames)]

        science_images.sort(["original", "reduced", "catalog"])

        log_text = ["General remarks:", ""]

        mask = ~science_images["invalid"].mask
        invalid_astrometry = ["astrometry" == value
                              for value in science_images["invalid"][mask]]
        if invalid_astrometry:
            log_text += [("Astrometry could not be determined for some of "
                          "the reduced images. Only photometry is provied "
                          "for those images."), ""]

        invalid_photometry = ["photometry" == value
                              for value in science_images["invalid"][mask]]
        if invalid_photometry:
            log_text += [("Valid photometry could not be extracted for some "
                          "of the reduced images. Those images contain no "
                          "'num_stars' and no 'fwhm' determination in the "
                          "list below."), ""]

        log_text += ["", "Reduced images:", ""]

        columns += ["num_stars", "fwhm"]
        log_text += science_images[columns].pformat(max_width=-1, max_lines=-1)

        filename = os.path.join(cat_path, "{}.txt".format(red_path))

        if os.path.isfile(filename):
            log_msg = "Adding data to an existing file '%s'"
            logger.warning(log_msg, filename)

        with open(filename, "a") as logf:
            logf.write("\n".join(log_text))


def create_permanent_logs(imdata, suffix):
    """
    Create permanent log file for imdata and link it to home path
    """

    log_name = None
    raw_images = imdata["filename"][~imdata["filename"].mask]

    default = Time.now().iso.split()[0].replace("-", "")
    destination = set(os.path.dirname(image) for image in raw_images)
    for path in destination:

        red_path = path.replace("rawdata", "reddata")

        try:
            os.makedirs(red_path)
        except OSError:
            pass

        today = os.path.basename(red_path)
        if today.isdigit():
            now = Time(int(today), format="jd").iso.split()[0].replace("-", "")
        else:
            now = default

        log_name = os.path.join(red_path, "{}_{}.txt".format(now, suffix))

        file_text = imdata.pformat(max_width=-1, max_lines=-1)
        with open(log_name, "a") as logf:
            logf.write("\n".join(file_text))
            logf.write("\n")

    link_name = os.path.join(os.path.expanduser("~"), "{}.txt".format(suffix))
    try:
        os.remove(link_name)
    except OSError:
        pass

    if log_name:
        os.symlink(log_name, link_name)


def read_input(input_paths, recursive, files_filter):
    """
    Transform input arguments into a list of images
    """

    input_images = set(system.path_flatten(input_paths, recursive,
                                           files_filter=files_filter))

    if not input_images:
        if files_filter:
            err_msg = "No input images were considered science images. "
            err_msg += "Filter '{}' ".format(files_filter)
            err_msg += "is being used to filter images."
        else:
            err_msg = "No images found to be reduced"

        raise MissingError(err_msg)

    return input_images


def find_destination(input_images, proposals):
    """
    Find invalid keywords or destination path for the input images
    belonging to proposals
    """

    try:
        valid = table.Table([proposals.keys()], names=["propcode"])
        propcode = input_images[~input_images["propcode"].mask]
        input_images = table.join(valid, propcode)
    except ValueError:
        if not valid:
            err_msg = "No valid proposal codes found in database"
        else:
            err_msg = "None of the input images have valid proposal codes"
        raise ppl.PipelineError(err_msg)

    if not input_images:
        err_msg = "No input images belonging to selected proposals could be "
        err_msg += "found in input path/s"
        raise MissingError(err_msg)

    options = ["reduced", "imgqs", "propcode", "date-obs or time-obs"]
    root_dest = config_option(__name__, "destination", ignore_missing=False)

    num_images = len(input_images)
    plural_images = "" if num_images == 1 else "s"
    num_proposals = len(proposals)
    plural_proposals = "" if num_proposals == 1 else "s"
    log_msg = "Selecting destination path for {} input image{} among "
    log_msg += "{} possible proposal{}..."

    logger.info(log_msg.format(num_images, plural_images,
                               num_proposals, plural_proposals))

    mask = []
    invalid_column = []
    destination_column = []
    for image in LogOrProgressBar(input_images.filled()):

        reduced = None if image["original_ids"] else False
        imgqs = match("[02][0-9][029][0-9]", image["imgqs"])
        propcode = image["propcode"]
        user_dir = proposals.get(propcode)
        try:
            real_date = Time(image["date-obs"]+" "+image["time-obs"])
            obs_date = real_date
            if obs_date.datetime.hour < 12:
                obs_date -= 1 * day
        except ValueError:
            obs_date = None

        parameters = (reduced, imgqs, user_dir, obs_date)
        invalid = ",".join(key for key, value in zip(options, parameters)
                           if value is None)

        if invalid:
            log_msg = "Image '{}' discarded due to invalid '{}'"
            logger.debug(log_msg.format(image["filename"], invalid))
            mask += [False]
            invalid_column += [invalid]
            destination_column += [""]

        else:
            log_msg = "Selecting '{}' for publication"
            logger.debug(log_msg.format(image["filename"]))

            prop_dir = "proposal_{}".format(propcode)
            date_dir = obs_date.datetime.strftime("%Y%m%d_rawdata")

            mask += [True]
            invalid_column += [""]
            destination_column += [os.path.join(root_dest, prop_dir,
                                                user_dir, date_dir)]

    input_images["invalid"] = invalid_column
    input_images["destination"] = destination_column

    input_images["invalid"].mask = mask
    input_images["destination"].mask = ~input_images["invalid"].mask

    valid_images = input_images[~input_images["destination"].mask]
    rejected_images = input_images[input_images["destination"].mask]

    num_input = len(input_images)
    num_valid = len(valid_images)
    num_rejected = len(rejected_images)
    if num_valid + num_rejected != num_input:
        err_msg = "Number of selected images ({}) ".format(num_valid)
        err_msg += "and number of discarded images ({}) ".format(num_rejected)
        err_msg += "do not match number of input images ({})".format(num_input)
        raise ppl.PipelineError(err_msg)

    if not valid_images:
        err_msg = "Input image was not selected for publication"
        if num_input > 1:
            err_msg = "None of the {} ".format(num_input)
            err_msg += "input images were selected for publication"

        raise MissingError(err_msg)

    return valid_images, rejected_images


def select_calibration(science, calib, verbose=False):
    """
    Retrieve calibration images matching each science image
    """

    num_science = len(science)
    num_calib = len(calib[~calib["obstype"].mask])
    plural_science = "" if num_science == 1 else "s"
    plural_calib = "" if num_calib == 1 else "s"
    log_msg = "Searching suitable calibration images for {} science image{} "
    log_msg += "among {} candidate{}..."
    logger.info(log_msg.format(num_science, plural_science,
                               num_calib, plural_calib))

    # Splitting calib in obstypes is a hack to speed up the process

    bias = calib[calib["obstype"] == "Bias"]
    darks = calib[calib["obstype"] == "Dark"]
    flats = calib[calib["obstype"] == "Flat"]
    for image in LogOrProgressBar(science):
        if not image["destination"]:
            continue    # Updating science table in place requires this hack

        select_bias(image, bias)
        select_darks(image, bias, darks)
        select_flats(image, bias, darks, flats)

    # Update invalid column and tables mask

    logger.debug("Merging selected calibration images and science images...")

    invalid_column = []
    for image in science:
        invalid = image["invalid"].split(",") if image["invalid"] else []

        for imtype in ("bias", "darks", "flats"):
            if len(image[imtype]) <= 2:
                invalid += [imtype]

        invalid_column += [",".join(invalid) if invalid else ""]

    science["invalid"] = invalid_column

    for col in "bias", "darks", "flats", "invalid":
        science[col].mask |= (science[col] == "||") | (science[col] == "")
        bias[col].mask |= (bias[col] == "||") | (bias[col] == "")
        darks[col].mask |= (darks[col] == "||") | (darks[col] == "")
        flats[col].mask |= (flats[col] == "||") | (flats[col] == "")

    # Selected calibration images are those having "bias", "darks" and "flats"

    has_bias_darks = ~science["bias"].mask & ~science["darks"].mask
    selected_science = science[has_bias_darks & ~science["flats"].mask]

    flats_ids = [int(num) for ids in selected_science["flats"]
                 for num in ids.strip("|").split("|")]
    selected_flats = flats[[num in set(flats_ids) for num in flats["id"]]]

    darks_ids = [int(num) for ids in selected_science["darks"]
                 for num in ids.strip("|").split("|")]
    darks_ids += [int(num) for ids in selected_flats["darks"]
                  for num in ids.strip("|").split("|")]
    selected_darks = darks[[num in set(darks_ids) for num in darks["id"]]]

    bias_ids = [int(num) for ids in selected_darks["bias"]
                for num in ids.strip("|").split("|")]
    bias_ids += [int(num) for ids in selected_flats["bias"]
                 for num in ids.strip("|").split("|")]
    bias_ids += [int(num) for ids in selected_science["bias"]
                 for num in ids.strip("|").split("|")]
    selected_bias = bias[[num in set(bias_ids) for num in bias["id"]]]

    return table.vstack([science, selected_bias,
                         selected_darks, selected_flats])


def copy_raw_images(images, verbose=False):
    """
    Copy original images to destination paths
    """

    copy_paths = set()
    ids = images["id"].astype(str)
    destinations = set(images["destination"][~images["destination"].mask])

    num_dest = len(destinations)
    plural_dest = "" if num_dest == 1 else "s"
    log_msg = "Grouping required calibration images for {} "
    log_msg += "destination path{}..."
    logger.info(log_msg.format(num_dest, plural_dest))

    for dest in LogOrProgressBar(destinations):

        # Select destination for science images

        mask = images["destination"] == [dest]
        science_ids = images["id"][mask].astype(str).tolist()

        # Select destination for calibration images required by science images

        has_bias = ~images["bias"].mask
        has_darks = ~images["darks"].mask
        has_flats = ~images["flats"].mask

        bias_ids = [num for calib_ids in images["bias"][mask & has_bias]
                    for num in calib_ids.strip("|").split("|")]
        dark_ids = [num for calib_ids in images["darks"][mask & has_darks]
                    for num in calib_ids.strip("|").split("|")]
        flat_ids = [num for calib_ids in images["flats"][mask & has_flats]
                    for num in calib_ids.strip("|").split("|")]

        # Select destination for calibration images required by flats

        flat_mask = [num in flat_ids for num in ids]

        bias_ids += [num for calib_ids in images["bias"][flat_mask]
                     for num in calib_ids.strip("|").split("|")]
        dark_ids += [num for calib_ids in images["darks"][flat_mask]
                     for num in calib_ids.strip("|").split("|")]

        # Select destination for calibration images required by darks

        dark_mask = [num in dark_ids for num in ids]

        bias_ids += [num for calib_ids in images["bias"][dark_mask]
                     for num in calib_ids.strip("|").split("|")]

        # Join everything

        requires = set(science_ids + flat_ids + dark_ids + bias_ids)
        filenames = images["filename"][[num in requires for num in ids]]
        copy_paths |= set((orig, dest) for orig in filenames)

        logger.debug("Selecting %s images for '%s'", len(filenames), dest)

    num_images = len(set(orig for orig, dest in copy_paths))
    num_copies = len(copy_paths)
    images_plural = "" if num_copies == 1 else "s"
    log_msg = "Copying {} image{} ({} different) to {} destination path{}..."
    logger.info(log_msg.format(num_copies, images_plural, num_images,
                               num_dest, plural_dest))

    copied_paths = set()
    copied_images = set()
    for orig, dest in LogOrProgressBar(copy_paths):
        compressed_image = os.path.join(dest, os.path.basename(orig))
        decompressed_image = os.path.splitext(compressed_image)[0]

        if not is_newer(orig, decompressed_image):
            err_msg = "Destination image '{}' already exists. Skipping..."
            logger.debug(err_msg.format(decompressed_image))
            copied_paths |= {dest}
            copied_images |= {orig}
            continue

        elif is_newer(orig, compressed_image):
            if os.path.isfile(compressed_image):
                os.chmod(compressed_image, 0644)
            system.new_copy(orig, compressed_image, makedirs=True,
                            overwrite=True)

        if os.path.isfile(decompressed_image):
            system.remove(decompressed_image)

        system.decompress(compressed_image)
        copied_paths |= {dest}
        copied_images |= {orig}

    num_done = len(copied_images)

    if num_done != num_images:
        err_msg = "Number of selected images ({}) does not match number of "
        err_msg += "copied images ({})"
        raise ppl.PipelineError(err_msg.format(num_images, num_done))

    else:
        log_msg = "Successfully copied {} different images to {} paths"
        logger.debug(log_msg.format(num_done, len(copied_paths)))


def clean_images(images, verbose=False):
    """
    Subtract bias, darks and flats to science images
    """

    invalid_images = images[~images["original_ids"].mask]
    if invalid_images:
        err_msg = "{} images in 'clean_images' cannot be processed because "
        err_msg += "they are not raw images"
        raise ppl.PipelineError(err_msg.format(len(invalid_images)))

    valid_science = (images["obstype"] == "Science") & images["invalid"].mask
    science_images = images[valid_science]

    if not science_images:
        num_images = len(images)
        err_msg = "No valid calibration images found to reduce "
        if num_images > 1:
            err_msg += "any of the {} science images".format(num_images)
        else:
            err_msg += "the science image"

        raise MissingError(err_msg)

    num_images = len(science_images)
    images_plural = "" if num_images == 1 else "s"
    log_msg = "Cleaning {} science image{}..."
    logger.info(log_msg.format(num_images, images_plural))

    # Ensure that reddata directories are empty before starting reduction

    for rawdir in set(science_images["destination"]):
        redpath = rawdir.replace("_rawdata", "_reddata")

        if os.path.isdir(str(redpath)):
            rmtree(redpath)

        try:
            os.makedirs(redpath)
        except OSError:
            pass

    raw_red = defaultdict(list)
    reduced_images = {}
    ids = images["id"].astype(str)
    for science in LogOrProgressBar(science_images):

        logger.debug("Cleaning image: '{}'".format(science["filename"]))

        origin = science["destination"]
        redpath = origin.replace("_rawdata", "_reddata")

        requires = science["bias"].strip("|").split("|")
        bias = images[[num in requires for num in ids]]

        master_bias = bias_combine(reduced_images, bias, origin, redpath)

        requires = science["darks"].strip("|").split("|")
        darks = images[[num in requires for num in ids]]

        master_dark = dark_combine(reduced_images, darks, images,
                                   origin, redpath)

        requires = science["flats"].strip("|").split("|")
        flats = images[[num in requires for num in ids]]
        master_flat = flat_combine(reduced_images, flats, images,
                                   origin, redpath)

        raw_basename = os.path.basename(science["filename"])
        if raw_basename.endswith(".gz"):
            raw_basename = os.path.splitext(raw_basename)[0]

        raw_name = os.path.join(origin, raw_basename)

        raw_data = ccdproc.CCDData.read(raw_name, unit="adu")
        bias_data = ccdproc.CCDData.read(master_bias, unit="adu")
        dark_data = ccdproc.CCDData.read(master_dark, unit="adu")
        flat_data = ccdproc.CCDData.read(master_flat, unit="adu")

        sub_data = ccdred(raw_data, master_bias=bias_data, add_keyword=False,
                          dark_frame=dark_data, master_flat=flat_data,
                          exposure_key="EXPTIME", exposure_unit=second)

        red_basename = raw_basename.replace("imr.fits", "iml.fits")
        red_name = os.path.join(redpath, red_basename)

        # Store clean image at 32 bits per pixel
        # Unfortunately, astropy.io.fits.ImageHDU.scale is not working
        # for float values in astropy 1.3

        sub_data.data = sub_data.data.astype("float32")
        sub_data.header.pop("BZERO")
        sub_data.write(red_name, overwrite=True)

        used_ids = set([science["id"].astype(str)])
        used_ids |= set(reduced_images[master_bias].strip("|").split("|"))
        used_ids |= set(reduced_images[master_dark].strip("|").split("|"))
        used_ids |= set(reduced_images[master_flat].strip("|").split("|"))

        reduced_images[red_name] = "|{}|".format("|".join(sorted(used_ids)))
        raw_red["reduced"] += [red_name]
        raw_red["filename"] += [science["filename"]]
        raw_red["master_bias"] += [master_bias]
        raw_red["master_dark"] += [master_dark]
        raw_red["master_flat"] += [master_flat]

    science_images.remove_column("original_ids")

    reduced = table.Table(raw_red)
    reduced_science = table.join(reduced, science_images, join_type="left")

    originals = [reduced_images.keys(), reduced_images.values()]
    original_ids = table.Table(originals, names=["reduced", "original_ids"])

    return table.join(original_ids, reduced_science, join_type="left")


def astrophot(images, verbose=False):
    """
    Compute astrometry and photometry on a set of images
    """

    invalid_images = images[~images["invalid"].mask]
    if invalid_images:
        err_msg = "{} images in 'astrophot' cannot be used due to 'invalid' "
        err_msg += "values"
        raise ppl.PipelineError(err_msg.format(len(invalid_images)))

    valid_science = images["obstype"].data == "Science"
    valid_science *= ~images["obstype"].mask
    valid_science *= ~images["reduced"].mask

    science_images = images[valid_science]

    num_images = len(science_images)
    images_plural = "" if num_images == 1 else "s"
    log_msg = "Computing astrometry for {} science image{}..."
    logger.info(log_msg.format(num_images, images_plural))

    astro_list = []
    invalid_list = []
    for science in LogOrProgressBar(science_images):

        filename = science["reduced"]
        astro_image = filename.replace("iml.fits", "imc.fits")

        logger.debug("Finding solution for: {}".format(filename))

        try:
            # Using imgqs configuration file changes nothing

            astrometry.solve(filename, ra=science["ra"], dec=science["dec"],
                             use_sextractor=True, new_fits=astro_image)

            astro_list += [astro_image]
            invalid_list += [""]

            os.remove(filename)

            log_msg = "Astrometry image '{}' created"
            logger.debug(log_msg.format(astro_image))

        except (astrometry.AstrometryError,
                astrometry.AstrometryWarning) as err:

            logger.debug(err)

            astro_list += [filename]
            invalid_list += ["astrometry"]

    science_images["astrometry"] = astro_list
    science_images["invalid"] = invalid_list

    log_msg = "Computing photometry for {} science image{}..."
    logger.info(log_msg.format(num_images, images_plural))

    # Ensure that catalog directories are empty before starting photometry

    for reddir in set(os.path.dirname(_image)
                      for _image in science_images["reduced"]):
        catpath = reddir.replace("_reddata", "_catalog")

        if os.path.isdir(str(catpath)):
            rmtree(catpath)

        try:
            os.makedirs(catpath)
        except OSError:
            pass

    cats = defaultdict(list)
    config_path = os.path.dirname(os.path.dirname(__file__))
    default_config = os.path.join(config_path, "config", "sextractor.conf")
    default_keys = os.path.join(config_path, "config", "sextractor.keys")
    stdout = stderr = None if verbose else open(os.devnull, "w")
    for science in LogOrProgressBar(science_images):
        config_file = default_config
        config_keys = default_keys

        astro_image = science["astrometry"]
        cat_dir = astro_image.replace("_reddata", "_catalog")
        phot_cat = cat_dir.replace(".fits", "_trl.dat")
        invalid = science["invalid"] if science["invalid"] else ""

        cats["reduced"] += [science["reduced"]]
        cats["astrometry"] += [astro_image]

        try:
            os.makedirs(os.path.dirname(phot_cat))
        except OSError:
            pass

        try:
            if 3600 < science["focuspos"] < 7000:
                config_file = config_file.replace("sextractor.conf",
                                                  "sextractor.defocused.conf")

            if "iml.fits" in astro_image:
                config_keys = config_keys.replace("sextractor.keys",
                                                  "sextractor.no_astro.keys")

            photometry.sextractor(astro_image, catalog_name=phot_cat,
                                  config=config_file, params=config_keys,
                                  stdout=stdout, stderr=stderr)

            statistics = imgqs_ppl.statistics(astro_image, phot_cat)

            try:
                fwhm = statistics["FWHM"].arcsec
            except KeyError:
                fwhm = 0

            nstars = statistics["NSTARS"]

            cats["catalog"] += [phot_cat]
            cats["invalid"] += [invalid or ""]
            cats["fwhm"] += [fwhm]
            cats["num_stars"] += [nstars]

        except Exception as err:

            logger.debug(err)

            try:
                os.remove(phot_cat)
            except OSError:
                pass

            cats["catalog"] += [""]
            cats["invalid"] += [invalid + ",photometry"
                                if invalid else "photometry"]
            cats["fwhm"] += [0]
            cats["num_stars"] += [0]

    catalogs = table.Table(cats, masked=True)
    catalogs["catalog"].mask = catalogs["catalog"] == ""
    catalogs["invalid"].mask = catalogs["invalid"] == ""
    catalogs["num_stars"].mask = catalogs["num_stars"] == 0
    catalogs["fwhm"].mask = catalogs["fwhm"] == 0

    images.remove_column("invalid")

    reduced_catalogs = table.join(images, catalogs, join_type="left")

    astro_mask = reduced_catalogs["astrometry"].mask
    merged_reduced = [row["reduced"] if mask else row["astrometry"]
                      for row, mask in zip(reduced_catalogs, astro_mask)]
    reduced_catalogs["reduced"] = merged_reduced

    reduced_catalogs.remove_column("astrometry")

    return reduced_catalogs


def copy_red_images(images, verbose=False):
    """
    Copy reduced images to permanent storage paths
    """

    # Use those images with an original image (with filename) to find
    # destination path and a reduced image to store (with original_ids)

    mask = ~images["filename"].mask & ~images["original_ids"].mask
    valid_images = images[mask]

    num_images = len(valid_images)
    images_plural = " to its" if num_images == 1 else "s to their"
    log_msg = "Copying {} science image{} permanent path..."
    logger.info(log_msg.format(num_images, images_plural))

    copied_paths = set()
    could_be_removed = set()
    for reduced in LogOrProgressBar(valid_images):
        raw_path = os.path.dirname(reduced["filename"])
        dirname = raw_path.replace("rawdata", "reddata")

        image_types = ("reduced", "master_bias", "master_dark", "master_flat",
                       "catalog")
        for imtype in image_types:
            if "master" in imtype:
                calibration = images[images["reduced"] == reduced[imtype]]
                original = calibration["reduced"][0]
                original_ids = calibration["original_ids"][0]
            else:
                original = reduced[imtype]
                original_ids = reduced["original_ids"]

            # Catalog could be missing due to an error in the photometry

            if not original:
                continue

            destination = os.path.join(dirname, os.path.basename(original))
            compressed = destination + ".gz"

            if not is_newer(original, compressed):
                err_msg = "Destination image '{}' already exists. Skipping..."
                logger.debug(err_msg.format(compressed))
                continue

            elif (imtype == "catalog") and not is_newer(original, destination):
                err_msg = "Destination catalog '{}' already exists. Skipped..."
                logger.debug(err_msg.format(destination))
                continue

            elif is_newer(original, destination):
                if os.path.isfile(destination):
                    imgqs_ppl.delete.run(destination)

                if os.path.isfile(compressed):
                    imgqs_ppl.delete.run(compressed)

                system.new_copy(original, destination, makedirs=True)

            if imtype == "catalog":
                catstore.run(destination)
            else:
                imgqs_ppl.store.run(destination, dependencies=original_ids)

            temporary_file = ""
            if "imc.fits" in destination:
                temporary_file = destination.replace("imc.fits", "imt.fits")
            elif "iml.fits" in destination:
                temporary_file = destination.replace("iml.fits", "imt.fits")

            if os.path.isfile(temporary_file):
                imgqs_ppl.delete.run(temporary_file)
            elif os.path.isfile(temporary_file+".gz"):
                imgqs_ppl.delete.run(temporary_file+".gz")

        if os.path.samefile(dirname, os.path.dirname(reduced["reduced"])):
            could_be_removed |= {dirname}

        copied_paths |= {dirname}

    if could_be_removed:
        err_msg = "In case of error, some of the reduced images could be "
        err_msg += "removed because reduction and storage paths are the same:"
        logger.warning("{} {}".format(err_msg, ",".join(could_be_removed)))

    log_msg = "Successfully copied {} reduced images to {} paths"
    logger.debug(log_msg.format(num_images, len(copied_paths)))


def create_user_logs(rawdata, reddata, sky_coverage):
    """
    Create log files to be submitted to observers
    """

    try:
        logger.info("Creating log files...")

        seeing = None
        coverage = None
        create_rawdata_log(rawdata, seeing, coverage)

        if reddata:
            create_reddata_log(reddata)

    except Exception as err:
        logger.exception(err)
        return mail_error()


def mail_statistics(rejected, published, reduced, istest):
    """
    Create mail message with statistics of the reduction
    """

    logger.info("Creating mail message...")

    try:
        # Create required tables from input data

        columns = ["filename", "invalid", "imgqs", "propcode", "object",
                   "filter", "date-obs", "time-obs"]

        empty_table = rejected[:0].copy()
        empty_table["reduced"] = []
        rejected = empty_table if rejected is None else rejected
        published = empty_table if published is None else published
        reduced = empty_table if reduced is None else reduced
        astrored = reduced[reduced["invalid"].mask] if reduced else empty_table

        rawdata = table.vstack([rejected, published])
        requested = rawdata["filename"][rawdata["obstype"] == "Science"]

        # Create internal reduction logs

        if not istest:
            create_permanent_logs(rawdata, "rawdata")
            create_permanent_logs(reduced, "reddata")

        # Prepare a table with the statistics of the reduction

        msg = ("A new CAT reduction was completed with the "
               "following results:\n\n{}\n\n")

        imtypes = ["Requested", "Published", "Reduced", "Valid catalog"]
        statistics = {"Image type": imtypes}

        published_mask = ~published["obstype"].mask
        published_mask *= published["obstype"] == "Science"
        reduced_mask = ~reduced["obstype"].mask
        reduced_mask *= reduced["obstype"] == "Science"
        astro_mask = ~astrored["obstype"].mask
        astro_mask *= astrored["obstype"] == "Science"

        imtables = (requested, published[published_mask],
                    reduced[reduced_mask], astrored[astro_mask])
        statistics["Science images"] = [len(step) for step in imtables]

        num_paths = [len(set(os.path.dirname(image) for image in requested))]
        num_paths += [len(set(published["destination"][published_mask]))]
        num_paths += [len(set(os.path.dirname(image)
                              for image in reduced["reduced"][reduced_mask]))]
        num_paths += [len(set(os.path.dirname(image)
                              for image in astrored["reduced"][astro_mask]))]

        statistics["Paths"] = num_paths

        statistics_table = table.Table(statistics)
        stat_cols = ("Image type", "Paths", "Science images")
        stats_table = statistics_table[stat_cols].pformat(max_width=-1,
                                                          max_lines=-1)

        mail_msg = msg.format("\n".join(stats_table))

        # Images dropped during the reduction process should also be reported

        non_reduced = published[~published["invalid"].mask]

        non_astrometry = reduced[~reduced["invalid"].mask]

        invalid = table.vstack([rejected[columns], non_reduced[columns],
                                non_astrometry[columns]])
        invalid.sort(columns)
        invalid_table = invalid.pformat(max_width=-1, max_lines=-1)

        host = gethostname()
        num_invalid = len(invalid)
        if num_invalid == 1:
            msg = ("However, the reduction could not be fully completed "
                   "for one science image:\n\n{}\n\nPlease, see files "
                   "'~/rawdata.txt' and '~/reddata.txt' in '{}' to retrieve "
                   "further information on the invalid image.")
            mail_msg += msg.format("\n".join(invalid_table), host)
        elif num_invalid > 1:
            msg = ("However, the reduction could not be fully completed "
                   "for {} science images:\n\n{}\n\nPlease, see files "
                   "'~/rawdata.txt' and '~/reddata.txt' in '{}' to retrieve "
                   "further information on the invalid images.")
            mail_msg += msg.format(num_invalid, "\n".join(invalid_table), host)
        else:
            path_plural = "" if num_paths[-1] == 1 else "s"
            root_dest = config_option(__name__, "destination")
            msg = "You can now proceed to review the log file{} in '{}'."
            mail_msg += msg.format(path_plural, root_dest)

        return base_msg.format(mail_msg)

    except Exception as err:
        logger.warning(err)
        mail_msg = ("Reduction finished successfully, but statistics mail "
                    "could not be created due to error:\n\n{}\n\n"
                    "Please, report this incident and proceed to review the "
                    "log file/s.")
        return base_msg.format(mail_msg.format(err))


def mail_error(error=None):
    """
    Create mail message reporting an error
    """

    if isinstance(error, MissingError):
        mail_msg = ("A new CAT reduction has been requested, but no images "
                    "were found to be reduced. ")
    else:
        mail_msg = "An error occurred when trying to reduce observations. "

    mail_msg += "The reported error is:\n\n{}\n\n".format(format_exc())
    mail_msg += "Please, see the full reduction details at the log file"

    try:
        log_name = config_option("Logging", "file_name")
        host = gethostname()
        mail_msg += " '{}' in '{}'.".format(log_name, host)
    except Exception:
        mail_msg += "."

    return base_msg.format(mail_msg)


def run(input_paths, recursive=False, files_filter=None, propcode=None,
        verbose=False, istest=False):
    """
    Pipeline execution function
    """

    title = mail_msg = ""
    input_images = rejected_images = selected_images = reduced_images = None

    try:

        input_images = read_input(input_paths, recursive, files_filter)

        db = TJODB(read_default_group="Database")

        images = db.get_images(input_images)

        active = db.proposals(propcode)
        valid_images, rejected_images = find_destination(images, active)

        calibration_images = db.find_calibration(valid_images)

        selected_images = select_calibration(valid_images, calibration_images)

        # coverage = None
        coverage = db.cloud_coverage(selected_images)

        copy_raw_images(selected_images)

#        previous_images = db.find_reduced(selected_images)

        reduced_images = clean_images(selected_images)

        reduced_images = astrophot(reduced_images, verbose)

        if istest:
            log_msg = "Execution run in test mode. "
            log_msg += "Images NOT stored in their permanent paths..."
            logger.warning(log_msg)
        else:
            copy_red_images(reduced_images)

        title = "Error creating log files when reducing CAT images"
        mail_msg = create_user_logs(selected_images, reduced_images, coverage)

        if not mail_msg:
            title = "New CAT images reduced"
            mail_msg = mail_statistics(rejected_images, selected_images,
                                       reduced_images, istest)

    except MissingError as err:
        logger.warning(err)

        if selected_images:
            # No images reduced
            title = "Error creating log files when reducing CAT images"
            mail_msg = create_user_logs(selected_images, reduced_images,
                                        coverage)

            if not mail_msg:
                title = "New CAT images selected for publication"
                mail_msg = mail_statistics(rejected_images, selected_images,
                                           reduced_images, istest)

        elif rejected_images:
            # No images copied
            title = "New CAT reduction requested"
            mail_msg = mail_statistics(rejected_images, selected_images,
                                       reduced_images, istest)
        else:
            # No images found for reduction
            title = "New CAT reduction requested"
            mail_msg = mail_error(err)

    except Exception as err:
        logger.exception(err)

        if selected_images and not istest:
            for dest_path in set(selected_images["destination"].filled("")):
                if not dest_path:
                    continue

                red_path = dest_path.replace("_rawdata", "_reddata")
                cat_path = dest_path.replace("_rawdata", "_catalog")
                for prop_path in (dest_path, red_path, cat_path):
                    if os.path.isdir(str(prop_path)):
                        rmtree(prop_path)

        elif selected_images:
            logger.debug("Copied and reduced images not removed!")

        title = "Error reducing CAT images"
        mail_msg = mail_error(err)

        raise

    finally:
        if title and mail_msg:
            send_mail(title=title, message=mail_msg)


def main():
    """
    Entry point
    """

    global logger

    try:
        origin_path = config_option(__name__, "origin")
        all_paths = [path for path in os.listdir(origin_path)
                     if path.isdigit()]
        default_path = os.path.join(origin_path, sorted(all_paths)[-1])
    except (TypeError, IndexError):
        default_path = os.path.expanduser("~")

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-r", "--recursive", action="store_true",
                        help="recursively enter in directories")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log debug messages")

    help_msg = "execute entire reduction, but do not store images "
    help_msg += "to their permanent paths"
    parser.add_argument("-t", "--test", action="store_true", help=help_msg)

    help_msg = "reduce only proposal/s having code PROPCODE "
    help_msg += "(it can be repeated several times for various proposals)"
    parser.add_argument("-p", "--propcode", action="append", help=help_msg)

    help_msg = "when path is a directory, reduce only images having FILTER "
    help_msg += "in name (default: '%(default)s')"
    parser.add_argument("-f", "--filter", default="imr.fits", help=help_msg)

    help_msg = "file/s or directory/ies with FITS image/s "
    help_msg += "(default: '%(default)s')"
    parser.add_argument("path", nargs="*", default=default_path, help=help_msg)

    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    run(args.path, args.recursive, args.filter,
        args.propcode, args.verbose, args.test)


if __name__ == "__main__":
    main()
