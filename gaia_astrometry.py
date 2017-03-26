"""
Pipeline used to reduce TJO images of the proposal entitled
'First astrometric measurements of stellar activity with spectroscopy'
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)
import os.path
import argparse
from shutil import rmtree
from tempfile import mkdtemp, NamedTemporaryFile

import numpy as np
from astropy import log as astrolog
from astropy.units import second
import ccdproc

import icat
from icat.tools.logs import get_logger
from icat.tools.photometry import sextractor
from icat.tools.ccdproc_database import create_database
from icat.tools.astrometry import solve, AstrometryError, AstrometryWarning


astrolog.setLevel("ERROR")
logger = get_logger(__name__)


def valid_values(image_list, median_range, std_max):
    """
    Ensure that images have values within a certain range
    """

    valid_images = []
    for name in image_list:
        image = ccdproc.CCDData.read(name, unit="adu")

        std = np.std(image.data)
        median = np.median(image.data)

        if median_range[0] < median < median_range[1] and std < std_max:
            valid_images += [name]

    if image_list and not valid_images:
        raise ValueError("No images left after checking for image values")

    return valid_images


def filter_bias(database, image):
    """
    Select bias images suitable to be used in a given image
    """

    ref_image = database[database["file"] == image]

    mask = database["obstype"] == "Bias"
    mask &= database["instrume"] == ref_image["instrume"]
    mask &= database["filter"] == "C"
    mask &= database["jd"].astype(int) == int(ref_image["jd"])
    mask &= database["naxis1"].astype(int) == int(ref_image["naxis1"])
    mask &= database["naxis2"].astype(int) == int(ref_image["naxis2"])
    mask &= database["exptime"].astype(int) == 0
    mask &= abs(database["camtemp"] - ref_image["camtemp"]) < 3

    return database["file"][mask].tolist()


def bias_combine(master_frames, bias_list, temp_path=None):
    """
    Create the master bias, if required
    """

    bias_list = valid_values(bias_list, median_range=(1700, 2200), std_max=150)

    biasname = [name for name, combolist in master_frames.iteritems()
                if bias_list == combolist]

    if biasname:
        master_bias = biasname[0]
        debug_string = "Using bias '{}'"
        logger.debug(debug_string.format(master_bias))

    else:
        with NamedTemporaryFile(suffix="_bias.fits", dir=temp_path,
                                delete=False) as biasfile:
            master_bias = biasfile.name

        log_msg = "Using {} bias to create master bias '{}'"
        logger.debug(log_msg.format(len(bias_list), master_bias))

        ccdproc.combine(bias_list, output_file=master_bias, mem_limit=1e9,
                        method="median", unit="adu", dtype="float32")

        master_frames[master_bias] = bias_list

    return master_bias


def filter_darks(database, image):
    """
    Select dark images suitable to be used in a given image
    """

    ref_image = database[database["file"] == image]

    mask = database["obstype"] == "Dark"
    mask &= database["instrume"] == ref_image["instrume"]
    mask &= database["filter"] == "C"
    mask &= database["jd"].astype(int) == int(ref_image["jd"])
    mask &= database["naxis1"].astype(int) == int(ref_image["naxis1"])
    mask &= database["naxis2"].astype(int) == int(ref_image["naxis2"])
    mask &= database["exptime"].astype(int) >= int(ref_image["exptime"])
    mask &= abs(database["camtemp"] - ref_image["camtemp"]) < 3

    return database["file"][mask].tolist()


def dark_combine(master_frames, dark_list, database, temp_path=None):
    """
    Create the master dark, if required
    """

    dark_list = valid_values(dark_list, median_range=(1700, 2200), std_max=150)

    darkname = [name for name, combolist in master_frames.iteritems()
                if dark_list == combolist]

    if darkname:
        master_dark = darkname[0]
        debug_string = "Using dark '{}'"
        logger.debug(debug_string.format(master_dark))

    else:
        with NamedTemporaryFile(suffix="_dark.fits", dir=temp_path,
                                delete=False) as darkfile:
            master_dark = darkfile.name

        log_msg = "Cleaning {} darks to create master dark '{}'"
        logger.debug(log_msg.format(len(dark_list), master_dark))

        scaling = []
        clean_images = []
        for image in dark_list:
            bias_list = filter_bias(database, image)

            master_bias = bias_combine(master_frames, bias_list, temp_path)

            raw_data = ccdproc.CCDData.read(image, unit="adu")
            master_data = ccdproc.CCDData.read(master_bias, unit="adu")

            sub_data = ccdproc.subtract_bias(raw_data, master_data,
                                             add_keyword=False)

            clean_images += [sub_data]

            read = ccdproc.CCDData.read
            scaling += [1.0 / read(image, unit="adu").header["EXPTIME"]]

        dark = ccdproc.combine(clean_images, mem_limit=1e9, scale=scaling,
                               method="median", unit="adu", dtype="float32")
        dark.header["EXPTIME"] = 1.0
        dark.write(master_dark)

        master_frames[master_dark] = dark_list

    return master_dark


def filter_flats(database, image):
    """
    Select flat images suitable to be used in a given image
    """

    ref_image = database[database["file"] == image]

    mask = database["obstype"] == "Flat"
    mask &= database["instrume"] == ref_image["instrume"]
    mask &= database["filter"] == ref_image["filter"]
    mask &= database["jd"].astype(int) == int(ref_image["jd"])
    mask &= database["naxis1"].astype(int) == int(ref_image["naxis1"])
    mask &= database["naxis2"].astype(int) == int(ref_image["naxis2"])
    mask &= abs(database["camtemp"] - ref_image["camtemp"]) < 3

    return database["file"][mask].tolist()


def flat_combine(master_frames, flat_list, database, temp_path=None):
    """
    Create the master flat, if required
    """

    flat_list = valid_values(flat_list, median_range=(20000, 40000),
                             std_max=1500)

    flatname = [name for name, combolist in master_frames.iteritems()
                if flat_list == combolist]

    if flatname:
        master_flat = flatname[0]
        debug_string = "Using flat '{}'"
        logger.debug(debug_string.format(master_flat))

    else:
        with NamedTemporaryFile(suffix="_flat.fits", dir=temp_path,
                                delete=False) as flatfile:
            master_flat = flatfile.name

        log_msg = "Cleaning {} flats to create master flat '{}'"
        logger.debug(log_msg.format(len(flat_list), master_flat))

        scaling = []
        clean_images = []
        for image in flat_list:
            bias_list = filter_bias(database, image)

            master_bias = bias_combine(master_frames, bias_list, temp_path)

            dark_list = filter_darks(database, image)

            master_dark = dark_combine(master_frames, dark_list,
                                       database, temp_path)

            raw_data = ccdproc.CCDData.read(image, unit="adu")
            bias_data = ccdproc.CCDData.read(master_bias, unit="adu")
            dark_data = ccdproc.CCDData.read(master_dark, unit="adu")

            sub_data = ccdproc.ccd_process(raw_data, master_bias=bias_data,
                                           dark_frame=dark_data,
                                           exposure_key="EXPTIME",
                                           exposure_unit=second,
                                           add_keyword=False)

            clean_images += [sub_data]

            scaling += [1.0 / sub_data.data.mean()]

        flat = ccdproc.combine(clean_images, mem_limit=1e9, scale=scaling,
                               method="median", unit="adu", dtype="float32")
        flat.write(master_flat)

        master_frames[master_flat] = flat_list

    return master_flat


def process(origin, suffix, temp_path):
    """
    Image processing pipeline
    """

    logger.debug("Creating database with all the images")

    database = create_database(origin, recursive=True)

    logger.info("Iterating over all the science images")

    master_frames = {}
    for image in database["file"]:
        if image.endswith(suffix):
            logger.debug("Reducing image: {}".format(image))

            bias_list = filter_bias(database, image)

            master_bias = bias_combine(master_frames, bias_list, temp_path)

            dark_list = filter_darks(database, image)

            master_dark = dark_combine(master_frames, dark_list,
                                       database, temp_path)

            flat_list = filter_flats(database, image)

            master_flat = flat_combine(master_frames, flat_list,
                                       database, temp_path)

            raw_data = ccdproc.CCDData.read(image, unit="adu")
            bias_data = ccdproc.CCDData.read(master_bias, unit="adu")
            dark_data = ccdproc.CCDData.read(master_dark, unit="adu")
            flat_data = ccdproc.CCDData.read(master_flat, unit="adu")

            sub_data = ccdproc.ccd_process(raw_data, master_bias=bias_data,
                                           dark_frame=dark_data,
                                           master_flat=flat_data,
                                           exposure_key="EXPTIME",
                                           exposure_unit=second,
                                           add_keyword=False)

            tmp_image = os.path.join(temp_path, os.path.basename(image))
            clean_image = tmp_image.replace(suffix, "imc.fits")

            sub_data.data = sub_data.data.astype("float32")
            sub_data.header.pop("BZERO")
            sub_data.write(clean_image)

            logger.debug("Doing astrometry and photometry")

            astro_image = clean_image.replace('imc.fits', 'ima.fits')
            phot_cat = clean_image.replace('imc.fits', 'ima.cat')
            try:
                solve(clean_image, new_fits=astro_image, use_sextractor=True,
                      ra=sub_data.header["RA"], dec=sub_data.header["DEC"])
                sextractor(astro_image, catalog_name=phot_cat)
            except (AstrometryError, AstrometryWarning):
                logger.exception("Astrometry failed, doing photometry only")
                sextractor(astro_image, catalog_name=phot_cat)


def run(origin, suffix, verbose=None):
    """
    Pipeline execution function
    """

    try:
        temp_path = mkdtemp(prefix=icat.__basename__+"_")
        logger.debug("Temporary directory created: {}".format(temp_path))

        process(origin, suffix, temp_path=temp_path)

    except Exception:
        logger.exception("")
        raise
    finally:
        if verbose <= 1:
            rmtree(temp_path)


def main():
    """
    Entry point
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="count",
                        help="log debug messages")
    parser.add_argument("-s", "--suffix", default="imr.fits",
                        help="Suffix of science images (default: imr.fits)")
    parser.add_argument("path", help="location of FITS images")
    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    run(args.path, args.suffix, args.verbose)


if __name__ == "__main__":
    main()
