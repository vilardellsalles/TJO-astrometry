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
import astropy.coordinates as coord
from astropy import table
from astropy.io import fits
from astropy import log as astrolog
from astropy.units import second, arcmin
import ccdproc
from astroquery.vizier import Vizier

import icat
import icat.astrometry as astro
from icat.logs import get_logger
from icat.photometry import sextractor
from icat.ccdproc_database import create_database

import bokeh.plotting as plt
from bokeh.palettes import Category20b as col


astrolog.setLevel("ERROR")
logger = get_logger(__name__)
Vizier.ROW_LIMIT = 100000


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


def get_zeropoint(catalog, imcenter, imfilter):
    """
    Return the zero point to transform intrumental photometry to
    standard photometry
    """

    obs_catalog = table.Table.read(catalog, format="ascii.sextractor")

    nomad = "I/297"
    apass = "II/336/apass9"
    ref_catalog = Vizier.query_region(imcenter, 10*arcmin, catalog=apass)[0]
    ref_catalog.rename_column("{}mag".format(imfilter), "MAG_STD")
    ref_catalog.rename_column("e_{}mag".format(imfilter), "MAGERR_STD")
    ref_catalog = ref_catalog[~ref_catalog["MAG_STD"].mask]

    ref = coord.SkyCoord(ref_catalog["RAJ2000"], ref_catalog["DEJ2000"])
    obs = coord.SkyCoord(obs_catalog["ALPHA_J2000"],
                         obs_catalog["DELTA_J2000"])

    matched = table.Table(obs.match_to_catalog_sky(ref), masked=True)
    matched["col0"].name = "index"
    matched["col1"].name = "distance"
    del matched["col2"]

    ins_mag = table.Table([obs_catalog["MAG_AUTO"]])
    ins_err = table.Table([obs_catalog["MAGERR_AUTO"]])
    ref_mag = table.Table([ref_catalog[matched["index"]]["MAG_STD"]])
    ref_err = table.Table([ref_catalog[matched["index"]]["MAGERR_STD"]])

    magnitudes = table.hstack([matched, ins_mag, ins_err, ref_mag, ref_err])
    magnitudes.sort(["index", "distance"])

    index = magnitudes["index"]
    good_match = [True]
    good_match += [star != oldstar for star, oldstar in zip(index, index[1:])]

    common = magnitudes[good_match]
    common["ZP"] = common["MAG_STD"] - common["MAG_AUTO"]
    common["ZP_ERR"] = 1.0/np.sqrt(common["MAGERR_STD"]**2+
                                   common["MAGERR_AUTO"]**2)
    print(catalog)
    for star in range(len(common)):
        print(star+1, np.average(common["ZP"][:star],
                                 weights=common["ZP_ERR"][:star]))

    return common


def ccdred(database, suffix, temp_path):
    """
    Clean raw FITS files
    """

    logger.info("Creating master frames and cleaning science images...")

    clean_images = []
    master_frames = {}
    for image in filter(lambda i: i.endswith(suffix), database["file"]):
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
        clean = tmp_image.replace(suffix, "imc.fits")

        sub_data.data = sub_data.data.astype("float32")
        sub_data.header.pop("BZERO")
        sub_data.write(clean)

        clean_images += [clean]

    return clean_images


def astrophot(database):
    """
    Compute astrometry and photometry for a list of images
    """

    logger.info("Computing astrometry and photometry for science images...")

    plot_color = 0
    plot = plt.figure()
    for image in filter(lambda i: i.endswith("imc.fits"), database["file"]):
        log_msg = "Computing astrometry and photometry on image '{}'"
        logger.debug(log_msg.format(image))

        header = database[database["file"] == image][0]

        astro_image = image.replace('imc.fits', 'ima.fits')
        phot_cat = image.replace('imc.fits', 'ima.cat')
        try:
            astro.solve(image, new_fits=astro_image, use_sextractor=True,
                        ra=header["ra"], dec=header["dec"])

            root_path = os.path.dirname(os.path.abspath(__file__))
            config_path = os.path.join(root_path, "config")
            config_file = os.path.join(config_path, "sextractor.conf")
            keys_file = os.path.join(config_path, "sextractor.keys")
            sextractor(astro_image, config=config_file, params=keys_file,
                       catalog_name=phot_cat)

            magnitudes = get_zeropoint(phot_cat, astro.center(astro_image),
                                       header["filter"])

            plot.circle(magnitudes["MAG_STD"],
                        magnitudes["MAG_STD"]-magnitudes["MAG_AUTO"],
                        color=col[20][plot_color],
                        size=magnitudes["distance"].to("arcsec"))
            plot_color += 1

        except (astro.AstrometryError, astro.AstrometryWarning):
            logger.exception("Astrometry failed, doing photometry only")
            sextractor(astro_image, catalog_name=phot_cat)

    plt.show(plot)


def run(origin, suffix, verbose=None):
    """
    Pipeline execution function
    """

    try:
        temp_path = mkdtemp(prefix=icat.__basename__+"_")
        logger.debug("Temporary directory created: {}".format(temp_path))

        logger.debug("Creating database with all the images")

        database = create_database(origin, recursive=True)

        clean_images = ccdred(database, suffix, temp_path=temp_path)

        if clean_images:
            clean_database = create_database(temp_path)
            database = table.vstack([database, clean_database],
                                    join_type="outer")

        astrophot(database)

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
