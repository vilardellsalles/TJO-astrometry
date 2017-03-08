"""
Pipeline used to reduce TJO images of the proposal entitled
'First astrometric measurements of stellar activity with spectroscopy'
"""

from __future__ import (unicode_literals, absolute_import,
                        division, print_function)
import argparse
from shutil import rmtree
from tempfile import mkdtemp, NamedTemporaryFile

from astropy import log as astrolog
import ccdproc

import icat
from icat.tools.logs import get_logger
from icat.tools.ccdproc_database import create_database


astrolog.setLevel("ERROR")
logger = get_logger(__name__)


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

        debug_string = "Using {} images to create master bias '{}'"
        logger.debug(debug_string.format(len(bias_list), master_bias))

        ccdproc.combine(bias_list, output_file=master_bias, mem_limit=1e9,
                        method="median", unit="adu")

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

        debug_string = "Cleaning {} images before creating master dark '{}'"
        logger.debug(debug_string.format(len(dark_list), master_dark))

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
                               method="median", unit="adu")
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


def run(origin, suffix, verbose=None):
    """
    Pipeline execution function
    """

    try:
        temp_path = mkdtemp(prefix=icat.__basename__+"_")
        logger.debug("Temporary directory created: {}".format(temp_path))

        process(origin, suffix, temp_path=temp_path)

    except Exception as err:
        logger.exception(err)
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
