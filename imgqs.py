"""
Verify the quality of the TJO images and store them
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os.path
import argparse
from re import search
from glob import glob
from shutil import rmtree
from socket import gethostname
from tempfile import mkdtemp
from tempfile import NamedTemporaryFile as NTF
from subprocess import call, CalledProcessError

import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization.mpl_normalize import simple_norm

import icat
import icat.pipelines as ppl
from icat import system
from icat import astrometry
from icat import photometry
from icat.logs import get_logger

store = ppl.import_pipeline("gaia-store")
delete = ppl.import_pipeline("gaia-remove")
catstat = ppl.import_pipeline("catstat")
logger = get_logger(__name__)
hostname = gethostname()

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config")


def validate(header, keyword, dtype=None, values=None, min_value=None,
             default=None):
    """
    Verify header values
    """

    if values is not None and min_value is not None:
        error_str = "Incompatible arguments used to validate keyword "
        error_str += "'{}': {} and {}"
        raise ppl.PipelineError(error_str.format(keyword,
                                                 "values", "min_value"))

    new_value = None
    header_value = header.get(keyword, default)
    try:
        if dtype is None:
            new_value = comp_value = header_value
        elif dtype == "str":
            new_value = comp_value = "{}".format(header_value or "")
        elif dtype == "int":
            new_value = comp_value = int(float(header_value))
        elif dtype == "float":
            new_value = comp_value = float(header_value)
        elif dtype == "datetime":
            comp_value = Time(header_value, scale="utc")
            new_value = comp_value.isot
        elif dtype == "date":
            comp_value = Time(header_value, scale="utc")
            new_value = comp_value.iso.split()[0]
        elif dtype == "time":
            comp_value = Angle("{} hours".format(header_value))
            new_value = comp_value.to_string(sep="::")
        elif dtype == "angle":
            comp_value = Angle("{} degrees".format(header_value))
            new_value = comp_value.to_string(sep="::")
        else:
            error_string = "Invalid type '{}' found in keyword '{}'"
            raise ppl.PipelineError(error_string.format(dtype, keyword))
    except (ValueError, TypeError):
        return [keyword]

    if new_value is not None:
        header[keyword] = new_value

    if values is not None and comp_value not in values:
        return [keyword]
    elif min_value is not None and comp_value < min_value:
        return [keyword]

    return []


def in_range(operation, obstype, value, min_value=None, max_value=None):
    """
    Report an error message if value is not between
    min_value and max_value
    """
    error_msg = ""
    error_string = "{} value ({}) for '{}' "

    if min_value is not None and value < min_value:
        error_msg = error_string.format(operation, str(value), obstype)
        error_msg += "below minimum value of {}".format(min_value)
        logger.warning(error_msg)

    if max_value is not None and value > max_value:
        error_msg = error_string.format(operation, str(value), obstype)
        error_msg += "over maximum value of {}".format(max_value)
        logger.warning(error_msg)

    return False if error_msg else True


def nameread(string):
    """
    Parse string, returning a dictionary with parsed values
    """

    logger.debug("Entering nameread...")

    # Option 1: string is a reduced filename

    result = {}
    if string.startswith("TJO"):
        return result

    # Option 2: string is an observe.sh filename

    regex = (r"(RAW(?P<PROPCODE>(?<=RAW)\d{5})|"
             r"D(?P<DEFOCUS>(?<=D)-?\d+)|"
             r".*_(?P<INSTRUME>[RUV])_|"
             r"(?P<OBSTYPE>cbs|cth|cf|imr)|"
             r"(?P<FILTER>[UBVRIN]).fits){,5}")

    result = search(regex, string).groupdict()

    # Option 3: string is an OCS filename

    if result["PROPCODE"] is None:
        regex = (r"(P(?P<PROPCODE>(?<=P)\d+)|"
                 r"S(?P<BLOCKID>(?<=S)\d+)|"
                 r".*T(?P<TARGETID>(?<=T)\d+)|"
                 r"O(?P<OBSERVID>(?<=O)\d+)|"
                 r"I(?P<INSTRUID>(?<=I)\d+)|"
                 r"D(?P<DEFOCUS>(?<=D)-?\d+)|"
                 r"E(?P<IMAGEID>(?<=E)\d+)|"
                 r"_(?P<INSTRUME>[RSTUV])_|"
                 r"(?P<OBSTYPE>cbs|cth|cf|cd|imr)|"
                 r"(?P<FILTER>[UBVRIN])_*\d*.fits){,10}")

        result = search(regex, string).groupdict()

    defocus = result["DEFOCUS"]
    result["DEFOCUS"] = defocus if defocus is None else int(defocus)

    instrume = result["INSTRUME"]
    options = {"R": "MEIA", "S": "ARES", "T": "SBIG",
               "U": "MEIA2", "V": "MEIA3"}
    result["INSTRUME"] = options.get(instrume, instrume)

    obstype = result["OBSTYPE"]
    options = {"cbs": "Bias", "cth": "Dark",
               "cf": "Sky Flat", "cd": "Dome Flat"}
    options["imr"] = "SST" if int(result["PROPCODE"] or 0) == 99 else "Science"
    result["OBSTYPE"] = options.get(obstype, obstype)

    header = {name: value for name, value in result.iteritems()
              if value is not None}

    if not header:
        log_msg = "File name '%s' does not match any of the known formats"
        logger.warning(log_msg, string)

    return header


def statistics(filename, catname, header=None):
    """
    Compute seeing and ellipticity for a given image with a catalog
    """

    try:
        if header is None:
            header = fits.open(filename)[0].header
    except Exception:
        error_string = "WCS information could not be read for image '{}'"
        raise ppl.PipelineError(error_string.format(filename))

    try:
        scales = wcs.utils.proj_plane_pixel_scales(wcs.WCS(header))
        pix_scale = Angle([Angle(scales[0], unit=header["CUNIT1"]),
                           Angle(scales[1], unit=header["CUNIT2"])]).mean()
    except Exception:
        err_msg = "Image resolution not found in image '{}'"
        logger.warning(err_msg.format(filename))

    result = {}
    stars = catstat.run(catname, show_spinner=False)
    try:
        fwhm = stars[stars["COLUMN"] == "FLUX_RADIUS"]["MEDIAN"][0]
        result["FWHM"] = 2 * fwhm * pix_scale

    except (IndexError, NameError) as err:
        logger.debug(err)

        try:
            fwhm = stars[stars["COLUMN"] == "FWHM_WORLD"]["MEDIAN"][0]
            result["FWHM"] = Angle(fwhm, unit="degree")

        except IndexError as err:
            logger.debug(err)

    ellipticity = stars[stars["COLUMN"] == "ELLIPTICITY"]["MEDIAN"][0]
    result["ELLIPTICITY"] = ellipticity
    result["NSTARS"] = stars["NSTARS"][0]

    return result


def imgqf(header, system_size=0):
    """
    Validate header and file size
    """

    logger.debug("Entering imgqf...")

    # Some keywords must simply be added (will never be invalid)

    header["ORIGIN"] = "OAdM"
    header["TELESCOP"] = "TJO"
    instrume = header.get("INSTRUME")
    if instrume == "MEIA":
        header["GAIN"] = 1.53
        header["RDNOISE"] = 8.0
    elif instrume == "MEIA2":
        header["GAIN"] = 1.04
        header["RDNOISE"] = 9.0
    elif instrume == "MEIA3":
        header["GAIN"] = 1.94
        header["RDNOISE"] = 8.3

    # OBSTYPE should be equal to IMAGETYP, if IMAGETYP exists
    # (introduced with INDI driver and OCS)

    obstype = header.get("OBSTYPE", "Manual")
    header["OBSTYPE"] = header.get("IMAGETYP", obstype)

    # Validate the rest of header keywords

    obstype = header.get("OBSTYPE", "Manual")
    if obstype == "Bias" or obstype == "Dark":
        valid_filters = ["C"]
        min_value = 0
    else:
        valid_filters = ["U", "B", "V", "R", "I", "N"]
        min_value = 0 if "flat" in obstype.lower() else 1

    valid_instruments = ["MEIA", "MEIA2", "MEIA3", "ARES", "SBIG"]

    invalid_keys = []
    invalid_keys += validate(header, "BITPIX", "int", values=[16])
    invalid_keys += validate(header, "NAXIS", "int", values=[2])
    invalid_keys += validate(header, "NAXIS1", "int", min_value=1)
    invalid_keys += validate(header, "NAXIS2", "int", min_value=1)
    invalid_keys += validate(header, "OFFSET1", "int", min_value=0)
    invalid_keys += validate(header, "OFFSET2", "int", min_value=0)
    invalid_keys += validate(header, "XFACTOR", "int", min_value=1)
    invalid_keys += validate(header, "YFACTOR", "int", min_value=1)
    invalid_keys += validate(header, "EXPTIME", "float", min_value=0.0)
    invalid_keys += validate(header, "OBSERVER", default="Talon")
    invalid_keys += validate(header, "OBJECT", min_value=chr(0), default="")
    invalid_keys += validate(header, "PRIORITY", "int", default=0)
    invalid_keys += validate(header, "INSTRUME", values=valid_instruments)
    invalid_keys += validate(header, "JD", "float")
    invalid_keys += validate(header, "TIME-OBS", "time")
    invalid_keys += validate(header, "PROPCODE", "int", min_value=1, default=0)
    invalid_keys += validate(header, "BLOCKID", "int", min_value=min_value,
                             default=0)
    invalid_keys += validate(header, "TARGETID", "int", min_value=min_value,
                             default=0)
    invalid_keys += validate(header, "OBSERVID", "int", min_value=min_value,
                             default=0)
    invalid_keys += validate(header, "INSTRUID", "int", min_value=1, default=0)
    invalid_keys += validate(header, "IMAGEID", "int", min_value=1, default=0)
    invalid_keys += validate(header, "DOMESTAT", min_value=chr(0), default="")
    invalid_keys += validate(header, "OBSTYPE", default="Manual")
    invalid_keys += validate(header, "LST", "time")
    invalid_keys += validate(header, "POSANGLE", "angle")
    invalid_keys += validate(header, "LATITUDE", "angle")
    invalid_keys += validate(header, "LONGITUD", "angle")
    invalid_keys += validate(header, "ELEVATIO", "angle")
    invalid_keys += validate(header, "AZIMUTH", "angle")
    invalid_keys += validate(header, "HA", "angle")
    invalid_keys += validate(header, "RAEOD", "angle")
    invalid_keys += validate(header, "DECEOD", "angle")
    invalid_keys += validate(header, "RA", "angle")
    invalid_keys += validate(header, "DEC", "angle")
    invalid_keys += validate(header, "OBJRA", "angle")
    invalid_keys += validate(header, "OBJDEC", "angle")
    invalid_keys += validate(header, "EQUINOX", "float", default=2000.0)
    invalid_keys += validate(header, "FILTER", values=valid_filters)
    invalid_keys += validate(header, "CAMTEMP", "float")
    invalid_keys += validate(header, "RAWHENC", "float")
    invalid_keys += validate(header, "RAWDENC", "float")
    invalid_keys += validate(header, "RAWOSTP", "float")
    invalid_keys += validate(header, "FOCUSPOS", "float")
    invalid_keys += validate(header, "DEFOCUS", "int")
    invalid_keys += validate(header, "RAWISTP", "float")
    invalid_keys += validate(header, "WXTEMP", "float")
    invalid_keys += validate(header, "WXPRES", "float")
    invalid_keys += validate(header, "WXWNDSPD", "int", min_value=0)
    invalid_keys += validate(header, "WXWNDDIR", "int", min_value=0)
    invalid_keys += validate(header, "WXHUMID", "int", min_value=0)

    if instrume == "ARES":
        invalid_keys += validate(header, "DATE-OBS", "datetime")
    else:
        invalid_keys += validate(header, "DATE-OBS", "date")

    observer = header.get("OBSERVER")
    target = header.get("OBJECT", "")
    if observer == "Indi":
        invalid_keys += validate(header, "MJD-OBS", "float")
        invalid_keys += validate(header, "IMAGETYP", default=obstype)
        invalid_keys += validate(header, "TARGET", min_value=chr(0),
                                 default=target)

    imgqs = 0
    if invalid_keys:
        invalid_msg = ", ".join(invalid_keys)
        logger.warning("Invalid keywords in IMGQF: {}".format(invalid_msg))
        imgqs += 2

    # Check that image has a size compatible with a FITS file

    hblocks = len(header) * fits.card.CARD_LENGTH
    dblocks = header["NAXIS1"] * header["NAXIS2"] * header["BITPIX"] / 8
    min_size = int(hblocks / fits.header.BLOCK_SIZE + 0.5)
    min_size += int(dblocks / fits.header.BLOCK_SIZE + 0.5)
    min_size *= fits.header.BLOCK_SIZE

    if system_size % fits.header.BLOCK_SIZE != 0 or system_size < min_size:
        warning_msg = "System file size ({}) is not equal ".format(system_size)
        warning_msg += "to the expected size ({})".format(min_size)
        logger.warning(warning_msg)
        imgqs += 1

    # Finally, update IMGQS

    header["IMGQS"] = "{}".format(imgqs)


def imgqv(image, obstype, instrume):
    """
    Compute statistical values for calibration images
    """

    logger.debug("Entering imgqv...")

    valid = True
    submin, submax = (image.shape * np.array([[0.3], [0.7]])).astype(int)
    subimage = image[submin[0]:submax[0], submin[1]:submax[1]]
    std = np.std(subimage)
    median = np.median(subimage)
    if obstype == "Bias" and instrume == "MEIA":
        valid = valid and in_range("Median", obstype, median, 1700, 2200)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=100)
    elif obstype == "Bias" and instrume == "MEIA2":
        valid = valid and in_range("Median", obstype, median, 100, 500)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=100)
    elif obstype == "Bias" and instrume == "MEIA3":
        valid = valid and in_range("Median", obstype, median, 800, 1200)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=100)
    elif obstype == "Dark" and instrume == "MEIA":
        valid = valid and in_range("Median", obstype, median, 1700, 2200)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=500)
    elif obstype == "Dark" and instrume == "MEIA2":
        valid = valid and in_range("Median", obstype, median, 100, 500)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=100)
    elif obstype == "Dark" and instrume == "MEIA3":
        valid = valid and in_range("Median", obstype, median, 800, 1200)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=100)
    elif obstype in ["Flat", "Sky Flat"] and instrume in ["MEIA", "MEIA2"]:
        valid = valid and in_range("Median", obstype, median, 20000, 40000)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=1500)
    elif obstype in ["Flat", "Sky Flat"] and instrume == "MEIA3":
        valid = valid and in_range("Median", obstype, median, 20000, 40000)
        valid = valid and in_range("Standard deviation", obstype, std,
                                   max_value=1500)
    else:
        return "9"

    imgqs = 0 if valid else 1

    return "{}".format(imgqs), median, std


def imgqa(filename, ra, dec, temp_path, verbose):
    """ Ensure that image center corresponds to header values """

    logger.debug("Entering imgqa...")

    try:
        config_file = os.path.join(CONFIG_PATH, "astrometry.cfg")
        wcsfile = NTF(dir=temp_path, delete=False)
        wcsfits = wcsfile.name

        wcsfits = astrometry.solve(filename, config_file, ra, dec, radius=1,
                                   use_sextractor=True, verbosity=verbose,
                                   new_fits=wcsfits)

        head_center = SkyCoord(ra + "hours", dec + "degrees")
        distance = astrometry.center(wcsfits).separation(head_center)

        good_sep = Angle("1 arcmin")
        max_sep = Angle("8 arcmin")
        if distance > max_sep:
            warning_str = "Center offset ({}) over maximum distance of {}"
            logger.warning(warning_str.format(distance.to_string(), max_sep))
            imgqs = 2
        elif distance > good_sep:
            warning_str = "Center offset ({}) over recommended value of {}"
            logger.warning(warning_str.format(distance.to_string(), good_sep))
            imgqs = 4
        else:
            logger.debug("Center offset: {}".format(distance.to_string()))
            imgqs = 0

    except astrometry.AstrometryWarning:
        warning_str = "No astrometric solution could be found for '{}'"
        logger.warning(warning_str.format(filename))
        imgqs = 1
        wcsfits = None

    except (astrometry.AstrometryError, KeyError, AttributeError) as err:
        logger.exception(err)
        imgqs = 9
        wcsfits = None

    return "{}".format(imgqs), wcsfits


def imgqp(filename, temp_path, defocus=0):
    """
    Compute quality of stellar like sources
    """

    logger.debug("Entering imgqp...")

    config_file = os.path.join(CONFIG_PATH, "sextractor.conf")
    config_keys = os.path.join(CONFIG_PATH, "sextractor.keys")

    try:
        with NTF(suffix=".cat", dir=temp_path, delete=False) as catfile:
            catname = catfile.name

        if abs(defocus) > 150:
            config_file = config_file.replace("sextractor.conf",
                                              "sextractor.defocused.conf")

        catname = photometry.sextractor(filename, config=config_file,
                                        params=config_keys,
                                        catalog_name=catname)
        logger.debug("Catalog file '{}' created".format(catname))

        properties = statistics(filename, catname)

        valid = True
        valid = valid and in_range("Median", "FWHM", properties["FWHM"],
                                   max_value=Angle("5 arcsec"))

        valid = valid and in_range("Median", "ELLIPTICITY",
                                   properties["ELLIPTICITY"], max_value=0.2)

        imgqs = 0 if valid else 2

    except Exception as err:
        logger.exception(err)
        imgqs = 9

    properties["IMGQS"] = "{}".format(imgqs)
    properties["CATALOG"] = catname
    return properties


def imgqc(filename=None):
    """ Check whether images conform to requirements in database """

    # pstoi_keys = ["PROPCODE", "BLOCKID", "TARGETID",
    #               "OBSERVID", "INSTRUID"]

    logger.debug("Entering imgqc...")

    return "9"


def namewrite(header, prefix="", sfx=None):
    """
    Create a file name from a header
    """

    logger.debug("Entering namewrite...")

    jd = header.get("JD", Time(header["DATE-OBS"], format="isot").jd)

    instrume = header.get("INSTRUME")
    if instrume == "MEIA":
        ins = "_R_"
    elif instrume == "ARES":
        ins = "_S_"
    elif instrume == "SBIG":
        ins = "_T_"
    elif instrume == "MEIA2":
        ins = "_U_"
    elif instrume == "MEIA3":
        ins = "_V_"
    else:
        return

    if sfx is None:
        obstype = header.get("OBSTYPE")
        if obstype == "Bias":
            sfx = "cbs"
        elif obstype == "Dark":
            sfx = "cth"
        elif obstype in ["Flat", "Sky Flat"]:
            sfx = "cf{}".format(header.get("FILTER", ""))
        elif obstype == "Science" or obstype == "SST":
            sfx = "imr"
        else:
            return

    outf_string = "{}{}/TJO{:.5f}{}{}.fits"
    return outf_string.format(prefix, int(jd), jd, ins, sfx)


def save_image(origin, destination, dependencies=None):
    """
    Store temporary images to destination path, with required dependencies
    """
    if destination:
        system.new_copy(origin, destination, makedirs=True)

        if dependencies is not None:
            dependencies = [dependencies]

        store.run([destination], dependencies=dependencies)

        # store.run usually compresses images, modifying their names

        compressed_name = destination + ".gz"
        if os.path.isfile(compressed_name):
            destination = compressed_name

    return destination


def plotfits(orig_file, image=None, header=None, title=""):
    """
    Create PNG image from orig_file
    """

    if image is None or header is None:
        hdu = fits.open(orig_file)[0]
        image = image or hdu.data
        header = header or hdu.header

    if orig_file.endswith(".gz"):
        orig_file = os.path.splitext(orig_file)[0]
    png_name = "{}.png".format(os.path.splitext(orig_file)[0])

    submin, submax = (image.shape * np.array([[0.3], [0.7]])).astype(int)
    subimage = image[submin[0]:submax[0], submin[1]:submax[1]]
    norm = simple_norm(subimage, percent=99)

    try:
        wcsheader = wcs.WCS(header)
        resolved = wcsheader.is_celestial
    except Exception:
        logger.warning("Header could not be passed to WCS")
        header = {}
        resolved = False

    projection = dict(projection=wcsheader) if resolved else None
    fig, norm_image = plt.subplots(subplot_kw=projection, figsize=(5, 5))
    norm_image.imshow(norm(image), cmap="Greys_r")

    if resolved:
        norm_image.coords.grid(color="yellow", ls="solid")

    norm_image.axis("off")
    norm_image.set_xticklabels([])
    norm_image.set_yticklabels([])

    fig.tight_layout(pad=0, rect=(-0.01, -0.01, 1.01, 1.01))

    long_title = "{0}{1}{0}".format(" "*50, title or orig_file)
    fig.text(0.5, 1, long_title, color="white",
             verticalalignment="top", horizontalalignment="center",
             weight="bold", size="x-large", backgroundcolor="black")

    footer = [header.get("OBJECT", ""), header.get("FILTER", ""),
              header.get("IMGQS", "")]
    if any(footer):
        empty_text = "I" + " " * 100
        fig.text(0.5, 0, empty_text, backgroundcolor="black", size="x-large",
                 verticalalignment="bottom", horizontalalignment="center")

    offset = 0
    for text, position in zip(footer, ["left", "center", "right"]):
        fig.text(offset, 0, "   {}   ".format(text), color="white",
                 verticalalignment="bottom", horizontalalignment=position,
                 weight="bold", size="x-large",  backgroundcolor="black")
        offset += 0.5

    fig.savefig(png_name)

    plt.close(fig)


def copy_png(orig_file):
    """
    plotfits returns nothing.
    Find PNG image resulting from orig_file and copy it to RUBIES
    """

    if orig_file.endswith(".gz"):
        orig_file = os.path.splitext(orig_file)[0]
    png_name = os.path.splitext(orig_file)[0]
    try:
        png_image = glob("{}*.png".format(png_name))[0]
    except IndexError:
        return

    logger.debug("PNG image '%s' created", png_image)

    if "estall" in hostname:
        cmd = "scp {} rubies:/tmp/meia.jpg"
        call(cmd.format(png_image).split())

    return png_image


def run(filename, temp_path=None, verbose=0, dry_run=False):
    """
    Pipeline execution function
    """

    # Make a temporary decompressed copy to work on it

    extension = os.path.splitext(filename)[-1]
    rawfname = system.make_temp(filename, suffix=extension, dir=temp_path)
    rawfname = system.decompress(rawfname)
    os.chmod(rawfname, 0644)

    logger.debug("Temporary copy created with name: {}".format(rawfname))

    with fits.open(rawfname) as rawim:
        header = rawim[0].header
        image = rawim[0].data

    # First step, update header with values in file name, if required

    header.update(nameread(os.path.basename(filename)))

    # Validate header values and file size

    imgqf(header, os.path.getsize(rawfname))

    # Update raw image header with valid values before continuing

    with fits.open(rawfname, mode="update") as rawim:
        rawim[0].header = header

    wcsfname = None
    obstype = header.get("OBSTYPE")
    instrument = header.get("INSTRUME")
    valid_obstypes = ["Science", "SST", "Acquisition"]
    if obstype in ["Bias", "Dark", "Flat", "Sky Flat"]:

        # Validate statistical values for calibration images

        imgqs = imgqv(image, obstype, instrument)[0]
        header["IMGQS"] = header.get("IMGQS", "") + imgqs

    elif "MEIA" in instrument and obstype in valid_obstypes:

        # Compute astrometric and photometric solution for science images

        imgqs, wcsfname = imgqa(rawfname, header["RA"], header["DEC"],
                                temp_path, verbose)
        if imgqs == "9":
            logger.warning("Repeating imgqa...")
            imgqs, wcsfname = imgqa(rawfname, header["RA"], header["DEC"],
                                    temp_path, 1)
        header["IMGQS"] = header.get("IMGQS", "") + imgqs
        logger.debug("WCS file '{}' created".format(wcsfname))

        if wcsfname:
            defocus = header.get("DEFOCUS", 0)
            properties = imgqp(wcsfname, temp_path, defocus=defocus)
            imgqs = properties["IMGQS"]
        else:
            imgqs = "9"

        header["IMGQS"] = header.get("IMGQS", "") + imgqs

    elif obstype in valid_obstypes:
        header["IMGQS"] = header.get("IMGQS", "") + "9"

    else:
        raise ppl.PipelineError("Invalid 'OBSTYPE': {}".format(obstype))

    # Ensure that image conforms to requirements in the database (TBD)

    header["IMGQS"] = header.get("IMGQS", "") + imgqc()

    # Update IMGQS in all headers

    with fits.open(rawfname, mode="update") as rawim:
        rawheader = rawim[0].header
        rawheader["IMGQS"] = header["IMGQS"]

    wcsheader = None
    if wcsfname:
        with fits.open(wcsfname, mode="update") as wcsim:
            wcsheader = wcsim[0].header
            wcsheader["IMGQS"] = header["IMGQS"]
            wcsheader["NSTARS"] = properties["NSTARS"]
            if wcsheader["NSTARS"]:
                wcsheader["SEEING"] = properties["FWHM"].arcsec
                wcsheader["ELLIPTIC"] = properties["ELLIPTICITY"]

    # Define output file name for images to be saved

    base_path = os.path.dirname(os.path.dirname(__file__))
    outfraw = namewrite(header, prefix=os.path.join(base_path, "rawdata"))

    if outfraw:
        logger.debug("Output file name for '{}': {}".format(filename, outfraw))

    outfwcs = None
    red_prefix = os.path.join(base_path, "reddata")
    if wcsheader:
        outfwcs = namewrite(wcsheader, prefix=red_prefix, sfx="imt")

        logger.debug("Output file name for '{}': {}".format(wcsfname, outfwcs))

    # If not a dry run, store changes permanently

    if dry_run:
        if outfraw and os.path.isfile(outfraw):
            raise OSError("Output file '{}' already exists.".format(outfraw))

        if outfwcs and os.path.isfile(outfwcs):
            raise OSError("Output file '{}' already exists.".format(outfwcs))

    else:
        try:
            outfraw = save_image(rawfname, outfraw)
            outfwcs = save_image(wcsfname, outfwcs, dependencies=outfraw)

            if not outfraw and not outfwcs:
                err_msg = "No files saved for '{}'"
                raise ppl.PipelineError(err_msg.format(filename))

        except Warning:
            pass

        except Exception:
            if outfwcs and os.path.isfile(outfwcs):
                delete.run(outfwcs)
            if outfraw and os.path.isfile(outfraw):
                delete.run(outfraw)
            raise

        # Create a PNG image from reduction for illustration purposes

        try:
            png_image = None
            if "MEIA" in instrument and png_image is None:
                valid_header = wcsheader or header
                plotfits(rawfname, image, header=valid_header,
                         title=os.path.basename(filename))
                copy_png(rawfname)

        except Exception:
            log_msg = "PNG image for '%s' could not be created"
            logger.warning(log_msg, filename, exc_info=True)

        # Tidy up

        try:
            system.remove(filename)
        except (OSError, CalledProcessError):
            warn_string = "Original image '{}' could not be deleted"
            logger.warning(warn_string.format(filename))
            raise RuntimeWarning(warn_string.format(filename))

    return outfraw, outfwcs


def main():
    """
    Entry point
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("filename", help="FITS image to verify and store")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="Test only, without modifying anything")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="count",
                           help="Log debug messages")
    verbosity.add_argument("-q", "--quiet", action="count",
                           help="Log errors only")
    args = parser.parse_args()

    logger = get_logger(__name__, verbose=bool(args.verbose)-bool(args.quiet))

    try:
        # Create temp_path here to delete it at the end

        temp_path = mkdtemp(prefix=icat.__basename__+"_")
        logger.debug("Temporary directory created: {}".format(temp_path))

        logger.info("Storing '{}'".format(args.filename))

        saved_files = run(args.filename, temp_path, args.verbose, args.dry_run)

        dry_msg = "would have been " if args.dry_run else ""
        if all(saved_files):
            outfiles = "' and '".join(saved_files)
            info_msg = "Two files {}successfully stored for '{}': '{}'"
            logger.info(info_msg.format(dry_msg, args.filename, outfiles))
        elif any(saved_files):
            outfile = filter(None, saved_files)[0]
            info_msg = "One file {}successfully stored for '{}': '{}'"
            logger.info(info_msg.format(dry_msg, args.filename, outfile))

    except Exception:
        if args.dry_run:
            error_string = "File '{}' would have NOT been stored."
        else:
            error_string = "File '{}' has NOT been stored."
        logger.exception(error_string.format(args.filename))
        raise

    finally:
        if args.verbose <= 1:
            rmtree(temp_path)


if __name__ == "__main__":
    main()
