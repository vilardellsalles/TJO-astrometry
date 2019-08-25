"""
Construct light curves from photometry stored in database
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse
import os.path

import numpy as np

# from bokeh.palettes import plasma
# from bokeh.models.annotations import Whisker, Span
# from bokeh.models.sources import ColumnDataSource
# from bokeh.resources import INLINE
# import bokeh.plotting as plt

# from astropy.modeling import models, fitting
# from astropy.stats import histogram
# from astropy.units import arcsec
from astropy.time import Time
import astropy.coordinates as coord
from astropy import table

from icat.logs import get_logger, LogOrProgressBar
from icat.mysql_database import MySQLDatabase
# from icat import system
import icat.pipelines as ppl


logger = get_logger(__name__)
reduce_ppl = ppl.import_pipeline("gaia-reduce")


def multiple_observations(tbl, keys):
    """
    Return whether light curve has more than one measurement
    """

    return len(tbl) > 1


def run(ref_path=".", max_sep=5, max_chi=9, proposal=None, name=None):
    """
    Pipeline execution function
    """

    log_msg = "Reading photometry"
    log_msg += " for object '{}'".format(name) if name else ""
    log_msg += " from proposal {}".format(proposal) if proposal else ""

    logger.info("%s...", log_msg)

    db = MySQLDatabase(option_groups=["Database"])

    phot_columns = ["x_image", "y_image", "alpha_j2000", "delta_j2000"]
    phot_columns += ["mag_auto", "magerr_auto", "fwhm_world", "flux_radius"]
    phot_columns += ["ellipticity", "theta_j2000", "theta_world", "flags"]
    im_columns = ["object", "jd", "filename", "elevation"]

    query = "SELECT {}, object, jd, filename, elevatio FROM ph3_photometry "
    query += "JOIN ph3_images ON ph3_images.id=ph3_photometry.image_id "

    if proposal and name:
        query += "WHERE propcode={} AND object='{}'".format(proposal, name)
    elif proposal:
        query += "WHERE propcode={}".format(proposal)
    elif name:
        query += "WHERE object='{}'".format(name)

    object_name = ""
    if name:
        object_name = "_" + name.lower().replace(" ", "_")

    result = db.execute(query.format(",".join(phot_columns)))

    cat_table = table.Table(rows=result or None, names=phot_columns+im_columns)

    if not cat_table:
        raise ppl.PipelineError("No photometry found!")

    # Compute airmass. Stars without astrometry, use header elevation

    elevation = coord.Angle(cat_table["elevation"], unit="degree")
    cat_table["airmass"] = 1 / np.sin(elevation)

    with_pos = (cat_table["alpha_j2000"] != None)
    with_pos &= (cat_table["delta_j2000"] != None)
    if any(with_pos):
        alpha = cat_table["alpha_j2000"][with_pos]
        delta = cat_table["delta_j2000"][with_pos]
        obstime = Time(cat_table["jd"][with_pos], format="jd")
        site = coord.AltAz(obstime=obstime, location=reduce_ppl.OAdM)
        cat_pos = coord.SkyCoord(alpha, delta, unit="degree, degree")
        cat_table["airmass"][with_pos] = cat_pos.transform_to(site).secz

    # Identify same star in all the images

    stars_table = cat_table.group_by("object")

    # Store results

    filename = "curves{}.html".format(object_name)
    result_name = os.path.join(ref_path, filename)
    logger.info("Photometry saved in '%s'", result_name)

    if os.path.isfile(result_name):
        os.remove(result_name)

    cat_table.write(result_name, format="jsviewer", max_lines=-1)


def main():
    """
    Entry point
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("path", default=".", nargs="*",
                        help="path where to store results")

    parser.add_argument("-p", "--proposal",
                        help="select photometry from a single proposal")

    parser.add_argument("-o", "--object",
                        help="select photometry from a single object")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log debug messages")

    help_msg = "maximum separation to match the same star in different "
    help_msg += "catalogs (in arc seconds)"
    parser.add_argument("-s", "--separation", default=5, help=help_msg)

    help_msg = "maximum chi2 to consider a star to be constant"
    parser.add_argument("-c", "--chi2", default=9, help=help_msg)

    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    try:
        run(args.path, args.separation, args.chi2, args.proposal, args.object)
    except Exception as err:
        logger.exception(err)
        raise


if __name__ == "__main__":
    main()
