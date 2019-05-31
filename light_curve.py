"""
Construct light curves from photometry stored in database
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse
import os.path

# import numpy as np

# from bokeh.palettes import plasma
# from bokeh.models.annotations import Whisker, Span
# from bokeh.models.sources import ColumnDataSource
# from bokeh.resources import INLINE
# import bokeh.plotting as plt

# from astropy.modeling import models, fitting
# from astropy.stats import histogram
# from astropy.units import arcsec
# from astropy.time import Time
import astropy.coordinates as coord
from astropy import table

from icat.logs import get_logger, LogOrProgressBar
from icat.mysql_database import MySQLDatabase
from icat.config import config_section
# from icat import system
import icat.pipelines as ppl


logger = get_logger(__name__)


def multiple_observations(tbl, keys):
    """
    Return whether light curve has more than one measurement
    """

    return len(tbl) > 1


def run(ref_path=".", max_sep=5, max_chi=9, proposal=None, name=None):
    """
    Pipeline execution function
    """

    log_msg = ""
    log_msg += " for object '{}'".format(name) if name else ""
    log_msg += " from proposal {}".format(proposal) if proposal else ""

    logger.info("Reading photometry%s...", log_msg)

    db = MySQLDatabase(option_groups=["Database"])

    phot_columns = ["x_image", "y_image", "alpha_j2000", "delta_j2000"]
    phot_columns += ["mag_auto", "magerr_auto", "fwhm_world", "flux_radius"]
    phot_columns += ["ellipticity", "theta_j2000", "theta_world", "flags"]
    im_columns = ["object", "jd", "filename"]

    query = "SELECT {}, object, jd, filename FROM ph3_photometry "
    query += "JOIN ph3_images ON ph3_images.id=ph3_photometry.image_id "
    if proposal and name:
        query += "WHERE propcode={} AND object='{}'".format(proposal, name)
    elif proposal:
        query += "WHERE propcode={}".format(proposal)
    elif name:
        query += "WHERE object='{}'".format(name)

    result = db.execute(query.format(",".join(phot_columns)))

    cat_table = table.Table(rows=result or None, names=phot_columns+im_columns)

    if not cat_table:
        raise ppl.PipelineError("No photometry found!")

    site = config_section("Location")
    OAdM = coord.EarthLocation(lat=site["latitude"], lon=site["longitude"],
                               height=site["elevation"])

    result_name = os.path.join(ref_path, "curves.html")
    logger.info("Light curves saved in '%s'", result_name)

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
