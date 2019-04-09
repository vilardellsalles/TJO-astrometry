"""
Add SExtractor photometry to ICAT database
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse

from astropy.table import Table

import icat.mysql_database as icatdb
from icat import system
from icat.logs import get_logger


logger = get_logger(__name__)


def run(input_catalog, recursive=False, files_filter="_trl.dat"):
    """
    Pipeline execution function
    """

    imdb = icatdb.ImagesDB(read_default_group="Database")
    photdb = icatdb.PhotometryDB(read_default_group="Database")

    num_stars = 0
    num_catalogs = 0
    for catalog in system.path_flatten(input_catalog, recursive,
                                       files_filter=files_filter):

        logger.debug("Adding '{}' to ICAT database".format(catalog))

        try:
            image = catalog.replace(files_filter, ".fits.gz")
            imid = imdb.image_value(image, value="id")[0]
        except icatdb.DBError:
            err_msg = "Image '{}' not found in database".format(image)
            raise icatdb.DBError(err_msg)

        catvalues = Table.read(catalog, format="ascii.sextractor")
        photdb.insert(catvalues, image_id=imid)

        num_stars += len(catvalues)
        num_catalogs += 1

    if num_stars:
        cat_plural = "" if num_catalogs == 1 else "s"
        star_plural = "" if num_stars == 1 else "s"

        log_msg = "Successfully added {} file{} with {} star{} to database"
        logger.info(log_msg.format(num_catalogs, cat_plural,
                                   num_stars, star_plural))
    else:
        warn_msg = "No stars in input files!"
        logger.warning(warn_msg)
        raise Warning(warn_msg)


def main():
    """
    Entry point
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-r", "--recursive", action="store_true",
                        help="recursively enter in directories")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log debug messages")
    help_string = "when path is a directory, store only files having FILTER "
    help_string += "in name (default: '%(default)s')"
    parser.add_argument("-f", "--filter", default="_trl.dat", help=help_string)
    parser.add_argument("path", nargs="+",
                        help="file or directory with photometry to add")
    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    try:
        run(args.path, args.recursive, args.filter)
    except Exception as err:
        logger.exception(err)
        raise


if __name__ == "__main__":
    main()
