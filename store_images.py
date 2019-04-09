"""
Add compressed images to ICAT database
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse
from os import chmod

from icat import system
from icat.logs import get_logger
from icat.mysql_database import ImagesDB


logger = get_logger(__name__)


def run(input_images, recursive=False, files_filter=".fits",
        dependencies=None):
    """
    Pipeline execution function
    """

    db = ImagesDB(read_default_group="Database")

    num_images = 0
    for image in system.path_flatten(input_images, recursive,
                                     files_filter=files_filter):
        num_images += 1

        logger.debug("Compressing '{}'".format(image))

        compressed_image = system.compress(image)

        logger.debug("Adding '{}' to ICAT database".format(compressed_image))

        try:
            db.insert(compressed_image, relations=dependencies)
        except Exception:
            if compressed_image != image:
                system.decompress(compressed_image)
            raise

        chmod(compressed_image, 0444)

    plural = "" if num_images == 1 else "s"
    log_string = "Successfully added {} image{} to database"
    logger.info(log_string.format(num_images, plural))


def main():
    """
    Entry point
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-r", "--recursive", action="store_true",
                        help="recursively enter in directories")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log debug messages")
    help_string = "when path is a directory, store only images having FILTER "
    help_string += "in name (default: '%(default)s')"
    parser.add_argument("-f", "--filter", default=".fits", help=help_string)
    parser.add_argument("path", nargs="+",
                        help="file or directory with FITS image to add")
    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    try:
        run(args.path, args.recursive, args.filter)
    except Exception as err:
        logger.exception(err)
        raise


if __name__ == "__main__":
    main()
