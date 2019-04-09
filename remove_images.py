"""
Remove images from ICAT database and from file system
"""

from __future__ import absolute_import, division, print_function
import argparse
from os import chmod
from os.path import isfile

from icat import system
from icat.logs import get_logger
from icat.mysql_database import ImagesDB, PhotometryDB, DBError


logger = get_logger(__name__)


def run(input_images, recursive=False, files_filter=".fits"):
    """
    Pipeline execution function
    """

    imdb = ImagesDB(read_default_group="Database")
    catdb = PhotometryDB(read_default_group="Database")

    num_images = 0
    for image in system.path_flatten(input_images, recursive,
                                     files_filter=files_filter):
        num_images += 1

        logger.debug("Removing '{}'".format(image))

        chmod(image, 0644)

        try:
            image_id = imdb.image_value(image, value="id")[0]
        except DBError:
            image_id = None

        try:
            if image_id:
                catdb.remove(image_id, commit=False)
                imdb.remove(image, commit=False)

            catalog = ""
            if image.endswith(".fits.gz"):
                catalog = image.replace(".fits.gz", "_trl.dat")
            elif image.endswith(".fits"):
                catalog = image.replace(".fits", "_trl.dat")

            if isfile(catalog):
                chmod(catalog, 0644)
                system.remove(catalog)
            system.remove(image)

            if image_id:
                catdb.conn.commit()
                imdb.conn.commit()

            logger.info("File '{}' successfully removed".format(image))
        except Exception as err:
            if image_id:
                imdb.conn.rollback()
                catdb.conn.rollback()

            if isfile(image):
                chmod(image, 0444)

            logger.exception(err)
            logger.error("File '{}' NOT removed".format(image))
            raise

    plural = "s" if num_images > 1 else ""
    log_string = "Successfully removed {} image{}"
    logger.info(log_string.format(num_images, plural))


def main():
    """
    Entry point
    """

    global logger

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-r", "--recursive", action="store_true",
                        help="recursively enter in directories")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log debug messages")
    help_string = "when path is a directory, remove only images having FILTER "
    help_string += "in name (default: '.fits')"
    parser.add_argument("-f", "--filter", default=".fits", help=help_string)
    parser.add_argument("path", nargs="+",
                        help="file or directory with FITS image to remove")
    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    run(args.path, args.recursive, args.filter)


if __name__ == "__main__":
    main()
