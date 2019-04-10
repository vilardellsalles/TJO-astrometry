"""
Update image keywords in ICAT database and in file system
"""

from __future__ import absolute_import, division, print_function
import argparse
from os import chmod
from shutil import move
from os.path import splitext

from astropy.io import fits

from icat import system
from icat.logs import get_logger
from icat.mysql_database import ImagesDB


logger = get_logger(__name__)


def run(input_images, recursive=False, **keywords):
    """
    Pipeline execution function
    """

    if not keywords:
        raise SyntaxError("No keywords defined")

    db = ImagesDB(read_default_group="Database")

    num_images = 0
    for filename in system.path_flatten(input_images, recursive):

        logger.debug("Updating '{}'".format(filename))

        # Make a temporary copy to be used in case of error

        rawfname = system.make_temp(filename, splitext(filename)[-1])

        chmod(filename, 0644)
        try:
            with fits.open(filename, mode="update") as image:
                header = image[0].header

                valid_keywords = fits.Header()
                for key in keywords:
                    try:
                        valid_keywords[key] = type(header[key])(keywords[key])
                    except KeyError:
                        log_msg = "Keyword '{}' not found in '{}'. Ignoring..."
                        logger.warning(log_msg.format(key, filename))

                header.update(valid_keywords)

            db.insert(filename, overwrite=True)

            system.remove(rawfname)
            num_images += 1

        except Exception as err:
            logger.exception(err)
            move(rawfname, filename)
            db.insert(filename, overwrite=True)
            raise

        finally:
            chmod(filename, 0444)

        logger.debug("File '{}' successfully updated".format(filename))

    plural = "s" if num_images > 1 else ""
    log_string = "Successfully updated {} image{}"
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
    parser.add_argument("--keyword", metavar="value",
                        help="Keyword to uptate (at least one required)")
    parser.add_argument("path", nargs="+",
                        help="file or directory with FITS image to update")

    # Adapted from: http://stackoverflow.com/questions/37367331/
    parsed, unknown = parser.parse_known_args()
    for key in unknown:
        if key.startswith(("-", "--")):
            parser.add_argument(key)

    args = parser.parse_args()

    keywords = {key: value for key, value in vars(args).iteritems()
                if key not in vars(parsed)}

    logger = get_logger(__name__, verbose=args.verbose)

    run(args.path, args.recursive, **keywords)


if __name__ == "__main__":
    main()
