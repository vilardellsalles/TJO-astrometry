from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import shutil
import os.path
import tempfile
import subprocess
from distutils.spawn import find_executable
import astropy.units as u
from astropy.io import fits
import astropy.coordinates as coord


ASTROMETRY_COMMAND = 'solve-field'

class AstrometryError(StandardError, subprocess.CalledProcessError):
    """ 
    Generic error to be kept in case that astrometry.net is no longer
    used
    """
    pass

class AstrometryNetNotInstalled(StandardError):
    """ Raised if Astrometry.net is not installed on the system """
    pass


class AstrometryNetError(subprocess.CalledProcessError):
    """ Raised if the execution of Astrometry.net fails """
    pass


class AstrometryNetUnsolvedField(subprocess.CalledProcessError):
    """ Raised if Astrometry.net could not solve the field """

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return self.path + ": could not solve field"


def astrometry_net(path, ra = None, dec = None, radius = 1, verbosity = 0):
                   
    """ Do astrometry on a FITS image using Astrometry.net.

    Use a local build of the amazing Astrometry.net software [1] in order to
    compute the astrometric solution of a FITS image. This software has many,
    many advantages over the well-respected SCAMP, but the most important one
    is that it is a blind astrometric calibration service. We do not need to
    know literally anything about the image, including approximate coordinates,
    scale and equinox. It just works, giving us a new FITS file containing the
    WCS header.

    In order for this function to work, you must have built and installed the
    Astrometry.net code in your machine [2]. The main high-level command-line
    user interface, 'solve-field', is expected to be available in your PATH;
    otherwise, the AstrometryNetNotInstalled exception is raised. Note that
    you also need to download the appropriate index files [3], which are 
    considerably heavy. 

    This code has been adapted from the Lemon repository [4], but stripping
    all the superfluous dependencies (except astropy).

    [1] http://astrometry.net/
    [2] http://astrometry.net/doc/build.html
    [3] http://astrometry.net/doc/readme.html#getting-index-files
    [4] https://github.com/vterron/lemon/blob/master/astrometry.py

    Keyword arguments:

    ra,
    dec,
    radius - restrict the Astrometry.net search to those indexes within
             'radius' degrees of the field center given by ('ra', 'dec').
             Both the right ascension and declination must be given in order
             for this feature to work. The three arguments must be expressed
             in degrees.
    verbosity - the verbosity level. The default value is zero, meaning that
                the function executes silently. A value of one makes both the
                standard output and standard error of Astrometry.net visible.
                Above that, the number of -v flags send to it equals the value
                of the argument minus one. For example: verbosity = 3 allows us
                to see stdout and stderr, and calls Astrometry.net with two -v
                flags. Most of the time, verbosities greater than one are only
                needed for debugging.
    """

    emsg = '{0} '.format(ASTROMETRY_COMMAND)
    emsg += 'not found in the current environment'
    if not find_executable(ASTROMETRY_COMMAND):
        raise AstrometryNetNotInstalled(emsg)

    root, ext = os.path.splitext(os.path.basename(path))
    # Place all output files in this directory
    kwargs = dict(prefix = root + '_', suffix = '_astrometry.net')
    output_dir = tempfile.mkdtemp(**kwargs)

    # Path to the temporary FITS file containing the WCS header
    kwargs = dict(prefix = root + '_astrometry_', suffix = ext)
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    # If the field solved, Astrometry.net creates a <base>.solved output file
    # that contains (binary) 1. That is: if this file does not exist, we know
    # that an astrometric solution could not be found.
    solved_file = os.path.join(output_dir, root + '.solved')

    # --dir: place all output files in the specified directory.
    # --no-plots: don't create any plots of the results.
    # --new-fits: the new FITS file containing the WCS header.
    # --no-fits2fits: don't sanitize FITS files; assume they're already valid.
    # --overwrite: overwrite output files if they already exist.

    args = [ASTROMETRY_COMMAND, path,
            '--dir', output_dir,
            '--no-plots',
            '--use-sextractor',
            '--new-fits', output_path,
    #Do not use this option, just in case we have bad FITS
    #        '--no-fits2fits',
            '--overwrite']

    # -3 / --ra <degrees or hh:mm:ss>: only search in indexes within 'radius'
    # of the field center given by 'ra' and 'dec'
    # -4 / --dec <degrees or [+-]dd:mm:ss>: only search in indexes within
    # 'radius' of the field center given by 'ra' and 'dec'
    # -5 / --radius <degrees>: only search in indexes within 'radius' of the
    # field center given by ('ra', 'dec')

    if ra is not None:
        args += ['--ra', '{0}'.format(
                 coord.Longitude(ra, unit=u.hourangle).degree)]

    if dec is not None:
        args += ['--dec', '{0}'.format(
                 coord.Latitude(dec, unit=u.degree).degree)]

    if radius is not None:
        args += ['--radius', '{0}'.format(radius)]

    # -v / --verbose: be more chatty -- repeat for even more verboseness. A
    # value of 'verbosity' equal to zero means that both the standard output
    # and error of Astrometry.net and redirected to the null device. Above
    # that, we send 'verbosity' minus one -v flags to Astrometry.net.

    kwargs = {}
    if not verbosity:
        # Needed when 'verbosity' is 0
        null_fd = open(os.devnull, 'w')
        kwargs['stdout'] = kwargs['stderr'] = null_fd
    elif verbosity > 1:
        args.append('-{0}'.format('v' * (verbosity - 1)))

    try:
        subprocess.check_call(args, **kwargs)

        # .solved file must exist and contain a binary one
        with open(solved_file, 'rb') as fd:
            if ord(fd.read()) != 1:
                raise AstrometryNetUnsolvedField(path)

        return output_path

    # Not a FITS image? Astrometry.net should raise the proper error
    except subprocess.CalledProcessError as e:
        raise AstrometryNetError(e.returncode, e.cmd)
    # If .solved file doesn't exist or contain one
    except (IOError, AstrometryNetUnsolvedField):
        raise AstrometryNetUnsolvedField(path)
    finally:
        shutil.rmtree(output_dir)
        if not verbosity:
            null_fd.close()


def solve(fitsname, ra = None, dec = None, radius = 1, verbosity = 0):
    """ Process a given image using Astrometry.net """

    try:
        return astrometry_net(fitsname, ra, dec, radius, verbosity)
    except (AstrometryNetError, AstrometryNetNotInstalled) as err:
        print(err.args)
        raise AstrometryError(err)
    except AstrometryNetUnsolvedField:
        return None


def center(fitsname):
    """
    Retrieve image center from a FITS image with WCS information in header
    """

    try:
        output = subprocess.check_output(['wcsinfo', fitsname], 
                                         stderr=subprocess.STDOUT)
        ra = output[output.rfind('ra_center_hms'):].split()[1]
        dec = output[output.rfind('dec_center_dms'):].split()[1]
        return coord.SkyCoord(ra + 'hours', dec + 'degrees')
    except subprocess.CalledProcessError as err:
        print(err.args)
        raise AstrometryError(err)

if __name__ == '__main__':
   filename = sys.argv[1]

   header = fits.open(filename)[0].header
   solve(filename, header['RA'], header['DEC'], verbosity=1)

