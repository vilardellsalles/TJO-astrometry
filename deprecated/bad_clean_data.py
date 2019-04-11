from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
import os
import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import ccdproc as ccd
from pyraf import iraf

def zerocombine(image_list):
    """ Combine all bias frames """

    #TODO example: Should check whether images have the proper median values
    data = [ccd.CCDData.read(image, unit='adu') for image in image_list]
    combiner = ccd.Combiner(data)
    return combiner.median_combine()


def darkcombine(image_list, bias=None):
    """ 
    Combine all dark frames by subtracting a bias frame and 
    scaling by the exposure time
    """

    data = [ccd.CCDData.read(image, unit='adu') for image in image_list]
    if bias is not None:
        data = [ccd.subtract_bias(image, bias) for image in data]
    combiner = ccd.Combiner(data)
    combiner.scaling = [1.0/image.header['exptime'] for image in data]
    dark = combiner.median_combine()
    dark.header['exptime'] = 1.0
    return dark


def flatcombine(image_list, bias=None, dark=None):
    """
    Combine all flat field images by subtracting master bias and dark frames 
    and scaling each flat by the mean value
    """

    data = [ccd.CCDData.read(image, unit='adu') for image in image_list]
    if bias is not None:
        data = [ccd.subtract_bias(image, bias) for image in data]
    if dark is not None:
        data = [ccd.subtract_dark(image, dark, exposure_time='exptime', 
                                  exposure_unit=u.second, scale=True) 
                                  for image in data]
    combiner = ccd.Combiner(data)
    combiner.scaling = lambda arr: 1/np.ma.average(arr)
    return combiner.median_combine()


if __name__ == '__main__':
    plot = False
    if len(sys.argv) == 2:
        if sys.argv[1] == 'plot':
            plot = True
        else:
            root = os.path.abspath(sys.argv[1]) + os.path.sep
    elif len(sys.argv) == 3:
         root = os.path.abspath(sys.argv[1]) + os.path.sep
         plot = True
    else:
        root = ''

    if plot:
        bias = iraf.imstat(root + "TJO*_R_cbs.fits", 
                          fields="image,npix,mean,stddev,midpt,mode", Stdout=1)
        jd = [float(stats.split()[0].replace(root + 'TJO','').replace(
              '_R_cbs.fits',''))-2457000 for stats in bias if 'TJO' in stats]
        mean = [float(stats.split()[2]) for stats in bias if 'TJO' in stats]
        stddev = [float(stats.split()[3]) for stats in bias if 'TJO' in stats]
        plt.errorbar(jd, mean, yerr=stddev, fmt='o')
    
        dark = iraf.imstat(root + "TJO*_R_cth.fits", 
                          fields="image,npix,mean,stddev,midpt,mode", Stdout=1)
        jd = [float(stats.split()[0].replace(root + 'TJO','').replace(
              '_R_cth.fits',''))-2457000 for stats in dark if 'TJO' in stats]
        mean = [float(stats.split()[2]) for stats in dark if 'TJO' in stats]
        stddev = [float(stats.split()[3]) for stats in dark if 'TJO' in stats]
        plt.errorbar(jd, mean, yerr=stddev, fmt='o')
        plt.show()
        raise SystemExit()

    print('Combining bias...')
    master_bias = zerocombine(glob(root + 'TJO*cbs.fits'))
    master_bias.write(root + 'bias.fits')
    mean = master_bias.data.mean()
    stddev = master_bias.data.std()
    print('Mean:', mean, '+/-', stddev, 'ADU')
    plt.imshow(master_bias, cmap='gray', 
               vmin=mean-3*stddev, vmax=mean+3*stddev)

    print('Combining darks...')
    master_dark = darkcombine(glob(root + 'TJO*cth.fits'), bias=master_bias)
    master_dark.write(root + 'dark.fits')
    mean = master_dark.data.mean()
    stddev = master_dark.data.std()
    print('Mean:', mean, '+/-', stddev, 'ADU/s')
    plt.imshow(master_dark, cmap='gray', vmin=mean-stddev, vmax=mean+stddev)

    master_bias = ccd.CCDData.read(root + 'bias.fits')
    master_dark = ccd.CCDData.read(root + 'dark.fits')
    master_flat = {}
    for filter_name in ['U', 'B', 'V', 'R', 'I']:
        images = glob(root + 'TJO*cf' + filter_name + '.fits')
        if images:
            print('Combining ' + filter_name + ' band flats...')
            master_flat[filter_name] = flatcombine(images, bias=master_bias, 
                                                           dark=master_dark)
            master_flat[filter_name].write(root + 'flat' + filter_name + '.fits')

    print('Cleaning images...')
    for image in glob(root + 'TJO*imr.fits'):
        rawdata = ccd.CCDData.read(image, unit='adu')
        bias_data = ccd.subtract_bias(rawdata, master_bias)
        dark_data = ccd.subtract_dark(bias_data, master_dark, 
                                      exposure_time='exptime', 
                                      exposure_unit=u.second, scale=True)
        filter_name = rawdata.header['filter']
        flat_data = ccd.flat_correct(dark_data, master_flat[filter_name], 
                                     min_value=0.1)
        outimage = image.replace('imr.fits', 'imc.fits')
        flat_data.write(outimage)

# Possible options for debugging
#    plt.imshow(bias)
#    plt.show()
#    plt.imshow(bias.uncertainty.array)
#    plt.show()
