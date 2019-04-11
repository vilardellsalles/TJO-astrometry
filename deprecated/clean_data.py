from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
from glob import glob
import numpy as np
import astropy.units as u
import ccdproc as ccd
#import matplotlib.pyplot as plt

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
    and scaling each flat by the median value
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
    flat = combiner.median_combine()
    return flat 


if __name__ == '__main__':
    print('Combining bias...')
    master_bias = zerocombine(glob('TJO*cbs.fits'))
    print('Combining darks...')
    master_dark = darkcombine(glob('TJO*cth.fits'), bias=master_bias)
    master_flat = {}
    for filter_name in ['U', 'B', 'V', 'R', 'I']:
        images = glob('TJO*cf' + filter_name + '.fits')
        if images:
            print('Combining ' + filter_name + ' band flats...')
            master_flat[filter_name] = flatcombine(images, bias=master_bias, 
                                                           dark=master_dark)
    print('Cleaning images...')
    for image in glob('TJO*imr.fits'):
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
    master_bias.write('bias.fits')
    master_dark.write('dark.fits')
    for filter_name, flat in master_flat.iteritems():
        flat.write('flat' + filter_name + '.fits')
