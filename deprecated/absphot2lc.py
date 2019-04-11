from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
import os.path
import sys
from glob import iglob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import astropy.io.ascii as ascii
import astropy.units as u
import astropy.coordinates as coord

def get_zp(ref_mags, mags, ref_errors, errors, image=''):
    """ Return zero-point in magnitudes and its uncertainty """
    orig_diffs = ref_mags - mags
    orig_errs = np.sqrt(ref_errors*ref_errors + errors*errors)

    niter = 0
    old_diffs = []
    diffs = orig_diffs
    while np.any(old_diffs != diffs) and niter<10:
        niter += 1
        old_diffs = diffs
        nstars = len(diffs)
        zp = np.median(diffs)
        zp_err = np.percentile(np.abs(diffs-zp),68.27) 
        sel_mags = ref_mags[np.where(np.abs(orig_diffs-zp) < 3*zp_err)]
        sel_errs = ref_errors[np.where(np.abs(orig_diffs-zp) < 3*zp_err)]
        diffs = orig_diffs[np.where(np.abs(orig_diffs-zp) < 3*zp_err)]
        errs = orig_errs[np.where(np.abs(orig_diffs-zp) < 3*zp_err)]
    zp_err /= np.sqrt(nstars)

    plt.title(image + ' ' + str(zp) + '+/-' + str(zp_err))
              
    plt.errorbar(ref_mags, orig_diffs, xerr=ref_errors, yerr=orig_errs, 
                 fmt='x')
    plt.errorbar(sel_mags, diffs, xerr=sel_errs, yerr=errs, 
                 fmt='o')
    min_mag = np.min(ref_mags)
    max_mag = np.max(ref_mags)
    plt.plot([min_mag,max_mag],[zp,zp], 'k-')
    plt.plot([min_mag,max_mag],[zp+zp_err,zp+zp_err], 'k--')
    plt.plot([min_mag,max_mag],[zp-zp_err,zp-zp_err], 'k--')
    plt.show()


    return zp, zp_err

if __name__ == '__main__':
    if len(sys.argv) == 2:
        root = os.path.abspath(sys.argv[1]) + os.path.sep
    else:
        root = ''

    ref_data = ascii.read(root + 'apass.csv')
    ref_coords = coord.SkyCoord(ref_data['radeg'] * u.degree, 
                                ref_data['decdeg'] * u.degree)

    lcs = np.array([[]])
    lcs_err = np.array([[]])
    times = np.array([])
    for image in iglob(root + 'TJO*.cat'):
        data = ascii.read(image)
        coords = coord.SkyCoord(data['ALPHA_J2000'] * u.degree, 
                                data['DELTA_J2000'] * u.degree)

        matched, sep, _ = ref_coords.match_to_catalog_sky(coords)  

        mags = np.array([])
        ref_mags = np.array([])
        mags_err = np.array([])
        ref_mags_err = np.array([])
        for i, idx in enumerate(matched):
            if sep[i] < 2 * u.arcsec:
                ref_mags = np.append(ref_mags, float(ref_data['Johnson_B'][i]))
                ref_mags_err = np.append(ref_mags_err, 
                                         float(ref_data['B_err'][i]))
                mags = np.append(mags, float(data['MAG_AUTO'][idx]))
                mags_err = np.append(mags_err, float(data['MAGERR_AUTO'][idx]))
            else:
                ref_mags = np.append(ref_mags, np.nan)
                ref_mags_err = np.append(ref_mags_err, np.nan)
                mags = np.append(mags, np.nan)
                mags_err = np.append(mags_err, np.nan)

        nstars = len(mags[~np.isnan(mags)]) 
        if nstars > 10:
            times = np.append(times, 
                      float(image.replace(root + 'TJO','').replace('.cat','')))
            zp, zp_err = get_zp(ref_mags[~np.isnan(ref_mags)], 
                                mags[~np.isnan(mags)], 
                                ref_mags_err[~np.isnan(ref_mags_err)],
                                mags_err[~np.isnan(mags_err)],
                                os.path.basename(image))

            mags += zp
            mags_err = np.sqrt(mags_err*mags_err + zp_err*zp_err)
            try:
                lcs = np.vstack([lcs, [mags]])
                lcs_err = np.vstack([lcs_err, [mags_err]])
            except ValueError:
                lcs = mags
                lcs_err = mags_err


    for idx, (mag, error) in enumerate(zip(lcs.transpose(), 
                                            lcs_err.transpose())):
        if len(mag[~np.isnan(mag)]) > 10 and mag[0] < 15:
            print(idx, ref_coords[idx].to_string('hmsdms'), end=' ')
            print(np.nanmean(mag), np.nanstd(mag), np.nanmean(error))
            ascii.write([times, mag, error], 'lc_' + str(idx) + '.dat',
                        format='no_header')
            plt.errorbar(times, mag, yerr=error, fmt='o')
    plt.show()
