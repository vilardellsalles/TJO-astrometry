from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
import sys
import os.path
from glob import iglob

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
from astropy.table import Table

VAR_RATIO = 5
MAX_SEP = 2 * u.arcsec


def get_lcs(path):
    """ Get light curves from a list of files in a path """

    data = {}
    max_len = 0
    best_image = None
    for image in iglob(os.path.join(path, "TJO*.cat")):
        photometry = ascii.read(image)

        nstars = len(photometry)
        photometry["INDEX"] = [-1] * nstars
        if nstars > max_len:
            max_len = nstars
            best_image = image

        data[image] = photometry

    ref_data = data[best_image]
    ref_data["INDEX"] = np.arange(max_len)
    ra = ref_data["ALPHA_J2000"] * u.deg
    dec = ref_data["DELTA_J2000"] * u.deg
    ref_coords = coord.SkyCoord(ra, dec)

    mags = {}
    for image, photometry in data.iteritems():
        if image != best_image:
            ra = photometry["ALPHA_J2000"] * u.deg
            dec = photometry["DELTA_J2000"] * u.deg
            coords = coord.SkyCoord(ra, dec)

            matched, sep, _ = ref_coords.match_to_catalog_sky(coords)

            for i, idx in enumerate(matched):
                if sep[i] < MAX_SEP:
                    photometry["INDEX"][idx] = i

                    if i not in mags:
                        mags[i] = [ref_data["MAG_AUTO"][i]]
                    mags[i] += [photometry["MAG_AUTO"][idx]]

    variables = []
    old_variables = None
    while old_variables is None or variables not in old_variables:
        print("New iteration with {} variables".format(len(variables)))
        if old_variables:
            old_variables += [variables]
        else:
            old_variables = [variables]

        median = {idx: np.median(curve) for idx, curve in mags.iteritems() 
                  if idx not in variables}
 
        light_curves = {}
        for image, photometry in data.iteritems():
            diffs = [row["MAG_AUTO"] - median[row["INDEX"]] 
                     for row in photometry
                     if row["INDEX"] >= 0 and row["INDEX"] in median]
            if diffs:
                zp = np.median(diffs)
                zp_err = np.percentile(np.abs(diffs-zp), 68.27) 
                zp_err /= np.sqrt(len(diffs))

                sigma = photometry["MAGERR_AUTO"]**2 + zp_err**2
                photometry["MAG_CORR"] = photometry["MAG_AUTO"] - zp
                photometry["MAGERR_CORR"] = np.sqrt(sigma)
 
                for star in photometry:
                    idx = star["INDEX"]
                    if idx >= 0:
                        filename = os.path.basename(image)
                        if idx not in light_curves:
                            light_curves[idx] = Table(star)
                            light_curves[idx]["IMAGE"] = filename
                        else:
                            new_row = Table(star)
                            new_row["IMAGE"] = filename
                            light_curves[idx].add_row(new_row[0])

        variables = []
        for idx, curve in light_curves.iteritems():
            mean_err = np.mean(curve["MAGERR_CORR"])
            dispersion = np.std(curve["MAG_CORR"])
            var_ratio = dispersion / mean_err
            if var_ratio > VAR_RATIO or var_ratio < 1.0/VAR_RATIO:
                variables += [idx]

    return light_curves


if __name__ == '__main__':
    if len(sys.argv) == 2:
        path = os.path.abspath(sys.argv[1])
    else:
        path = ""

    lcs = get_lcs(path)

    columns = ["INDEX", "RA", "DEC", "MAG_CORR", "MAGERR_CORR",
               "MAGSTD_CORR", "VAR_RATIO", "NPOINTS"]
    stars = Table(names=columns, dtype=["object"]*len(columns))
    for idx, curve in lcs.iteritems():
        ra = coord.Angle(curve["ALPHA_J2000"][0] * u.deg)
        dec = coord.Angle(curve["DELTA_J2000"][0] * u.deg)
        mean_err = np.mean(curve["MAGERR_CORR"])
        dispersion = np.std(curve["MAG_CORR"])
        var_ratio = dispersion / mean_err
        row = [idx, ra.to_string(unit=u.hour), dec.to_string(sep="dms"), 
               np.mean(curve["MAG_CORR"]), mean_err, dispersion, var_ratio,
               len(curve)]

        stars.add_row(row)
                                  
        filename = os.path.join(path, "lc_{}.dat".format(idx))
        curve.write(filename, format="ascii.commented_header")

    filename = os.path.join(path, "lcs.dat")
    stars.write(filename, format="ascii.commented_header")
    stars.show_in_browser(jsviewer=True)
