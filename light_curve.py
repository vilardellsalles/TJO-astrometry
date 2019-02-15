"""
Plot light curves for all stars in a list of SExtractor files
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from re import findall
import argparse
import os.path

import numpy as np

from bokeh.palettes import plasma
from bokeh.models.annotations import Whisker, Span
from bokeh.models.sources import ColumnDataSource
from bokeh.resources import INLINE
import bokeh.plotting as plt

from astropy.modeling import models, fitting
from astropy.stats import histogram
from astropy.units import arcsec
from astropy.time import Time
from astropy import table
import astropy.coordinates as coord

from icat.logs import get_logger, LogOrProgressBar
from icat.config import config_section
from icat import system
import icat.pipelines as ppl


logger = get_logger(__name__)


def multiple_observations(tbl, keys):
    """
    Return whether light curve has more than one measurement
    """

    return len(tbl) > 1


def run(input_catalogs, recursive=False, files_filter="trl.dat", max_sep=5,
        max_chi=9):
    """
    Pipeline execution function
    """

    cat_files = list(system.path_flatten(input_catalogs, recursive,
                                         files_filter=files_filter))

    logger.info("Reading {} catalogs...".format(len(cat_files)))

    max_len = 0
    catalogs = []
    ref_idx = None
    for idx, cat in enumerate(LogOrProgressBar(cat_files)):
        time = float(findall("\d+\.\d+", cat)[0])
        cat_table = table.Table.read(cat, format="ascii.sextractor")
        cat_table["JD"] = time
        cat_table["Name"] = cat
        catalogs += [cat_table]
        cat_len = len(cat_table)
        if cat_len > max_len:
            max_len = cat_len
            ref_idx = idx

    if ref_idx is None:
        raise ppl.PipelineError("No valid catalogs found!")

    site = config_section("Location")
    OAdM = coord.EarthLocation(lat=site["latitude"], lon=site["longitude"],
                               height=site["elevation"])

    ref_file = cat_files.pop(ref_idx)
    ref_path = os.path.dirname(ref_file)
    ref_catalog = catalogs.pop(ref_idx)
    ref_time = int(ref_catalog["JD"][0])
    ref_pos = coord.SkyCoord(ref_catalog["ALPHA_J2000"],
                             ref_catalog["DELTA_J2000"])
    jd = Time(ref_catalog["JD"], format="jd")
    ref_altaz = ref_pos.transform_to(coord.AltAz(obstime=jd, location=OAdM))

    log_msg = "Matching {} catalogs with '{}'..."
    logger.info(log_msg.format(len(catalogs), ref_file))

    num_stars = len(ref_catalog)
    photometry = ref_catalog.copy()
    photometry["Star"] = list(range(num_stars))
    photometry["JD"] = photometry["JD"]
    photometry["Airmass"] = ref_altaz.secz

    all_seps = []
    for cat in LogOrProgressBar(catalogs):
        cat_pos = coord.SkyCoord(cat["ALPHA_J2000"], cat["DELTA_J2000"])
        jd = Time(cat["JD"], format="jd")
        altaz = cat_pos.transform_to(coord.AltAz(obstime=jd, location=OAdM))
        idx, sep, _ = cat_pos.match_to_catalog_sky(ref_pos)
        cat["Star"] = idx
        cat["Airmass"] = altaz.secz
        photometry = table.vstack([photometry, cat[sep < max_sep*arcsec]])
        all_seps += sep.arcsec.tolist()

    y, x = histogram(all_seps, bins="knuth")
    hist = plt.figure()
    hist.quad(top=y, bottom=0, left=x[:-1], right=x[1:])
    hist.add_layout(Span(location=max_sep, dimension="height"))
    hist_name = os.path.join(ref_path, "separation.html")
    plt.save(hist, filename=hist_name, title="Matched stars separation",
             resources=INLINE)

    logger.info("Separation histogram saved in '%s'", hist_name)

    light_curves = photometry.group_by("Star")
    light_curves = light_curves.groups.filter(multiple_observations)

    log_msg = "Selecting constant stars among %s light curves..."
    logger.info(log_msg, len(light_curves.groups))

    min_mag = light_curves["MAG_AUTO"] - light_curves["MAGERR_AUTO"]
    max_mag = light_curves["MAG_AUTO"] + light_curves["MAGERR_AUTO"]
    errors = dict(base=light_curves["Airmass"], lower=min_mag, upper=max_mag)
    error_bars = ColumnDataSource(data=errors)

    y_axis = (max(light_curves["MAG_AUTO"]) + 0.1,
              min(light_curves["MAG_AUTO"]) - 0.1)

    chi_plot = plt.figure(tools="hover,box_zoom,undo,redo,reset,save",
                          y_range=y_axis)
    chi_plot.add_layout(Whisker(source=error_bars, base="base",
                                lower="lower", upper="upper"))
    palette = plasma(num_stars)

    ref_stars = {}
    for star in LogOrProgressBar(light_curves.groups):
        color = star["Star"][0]

        chi_plot.line(star["Airmass"], star["MAG_AUTO"],
                      line_color=palette[color])

        linear_model = models.Linear1D(bounds={"slope": (0, 1)})
        linear_fitter = fitting.LevMarLSQFitter()
        weights = 1 / star["MAGERR_AUTO"]**2
        linear_fit = linear_fitter(linear_model, star["Airmass"],
                                   star["MAG_AUTO"], weights=weights)

        chi_plot.line(star["Airmass"], linear_fit(star["Airmass"]),
                      line_color=palette[color], line_dash="dashed")

        residuals = (star["MAG_AUTO"] - linear_fit(star["Airmass"]))**2
        num_obs = len(residuals)
        chi2 = sum(weights * residuals) / (num_obs - 2) if num_obs > 2 else 0

        logger.debug("Star %s (%s measurements): %s +/- %s * X (Chi2: %s)",
                     color, num_obs, linear_fit.intercept.value,
                     linear_fit.slope.value, chi2)

        if 0 < chi2 < max_chi:
            ref_stars[color] = num_obs

    airmass_name = os.path.join(ref_path, "airmass.html")
    plt.save(chi_plot, filename=airmass_name, title="Airmass",
             resources=INLINE)

    logger.info("Airmass curves saved in '%s'", airmass_name)

    max_obs = max(ref_stars.values())
    ref_stars = [key for key, value in ref_stars.iteritems()
                 if value == max_obs]
    num_standards = len(ref_stars)

    logger.info("Computing differential magnitudes with %s standard stars",
                len(ref_stars))

    light_curves["Standard"] = [star in ref_stars
                                for star in light_curves["Star"]]
    light_curves["FLUX_AUTO"] = pow(10, -0.4*light_curves["MAG_AUTO"])

    ref_mag = {}
    ref_err = {}
    grouped_cats = light_curves.group_by("Name")
    for cat in LogOrProgressBar(grouped_cats.groups):
        if np.count_nonzero(cat["Standard"]) != num_standards:
            continue
        ref_flux = cat["FLUX_AUTO"][cat["Standard"]]
        ref_mag[cat["Name"][0]] = -2.5 * np.log10(sum(ref_flux))
        ref_errflux = ref_flux**2 * cat["MAGERR_AUTO"][cat["Standard"]]**2
        ref_err[cat["Name"][0]] = np.sqrt(sum(ref_errflux)) / sum(ref_flux)

    mag_corr = []
    magerr_corr = []
    for row in light_curves:
        mag_corr += [row["MAG_AUTO"] - ref_mag.get(row["Name"], 0)]
        magerr_corr += [np.sqrt(row["MAGERR_AUTO"]**2 +
                                ref_err.get(row["Name"], 0)**2)]

    light_curves["MAG_CORR"] = mag_corr
    light_curves["MAGERR_CORR"] = magerr_corr
    light_curves["JD_CORR"] = light_curves["JD"] - ref_time

    light_curves.sort("JD_CORR")
    mask = [name in ref_mag.keys() for name in light_curves["Name"]]
    light_curves = light_curves[mask].group_by("Star")

    min_mag = light_curves["MAG_CORR"] - light_curves["MAGERR_CORR"]
    max_mag = light_curves["MAG_CORR"] + light_curves["MAGERR_CORR"]
    errors = dict(base=light_curves["JD_CORR"], lower=min_mag, upper=max_mag)
    error_bars = ColumnDataSource(data=errors)

    y_axis = (max(light_curves["MAG_CORR"]) + 0.1,
              min(light_curves["MAG_CORR"]) - 0.1)

    light_plot = plt.figure(tools="hover,box_zoom,undo,redo,reset,save",
                            y_range=y_axis)
    light_plot.add_layout(Whisker(source=error_bars, base="base",
                                  lower="lower", upper="upper"))

    for star in light_curves.groups:
        color = star["Star"][0]

        light_plot.line(star["JD_CORR"], star["MAG_CORR"],
                        line_color=palette[color])

    curves_name = os.path.join(ref_path, "light_curves.html")
    plt.save(light_plot, filename=curves_name, title="Light curves",
             resources=INLINE)

    logger.info("Light curves saved in '%s'", curves_name)

    result_name = os.path.join(ref_path, "curves.html")
    if os.path.isfile(result_name):
        os.remove(result_name)
    light_curves.write(result_name, format="jsviewer")


def main():
    """
    Entry point
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-r", "--recursive", action="store_true",
                        help="recursively enter in directories")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log debug messages")

    help_msg = "when path is a directory, read only files having FILTER "
    help_msg += "in name (default: '%(default)s')"
    parser.add_argument("-f", "--filter", default="trl.dat", help=help_msg)

    help_msg = "maximum separation to match the same star in different "
    help_msg += "catalogs (in arc seconds)"
    parser.add_argument("-s", "--separation", default=5, help=help_msg)

    help_msg = "maximum chi2 to consider a star to be constant"
    parser.add_argument("-c", "--chi2", default=9, help=help_msg)

    parser.add_argument("path", nargs="+",
                        help="file or directory with SExtractor files to use")

    args = parser.parse_args()

    logger = get_logger(__name__, verbose=args.verbose)

    try:
        run(args.path, args.recursive, args.filter, args.separation, args.chi2)
    except Exception as err:
        logger.exception(err)
        raise


if __name__ == "__main__":
    main()
