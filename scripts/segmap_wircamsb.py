#!/usr/bin/env python
# encoding: utf-8
"""
Compute WIRCam synthetic surface brightnesses within regions of a segmentation
map.

Accepts segmentation maps and pixel tables made by, e.g. andromass.

2014-11-18 - Created by Jonathan Sick
"""

import argparse
# from collections import defaultdict
# import math

import numpy as np
from astropy import log
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.table import Table

# from sqlalchemy import func
from sqlalchemy.orm import aliased
from starplex.database import connect_to_server, Session
from starplex.database import Catalog, Bandpass, CatalogStar, Observation

from androphotsys import wircam_vega_to_ab

# from starplex.utils.timer import Timer

from synthsb.directsb import NoDataError
# from synthsb.directsb import compute_sb


def main():
    log.setLevel("INFO")
    args = parse_args()

    segmap_fits = fits.open(args.seg_path)
    segmap = segmap_fits[0].data
    wcs = WCS(segmap_fits[0].header)
    pixel_table = Table.read(args.pix_table_path,
                             format='ascii.commented_header')
    fluxsum_J = np.full(len(pixel_table), 0, dtype=np.float)
    varsum_J = np.full(len(pixel_table), 0, dtype=np.float)
    fluxsum_Ks = np.full(len(pixel_table), 0, dtype=np.float)
    varsum_Ks = np.full(len(pixel_table), 0, dtype=np.float)
    star_count = np.zeros(len(pixel_table), dtype=np.int)

    fields = ["M31-{0:d}".format(i) for i in range(1, 28)] + \
             ["M31-{0:d}".format(i) for i in range(47, 72)]
    # fields = ['M31-1']
    for field in fields:
        print "Processing", field
        data = load_photometry(field)
        x, y = wcs.wcs_world2pix(data['ra'], data['dec'], 0)
        # Round down to pixel indices
        x = x.astype(np.int)
        y = y.astype(np.int)
        # Filter out stars contained inside the image footprint
        ny, nx = segmap.shape
        s = np.where((x >= 0) & (y >= 0) &
                     (x < nx) & (y < ny) &
                     np.isfinite(data['J']) & np.isfinite(data['Ks']) &
                     np.isfinite(data['J_err']) & np.isfinite(data['Ks_err']) &
                     (data['cfrac'] > 0.))[0]
        data = data[s]
        n_stars = data.shape[0]
        flux_J, flux_var_J = mag_to_mjy(data['J'], data['J_err'])
        flux_Ks, flux_var_Ks = mag_to_mjy(data['Ks'], data['Ks_err'])
        for i in xrange(n_stars):
            bin_id = segmap[y[i], x[i]]
            if bin_id >= 0:
                # add light to bin
                fluxsum_J[bin_id] += flux_J[i] / data['cfrac'][i]
                fluxsum_Ks[bin_id] += flux_Ks[i] / data['cfrac'][i]
                varsum_J[bin_id] += flux_var_J[i]
                varsum_Ks[bin_id] += flux_var_Ks[i]
                star_count[bin_id] += 1

    empty = np.where(star_count == 0)[0]
    fluxsum_J[empty] = np.nan
    fluxsum_Ks[empty] = np.nan
    varsum_J[empty] = np.nan
    varsum_Ks[empty] = np.nan

    flux_err_J = np.sqrt(varsum_J)
    flux_err_Ks = np.sqrt(varsum_Ks)
    pixel_table['n_stars'] = star_count
    pixel_table['synth_J'] = fluxsum_J
    pixel_table['synth_J_err'] = flux_err_J
    pixel_table['synth_Ks_err'] = flux_err_Ks
    pixel_table.write(args.output_path,
                      format='ascii.commented_header')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('pix_table_path')
    parser.add_argument('seg_path')
    parser.add_argument('output_path')
    return parser.parse_args()


def load_photometry(fieldname,
                    use_vega=False, apply_intercal=False,
                    server='marvin'):
    """Load WIRCam photometry from Starplex, converted to AB mag.

    Filter out MW stars with a rudimentary J-Ks > 0.9 (Vega) color cut.
    """
    instrument = "wircam"

    connect_to_server(server)
    session = Session()

    mag1obs = aliased(Observation)
    mag2obs = aliased(Observation)
    bp1 = aliased(Bandpass)
    bp2 = aliased(Bandpass)

    catalog = session.query(Catalog).\
        filter(Catalog.name == fieldname).\
        filter(Catalog.instrument == instrument).\
        one()
    q = session.query(CatalogStar.cfrac, CatalogStar.ra, CatalogStar.dec,
                      mag1obs.mag, mag1obs.mag_err,
                      mag2obs.mag, mag2obs.mag_err).\
        join(mag1obs, CatalogStar.observations).\
        join(mag2obs, CatalogStar.observations).\
        join(Catalog).\
        filter(Catalog.name == fieldname).\
        filter(Catalog.instrument == instrument).\
        join(bp1, mag1obs.bandpass).\
        filter(bp1.name == "J").\
        join(bp2, mag2obs.bandpass).\
        filter(bp2.name == "Ks")
    dt = [('cfrac', np.float), ('ra', np.float), ('dec', np.float),
          ('J', np.float), ('J_err', np.float),
          ('Ks', np.float), ('Ks_err', np.float)]
    data = np.array(q.all(), dtype=np.dtype(dt))

    # Filter out MW stars
    # FIXME rudimentary
    # Using Vega Mag here!
    sel = np.where((data['J'] - data['Ks']) > 0.9)[0]
    data = data[sel]

    # Apply the intercal ZP correction
    if apply_intercal:
        if 'intercal' in catalog.meta:
            for band in ['J', 'Ks']:
                if band in catalog.meta['intercal']:
                    data[band] += catalog.meta['intercal'][band]['zp']

    # Convert to AB
    if not use_vega:
        data['J'] = wircam_vega_to_ab(data['J'], "J")
        data['Ks'] = wircam_vega_to_ab(data['Ks'], "Ks")

    log.info("Field {0} {2} has {1:d} stars".
             format(fieldname, data.shape[0], instrument))
    session.close()
    if len(data) == 0:
        raise NoDataError
    return data


def mag_to_mjy(mag, mag_err):
    MICROJY_ZP = 10. ** 6. * 10. ** 23. * 10. ** (-48.6 / 2.5)
    mjy = MICROJY_ZP * np.power(10., -mag / 2.5)
    mjy_err = (mjy * mag_err) / 1.0875
    return mjy, mjy_err


if __name__ == '__main__':
    main()
