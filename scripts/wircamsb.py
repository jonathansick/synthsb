#!/usr/bin/env python
# encoding: utf-8
"""
SB computations from WIRCam fields.

To process all 07B fields::

    python synthsb/synthsb/wircamsb.py M31-{1..27}

To process all fields::

    python synthsb/synthsb/wircamsb.py M31-{1..27} M31-{47..71} \
        M31-sky-{28..31} M31-skyr-{01,13,27,39} --no-intercal
"""

import argparse
from collections import defaultdict
import math

import numpy as np
from astropy import log
import astropy.io.fits as fits
import astropy.wcs
from astropy.table import Table

from sqlalchemy import func
from sqlalchemy.orm import aliased
from starplex.database import connect_to_server, Session
from starplex.database import Catalog, Bandpass, CatalogStar, Observation

from androphotsys import wircam_vega_to_ab

# from starplex.utils.timer import Timer

from synthsb.directsb import compute_sb, NoDataError


def main():
    log.setLevel("INFO")

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fields",
        nargs='*',
        help="WIRCam field names(s)")
    parser.add_argument(
        '-n',
        type=int, default=8,
        help="Number of segments on each side of a WIRCam field")
    parser.add_argument(
        '--vega',
        action='store_true', default=False,
        help="Present SB in VEGAMAG, rather than ABMAG.")
    parser.add_argument(
        '--no-intercal',
        dest='disable_intercal',
        action='store_true', default=False,
        help="Disable starplex.intercal ZP corrections")

    args = parser.parse_args()

    for field in args.fields:
        process_wircam_field(field, args.n, args.vega, args.disable_intercal)


def process_wircam_field(fieldname, n_seg, use_vega, disable_intercal):
    """Process a single WIRCam field."""
    log.info("Processing {0}".format(fieldname))

    # Get path to J-band reference image
    refpath = "/Volumes/Zaphod/m31/pipe/blockphot/wircam/"\
        "blocks/{}_J_nightset.fits".format(fieldname)
    # refwpath = "/Volumes/Zaphod/m31/pipe/blockphot/blocks/{}_J_nightset."\
    #     "weight.fits".format(fieldname)

    cols = defaultdict(list)
    with fits.open(refpath) as ref_fits:
        header = ref_fits[0].header
        pix_scale = np.sqrt(
            header['CD1_1'] ** 2. + header['CD1_2'] ** 2.) * 3600.
        img = ref_fits[0].data
        for i, (radec_seg, yx_seg) in enumerate(image_segments(header, n_seg)):
            # process each RA, Dec subsection
            area = compute_area(img, yx_seg, pix_scale)  # sq arcsec
            try:
                phot = load_photometry(
                    fieldname, i, radec_seg, use_vega, disable_intercal)
            except NoDataError:
                continue
            ra0, dec0 = segment_center(radec_seg)
            cols['field'].append(fieldname)
            cols['tile'].append(i)
            cols['nstars'].append(len(phot))
            cols['ra'].append(ra0)
            cols['dec'].append(dec0)
            for band in ['J', 'Ks']:
                _sb, _err = compute_sb(phot[band],
                                       phot["{0}_err".format(band)],
                                       phot['cfrac'],
                                       area)
                cols[band].append(_sb)
                cols["{}_err".format(band)].append(_err)

    # Write out table of all segments for this field
    names = ['field', 'tile', 'nstars', 'ra', 'dec',
             'J', 'J_err', 'Ks', 'Ks_err']
    collist = [cols[k] for k in names]
    tbl = Table(collist, names=names)
    tbl = tbl[tbl['nstars'] >= 100]  # cut off low-density regions
    path_root = "wircam_sb_{0}".format(fieldname)
    if use_vega:
        path_root += "_vega"
    else:
        path_root += "_ab"
    if disable_intercal:
        path_root += "_nointercal"
    else:
        path_root += "_intercal"
    tbl.write(path_root + ".txt", format='ascii.commented_header')


def segment_center(poly):
    """docstring for segment_center"""
    ra0 = np.array([v[0] for v in poly]).mean()
    dec0 = np.array([v[1] for v in poly]).mean()
    return ra0, dec0


def compute_area(img, yxpoly, pix_scale):
    """Compute unmasked area in sq. arcseconds."""
    x = [v[1] for v in yxpoly]
    y = [v[0] for v in yxpoly]
    x1 = min(x)
    x2 = max(x)
    y1 = min(y)
    y2 = max(y)
    seg = img[y1:y2, x1:x2]
    badpix = np.where(np.isfinite(seg) == False)
    npix = seg.shape[0] * seg.shape[1] - len(badpix)
    return npix * pix_scale ** 2.


def load_photometry(fieldname, tile, radec_poly,
                    use_vega, disable_intercal,
                    server='marvin'):
    """Load WIRCam photometry from Starplex, converted to AB mag."""
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
    q = session.query(CatalogStar.cfrac, mag1obs.mag, mag2obs.mag_err,
                      mag2obs.mag, mag2obs.mag_err).\
        join(mag1obs, CatalogStar.observations).\
        join(mag2obs, CatalogStar.observations).\
        join(Catalog).\
        filter(Catalog.name == fieldname).\
        filter(Catalog.instrument == instrument).\
        join(bp1, mag1obs.bandpass).\
        filter(bp1.name == "J").\
        join(bp2, mag2obs.bandpass).\
        filter(bp2.name == "Ks").\
        filter(func.q3c_poly_query(CatalogStar.ra,
                                   CatalogStar.dec,
                                   np.array(radec_poly).flatten().tolist()))
    dt = [('cfrac', np.float), ('J', np.float), ('J_err', np.float),
            ('Ks', np.float), ('Ks_err', np.float)]
    data = np.array(q.all(), dtype=np.dtype(dt))

    # Apply the intercal ZP correction
    if not disable_intercal:
        if 'intercal' in catalog.meta:
            for band in ['J', 'Ks']:
                if band in catalog.meta['intercal']:
                    data[band] += catalog.meta['intercal'][band]['zp']

    # Convert to AB
    if not use_vega:
        data['J'] = wircam_vega_to_ab(data['J'], "J")
        data['Ks'] = wircam_vega_to_ab(data['Ks'], "Ks")

    log.info("Field {0} {2} {3:d} has {1:d} stars".
        format(fieldname, data.shape[0], instrument, tile))
    session.close()
    if len(data) == 0:
        raise NoDataError
    return data


def image_segments(header, n):
    """Generates RA,Dec and yx polygon segments if the image is divided into
    n-by-n square segments.
    """
    wcs = astropy.wcs.WCS(header)
    nx = header['NAXIS1']
    ny = header['NAXIS2']
    dx = int(math.floor(float(nx) / n))
    dy = int(math.floor(float(ny) / n))
    x0, y0 = 0, 0
    for i in xrange(n):
        x0 = 0

        for j in xrange(n):
            x = [x0, x0, x0 + dx, x0 + dx]
            y = [y0, y0 + dx, y0 + dy, y0]
            yx_poly = zip(y, x)
            ra, dec = wcs.wcs_pix2world(np.array(x), np.array(y), 0)
            radec_poly = zip(ra, dec)
            yield radec_poly, yx_poly
            
            x0 += dx

        y0 += dx

if __name__ == '__main__':
    main()
