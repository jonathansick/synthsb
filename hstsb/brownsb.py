#!/usr/bin/env python
# encoding: utf-8
"""
SB computation for the Brown data set.
"""

from collections import defaultdict

from astropy import log
import astropy.io.fits as fits
from astropy.table import Table

from m31hst import brown_image_path, brown_phot_path
from androphotsys import ACS_606_814_to_VRI, VRI_to_gri

from hstsb.directsb import load_photometry, compute_field_coord, \
    compute_sb, compute_area


def main():
    log.setLevel("INFO")
    fields = ('halo11', 'stream', 'disk', 'halo21', 'halo35a', 'halo35b')
    cols = defaultdict(list)
    for fieldname in fields:
        result = process_field(fieldname)
        for k, v in result.iteritems():
            cols[k].append(v)
    names = ['name', 'ra', 'dec', 'radius', 'f606w', 'f606w_err',
        'f814w', 'f814w_err', 'g', 'g_err', 'r', 'r_err', 'i', 'i_err']
    collist = [cols[k] for k in names]
    tbl = Table(collist, names=names)
    tbl.write("brown_sb.txt", format='ascii.commented_header')


def process_field(fieldname):
    """Compute SB for a single Brown HST field."""
    data = load_photometry(fieldname, 'brown', 'f606w', 'f814w')

    image_path = brown_image_path(fieldname, "f606w")
    msk_path = brown_phot_path(fieldname, kind="msk")
    header = fits.getheader(image_path)
    A = compute_area(fits.getdata(msk_path), header)

    sb606, e606 = compute_sb(data['m1'], data['e1'], data['cfrac'], A)
    sb814, e814 = compute_sb(data['m2'], data['e2'], data['cfrac'], A)
    V, Ve, R, Re, I, Ie = ACS_606_814_to_VRI(sb606, e606, sb814, e814,
        zp='STMAG')
    g, g_e, r, r_e, i, i_e = VRI_to_gri(V, Ve, R, Re, I, Ie)
    ra0, dec0, rad = compute_field_coord(header)
    log.info(fieldname)
    log.info("R {:.4f} kpc".format(rad.kpc))
    log.info("mu_606: {:.6f}".format(sb606))
    log.info("mu_814: {:.6f}".format(sb814))
    log.info("mu_V: {:.6f}".format(V))
    log.info("mu_R: {:.6f}".format(R))
    log.info("mu_I: {:.6f}".format(I))
    log.info("mu_g: {:.6f}".format(g))
    log.info("mu_r: {:.6f}".format(r))
    log.info("mu_i: {:.6f}".format(i))
    return {"name": fieldname, "ra": ra0, "dec": dec0, "radius": rad.kpc,
            "f606w": sb606, "f814w": sb814, "g": g, "r": r, "i": i,
            "f606w_err": e606, "f814w_err": e814, "g_err": g_e,
            "r_err": r_e, "i_err": i_e}


if __name__ == '__main__':
    main()
