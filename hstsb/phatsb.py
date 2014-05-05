#!/usr/bin/env python
# encoding: utf-8
"""
SB computation of PHAT fields.
"""

import argparse
from collections import defaultdict

from astropy import log
import astropy.io.fits as fits
from astropy.table import Table

from m31hst import phat_field_path
from androphotsys import ACS_475_814_to_BVRI, BVRI_to_ugri

from hstsb.directsb import load_photometry, compute_field_coord, \
    compute_sb, compute_area_from_weights


def main():
    log.setLevel("INFO")

    parser = argparse.ArgumentParser()
    parser.add_argument("brick", type=int, help="Brick number",
        choices=range(1, 24))
    parser.add_argument("instrument", help="PHAT instrument",
        choices=['phat_acs', 'phat_ir', 'phat_uv'])
    args = parser.parse_args()

    if args.instrument == 'phat_acs':
        process_acs(args.brick)
    elif args.instrument == 'phat_ir':
        pass
    elif args.instrument == 'phat_uv':
        pass


def process_acs(brick):
    """Process ACS observations of this brick."""
    cols = defaultdict(list)
    for fieldnum in xrange(1, 19):
        result = process_acs_field(brick, fieldnum)
        for k, v in result.iteritems():
            cols[k].append(v)
    names = ['name', 'brick', 'field', 'ra', 'dec', 'radius',
            'f475w', 'f475w_err',
            'f814w', 'f814w_err', 'g', 'g_err', 'r', 'r_err', 'i', 'i_err']
    collist = [cols[k] for k in names]
    tbl = Table(collist, names=names)
    tbl.write("phat_acs_sb_{0:02d}.txt".format(brick),
        format='ascii.commented_header')


def process_acs_field(brick, fieldnum):
    """Compute ACS SB in this brick+field."""
    fieldname = "phat_b{0:02d}_f{1:02d}".format(brick, fieldnum)
    data = load_photometry(fieldname, 'phat_acs', 'f475w', 'f814w')

    image_path = phat_field_path(brick, fieldnum, "f814w")
    header = fits.getheader(image_path, 'SCI')
    whtdata = fits.getdata(image_path, 'WHT')
    A = compute_area_from_weights(whtdata, header)
    
    sb475, e475 = compute_sb(data['m1'], data['e1'], data['cfrac'], A)
    sb814, e814 = compute_sb(data['m2'], data['e2'], data['cfrac'], A)
    B, Be, V, Ve, R, Re, I, Ie = ACS_475_814_to_BVRI(sb475, e475,
        sb814, e814, zp='VEGAMAG')
    u, u_e, g, g_e, r, r_e, i, i_e = BVRI_to_ugri(B, Be, V, Ve, R, Re, I, Ie)
    ra0, dec0, rad = compute_field_coord(header)
    log.info(fieldname)
    return {"name": fieldname, "brick": brick, "field": fieldnum,
            "ra": ra0, "dec": dec0, "radius": rad.kpc,
            "f475w": sb475, "f814w": sb814, "g": g, "r": r, "i": i,
            "f475w_err": e475, "f814w_err": e814, "g_err": g_e, "r_err": r_e,
            "i_err": i_e, "u": u, "u_err": u_e}


if __name__ == '__main__':
    main()
