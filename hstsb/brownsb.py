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

from hstsb.phottrans import VRI_from_ACS, gri_from_VRI
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
    names = ['name', 'ra', 'dec', 'radius', '606', '814', 'g', 'r', 'i']
    collist = [cols[k] for k in names]
    tbl = Table(collist, names=names)
    tbl.write("direct_brown_sb.txt", format='ascii.commented_header')


def process_field(fieldname):
    """Compute SB for a single Brown HST field."""
    data = load_photometry(fieldname, 'brown', 'f606w', 'f814w')

    image_path = brown_image_path(fieldname, "f606w")
    msk_path = brown_phot_path(fieldname, kind="msk")
    header = fits.getheader(image_path)
    A = compute_area(fits.getdata(msk_path), header)

    sb606 = compute_sb(data['cfrac'], data['m1'], A)
    sb814 = compute_sb(data['cfrac'], data['m2'], A)
    V, R, I = VRI_from_ACS(sb606, sb814)
    g, r, i = gri_from_VRI(V, R, I)
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
            "606": sb606, "814": sb814, "g": g, "r": r, "i": i}


if __name__ == '__main__':
    main()
