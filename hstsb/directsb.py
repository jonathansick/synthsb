#!/usr/bin/env python
# encoding: utf-8
"""
Module for directly estimating SB by adding up the light of stars in the Brown
catalog and transforming that total light into an SDSS magnitude.
"""

import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import ICRS
from astropy import units as u
from astropy import log
import astropy.io.fits as fits

from sqlalchemy.orm import aliased
from starplex.database import connect_to_server, Session

from starplex.database import Catalog, Bandpass, CatalogStar, Observation

from m31hst import brown_image_path, brown_phot_path

from hstsb.galradii import correct_rgc


def main():
    log.setLevel("INFO")
    fieldname = "disk"
    data = load_photometry(fieldname)
    A = compute_area(fieldname)
    sb606 = compute_sb(data['cfrac'], data['m606'], A)
    sb814 = compute_sb(data['cfrac'], data['m814'], A)
    R = compute_gal_radius(fieldname)
    log.info("R {:.4f} kpc".format(R.kpc))
    log.info("mu_606: {:.6f}".format(sb606))
    log.info("mu_814: {:.6f}".format(sb814))


def load_photometry(fieldname):
    """Get photometry from starplex."""
    connect_to_server('marvin', echo=True)
    session = Session()
    mag606obs = aliased(Observation)
    mag814obs = aliased(Observation)
    bp606 = aliased(Bandpass)
    bp814 = aliased(Bandpass)
    q = session.query(CatalogStar.cfrac, mag606obs.mag, mag814obs.mag)\
            .join(mag606obs, CatalogStar.observations)\
            .join(mag814obs, CatalogStar.observations)\
            .join(Catalog)\
            .filter(Catalog.name == fieldname)\
            .join(bp606, mag606obs.bandpass)\
            .filter(bp606.name == "f606w")\
            .join(bp814, mag814obs.bandpass)\
            .filter(bp814.name == "f814w")
    dt = [('cfrac', np.float), ('m606', np.float), ('m814', np.float)]
    data = np.array(q.all(), dtype=np.dtype(dt))
    log.info("Field {0} has {1:d} stars".format(fieldname, data.shape[0]))
    session.close()
    return data


def compute_area(fieldname):
    """Get the unmasked area for this field from the MSK image."""
    image_path = brown_image_path(fieldname, "f606w")
    msk_path = brown_phot_path(fieldname, kind="msk")
    header = fits.getheader(image_path)
    msk_pixels = fits.getdata(msk_path)
    pix_scale = np.sqrt(header['CD1_1'] ** 2. + header['CD1_2'] ** 2.) * 3600.
    # Masked pixels have values of 1, so nx*ny - sum(msk) gives N unmasked pix
    npix = msk_pixels.shape[0] * msk_pixels.shape[1] - msk_pixels.sum()
    return npix * pix_scale


def compute_sb(cfrac, mag, A):
    """Compute a surface brightness for a single bandpass from the sum of
    fluxes of individual stars, given a completeness estimate.

    Returns
    -------
    sb : float
        mag per square arcsecond
    """
    s = np.where(cfrac > 0.)[0]
    sb = -2.5 * np.log10(np.sum(10. ** (-0.4 * mag[s]) / cfrac[s]) / A)
    return sb


def compute_gal_radius(fieldname):
    """Compute galactocentric radius for this ACS field."""
    image_path = brown_image_path(fieldname, "f606w")
    header = fits.getheader(image_path)
    wcs = WCS(header)
    footprint = wcs.calcFootprint()
    ra0 = np.array([v[0] for v in footprint]).mean()
    dec0 = np.array([v[1] for v in footprint]).mean()
    coord = ICRS(ra=ra0, dec=dec0, unit=(u.degree, u.degree))
    return correct_rgc(coord)


if __name__ == '__main__':
    main()
