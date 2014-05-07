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

from sqlalchemy.orm import aliased
from starplex.database import connect_to_server, Session
from starplex.database import Catalog, Bandpass, CatalogStar, Observation

from hstsb.galradii import correct_rgc


class NoDataError(BaseException):
    pass


def load_photometry(fieldname, instrument, band1, band2, server="marvin"):
    """Get photometry from starplex."""
    connect_to_server(server, echo=True)
    session = Session()
    # potentially a factor 10 speed gain with this
    session.execute('SET ENABLE_NESTLOOP TO FALSE')
    mag1obs = aliased(Observation)
    mag2obs = aliased(Observation)
    bp1 = aliased(Bandpass)
    bp2 = aliased(Bandpass)
    q = session.query(CatalogStar.cfrac, mag1obs.mag, mag2obs.mag_err,
                mag2obs.mag, mag2obs.mag_err)\
            .join(mag1obs, CatalogStar.observations)\
            .join(mag2obs, CatalogStar.observations)\
            .join(Catalog)\
            .filter(Catalog.name == fieldname)\
            .filter(Catalog.instrument == instrument)\
            .join(bp1, mag1obs.bandpass)\
            .filter(bp1.name == band1)\
            .join(bp2, mag2obs.bandpass)\
            .filter(bp2.name == band2)
    dt = [('cfrac', np.float), ('m1', np.float), ('e1', np.float),
            ('m2', np.float), ('e2', np.float)]
    data = np.array(q.all(), dtype=np.dtype(dt))
    log.info("Field {0} {2} has {1:d} stars".
        format(fieldname, data.shape[0], instrument))
    session.execute('SET ENABLE_NESTLOOP TO TRUE')
    session.close()
    if len(data) == 0:
        raise NoDataError
    return data


def compute_area(msk_image, header):
    """Get the unmasked area for this field from the MSK image."""
    pix_scale = np.sqrt(header['CD1_1'] ** 2. + header['CD1_2'] ** 2.) * 3600.
    # Masked pixels have values of 1, so nx*ny - sum(msk) gives N unmasked pix
    npix = msk_image.shape[0] * msk_image.shape[1] - msk_image.sum()
    return npix * pix_scale ** 2.


def compute_area_from_weights(weight_image, header):
    """Get area of pixels with weight > 0."""
    pix_scale = np.sqrt(header['CD1_1'] ** 2. + header['CD1_2'] ** 2.) * 3600.
    badpix = np.where(weight_image <= 0.)[0]
    npix = weight_image.shape[0] * weight_image.shape[1] - len(badpix)
    return npix * pix_scale ** 2.


def compute_sb(mag, mag_err, cfrac, A):
    """Compute a surface brightness for a single bandpass from the sum of
    fluxes of individual stars, given a completeness estimate.

    Returns
    -------
    sb : float
        mag per square arcsecond
    sb_err : float
        mag uncertainty per square arcsecond
    """
    # Compute surface brightness given stars with a real completeness
    s = np.where(cfrac > 0.)[0]
    sb = -2.5 * np.log10(np.sum(10. ** (-0.4 * mag[s]) / cfrac[s]) / A)

    # Compute uncertainty given mag_err
    F = 10. ** (-0.4 * mag[s]) / cfrac[s]
    # gradient of d mu / d m_i
    grad = F / F.sum()
    sb_err = np.sqrt(np.sum((grad * mag_err[s]) ** 2.))
    return sb, sb_err


def compute_field_coord(header):
    """Compute galactocentric radius for this ACS field as well as the
    central RA, Dec.
    """
    wcs = WCS(header)
    footprint = wcs.calcFootprint()
    ra0 = np.array([v[0] for v in footprint]).mean()
    dec0 = np.array([v[1] for v in footprint]).mean()
    coord = ICRS(ra=ra0, dec=dec0, unit=(u.degree, u.degree))
    return ra0, dec0, correct_rgc(coord)
