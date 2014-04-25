#!/usr/bin/env python
# encoding: utf-8
"""
Module for directly estimating SB by adding up the light of stars in the Brown
catalog and transforming that total light into an SDSS magnitude.
"""

from astropy import log
import numpy as np

from sqlalchemy.orm import aliased
from starplex.database import connect_to_server, Session
from starplex.database import Catalog, Bandpass, CatalogStar, Observation


def main():
    log.setLevel("INFO")
    fieldname = "halo35b"
    data = load_photometry(fieldname)
    A = compute_area(fieldname)
    sb606 = compute_sb(data['cfrac'], data['m606'], A)
    sb814 = compute_sb(data['cfrac'], data['m814'], A)
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
    pass


def compute_sb(cfrac, mag, A):
    """Compute a surface brightness for a single bandpass."""
    pass


if __name__ == '__main__':
    main()
