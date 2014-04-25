#!/usr/bin/env python
# encoding: utf-8
"""
Photometric transformations between ACS magnitudes, Johnson/Cousins and
the SDSS system.
"""


def VRI_from_ACS(m606_stmag, m814_stmag):
    """Uses the transformations presented by Sirianni 2005 Table 18.
    Assumes that ACS magnitudes are in the STMAG system.
    
    See equation solutions at http://www.sagenb.org/home/pub/5048
    """
    c0v = 26.325  # pm 0.057
    c1v = 0.236  # pm 0.058
    # c2v = 0.
    c0i = 25.495  # pm 0.015
    c1i = -0.002  # pm 0.017
    # c2i = 0.
    c0ri = 25.492  # pm 0.013
    c1ri = 0.002  # pm 0.003

    # Subtract the STMAG zeropoints to get OBMAG (Sirianni 2005 Table 10)
    m606 = m606_stmag - 26.655
    m814 = m814_stmag - 26.776

    V = (m606 * (c1i + 1) + c0v * (c1i + 1) - (m814 + c0i) * c1v) \
            / (c1i - c1v + 1)
    I = (m606 * c1i + c0v * c1i - (m814 + c0i) * c1v + m814 + c0i) \
            / (c1i - c1v + 1)
    R = (c0v * c1i * (c1ri + 1) + c1i * (c1ri + 1) * m606 - (m814 + c0ri) * c1i
        + c0i * (c1ri + 1)
        - (c0i * (c1ri + 1) + (c1ri + 1) * m814 - m814 - c0ri) * c1v
        + (c1ri + 1) * m814 - m814 - c0ri) / (c1i * c1ri - c1ri * c1v + c1ri)
    return V, R, I


def gri_from_VRI(V, R, I):
    """Transformation of Johnson Cousins VRI magnitudes to SDSS AB mags.
    
    Based on Lupton 2005
    http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php#Lupton2005

    See equation solutions at http://www.sagenb.org/home/pub/5049
    """
    g = -5784. / 6053. * R + 11837. / 6053. * V - 2583229. / 30265000.
    r = 4216. / 6053. * R + 1837. / 6053. * V + 2081771. / 30265000.
    i = 2500. / 3111. * I + 151528. / 1107699. * R \
        + 1122407. / 18830883. * V + 30175037081. / 94154415000.
    return g, r, i
