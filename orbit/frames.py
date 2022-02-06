import numpy as np
from numpy import sin, cos


def perifocal_to_geocentric_equatorial(RA, i, AP, r_p):
    ''' Transform a vector's frame from the perifocal frame
    to the geocentric equatorial frame. Assumes no correction
    for Coriolis Effect necessary.

    :param RA: `float` (rad) right ascension
    :param i: `float` (rad) inclination
    :param AP: `float` (rad) argument of periapsis
    :param r_p: `np.array, 1x3` (--) vector in perifocal frame
    '''
    Q = np.array([
        [-sin(RA)*cos(i)*sin(AP) + cos(RA)*cos(AP), -sin(RA)*cos(i)*cos(AP) - cos(RA)*sin(AP), sin(RA)*sin(i)],
        [cos(RA)*cos(i)*sin(AP) + sin(RA)*cos(AP), cos(RA)*cos(i)*cos(AP) - sin(RA)*sin(AP), -cos(RA)*sin(i)],
        [sin(i)*sin(AP), sin(i)*cos(AP), cos(i)]
    ])

    return np.dot(Q, r_p)


def geocentric_equatorial_to_perifocal(RA, i, AP, r_ge):
    ''' Transform a vector's frame from the geocentric equatorial
    to the perifocal frame. Assumes no correction for the Coriolis
    Effect necessary.

    :param RA: `float` (rad) right ascension
    :param i: `float` (rad) inclination
    :param AP: `float` (rad) argument of periapsis
    :param r_ge: `np.array, 1x3` (--) vector in geocentric equatorial frame
    '''
    Q = np.array([
        [-sin(RA)*cos(i)*sin(AP) + cos(RA)*cos(AP), cos(RA)*cos(i)*sin(AP) + sin(RA)*cos(AP), sin(i)*sin(AP)],
        [-sin(RA)*cos(i)*cos(AP) - cos(RA)*sin(AP), cos(RA)*cos(i)*cos(AP) - sin(RA)*sin(AP), sin(i)*cos(AP)],
        [sin(RA)*sin(i), -cos(RA)*sin(i), cos(i)]
    ])

    return np.dot(Q, r_ge)