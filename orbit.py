import numpy as np


def calc_T(mu, a):
    ''' Calculate orbital period given gravitational parameter
    and semi-major axis. '''
    return 2*np.pi*a**1.5/np.sqrt(mu)


def calc_RA_dot(mu, J2, R, e, a, i):
    ''' Calculate time rate of change of Right Ascension (RA)
    of an orbit given gravitational parameter, J2, equatorial
    radius of primary body, eccentricity, semi-major axis,
    and inclination. '''
    return -1.5*(np.sqrt(mu)*J2*R**2)/((1 - e**2)**2*a**3.5)*np.cos(i)


def calc_AP_dot_1(mu, J2, R, e, a, i):
    ''' Calculate time rate of change of Argument of Periapsis (AP)
    given gravitational parameter, J2, equatorial radius of primary
    body, eccentricity, semi-major axis, and inclination. '''
    return -1.5*(np.sqrt(mu)*J2*R**2)/((1 - e**2)**2*a**3.5)*(2.5*(np.sin(i))**2 - 2)


def calc_AP_dot_2(RAdot, i):
    ''' Calculate time rate of change of Argument of Periapsis (AP)
    given time rate of change of Right Ascension (RAdot) and
    inclination. '''
    return RAdot*(2.5*(np.sin(i))**2 - 2)/(np.cos(i))