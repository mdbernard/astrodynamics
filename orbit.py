import numpy as np


def calc_a()


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


def find_orbital_elements(r_, v_, mu):
    ''' Calculate all six orbital elements given a position and
    velocity vector at some point in the orbit, and the gravitational
    parameter of the primary body. A variable with a trailing
    underscore denotes a vector: x_ is a vector, x is the magnitude
    of x_. Applies Algorithm 4.2 from Orbital Mechanics for Engineering
    Students, 4 ed, Curtis. Returns angles in radians. '''
    r = np.linalg.norm(r_)
    v = np.linalg.norm(v_)
    vr = np.dot(r, v)/r

    # Specific Angular Momentum
    h_ = np.cross(r, v)
    h = np.linalg.norm(h_)

    # Inclination
    i = np.arccos(h_[2]/h)

    # Right Ascension
    Khat_ = np.array([0.0, 0.0, 1.0])
    N_ = np.cross(Khat_, h_)
    N = np.linalg.norm(N_)

    RA = np.arccos(N_[2]/N)
    RA = RA if N_[1] >= 0 else 2*np.pi - RA

    # Eccentricity
    e_ = (1/mu)*(np.cross(v_, h_) - mu*r_/r)
    e = np.linalg.norm(e_)

    # Argument of Periapsis
    AP = np.arccos(np.dot(N_/N, e_/e))
    AP = AP if e_[2] >= 0 else 2*np.pi - AP

    # True Anomaly
    TA = np.arccos(np.dot(e_/e, r_/r))
    TA = TA if vr >= 0 else 2*np.pi - TA

    return h, i, RA, e, AP, TA
