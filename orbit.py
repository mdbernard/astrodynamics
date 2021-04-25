import numpy as np


def calc_v(mu, r, a):
    ''' Calculate velocity at a point in an orbit given
    the gravitational parameter of the primary body, radius
    at the point of interest, and the semi-major axis of
    the orbit. '''
    return np.sqrt(2*(mu/r - mu/2/a))

def calc_a(h, mu, e):
    ''' Calculate semi-major axis given specific angular momentum,
    gravitational parameter of primary body, and eccentricity. '''
    return (h**2/mu)*(1/(1 - e**2))


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
    h_ = np.cross(r_, v_)
    h = np.linalg.norm(h_)

    # Inclination
    i = np.arccos(h_[2]/h)

    # Right Ascension
    Khat_ = np.array([0.0, 0.0, 1.0])
    N_ = np.cross(Khat_, h_)
    N = np.linalg.norm(N_)

    if N != 0:
        RA = np.arccos(N_[0]/N)
        RA = RA if N_[1] >= 0 else 2*np.pi - RA
    else:
        RA = 0.0

    # Eccentricity
    e_ = (1/mu)*(np.cross(v_, h_) - mu*r_/r)
    e = np.linalg.norm(e_)

    # Argument of Periapsis
    if N != 0:
        if e > 1e-7:
            AP = np.arccos(np.dot(N_/N, e_/e))
            AP = AP if e_[2] >= 0 else 2*np.pi - AP        
        else:
            AP = 0.0
    else:
        AP = 0.0

    # True Anomaly
    if e > 1e-7:
        TA = np.arccos(np.dot(e_/e, r_/r))
        TA = TA if vr >= 0 else 2*np.pi - TA
    else:
        cp = np.cross(N_, r_)
        TA = np.arccos(np.dot(N_, r_)/N/r)
        TA = TA if cp[2] >= 0 else 2*np.pi - TA

    return h, i, e, RA, AP, TA


def calc_v_infinity(mu, R1, R2):
    ''' Calculate the hyperbolic excess speed for an escape trajectory
    from planet 1 with circular orbit radius R1 about the primary with
    gravitational parameter mu headed toward planet 2 with radius R2
    about the primary.
    :param R1: `float` (km) circular orbital radius, planet 1
    :param R2: `float` (km) circular orbital radius, planet 2
    :return: `float` (km/s) hyperbolic excess speed of escape trajectory
    '''
    return np.sqrt((mu/R1)*(np.sqrt((2*R2)/(R1 + R2)) - 1))


def calc_dv_hohmann_common_apseline(mu, ra1, rp1, ra2, rp2):
    ''' Calculate total delta-V (m/s) required for a Hohmann transfer between
    two coplanar orbits with a common apse line.
    :param mu: `float` (km**3/s**2) gravitational parameter of primary body
    :param ra1: `float` (km) apoapsis radius, initial orbit
    :param rp1: `float` (km) periapsis radius, initial orbit
    :param ra2: `float` (km) apoapsis radius, final orbit
    :param rp2: `float` (km) periapsis radius, final orbit
    :return: `float` (m/s) total delta-V required for maneuver
    '''
    a1 = 0.5*(ra1 + rp1)
    ai = 0.5*(ra2 + rp1)
    a2 = 0.5*(ra2 + rp2)

    vp1 = calc_v(mu, rp1, a1)
    vpi = calc_v(mu, rp1, ai)
    vai = calc_v(mu, ra2, ai)
    va2 = calc_v(mu, ra2, a2)

    dv1 = vpi - vp1  # (km/s) delta-V of first burn
    dv2 = va2 - vai  # (km/s) delta-V of second burn

    return (np.abs(dv1) + np.abs(dv2))  # (km/s) total maneuver delta-V
