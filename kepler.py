import numpy as np
from stumpff import C, S
from CelestialBody import BODIES
from numerical import newton, laguerre
from lagrange import calc_f, calc_fd, calc_g, calc_gd


def kepler_chi(chi, alpha, r0, vr0, mu, dt):
    ''' Kepler's Equation of the universal anomaly, modified
    for use in numerical solvers. '''
    z = alpha*chi**2
    return (r0*vr0/np.sqrt(mu))*chi**2*C(z) + \
           (1 - alpha*r0)*chi**3*S(z) + \
           r0*chi - np.sqrt(mu)*dt


def dkepler_dchi(chi, alpha, r0, vr0, mu, dt):
    ''' Derivative of Kepler's Equation of the universal anomaly,
    modified for use in numerical solvers. '''
    z = alpha*chi**2
    return (r0*vr0/np.sqrt(mu))*chi*(1 - alpha*chi**2*S(z)) + \
           (1 - alpha*r0)*chi**2*C(z) + r0


def d2kepler_dchi2(chi, alpha, r0, vr0, mu, dt):
    ''' Second derivative of Kepler's Equation of the universal
    anomaly, modified for use in numerical solvers. '''
    z = alpha*chi**2
    S_ = S(z)
    return (r0*vr0/np.sqrt(mu))*(1 - 3*z*S_ + z*(C(z) - 3*S_)) + \
           chi*(1 - z*S_)*(1 - alpha*r0)


def solve_kepler_chi(r_0, v_0, dt, body=BODIES['Earth'], method='laguerre', tol=1e-7, max_iters=100):
    ''' Solve Kepler's Equation of the universal anomaly chi using the specified
    numerical method. Applies Algorithm 3.4 from Orbital Mechanics for Engineering
    Students, 4 ed, Curtis.

    :param r_0: `iterable` (km) initial position 3-vector
    :param v_0: `iterable` (km/s) initial velocity 3-vector
    :param dt: `float` (s) time after initial state to solve for r, v as 3-vectors
    :param body: `CelestialBody` (--) the celestial body to use for orbital parameters
    :param method: `str` (--) which numerical method to use to solve Kepler's Equation
    :param tol: `float` (--) decimal tolerance for numerical method (default 1e-7 is IEEE 745 single precision)
    :param max_iters: `int` (--) maximum number of iterations in numerical method before breaking
    :return: (km) final position 3-vector, (km/s) final velocity 3-vector
    '''
    VALID_METHODS = ('laguerre', 'newton')

    mu = body.mu  # (km**3/s**2) gravitational parameter of the specified primary body

    r0 = np.linalg.norm(r_0)  # (km) initial position magnitude
    v0 = np.linalg.norm(v_0)  # (km/s) initial velocity magnitude
    vr0 = np.dot(v_0, r_0)/r0  # (km/s) initial radial velocity magnitude
    alpha = 2/r0 - v0**2/mu  # (1/km) inverse of semi-major axis

    chi0 = np.sqrt(mu)*np.abs(alpha)*dt

    if method not in VALID_METHODS:
        print(f'Method \'{method}\' is not valid, must be one of {VALID_METHODS}.\nDefaulting to laguerre method.')
        chi, _, _ = laguerre(chi0, kepler_chi, dkepler_dchi, d2kepler_dchi2, alpha, r0, vr0, mu, dt)
    elif method == 'newton':
        chi, _, _ = newton(chi0, kepler_chi, dkepler_dchi, alpha, r0, vr0, mu, dt)
    else:  # method == 'laguerre'
        chi, _, _ = laguerre(chi0, kepler_chi, dkepler_dchi, d2kepler_dchi2, alpha, r0, vr0, mu, dt)

    f = calc_f(chi, r0, alpha)
    g = calc_g(dt, mu, chi, alpha)
    r_1 = f*r_0 + g*v_0
    r1 = np.linalg.norm(r_1)

    fd = calc_fd(mu, r1, r0, alpha, chi)
    gd = calc_gd(chi, r1, alpha)
    v_1 = fd*r_0 + gd*v_0

    return r_1, v_1


def solve_kepler_E(e, Me, tol=1e-7, max_iters=100):
    ''' Solve Kepler's Equation in the form containing Eccentric Anomaly (E),
    eccentricity (e), and Mean Anomaly of Ellipse (Me). Uses Algorithm 3.1 from Orbital
    Mechanics for Engineering Students, 4 ed, Curtis. '''

    # TODO: have this function make use of one of the numerical methods in numerical.py

    def f(E, e, Me):
        return E - e*np.sin(E) - Me

    def fp(E, e):
        return 1 - e*np.cos(E)

    E = Me + e/2 if Me < np.pi else Me - e/2
    ratio = f(E, e, Me)/fp(E, e)
    iters = 0
    while abs(ratio) > tol and iters < max_iters:
        E -= ratio
        ratio = f(E, e, Me)/fp(E, e)
        iters += 1

    E -= ratio
    converged = np.abs(ratio) <= tol

    return E, iters, converged


def test():
    ''' Test the functionality of solve_kepler_chi 
    and solve_kepler_laguerre using Problem 3.20 from
    Orbital Mechanics for Engineering Students, 4 ed, Curtis.
    '''
    # given starting information
    Earth = BODIES['Earth']  # `CelestialBody` (--) Earth and all the Earth things
    r_0 = np.array([20000, -105000, -19000])  # (km) initial position vector
    v_0 = np.array([0.9, -3.4, -1.5])  # (km/s) initial velocity vector
    dt = 2*60*60  # (s) time of interest after initial time

    # given correct answer from textbook
    correct_r_1 = np.array([26338, -128750, -29656])  # (km) final position vector
    correct_v_1 = np.array([0.86280, -3.2116, -1.4613])  # (km/s) final velocity vector

    # solve using above methods
    r_n, v_n = solve_kepler_chi(r_0, v_0, dt, Earth, method='newton')
    r_l, v_l = solve_kepler_chi(r_0, v_0, dt, Earth, method='laguerre')

    # check correctness
    # tolerance based on significant figures of given answers
    newton_valid = np.allclose(r_n, correct_r_1, atol=1) and np.allclose(v_n, correct_v_1, atol=1e-4)
    laguerre_valid = np.allclose(r_l, correct_r_1, atol=1) and np.allclose(v_l, correct_v_1, atol=1e-4)

    return all([newton_valid, laguerre_valid])


if __name__ == '__main__':
    print(test())
