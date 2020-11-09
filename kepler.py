import numpy as np
from lagrange import calc_f, calc_fd, calc_g, calc_gd
from stumpff import C, S, dCdz, dSdz


def f(chi_i, alpha, r_0, v_r_0, mu, dt):
    z_i = alpha*chi_i**2
    return (r_0*v_r_0/np.sqrt(mu))*chi_i**2*C(z_i) + \
           (1 - alpha*r_0)*chi_i**3*S(z_i) + \
           r_0*chi_i - np.sqrt(mu)*dt


def fp(chi_i, alpha, r_0, v_r_0, mu):
    z_i = alpha*chi_i**2
    return (r_0*v_r_0/np.sqrt(mu))*chi_i*(1 - alpha*chi_i**2*S(z_i)) + \
           (1 - alpha*r_0)*chi_i**2*C(z_i) + r_0


def fpp(chi, alpha, r_0, v_r_0, mu):
    z = alpha*chi**2
    S_ = S(z)
    return (r_0*v_r_0/np.sqrt(mu))*(1 - 3*z*S_ + z*(C(z) - 3*S_)) +  \
           chi*(1 - z*S_)*(1 - alpha*r_0)


def solve_kepler_newton(chi_0, r_0, v_r_0, mu, dt, alpha, tol=1e-7, max_iters=100):
    ''' Iteratively solve Kepler's Equation using Newton's method, assuming
    a two-body problem.

    $$\Chi_{i+1} = \Chi_{i} - f(\alpha \Chi_i^2)/f'(\alpha \Chi_i^2)$$

    :param chi_0: (rad) initial estimate for Chi
    :param r_0: (km) initial position magnitude
    :param v_r_0: (km/s) initial radial velocity magnitude
    :param mu: (km**3/s**2) gravitational parameter of primary body
    :param dt: (s) delta time from initial conditions at which to find Chi
    :param alpha: (1/km) inverse of semi-major axis of orbit
    :param tol: (--) decimal precision of returned Chi (1e-7 is IEEE 754 single precision)
    :param max_iters: (--) maximum number of iterations to attempt without convergence

    :return: (rad) Chi at dt sections after initial conditions
    '''
    chi_i = chi_0
    iters = 0
    ratio = f(chi_i, alpha, r_0, v_r_0, mu, dt) / \
            fp(chi_i, alpha, r_0, v_r_0, mu)
    
    while np.abs(ratio) > tol and iters < max_iters:
        chi_i -= ratio
        ratio = f(chi_i, alpha, r_0, v_r_0, mu, dt) / \
                fp(chi_i, alpha, r_0, v_r_0, mu)
        iters += 1

    chi_i -= ratio
    converged = np.abs(ratio) <= tol

    return chi_i, iters, converged


def solve_kepler_laguerre(chi_0, r_0, v_r_0, mu, dt, alpha, n=5, tol=1e-7, max_iters=100):
    ''' Iteratively solve Kepler's Equation using Laguerres's method, assuming
    a two-body problem.

    $$\Chi_{i+1} = \Chi_{i} - f(\alpha \Chi_i^2)/f'(\alpha \Chi_i^2)$$

    :ref: "An Improved Algorithm Due to Laguerre for the Solution of Kepler's Equation", Conway (1985)

    :param chi_0: (rad) initial estimate for Chi
    :param n: (--) degree of "polynomial" (see reference, default 5)
    :param r_0: (km) initial position magnitude
    :param v_r_0: (km/s) initial radial velocity magnitude
    :param mu: (km**3/s**2) gravitational parameter of primary body
    :param dt: (s) delta time from initial conditions at which to find Chi
    :param alpha: (1/km) inverse of semi-major axis of orbit
    :param tol: (--) decimal precision of returned Chi (1e-7 is IEEE 754 single precision)
    :param max_iters: (--) maximum number of iterations to attempt without convergence

    :return: (rad) Chi at dt sections after initial conditions
    '''
    chi_i = chi_0
    iters = 0

    f_ = f(chi_i, alpha, r_0, v_r_0, mu, dt)
    fp_ = fp(chi_i, alpha, r_0, v_r_0, mu)
    fpp_ = fpp(chi_i, alpha, r_0, v_r_0, mu)
    num = n*f_
    rad = np.sqrt(np.abs((n - 1)**2*fp_**2 - n*(n - 1)*f_*fpp_))
    den = np.amax([fp_ + rad, fp_ - rad])
    ratio = num/den

    while np.abs(ratio) > tol and iters < max_iters:
        chi_i -= ratio
        f_ = f(chi_i, alpha, r_0, v_r_0, mu, dt)
        fp_ = fp(chi_i, alpha, r_0, v_r_0, mu)
        fpp_ = fpp(chi_i, alpha, r_0, v_r_0, mu)
        num = n*f_
        rad = np.sqrt(np.abs((n - 1)**2*fp_**2 - n*(n - 1)*f_*fpp_))
        den = np.amax([fp_ + rad, fp_ - rad])
        ratio = num/den
        iters += 1

    chi_i -= ratio
    converged = np.abs(ratio) <= tol

    return chi_i, iters, converged


def solve_kepler(r_0, v_0, dt, mu=398600, method='laguerre', tol=1e-7, max_iters=100):
    ''' Wrapper for either solve_kepler_newton or solve_kepler_laguerre,
    whicever is specified. Applies Algorithm 3.4 from Orbital Mechanics for
    Engineering Students, 4 ed, Curtis.

    :param r_0: `iterable` (km) initial position 3-vector
    :param v_0: `iterable` (km/s) initial velocity 3-vector
    :param dt: `float` (s) time after initial state to sovle for r, v
    :param mu: `float` (km**3/s**2) gravitational parameter, default Earth
    :param method: `str` (--) which numerical method to use to solve Kepler's Equation
    :param tol: `float` (--) decimal tolerance for numerical method (default 1e-7 is IEEE 745 single precision)
    :param max_iters: `int` (--) maximum number of iterations in numerical method before breaking
    :return: (km) final position 3-vector, (km/s) final velocity 3-vector
    '''
    VALID_METHODS = ('laguerre', 'newton')

    r0 = np.linalg.norm(r_0)  # (km) initial position magnitude
    v0 = np.linalg.norm(v_0)  # (km/s) initial velocity magnitude
    vr0 = np.dot(v_0, r_0)/r0  # (km/s) initial radial velocity magnitude
    alpha = 2/r0 - v0**2/mu  # (1/km) inverse of semi-major axis

    chi0 = np.sqrt(mu)*np.abs(alpha)*dt

    if method not in VALID_METHODS:
        print(f'Method {method} is not valid, must be one of {VALID_METHODS}.\nDefaulting to laguerre method.')
        chi, _, _ = solve_kepler_laguerre(chi0, r0, vr0, mu, dt, alpha)
    elif method == 'newton':
        chi, _, _ = solve_kepler_newton(chi0, r0, vr0, mu, dt, alpha)
    else:  # method == 'laguerre'
        chi, _, _ = solve_kepler_laguerre(chi0, r0, vr0, mu, dt, alpha)

    f = calc_f(chi, r0, alpha)
    g = calc_g(dt, mu, chi, alpha)
    r_1 = f*r_0 + g*v_0
    r1 = np.linalg.norm(r_1)

    fd = calc_fd(mu, r1, r0, alpha, chi)
    gd = calc_gd(chi, r1, alpha)
    v_1 = fd*r_0 + gd*v_0

    return r_1, v_1


def test():
    ''' Test the functionality of solve_kepler_newton 
    and solve_kepler_laguerre using Problem 3.20 from
    Orbital Mechanics for Engineering Students, 4 ed, Curtis.
    '''
    # given starting information
    mu = 398600  # (km**3/s**2) gravitational parameter, Earth
    r_0 = np.array([20000, -105000, -19000])  # (km) initial position vector
    v_0 = np.array([0.9, -3.4, -1.5])  # (km/s) initial velocity vector
    dt = 2*60*60  # (s) time of interest after initial time

    # given correct answer from textbook
    correct_r_1 = np.array([26338, -128750, -29656])  # (km) final position vector
    correct_v_1 = np.array([0.86280, -3.2116, -1.4613])  # (km/s) final velocity vector

    # solve using above methods
    r_n, v_n = solve_kepler(r_0, v_0, dt, method='newton')
    r_l, v_l = solve_kepler(r_0, v_0, dt, method='laguerre')

    # check correctness
    # tolerance based on significant figures of given answers
    newton_valid = np.allclose(r_n, correct_r_1, atol=1) and np.allclose(v_n, correct_v_1, atol=1e-4)
    laguerre_valid = np.allclose(r_l, correct_r_1, atol=1) and np.allclose(v_l, correct_v_1, atol=1e-4)

    return all([newton_valid, laguerre_valid])


if __name__ == '__main__':
    print(test())
