import numpy as np
from stumpff import C, S


def calc_f(chi, r_0, alpha):
    return 1 - chi**2/r_0*C(alpha*chi**2)


def calc_g(dt, mu, chi, alpha):
    return dt - 1/np.sqrt(mu)*chi**3*S(alpha*chi**2)


def calc_fd(mu, r, r_0, alpha, chi):
    return np.sqrt(mu)/r/r_0*(alpha*chi**3*S(alpha*chi**2) - chi)


def calc_gd(chi, r, alpha):
    return 1 - chi**2/r*C(alpha*chi**2)

