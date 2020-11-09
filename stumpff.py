import numpy as np


def C(z):
    if z > 0:
        ret = (1 - np.cos(np.sqrt(z)))/z
    elif z < 0:
        ret = (np.cosh(np.sqrt(-z)) - 1)/(-z)
    else:
        ret = 0.5
    return ret


def S(z):
    if z > 0:
        ret = (np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3
    elif z < 0:
        ret = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3
    else:
        ret = 1/6
    return ret


def dCdz(z):
    return (1/2/z)*(1 - z*S(z) - 2*C(z))


def dSdz(z):
    return (1/2/z)*(C(z) - 3*S(z))