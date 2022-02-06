import numpy as np


def newton(x0, f, fp, *args, tol=1e-7, max_iters=100):
    ''' Solve f(x) for x numerically using Newton's Method. '''
    ratio = f(x0, *args)/fp(x0, *args)
    iters = 0
    while np.abs(ratio) > tol and iters < max_iters:
        x0 -= ratio
        ratio = f(x0, *args) / fp(x0, *args)
        iters += 1
    x0 -= ratio
    converged = np.abs(ratio) <= tol
    return x0, iters, converged


def laguerre(x0, f, fp, fpp, *args, n=5, tol=1e-7, max_iters=100):
    ''' Solve f(x) for x numerically using Laguerre's Method. '''

    def calc_ratio(x0, *args):
        f_ = f(x0, *args)
        fp_ = fp(x0, *args)
        fpp_ = fpp(x0, *args)
        num = n*f_
        rad = np.sqrt(np.abs((n - 1)**2*fp_**2 - n*(n - 1)*f_*fpp_))
        den = np.amax([fp_ + rad, fp_ - rad])
        return num/den

    ratio = calc_ratio(x0, *args)
    iters = 0
    while np.abs(ratio) > tol and iters < max_iters:
        x0 -= ratio
        ratio = calc_ratio(x0, *args)
        iters += 1
    x0 -= ratio
    converged = np.abs(ratio) <= tol
    return x0, iters, converged
