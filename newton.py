from thermostate import Q_
import numpy as np
from constants import G

def F_g(m1, r1, m2, r2):
    '''
    Calculates the gravitational force exerted on body 2 by body 1.

    :param m1: `float` (kg) mass of body 1
    :param m2: `float` (kg) mass of body 2
    :param r1: `array` ([m, m, m]) vector from origin to body 1
    :param r2: `array` ([m, m, m]) vector from origin to body 2

    :return: `array` ([N, N, N]) gravitational force exterted on body 2
    '''
    return -(G*m1*m2/(np.linalg.norm(r2-r1))**3)*(r2-r1)


def rdd(i, m_array, r_array):
    '''
    Computes the vector acceleration of body i given a list of body
    masses and positions. Assumes more than 1 body is given.

    :param i: `int` (--) the index of the body to find acceleration of
    :param m_array: `array` ([kg, kg, ...]) masses of all bodies
    :param r_array: `array` ([m, m, ...]) positions of all bodies

    :return: `array` ([m/s**2, m/s**2, m/s**2]) vector acceleration of body i
    '''
    mi = m_array[i]
    ri = r_array[i]
    m_array = np.array([m_array[j] for j in range(len(m_array)) if j != i])
    r_array = np.array([r_array[j] for j in range(len(r_array)) if j != i])

    r = np.zeros(3)
    for m, r in zip(m_array, r_array):
        r -= (G*m/(np.linalg.norm(ri-r))**3)*(ri-r)
    
    return r
