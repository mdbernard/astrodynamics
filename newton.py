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

    Assumes mass of body i is not changing, and all forces acting on
    body i are gravitational.

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


def rdd_rel(i, j, m_array, r_array):
    '''
    Computes the relative acceleration vector of body j with respect
    to body i. Assumes i and j are valid indices in m_array and r_array.

    Assumes all bodies are point masses of constant magnitude. Assumes all
    forces are gravitational between bodies.

    :param i: `int` (--) index of body i in m_array and r_array
    :param j: `int` (--) index of body j in m_array and r_array
    :param m_array: `1x3 array` (kg) masses of all bodies
    :param r_array: `1x3 array` (m) positions of all bodies relative to
                    origin of inertial frame
    
    :return: `1x3 array` acceleration vector of body j relative to body i
    '''
    mi, mj = m_array[i], m_array[j]
    ri, rj = r_array[i], r_array[j]

    m_array = np.array([m_array[m] for m in range(len(m_array)) if m != i and m != j])    
    r_array = np.array([r_array[r] for r in range(len(r_array)) if r != i and r != j])

    rdd_ij = -G*(mi+mj)/(np.linalg.norm(rj-ri)**3)*(rj-ri)
    for m, r in zip(m_array, r_array):
        rdd_ij -= G*m*((rj-r)/np.linalg.norm(rj-r)**3 - (ri-r)/np.linalg.norm(ri-r)**3)
    
    return rdd_ij
