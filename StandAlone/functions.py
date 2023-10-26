#!/usr/bin/env python
#
# Description: 
# 
import numpy as np
from math import acos, pi, cos, sin

def a_from_vectors(r_vector: np.array,
                   v_vector: np.array,
                   mu_val: float) -> float:
    """
    Calculates semi-major axis using rearranged vis-viva and vector magnitudes

    """
    r_mag = np.linalg.norm(r_vector)
    v_mag = np.linalg.norm(v_vector)

    a = 1 / ((2 / r_mag) - (v_mag ** 2 / mu_val))

    return a


def eccentricity_from_vectors(r_vector: np.array,
                              v_vector: np.array,
                              mu_val: float) -> float:
    """
    Calculates the eccentriciy of the orbit given r and v using formula from Lecture 5, slide 11

    """

    r_v_dot = np.dot(r_vector, v_vector)
    r_mag = np.linalg.norm(r_vector)
    v_mag = np.linalg.norm(v_vector)
    
    ecc_vector = (((v_mag ** 2) / mu_val) - (1 / r_mag)) * r_vector \
                 - (1 / mu_val) * r_v_dot * v_vector

    return ecc_vector


def inclination(r_vector: np.array,
                v_vector: np.array,) -> tuple:
    """
    Return inclination angle of orbit
    """
    h_vector = np.cross(r_vector, v_vector)
    h_mag = np.linalg.norm(h_vector)
    h_norm = h_vector / h_mag
    
    k = np.array([0, 0, 1])

    cos_i = np.dot(h_norm, k)

    # math.acos returns values in radians
    i_radians = acos(cos_i)

    return i_radians * 180/pi


def ra_o_an(r_vector, v_vector) -> tuple:
    """
    Return Right Angle of Ascending Node of orbit

    """

    h_vector = np.cross(r_vector, v_vector)
    h_mag = np.linalg.norm(h_vector)
    h_norm = h_vector / h_mag

    n_vector = np.cross([0,0,1], h_norm)
    n_mag = np.linalg.norm(n_vector)
    n_norm = n_vector / n_mag

    cos_omega = np.dot(n_norm, [1, 0, 0])
    
    # math.acos returns values in radians
    ra_o_an = acos(cos_omega)

    if np.dot(n_vector, [0, 1, 0]) < 0:
        ra_o_an = 2*pi - ra_o_an

    return ra_o_an  * 180/pi


def arg_of_periapse(r_vector, v_vector, mu) -> tuple:
    """
    Return the argument of periapce of orbit

    """
    e_vector = eccentricity_from_vectors(r_vector, v_vector, mu)

    h_vector = np.cross(r_vector, v_vector)
    h_mag = np.linalg.norm(h_vector)
    h_norm = h_vector / h_mag

    n_vector = np.cross([0,0,1], h_norm)
    n_mag = np.linalg.norm(n_vector)

    cos_w = np.dot(n_vector, e_vector) / (n_mag * np.linalg.norm(e_vector))

    # math.acos returns values in radians
    arg_of_p = acos(cos_w)

    if np.dot(e_vector, [0, 0, 1]) < 0:
        arg_of_p = 2*pi - arg_of_p
    
    return arg_of_p * 180/pi


def true_anom_f(r_vector, v_vector, mu) -> tuple:
    """
    Returns true anomaly of position in orbit

    """

    e_vector = eccentricity_from_vectors(r_vector, v_vector, mu)
    e_mag = np.linalg.norm(e_vector)
    r_mag = np.linalg.norm(r_vector)

    cos_f = np.dot(e_vector, r_vector) / (e_mag * r_mag)

    # math.acos returns values in radians
    true_anom_f = acos(cos_f)

    if np.dot(r_vector, v_vector) < 0:
        true_anom_f = 2*pi - true_anom_f
    
    return true_anom_f * 180/pi

def orbital_elems_from_vectors(r, v, mu):

    elems_dict = {}

    elems_dict['a'] = a_from_vectors(r, v, mu)
    elems_dict['e'] = np.linalg.norm(eccentricity_from_vectors(r, v, mu))
    elems_dict['i'] = inclination(r, v)
    elems_dict['cap_ohm'] = ra_o_an(r, v)
    elems_dict['low_ohm'] = arg_of_periapse(r, v, mu)
    elems_dict['f'] = true_anom_f(r, v, mu)

    return elems_dict

def radial_vect_from_orbital_elems(a, e, inc, cap_ohm, low_ohm, f):

    i, j, k = np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])

    # work out radial distance using standard eclipse formula
    r = (a * (1 - e**2) ) / \
        ( 1 + e*cos(f * pi/180))

    # convert angles into radians
    inc = inc * (pi/180)
    low_ohm = low_ohm * (pi/180)
    cap_ohm = cap_ohm * (pi/180)
    theta = low_ohm + f * (pi/180)

    # euler rotation vector transformation terms
    i_term = cos(theta)*cos(cap_ohm) - cos(inc)*sin(cap_ohm)*sin(theta)
    j_term = cos(theta)*sin(cap_ohm) + cos(inc)*cos(cap_ohm)*sin(theta)
    k_term = sin(inc)*sin(theta)

    r = r * (i_term*i + j_term*j + k_term*k)

    return r

def main():

    mu = float(1.32712440018E+11)
    r = np.array([1.269335017088585E+08, 7.759541255851591E+07, -5.128282554760575E+03])
    v = np.array([-1.602222110677506E+01, 2.530395196294028E+01, -1.292576889628805E-03])

    elems = orbital_elems_from_vectors(r, v, mu)

    print(elems)

    r_vector_check = radial_vect_from_orbital_elems(elems['a'], elems['e'], elems['i'], elems['cap_ohm'], elems['low_ohm'], elems['f'])
    print(r_vector_check)

    pass


if __name__ == '__main__':
    main()