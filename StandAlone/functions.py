#!/usr/bin/env python
#
# Description: 
# 
import numpy as np
from math import acos, pi, cos, sin, sqrt
from math import acos, pi, cos, sin, atan2

def v_vis_viva(r_t: float,
               a: float,
               mu_val: float) -> float:
    """
    Returns velocity magnitude from standard vis-viva eqaution for elliptical orbits

    """
    v = sqrt(mu_val * ( (2/r_t) - (1/a) ))

    return v


def a_from_vectors(r_vector: np.array,
                   v_vector: np.array,
                   mu_val:   float) -> float:
    """
    Calculates semi-major axis using rearranged vis-viva and vector magnitudes
    :param r_vector: position vector as dictionary of {vector, magnitude}
    :param v_vector: velocity vector as dictionary of {vector, magnitude}
    :param mu_val:   mu value for the central body
    :return: tuple of (vector, mag)
    """

    r_mag = np.linalg.norm(r_vector)
    v_mag = np.linalg.norm(v_vector)

    a = 1 / ((2 / r_mag) - (v_mag ** 2 / mu_val))

    return a


def eccentricity_from_vectors(r_vector: np.array,
                              v_vector: np.array,
                              mu_val:   float) -> tuple:
    """
    Calculates the eccentriciy of the orbit given r and v using formula from Lecture 5, slide 11

    :param r_vector: position vector as dictionary of {vector, magnitude}
    :param v_vector: velocity vector as dictionary of {vector, magnitude}
    :param mu_val:   mu value for the central body
    :return: tuple of (vector, mag)
    """
        
    r_mag = np.linalg.norm(r_vector)
    v_mag = np.linalg.norm(v_vector)
    r_v_dot = np.dot(r_vector, v_vector)

    ecc_vector = (((v_mag ** 2) / mu_val) - (1 / r_mag)) * r_vector \
                 - (1 / mu_val) * r_v_dot * v_vector

    ecc_mag = np.linalg.norm(ecc_vector)

    return ecc_vector, ecc_mag

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

    # acos returns values in radians
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
    
    # acos returns values in radians
    ra_o_an = acos(cos_omega)

    if np.dot(n_vector, [0, 1, 0]) < 0:
        ra_o_an = 2*pi - ra_o_an

    return ra_o_an  * 180/pi


def arg_of_periapse(r_vector, v_vector, mu) -> tuple:
    """
    Return the argument of periapce of orbit

    """
    e_vector = eccentricity_from_vectors(r_vector, v_vector, mu)[0]

    h_vector = np.cross(r_vector, v_vector)
    h_mag = np.linalg.norm(h_vector)
    h_norm = h_vector / h_mag

    n_vector = np.cross([0,0,1], h_norm)
    n_mag = np.linalg.norm(n_vector)

    cos_w = np.dot(n_vector, e_vector) / (n_mag * np.linalg.norm(e_vector))

    # acos returns values in radians
    arg_of_p = acos(cos_w)

    if np.dot(e_vector, [0, 0, 1]) < 0:
        arg_of_p = 2*pi - arg_of_p
    
    return arg_of_p * 180/pi


def true_anom_f(r_vector, v_vector, mu) -> tuple:
    """
    Returns true anomaly of position in orbit

    """

    e_vector = eccentricity_from_vectors(r_vector, v_vector, mu)[0]
    e_mag = np.linalg.norm(e_vector)
    r_mag = np.linalg.norm(r_vector)

    cos_f = np.dot(e_vector, r_vector) / (e_mag * r_mag)

    # acos returns values in radians
    true_anom_f = acos(cos_f)

    if np.dot(r_vector, v_vector) < 0:
        true_anom_f = 2*pi - true_anom_f
    
    return true_anom_f * 180/pi

def orbital_elems_from_vectors(r, v, mu):

    elems_dict = {}

    elems_dict['a'] = a_from_vectors(r, v, mu)
    elems_dict['e'] = eccentricity_from_vectors(r, v, mu)[1]
    elems_dict['i'] = inclination(r, v)
    elems_dict['cap_ohm'] = ra_o_an(r, v)
    elems_dict['low_ohm'] = arg_of_periapse(r, v, mu)
    elems_dict['f'] = true_anom_f(r, v, mu)

    return elems_dict

def initial_eccentric_anom(r_vector: np.array,
                           v_vector: np.array,
                           mu_val: float,
                           semi_maj_ax: float) -> tuple:
    """
    Return an initial eccentric anomaly value from vectors, mu and a

    :param r_vector:        position vector as dictionary of {vector, magnitude}
    :param v_vector:        velocity vector as dictionary of {vector, magnitude}
    :param mu_val:          mu value for the central body
    :param semi_maj_ax:     a value for orbit
    :return:                E0 value as float
    """
    tan_E0_num = np.dot(r_vector, v_vector)
    r_mag = np.linalg.norm(r_vector)
    tan_E0_denom = ((semi_maj_ax * mu_val) ** 0.5) * (1 - (r_mag / semi_maj_ax))

    return atan2(tan_E0_num, tan_E0_denom)


def kepler_E_solution_iteration(eccentricity, n, delta_t, E_0):
    """
    Kepler iteration function to give E, with an E0 not at perigee

    :param eccentricity: eccentricity
    :param n: mean angular rate
    :param delta_t: timestep for new E since E0
    :param E_0: value of initial eccentric anomaly
    :return: eccentric anomaly in radians
    """

    # Because might not be at perigee, need to use modified Kepler's with E0 as extra constant:
    adj_constant = E_0 - eccentricity * sin(E_0)

    E_iter = n * delta_t
    g = 1
    i = 0

    # Actual iteration loop
    while abs(g) > 1e-13:
        g = (E_iter - eccentricity * sin(E_iter)) - (n * delta_t) - adj_constant
        dgdE = 1 - eccentricity * cos(E_iter)

        E_1 = E_iter - g / dgdE

        # updates
        E_iter = E_1
        i += 1

    return E_iter


def f_value_for_r_t1(semi_maj_ax: float,
                     r_vector: dict,
                     initial_e0: float,
                     final_e: float) -> float:
    
    r_mag = np.linalg.norm(r_vector)
    F_value = 1 - (semi_maj_ax / r_mag) * (1 - cos(final_e - initial_e0))

    return F_value


def g_value_for_r_t1(semi_maj_ax: float,
                     mu: float,
                     delta_t: float,
                     initial_e0: float,
                     final_e: float) -> float:

    G_value = (delta_t) - ((semi_maj_ax ** 3 / mu) ** 0.5) * \
              ((final_e - initial_e0) - sin((final_e - initial_e0)))

    return G_value


def f_dot_value_for_v_t1(semi_maj_ax: float,
                         mu: float,
                         r0_vector: np.array,
                         r1_vector: np.array,
                         initial_e0: float,
                         final_e: float) -> float:

    r0_mag = np.linalg.norm(r0_vector)
    r1_mag = np.linalg.norm(r1_vector)
    
    F_dot_value = - ((mu * semi_maj_ax) ** 0.5 / (r0_mag * r1_mag)) \
                  * sin(final_e - initial_e0)

    return F_dot_value


def g_dot_value_for_v_t1(semi_maj_ax: float,
                         mu: float,
                         r1_vector: np.array,
                         delta_t: float,
                         initial_e0: float,
                         final_e: float) -> float:
    
    r1_mag = np.linalg.norm(r1_vector)
    G_dot_value = 1 - (semi_maj_ax / r1_mag) * (1 - cos(final_e - initial_e0))

    return G_dot_value


def r_as_function_of_t(mu, e, a, E0, r0_vector, v0_vector, dt):
    """
    Return radius vector for the given input values - Calculates new F and G

    :param mu: gravitational constant for central body
    :param e: eccentricity
    :param a: semi-major axis
    :param E0: value of initial eccentric anomaly
    :param r0_vector: initial state of r
    :param v0_vector: initial state of v
    :param dt: timestep for new E since E0
    :return: new r vector from calculation involing F and G
    """
    n = (mu / (a ** 3)) ** 0.5

    E_t = kepler_E_solution_iteration(e, n, dt, E0)

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    return r_t1

def r_and_v_as_function_of_t(mu, r0_vector, v0_vector, dt):
    """
    Return radius vector for the given input values - Calculates new F and G

    :param mu: gravitational constant for central body
    :param e: eccentricity
    :param a: semi-major axis
    :param E0: value of initial eccentric anomaly
    :param r0_vector: initial state of r
    :param v0_vector: initial state of v
    :param dt: timestep for new E since E0
    :return: new r vector from calculation involing F and G
    """

    # determine semi-major axis
    a = a_from_vectors(r0_vector, v0_vector, mu)

    # determine eccentricity
    e = eccentricity_from_vectors(r0_vector, v0_vector, mu)[0]
    e = np.linalg.norm(e)

    # h value (for check later on)
    # h_mag = h_value_from_elements(a, e[1], mu)

    # initial eccentric anom
    E0 = initial_eccentric_anom(r0_vector, v0_vector, mu, a)
    
    n = (mu / (a ** 3)) ** 0.5

    E_t = kepler_E_solution_iteration(e, n, dt, E0)

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    F_dot = f_dot_value_for_v_t1(a, mu, r0_vector, r_t1, E0, E_t)
    G_dot = g_dot_value_for_v_t1(a, mu, r_t1, dt, E0, E_t)

    v_t1 = F_dot * r0_vector + G_dot * v0_vector

    return r_t1, v_t1


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

def kepler_E_solution_iteration(eccentricity, n, delta_t, E_0):
    """
    Kepler iteration function to give E, with an E0 not at perigee

    :param eccentricity: eccentricity
    :param n: mean angular rate
    :param delta_t: timestep for new E since E0
    :param E_0: value of initial eccentric anomaly
    :return: eccentric anomaly in radians
    """

    # Because might not be at perigee, need to use modified Kepler's with E0 as extra constant:
    adj_constant = E_0 - eccentricity * sin(E_0)

    E_iter = n * delta_t
    g = 1
    i = 0

    # Actual iteration loop
    while abs(g) > 1e-13:
        g = (E_iter - eccentricity * sin(E_iter)) - (n * delta_t) - adj_constant
        dgdE = 1 - eccentricity * cos(E_iter)

        E_1 = E_iter - g / dgdE

        # updates
        E_iter = E_1
        i += 1

    return E_iter


def f_value_for_r_t1(semi_maj_ax: float,
                     r_vector: dict,
                     initial_e0: float,
                     final_e: float) -> float:
    
    r_mag = np.linalg.norm(r_vector)
    F_value = 1 - (semi_maj_ax / r_mag) * (1 - cos(final_e - initial_e0))

    return F_value


def g_value_for_r_t1(semi_maj_ax: float,
                     mu: float,
                     delta_t: float,
                     initial_e0: float,
                     final_e: float) -> float:

    G_value = (delta_t) - ((semi_maj_ax ** 3 / mu) ** 0.5) * \
              ((final_e - initial_e0) - sin((final_e - initial_e0)))

    return G_value


def f_dot_value_for_v_t1(semi_maj_ax: float,
                         mu: float,
                         r0_vector: np.array,
                         r1_vector: np.array,
                         initial_e0: float,
                         final_e: float) -> float:

    r0_mag = np.linalg.norm(r0_vector)
    r1_mag = np.linalg.norm(r1_vector)
    
    F_dot_value = - ((mu * semi_maj_ax) ** 0.5 / (r0_mag * r1_mag)) \
                  * sin(final_e - initial_e0)

    return F_dot_value


def g_dot_value_for_v_t1(semi_maj_ax: float,
                         mu: float,
                         r1_vector: np.array,
                         delta_t: float,
                         initial_e0: float,
                         final_e: float) -> float:
    
    r1_mag = np.linalg.norm(r1_vector)
    G_dot_value = 1 - (semi_maj_ax / r1_mag) * (1 - cos(final_e - initial_e0))

    return G_dot_value


def r_as_function_of_t(mu, e, a, E0, r0_vector, v0_vector, dt):
    """
    Return radius vector for the given input values - Calculates new F and G

    :param mu: gravitational constant for central body
    :param e: eccentricity
    :param a: semi-major axis
    :param E0: value of initial eccentric anomaly
    :param r0_vector: initial state of r
    :param v0_vector: initial state of v
    :param dt: timestep for new E since E0
    :return: new r vector from calculation involing F and G
    """
    n = (mu / (a ** 3)) ** 0.5

    E_t = kepler_E_solution_iteration(e, n, dt, E0)

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    return r_t1

def initial_eccentric_anom(r_vector: np.array,
                           v_vector: np.array,
                           mu_val: float,
                           semi_maj_ax: float) -> tuple:
    """
    Return an initial eccentric anomaly value from vectors, mu and a

    :param r_vector:        position vector as dictionary of {vector, magnitude}
    :param v_vector:        velocity vector as dictionary of {vector, magnitude}
    :param mu_val:          mu value for the central body
    :param semi_maj_ax:     a value for orbit
    :return:                E0 value as float
    """
    tan_E0_num = np.dot(r_vector, v_vector)
    r_mag = np.linalg.norm(r_vector)
    tan_E0_denom = ((semi_maj_ax * mu_val) ** 0.5) * (1 - (r_mag / semi_maj_ax))

    return atan2(tan_E0_num, tan_E0_denom)

def r_and_v_as_function_of_t(mu, r0_vector, v0_vector, dt):
    """
    Return radius vector for the given input values - Calculates new F and G

    :param mu: gravitational constant for central body
    :param e: eccentricity
    :param a: semi-major axis
    :param E0: value of initial eccentric anomaly
    :param r0_vector: initial state of r
    :param v0_vector: initial state of v
    :param dt: timestep for new E since E0
    :return: new r vector from calculation involing F and G
    """

    # determine semi-major axis
    a = a_from_vectors(r0_vector, v0_vector, mu)

    # determine eccentricity
    e = eccentricity_from_vectors(r0_vector, v0_vector, mu)[0]
    e = np.linalg.norm(e)

    # h value (for check later on)
    # h_mag = h_value_from_elements(a, e[1], mu)

    # initial eccentric anom
    E0 = initial_eccentric_anom(r0_vector, v0_vector, mu, a)
    
    n = (mu / (a ** 3)) ** 0.5

    E_t = kepler_E_solution_iteration(e, n, dt, E0)

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    F_dot = f_dot_value_for_v_t1(a, mu, r0_vector, r_t1, E0, E_t)
    G_dot = g_dot_value_for_v_t1(a, mu, r_t1, dt, E0, E_t)

    v_t1 = F_dot * r0_vector + G_dot * v0_vector

    return r_t1, v_t1

def main():

    mu = float(1.32712440018E+11)
    r = np.array([1.269335017088585E+08, 7.759541255851591E+07, -5.128282554760575E+03])
    v = np.array([-1.602222110677506E+01, 2.530395196294028E+01, -1.292576889628805E-03])

    elems = orbital_elems_from_vectors(r, v, mu)

    # print(elems)

    r_vector_check = radial_vect_from_orbital_elems(elems['a'], elems['e'], elems['i'], elems['cap_ohm'], elems['low_ohm'], elems['f'])
    # print(r_vector_check)

    print(v_vis_viva(1.515e+8, 2.457e+8, mu))

    pass


if __name__ == '__main__':
    main()