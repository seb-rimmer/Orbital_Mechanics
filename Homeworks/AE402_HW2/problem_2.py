"""
    problem_2.py created by Seb at 10:12 18/10/2021
    Path: HW2

    Write a computer program to convert orbital elements to position and velocity.
    Run the program for the following cases and upload your numbers and code into
    Gradescope.

    (a) Case 1: 𝑎 = 8000 km, 𝑒 = 0.125, 𝑖 = 10°, Ω = 45°, 𝜔 = 10°, 𝑀 = 170°
    (b) Case 2: 𝑟_p = 6611 km, 𝑒 = 0.01, 𝑖 = 90°, Ω = 0°, 𝜔 = 0°, 𝑀 = 355°
    
"""
from math import pi, sin, cos, tan, atan
import numpy as np
from problem_1 import sf_vector


def a_from_radii_and_e(orb_elems_dict: dict) -> float:
    """
    Returning a value for semi-maj axis from orbital elems dictionary, based on
    e and rp or ra values
    :param orb_elems_dict: dictionary of orbital elements
    :return:               a value as float
    """
    if orb_elems_dict['r_p']:
        a = orb_elems_dict['r_p'] / (1 - orb_elems_dict['e'])
    elif orb_elems_dict['r_a']:
        a = orb_elems_dict['r_a'] / (1 + orb_elems_dict['e'])
    else:
        print("Not enough information to calculate semi-maj axis from elements")

    return a


def kepler_E_solution_iteration(eccentricity, n, delta_t, E_0, mean_anom):
    """
    Kepler iteration function to give E, with an E0 not at perigee

    :param eccentricity: eccentricity
    :param n: mean angular rate
    :param delta_t: timestep for new E since E0
    :param E_0: value of initial eccentric anomaly
    :return: eccentric anomaly in radians
    """
    # M mean anomaly
    if mean_anom == None:
        mean_anom = n * delta_t

    # Because not at perigee, might need to use modified Kepler's with E0 as extra constant:
    # adj_constant = E_0 - sin(E_0)

    E_iter = mean_anom
    g = 1
    i = 0
    while abs(g) > 1e-13:

        g = (E_iter - eccentricity * sin(E_iter)) - mean_anom  # - adj_constant
        dgdE = 1 - eccentricity * cos(E_iter)

        E_1 = E_iter - g / dgdE

        # updates
        E_iter = E_1
        i += 1

    return E_iter, E_iter * 180 / pi


def true_anom(orb_elems: dict) -> float:
    """

    :param orb_elems:  dictionary of orbital elements
    :return:                Using tan formula to calcualte f from E and e
    """
    ecc = orb_elems['e']
    big_e = orb_elems['E'][0]

    tan_val = (((1 + ecc) / (1 - ecc)) ** 0.5) * tan(big_e / 2)

    f = atan(tan_val) * 2

    if f < 0:
        f += 2 * pi

    return f


def radius_magnitude(orb_elems_dict: dict) -> float:
    """
    Using Vis-Viva equation to calculate magnitude of r for the vector

    :param orb_elems_dict: dictionary of orbital elements
    :return:               radius magnitude as float
    """
    numer = orb_elems_dict['a'] * (1 - orb_elems_dict['e'] ** 2)
    denom = 1 + orb_elems_dict['e'] * cos(orb_elems_dict['f'])

    r_mag = numer / denom

    return r_mag


def h_value_from_elements(orb_elems_dict: dict, mu) -> float:
    """
    Return magnitude of h paramter from mu, e and a

    :param orb_elems_dict:  dictionary of orbital elements
    :param mu:              Mu value of central body
    :return:                h value, float
    """

    h = (mu * orb_elems_dict['a'] * (1 - orb_elems_dict['e'] ** 2)) ** 0.5

    return h


def radius_vector(case_dict: dict, theta_tup:float) -> tuple:
    """
    Compute full radius vector from orb elems dict
    :param orb_elems_dict: dictionary of orbital elements
    :return:               radius vector as tuple (vector, mag)
    """
    theta = theta_tup

    omega = case_dict['omega'] * pi/180
    i = case_dict['i'] * pi/180

    i_vector = cos(theta)*cos(omega) - \
               cos(i)*sin(omega)*sin(theta)

    j_vector = cos(theta)*sin(omega) + \
               cos(i)*cos(omega)*sin(theta)

    k_vector = sin(i)*sin(theta)

    unit_vectors = [i_vector, j_vector, k_vector]
    # print("Unit vectors: ", unit_vectors)

    r_vector = np.array(unit_vectors)

    r_vector = r_vector * case_dict['r_mag']

    return r_vector, np.linalg.norm(r_vector)


def velocity_vector(case_dict: dict, mu) -> tuple:

    theta = case_dict['theta'][0]
    mu_h_ratio = mu / case_dict['h_mag']

    w = case_dict['w'] * pi/180
    omega = case_dict['omega'] * pi/180
    i = case_dict['i'] * pi/180
    e = case_dict['e']

    i_vector = cos(omega) * (sin(theta) + e * sin(w)) + \
               sin(omega) * (cos(theta) + e * cos(w)) * cos(i)

    j_vector = sin(omega) * (sin(theta) + e * sin(w)) - \
               cos(omega) * (cos(theta) + e * cos(w)) * cos(i)

    k_vector = sin(i) * (cos(theta) + e * cos(w))

    scaled_vectors = [-mu_h_ratio * i_vector,
                      -mu_h_ratio * j_vector,
                      mu_h_ratio  * k_vector]

    v_vector = np.array(scaled_vectors)

    return v_vector, np.linalg.norm(v_vector)


def check_r_and_v_with_h(case_dict: dict) -> bool:
    """
    The magnitude of cross product of r and v should equal the mag of h

    :param case_dict:
    :return:
    """

    cross_product = np.cross(case_dict['r_vector'][0], case_dict['v_vector'][0])
    mag = np.linalg.norm(cross_product)
    # print(mag, case_dict['h_mag'])

    if (mag - case_dict['h_mag']) < 0.1:

        return True
    else:
        return False


def main():
    # variables
    mu = 398600
    earth_radius = 6378

    # defining cases as a dictionary
    cases_dict = {"Case 1":
                      {"a": 0,
                       "r_p": 6611,
                       "r_a": 0,
                       "e": 0.01,
                       "i": 90,
                       "omega": 0,
                       "w": 0,
                       "f": 0,
                       "M": 355
                       },
                  # "Case 2":
                  #     {"a": 0,
                  #      "r_p": 6611,
                  #      "r_a": 0,
                  #      "e": 0.01,
                  #      "i": 90,
                  #      "omega": 0,
                  #      "w": 0,
                  #      "M": 355
                  #      }
                  }

    # updating case 2 to inlcude a value of a in the orbital elements
    cases_dict['Case 1']['a'] = a_from_radii_and_e(cases_dict['Case 1'])

    intro_string = f"HW2 Problem 2 - r and v vectors from orbital elements.\n" \
                   f"------------------------------------------------------\n"

    with open('output/problem_2_output.txt', 'w') as output:
        output.write(intro_string)

    # Going through each case and the steps to find r and v vector for each
    for case in cases_dict:
        cases_dict[case]['E'] = kepler_E_solution_iteration(eccentricity=cases_dict[case]['e'],
                                                            n=None,
                                                            delta_t=None,
                                                            mean_anom=cases_dict[case]['M'],
                                                            E_0=0)

        if 'f' not in cases_dict[case]:
            cases_dict[case]['f'] = true_anom(cases_dict[case])

        cases_dict[case]['theta'] = (cases_dict[case]['w'] * pi/180 + cases_dict[case]['f'],
                                     cases_dict[case]['w'] + cases_dict[case]['f'] * 180/pi)
        # cases_dict[case]['theta'] = [pi/180 * (cases_dict[case]['w'] + cases_dict[case]['f'])]

        cases_dict[case]['r_mag'] = radius_magnitude(cases_dict[case])
        cases_dict[case]['r_vector'] = radius_vector(cases_dict[case], cases_dict[case]['theta'][0])
        cases_dict[case]['h_mag'] = h_value_from_elements(cases_dict[case], mu)

        cases_dict[case]['v_vector'] = velocity_vector(cases_dict[case], mu)

        output_string = f"Orbital elements for {case}:\n" \
                        f"a :                          {cases_dict[case]['a']:.2f} km\n" \
                        f"r_periapse :                 {cases_dict[case]['r_p']}\n" \
                        f"e :                          {cases_dict[case]['e']}\n" \
                        f"i (deg):                     {cases_dict[case]['i']}\n" \
                        f"RAoAN (deg):                 {cases_dict[case]['omega']}\n" \
                        f"Longit. of AN (deg):         {cases_dict[case]['w']}\n" \
                        f"Mean Anom (deg) :            {cases_dict[case]['M']}\n" \
                        f"Eccentric Anom (rad, deg):   {cases_dict[case]['E'][0]:.3f}, {cases_dict[case]['E'][1]:.3f}\n" \
                        f"Theta (rad):                 {cases_dict[case]['theta'][0]:.3f}\n" \
                        f"h vector :                   {cases_dict[case]['h_mag']:.3f}\n" \
                        f"True Anom:                   {cases_dict[case]['f']:.3f} radians\n" \
                        f"Radius Mag:                  {cases_dict[case]['r_mag']:.3f} km\n" \
                        f"Radius vector:               {sf_vector(cases_dict[case]['r_vector'][0], 3)} km\n" \
                        f"Radius mag (2):              {cases_dict[case]['r_vector'][1]:.3f} km\n" \
                        f"Velocity vector:             {sf_vector(cases_dict[case]['v_vector'][0], 3)} km/s\n" \
                        f"Velocity mag (2):            {cases_dict[case]['v_vector'][1]:.3f} km/s\n" \
                        f"Return of check function:    {check_r_and_v_with_h(cases_dict[case])}\n"\
                        f"---------------------------------------------\n"


        with open('output/problem_2_output.txt', 'a') as output:
            output.write(output_string)

    # calcaulte E from M
    # calcualte f from E using tan formula
    # calculate r mag from elems and f
    # calcualte r vector
    # calculate h value
    # calcualte v vector

    return 0


if __name__ == '__main__':
    main()
