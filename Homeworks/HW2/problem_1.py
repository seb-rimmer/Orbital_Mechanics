"""
    problem_1.py created by Seb at 13:55 15/10/2021
    Path: HW2
    -----------------------------------------------------------------------------
    Write a computer program to convert position and velocity to orbital elements.
    Run the program for the following cases and upload your numbers and code into
    Gradescope. For each case determine if the orbit is a low-Earth-orbit (LEO),
    medium-Earth-orbit (MEO), Geostationary orbit (GEO) or aGeostationary-transfer-
    orbit (GTO).

    mu = 398600 km^3/s^2

    Case 1
    r0 = [−6115.75, −6586.18, −58.65] km
    v0 = [4.42, −4.26, −1.08] km/s

    Case 2
    r0 = [6590, 0,0] km
    v0 = [0,10.153, 1.247] km/s

"""
import math
import numpy as np


def a_from_vectors(r_vector_mag: dict,
                   v_vector_mag: dict,
                   mu_val: float) -> float:
    """
    Calculates semi-major axis using rearranged vis-viva and vector magnitudes
    :param r_vector_mag: position vector as dictionary of {vector, magnitude}
    :param v_vector_mag: velocity vector as dictionary of {vector, magnitude}
    :param mu_val:   mu value for the central body
    :return: tuple of (vector, mag)
    """
    a = 1 / ((2 / r_vector_mag['mag']) - (v_vector_mag['mag'] ** 2 / mu_val))

    return a


def eccentricity_from_vectors(r_vector: dict,
                              v_vector: dict,
                              mu_val: float) -> tuple:
    """
    Calculates the eccentriciy of the orbit given r and v using formula from Lecture 5, slide 11

    :param r_vector: position vector as dictionary of {vector, magnitude}
    :param v_vector: velocity vector as dictionary of {vector, magnitude}
    :param mu_val:   mu value for the central body
    :return: tuple of (vector, mag)
    """

    r_v_dot = np.dot(r_vector['vector'], v_vector['vector'])
    # print(F"DOT product is {r_v_dot}")

    ecc_vector = (((v_vector['mag'] ** 2) / mu_val) - (1 / r_vector['mag'])) * r_vector['vector'] \
                 - (1 / mu_val) * r_v_dot * v_vector['vector']

    ecc_mag = np.linalg.norm(ecc_vector)

    return ecc_vector, ecc_mag


def h_value_from_vectors(r_vector: dict,
                         v_vector: dict) -> tuple:
    """
    Calculates the h value of the orbit given r and v

    :param r_vector: position vector as dictionary of {vector, magnitude}
    :param v_vector: velocity vector as dictionary of {vector, magnitude}
    :return: tuple of (vector, mag)
    """
    h_vector = np.cross(r_vector['vector'], v_vector['vector'])
    h_mag = np.linalg.norm(h_vector)

    return h_vector, h_mag


def n_value(h_value: tuple, k_unit_vector: np.ndarray) -> tuple:
    """
    Returns the n value of an orbit from h and k vectors as tuple of (vector, mag)
    :param h_value:         h_value of orbit as tuple (vector, magnitude)
    :param k_unit_vector:   unit vector in k plane
    :return:                tuple of (vector, mag)
    """

    h_norm = h_value[0] / np.linalg.norm(h_value[1])
    n_vector = np.cross(k_unit_vector, h_norm)

    return n_vector, np.linalg.norm(n_vector)


def inclination(h_value: tuple, k_unit_vector: np.ndarray) -> tuple:
    """
    Return inclination angle of orbit as tuple of (degrees, radians)

    :param h_value:         h_value of orbit as tuple (vector, magnitude)
    :param k_unit_vector:   unit vector in k plane
    :return:                inclination as tuple of (degrees, radians)
    """
    h_norm = h_value[0] / np.linalg.norm(h_value[1])

    cos_i = np.dot(h_norm, k_unit_vector)

    # math.acos returns values in radians
    i_radians = math.acos(cos_i)
    i_degrees = i_radians * 180 / math.pi

    return i_degrees, i_radians


def ra_o_an(n_value: tuple, i_unit_vector: np.ndarray) -> tuple:
    """
    Return Right Angle of Ascending Node of orbit as tuple of (degrees, radians)

    :param n_value:         previously calculated n value for orbit
    :param i_unit_vector:   unit vector in i plane
    :return:                tuple of (degrees, radians)
    """

    cos_omega = np.dot(n_value[0], i_unit_vector) / n_value[1]

    # math.acos returns values in radians
    ra_o_an_radians = math.acos(cos_omega)
    ra_o_an_degrees = ra_o_an_radians * 180 / math.pi

    return ra_o_an_degrees, ra_o_an_radians


def arg_of_periapse(eccentricity_value: tuple, n_value: tuple) -> tuple:
    """
    Return the argument of periapce of orbit as tuple of (degrees, radians)

    :param eccentricity_value:     eccentricity vector as tuple(vector, mag)
    :param n_value:                previously calculated n value for orbit
    :return:                       tuple of (degrees, radians)
    """

    cos_w = np.dot(n_value[0], eccentricity_value[0]) / (n_value[1] * eccentricity_value[1])

    # math.acos returns values in radians
    arg_of_p_radians = math.acos(cos_w)
    arg_of_p_degrees = arg_of_p_radians * 180 / math.pi

    return arg_of_p_degrees, arg_of_p_radians


def true_anom_f(r_vector_mag: dict, eccentricity_value: tuple) -> tuple:
    """
    Returns true anomaly of orbit given radius and orbit eccentricity

    :param r_vector_mag:           position vector as dictionary of {vector, magnitude}
    :param eccentricity_value:     eccentricity vector as tuple(vector, mag)
    :return:                       true anomaly f as tuple of (degrees, radians)
    """

    cos_f = np.dot(r_vector_mag['vector'], eccentricity_value[0]) / (r_vector_mag['mag'] * eccentricity_value[1])

    # math.acos returns values in radians
    true_anom_f_radians = math.acos(cos_f)
    true_anom_f_degrees = true_anom_f_radians * 180 / math.pi

    return true_anom_f_degrees, true_anom_f_radians

    return 0


def sf_vector(vector_arr: np.ndarray, num_sig_fig) -> list:
    """
    Simple function that just returns vector array with values at a specfied number of sig figs
    :param vector_arr:      np vector array with float values to lots of sf
    :param num_sig_fig:     num sig fig to return values with
    :return:                tuple of prettied up nd.np array
    """
    # conver to list
    vector_arr = list(vector_arr)

    vector_arr = [round(val, num_sig_fig) for val in vector_arr]

    return tuple(vector_arr)


def main():
    # variables
    mu = 398600
    earth_radius = 6378

    # Vectors as numpy arrays
    r0_vector_1, v0_vector_1 = np.array([-6115.75, -6586.18, -58.65]), np.array([4.42, -4.26, -1.08])
    r0_vector_2, v0_vector_2 = np.array([6590, 0, 0]), np.array([0, 10.153, 1.247])

    # Checking data from a HW4 Problem 1 question
    # r0_vector_2, v0_vector_2 = np.array([-3682.354, 5198.614, 1907.142]), np.array([-4.854, -4.884, 3.941])

    # Cases as dictionaries, with mag, to reduce computation in other functions
    r0_case_1 = {"vector": r0_vector_1, "mag": np.linalg.norm(r0_vector_1)}
    v0_case_1 = {"vector": v0_vector_1, "mag": np.linalg.norm(v0_vector_1)}

    r0_case_2 = {"vector": r0_vector_2, "mag": np.linalg.norm(r0_vector_2)}
    v0_case_2 = {"vector": v0_vector_2, "mag": np.linalg.norm(v0_vector_2)}

    case_1 = {'radius': r0_case_1, 'velocity': v0_case_1}
    case_2 = {'radius': r0_case_2, 'velocity': v0_case_2}

    # i, j, k vectors
    unit_vectors = {"i": np.array([1, 0, 0]),
                    "j": np.array([0, 1, 0]),
                    "k": np.array([0, 0, 1])}

    intro_string = f"HW2 Problem 1 - orbital elements from r and v vectors.\n" \
                   f"------------------------------------------------------\n" \
                   f"Case 1:  radius vector {list(case_1['radius']['vector'])} velocity vector {list(case_1['velocity']['vector'])}\n" \
                   f"Case 2:  radius vector {list(case_2['radius']['vector'])} velocity vector {list(case_2['velocity']['vector'])}\n\n"

    print(intro_string)

    with open('output/problem_1_output.txt', 'w') as output:
        output.write(intro_string)

    for i, case in enumerate([case_1, case_2]):
        e = eccentricity_from_vectors(case['radius'], case['velocity'], mu)
        a = a_from_vectors(case['radius'], case['velocity'], mu)
        h = h_value_from_vectors(case['radius'], case['velocity'])
        n = n_value(h, unit_vectors['k'])
        inc = inclination(h, unit_vectors['k'])
        asc_nod = ra_o_an(n, unit_vectors['i'])
        arg_p = arg_of_periapse(e, n)
        f = true_anom_f(case['radius'], e)
        periapse_alt = a * (1 - e[1]) - earth_radius
        apoapse_alt = a * (1 + e[1]) - earth_radius

        num_dp = 3
        output_string = f"Case {i+1} values:\n" \
                        f"---------------------------\n" \
                        f"Radius vector (mag):                 {np.linalg.norm(case['radius']['vector']):.{num_dp}f}\n" \
                        f"Velocity vector (mag):               {np.linalg.norm(case['velocity']['vector']):.{num_dp}f}\n" \
                        f"Eccentricity (mag, vector):          {e[1]:.{num_dp}f}, {sf_vector(e[0], 3)}\n" \
                        f"Semi-maj axis (mag):                 {a:.{num_dp}f} km\n" \
                        f"h value (mag, vector):               {h[1]:.{num_dp}f}, {sf_vector(h[0], 3)}\n" \
                        f"n value (mag, vector):               {n[1]:.{num_dp}f}, {sf_vector(n[0], 3)}\n" \
                        f"Inclination (degrees, radians):      {inc[0]:.{num_dp}f}, {inc[1]:.{num_dp}f}\n" \
                        f"RAoAN (degrees, radians):            {asc_nod[0]:.{num_dp}f}, {asc_nod[1]:.{num_dp}f}\n" \
                        f"Arg of periapse (degrees, radians):  {arg_p[0]:.{num_dp}f}, {arg_p[1]:.{num_dp}f}\n" \
                        f"True Anomaly f (degrees, radians):   {f[0]:.{num_dp}f}, {f[1]:.{num_dp}f}\n" \
                        f"Altitude at periapse (km):           {periapse_alt:.{num_dp}f}\n" \
                        f"Altitude at apoapse (km):            {apoapse_alt:.{num_dp}f}\n\n"
        print(output_string)

        with open('output/problem_1_output.txt', 'a') as output:
            output.write(output_string)

    conclusion_string = "Case 1 is in a LEO because the mean orbit altitude is between 200 and 2000km\n" \
                        "----------------------------------------------------------------------------------\n" \
                        "Case 2 is in a GTO because it is on a highly elliptical orbit, with the apoapse\n" \
                        "coinciding with the GEO altitude and periapse coinciding with a LEO altitude; it\n" \
                        "is likely it will perform some manoeuvre at apoapse to permanently move to GEO."
    print(conclusion_string)

    with open('output/problem_1_output.txt', 'a') as output:
        output.write(conclusion_string)

    return 0


if __name__ == '__main__':
    main()
