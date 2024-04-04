"""
    Write a computer program to compute the position 𝒓(𝑡) and velocity 𝒗(𝑡) for the two-body
    problem at arbitrary time, given 𝒓 (𝑡) and 𝒗 (𝑡).

    Run the program for the following initial conditions and determine position 𝒓(𝑡) and
    velocity = 𝒗(𝑡) after 𝑡 = 3 hours. Compute 𝒓(𝑡) and 𝒗(𝑡) for 100 equally spaced times
    between 𝑡 = 0 and 𝑡 = 3 hours. Plot the resulting orbit.

    𝒓 (𝑡) = [−8903.833 1208.356 213.066] km
    𝒗 (𝑡) = [−0.971 −6.065 −1.069] km/s

    Answers:
    r_t1_vector = [6974, -620, -109]  # km
    v_t1_vector = [0.64, 7.85, 1.38]  # km/s

    Steps:
    1. Find a
    2. Find e
    3. Find E0 and then E
    4. Determine F and G constants for r_t1
    5. Calculate r_t1
    6. Determine F* and G* constants for v_t1
    7. Calculate v_t1

"""
import numpy as np
import math
import matplotlib.pyplot as plt
import os

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

    return math.atan2(tan_E0_num, tan_E0_denom)


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
    adj_constant = E_0 - eccentricity * math.sin(E_0)

    E_iter = n * delta_t
    g = 1
    i = 0

    # Actual iteration loop
    while abs(g) > 1e-13:
        g = (E_iter - eccentricity * math.sin(E_iter)) - (n * delta_t) - adj_constant
        dgdE = 1 - eccentricity * math.cos(E_iter)

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
    F_value = 1 - (semi_maj_ax / r_mag) * (1 - math.cos(final_e - initial_e0))

    return F_value


def g_value_for_r_t1(semi_maj_ax: float,
                     mu: float,
                     delta_t: float,
                     initial_e0: float,
                     final_e: float) -> float:

    G_value = (delta_t) - ((semi_maj_ax ** 3 / mu) ** 0.5) * \
              ((final_e - initial_e0) - math.sin((final_e - initial_e0)))

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
                  * math.sin(final_e - initial_e0)

    return F_dot_value


def g_dot_value_for_v_t1(semi_maj_ax: float,
                         mu: float,
                         r1_vector: np.array,
                         delta_t: float,
                         initial_e0: float,
                         final_e: float) -> float:
    
    r1_mag = np.linalg.norm(r1_vector)
    G_dot_value = 1 - (semi_maj_ax / r1_mag) * (1 - math.cos(final_e - initial_e0))

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

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - math.cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - math.sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    return r_t1

def r_as_function_of_t_2(mu, r0_vector, v0_vector, dt):
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

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - math.cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - math.sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    return r_t1


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

    return list(vector_arr)


def h_value_from_elements(semi_maj_axis: float,
                          eccentricity: float,
                          mu: float) -> float:
    """
    Return magnitude of h paramter from mu, e and a

    :param orb_elems_dict:  dictionary of orbital elements
    :param mu:              Mu value of central body
    :return:                h value, float
    """

    h = (mu * semi_maj_axis * (1 - eccentricity ** 2)) ** 0.5

    return h


def check_r_and_v_with_h(h_mag_value: float,
                         r_vector: np.ndarray,
                         v_vector: np.ndarray) -> bool:
    """
    The magnitude of cross product of r and v should equal the mag of h

    :param case_dict:
    :return:
    """

    cross_product = np.cross(r_vector, v_vector)
    mag = np.linalg.norm(cross_product)

    if (mag - h_mag_value) < 0.1:

        return True
    else:
        return False


def main():

    # Some test vectors for a lecture example (Lec 5 slide 6)
    test_r_t0_vector = np.array([-4743, 4743, 0])  # km
    test_v_t0_vector = np.array([-5.879, -4.223, 0])  # km/s
    test_delta_t = 600  # seconds

    # Answers:
    # r_t1_vector = [-7027.094, 1382.68, 0]  # km
    # v_t1_vector = [1.603, -6.503, 0]       # km/s

    # Data from question
    hw_problem_r_t0_vector = np.array([-8903.833, 1208.356, 213.066])  # km
    hw_problem_v_t0_vector = np.array([-0.971, -6.065, -1.069])  # km/s
    hw_problem_delta_t = 3600 * 3  # seconds

    # Answers:
    # r_t1_vector = [6974, -620, -109]  # km
    # v_t1_vector = [0.64, 7.85, 1.38]  # km/s

    # Constants
    mu_Earth = 398600

    cases = {"Lec 6 Test Case":
                 {"r_t0":
                      {"vector": np.array(test_r_t0_vector),  # km
                       "mag": np.linalg.norm(test_r_t0_vector)
                       },
                  "v_t0":
                      {"vector": np.array(test_v_t0_vector),  # km
                       "mag": np.linalg.norm(test_v_t0_vector)
                       },
                  "delta_t": test_delta_t},

             "HW1 Problem 3 Case":
                 {"r_t0":
                      {"vector": np.array(hw_problem_r_t0_vector),  # km
                       "mag": np.linalg.norm(hw_problem_r_t0_vector)
                       },
                  "v_t0":
                      {"vector": np.array(hw_problem_v_t0_vector),  # km
                       "mag": np.linalg.norm(hw_problem_v_t0_vector)
                       },
                  "delta_t": hw_problem_delta_t},
             }

    header_string = f"HW1 Problem 3 - r_t1/v_t1 from r_t0/v_t0 .\n" \
                    f"------------------------------------------------------\n"
    pwd = os.getcwd()
    with open('output/problem_3_output_new.txt', 'w') as output:
        output.write(header_string)

    for case in cases:

        intro_string = f"Case: {case}\n" \
                       f"Radius vector   {list(cases[case]['r_t0']['vector'])} km\n" \
                       f"Velocity vector {list(cases[case]['v_t0']['vector'])} km/s\n" \
                       f"Delta t:        {cases[case]['delta_t']}\n" \
                       f"------------------------------------------------------\n"

        r_t0 = cases[case]['r_t0']['vector']
        v_t0 = cases[case]['v_t0']['vector']
        delta_t = cases[case]['delta_t']

        # determine semi-major axis
        a = a_from_vectors(r_t0, v_t0, mu_Earth)

        # determine eccentricity
        e = eccentricity_from_vectors(r_t0, v_t0, mu_Earth)

        # h value (for check later on)
        h_mag = h_value_from_elements(a, e[1], mu_Earth)

        # initial eccentric anom
        E0 = initial_eccentric_anom(r_t0, v_t0, mu_Earth, a)

        # mean angular rate and delta_t to work out Mean Anomaly M for finding E
        n = (mu_Earth / (a ** 3)) ** 0.5
        E = kepler_E_solution_iteration(e[1], n, delta_t, E0)

        # Constants of F and G values for r(t)
        F = f_value_for_r_t1(a, r_t0, E0, E)
        G = g_value_for_r_t1(a, mu_Earth, delta_t, E0, E)

        r_t1_initial_calc = F * r_t0 + G * v_t0

        r_t1 = {"vector": np.array(r_t1_initial_calc),  # km
                "mag": np.linalg.norm(r_t1_initial_calc)}

        # Constants of F_dot and G_dot values for v(t)
        F_dot = f_dot_value_for_v_t1(a, mu_Earth, r_t0, r_t1['vector'], E0, E)
        G_dot = g_dot_value_for_v_t1(a, mu_Earth, r_t1['vector'], delta_t, E0, E)

        v_t1_initial_calc = F_dot * r_t0 + G_dot * v_t0

        v_t1 = {"vector": np.array(v_t1_initial_calc),  # km
                "mag": np.linalg.norm(v_t1_initial_calc)}

        # Now performing h check
        h_check = check_r_and_v_with_h(h_mag, r_t1['vector'] ,v_t1['vector'])

        # Checking step function
        r_new = r_as_function_of_t_2(mu_Earth, r_t0, v_t0, cases[case]['delta_t'])

        # stating the number of decimal places to put things to
        num_dp = 3
        calculated_parameters = f"Eccentricity (mag, vector):          {e[1]:.{num_dp}f}\n" \
                                f"Semi-maj axis (mag):                 {a:.{num_dp}f} km\n" \
                                f"Initial E0 value: (rad):             {E0:.{num_dp}f} rad\n" \
                                f"E value for t1: (rad):               {E:.{num_dp}f} rad\n" \
                                f"F constant for r_t1:                 {F:.{num_dp}f} \n" \
                                f"G constant for r_t1:                 {G:.{num_dp}f} \n" \
                                f"r_t1 vector:                         [{r_t1['vector'][0]:.{num_dp}f}, {r_t1['vector'][1]:.{num_dp}f}, 0] km \n" \
                                f"F_dot constant for r_t1:             {F_dot:.{num_dp}f} \n" \
                                f"G_dot constant for r_t1:             {G_dot:.{num_dp}f} \n" \
                                f"v_t1 vector:                         [{v_t1['vector'][0]:.{num_dp}f}, {v_t1['vector'][1]:.{num_dp}f}, 0,] km/s \n" \
                                f"h_mag check result:                  {h_check} \n" \
                                f"------------------------------------------------------\n"

        print(intro_string + calculated_parameters)

        with open('output/problem_3_output_new.txt', 'a') as output:
            output.write(intro_string + calculated_parameters)

    t_step = np.linspace(0, 3 * 3600, 100)
    r_steps = [r_as_function_of_t(mu_Earth, e[1], a, E0, hw_problem_r_t0_vector, hw_problem_v_t0_vector, i) for i in t_step]

    x, y, z = [r[0] for r in r_steps], [r[1] for r in r_steps], [r[2] for r in r_steps]

    fig = plt.figure()
    ax1 = plt.axes(projection='3d')

    ax1.scatter3D(x, y, z, s=10, label="Plot of r(t) vector between t=0 and t=3 hours")
    ax1.scatter3D([0], [0], [0], c='g', s=6378*2)
    ax1.scatter3D([0], [0], [0], c='g', s=10, label="Representation of Earth (radius=6378km)")

    ax1.set_xlabel("i vector direction (x, km)")
    ax1.set_ylabel("j vector direction (y, km)")
    ax1.set_zlabel("k vector direction (z, km)")

    ax1.legend()
    fig.savefig("plots/rt_vector_plot_in_space.png")

    plt.show()

    return 0


if __name__ == '__main__':
    main()
