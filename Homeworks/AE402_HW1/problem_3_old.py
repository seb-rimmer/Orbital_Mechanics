"""
    Write a computer program to compute the position ð’“(ð‘¡) and velocity ð’—(ð‘¡) for the two-body
    problem at arbitrary time, given ð’“ (ð‘¡) and ð’— (ð‘¡).

    Run the program for the following initial conditions and determine position ð’“(ð‘¡) and
    velocity = ð’—(ð‘¡) after ð‘¡ = 3 hours. Compute ð’“(ð‘¡) and ð’—(ð‘¡) for 100 equally spaced times
    between ð‘¡ = 0 and ð‘¡ = 3 hours. Plot the resulting orbit.

    ð’“ (ð‘¡) = [âˆ’8903.833 1208.356 213.066] km
    ð’— (ð‘¡) = [âˆ’0.971 âˆ’6.065 âˆ’1.069] km/s

"""
import numpy as np
import math
import matplotlib.pyplot as plt

# Constants
mu_Earth = 398600


def kepler_E_solution_iteration(eccentricity, n, delta_t, E_0):
    """
    Kepler iteration function to give E, with an E0 not at perigee

    :param eccentricity: eccentricity
    :param n: mean angular rate
    :param delta_t: timestep for new E since E0
    :param E_0: value of initial eccentric anomaly
    :return: eccentric anomaly in radians
    """
    # M mean anomaly
    mean_anom = n * delta_t

    # Because not at perigee, need to use modified Kepler's with E0 as extra constant:
    adj_constant = E_0 - math.sin(E_0)

    E_iter = mean_anom
    g = 1
    i = 0
    while abs(g) > 1e-13:
        g = (E_iter - eccentricity * math.sin(E_iter)) - mean_anom - adj_constant
        dgdE = 1 - eccentricity * math.cos(E_iter)
        E_1 = E_iter - g / dgdE
        # updates
        E_iter = E_1
        i+=1
    # print(E_iter)
    return E_iter # 3.37


def r_as_function_of_t(e, a, E0, r0_vector, v0_vector, dt):
    """
    Return radius vector for the given input values.
    Calculates new F and G
    :param e: eccentricity
    :param a: semi-major axis
    :param E0: value of initial eccentric anomaly
    :param r0_vector: initial state of r
    :param v0_vector: initial state of v
    :param dt: timestep for new E since E0
    :return: new r vector from calculation involing F and G
    """
    n = (mu_Earth / (a ** 3)) ** 0.5

    E_t = kepler_E_solution_iteration(e, n, dt, E0)

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - math.cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu_Earth) ** 0.5) * ((E_t - E0) - math.sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    # print(E_t, dt, r_t1)
    return r_t1


def main():

    # Data from question
    r0_t = np.array([-8903.833, 1208.356, 213.066])  # km
    v0_t = np.array([-0.971, -6.065, -1.069])  # km/s

    r0_t_mag = np.linalg.norm(r0_t)
    v0_t_mag = np.linalg.norm(v0_t)

    # h as magnitude of cross product of r0_t and v0_t
    h = np.cross(r0_t, v0_t)
    h_mag = np.linalg.norm(h)

    # semi-major axis from vis-viva equation (km)
    a = (2 / r0_t_mag - (v0_t_mag ** 2) / mu_Earth) ** -1

    # time_period of orbit (s)
    T = 2 * math.pi * (a ** 3 / mu_Earth) ** 0.5

    # eccentricity from rearranged h equation
    e = (1 - (h_mag ** 2) / (mu_Earth * a)) ** 0.5
    e_mag = np.linalg.norm(e)

    # Testing kepler iteration function with example from Lecture 4 slide 5
    # ---------------------------------------------------------------------
    # example_ecc = 0.625
    # example_n = (mu_Earth / (6378 * 4) ** 3) ** 0.5  # rad/s^-1
    # example_delta_t = 4 * 3600   # seconds
    # example_E = kepler_E_solution_iteration(example_ecc,
    #                                         example_n,
    #                                         example_delta_t,
    #                                         0)
    # print(f"Testing the modified E: {example_E}")

    # Calculation of E0 - 3 methods
    # Method 1
    num = np.dot(r0_t, v0_t)
    denom = (a * mu_Earth) ** 0.5 * (1 - (r0_t_mag / a))
    E0 = math.atan2(float(num), float(denom))

    # mean angular rate and delta_t to work out Mean Anomaly M for finding E
    n = (mu_Earth / (a ** 3)) ** 0.5
    delta_t = 3600 * 3
    E = kepler_E_solution_iteration(e_mag, n, delta_t - T, E0)

    # Equations for F and G values for r(t)
    F_value = 1 - (a / r0_t_mag * (1 - math.cos(E - E0)))
    G_value = (delta_t - T) - ((a ** 3 / mu_Earth) ** 0.5 * (E - E0 - math.sin(E - E0)))

    r_t1 = F_value * r0_t + G_value * v0_t

    # Equations for F dot and G dot values for v(t)
    F_dot_val = -((mu_Earth * a) ** 0.5) / (np.linalg.norm(r_t1) * np.linalg.norm(r0_t)) * math.sin(E - E0)
    G_dot_val = 1 - a / np.linalg.norm(r_t1) * (1 - math.cos(E - E0))

    v_t1 = F_dot_val * r0_t + G_dot_val * v0_t

    # Output to console
    print(f"v mag: {v0_t_mag}, r mag: {r0_t_mag}")
    print(f"Semi major axis : {a:.3f}")
    # print(f"Eccentricity : {e}")
    print(f"E mag: {e_mag:.5f}")

    print(f"Orbit time period is: {T:.2f} seconds")

    print(f"h value: {h_mag:.3f}")
    print(f"E0 value  {E0}")
    print(f"E value {E}")

    print(f"F value:  {F_value:.2f}, G value:  {G_value:.2f}") # 0.09 , 3677
    print(f"r_1 value could be: {r_t1}, with a magnitude of {np.linalg.norm(r_t1):.2f}")
    print(f"v_1 value could be: {v_t1}, with a magnitude of {np.linalg.norm(v_t1):.2f}")

    # t_step = [i for i in range(0, 3*3600, 1)]
    t_step = np.linspace(0, 3*3600, 100)
    r_steps = [r_as_function_of_t(e, a, E0, r0_t, v0_t, i) for i in t_step]

    x, y, z = [r[0] for r in r_steps], [r[1] for r in r_steps], [r[2] for r in r_steps]

    print(f"last rt is : {r_steps[-1]}")

    fig = plt.figure()
    ax1 = plt.axes(projection='3d')

    ax1.scatter3D(x, y, z, s=10, label="Plot of r(t) vector between t=0 and t=3 hours")
    ax1.scatter3D([0], [0], [0], c='g', s=6378*2)
    ax1.scatter3D([0], [0], [0], c='g', s=10, label="Representation of Earth (radius=6378km)")

    ax1.set_xlabel("i vector direction (x, km)")
    ax1.set_ylabel("j vector direction (y, km)")
    ax1.set_zlabel("k vector direction (z, km)")

    ax1.legend()
    fig.savefig("rt_vector_plot_in_space.png")

    plt.show()

    return 0


if __name__ == '__main__':
    main()
