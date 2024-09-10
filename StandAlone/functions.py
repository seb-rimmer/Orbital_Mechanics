#!/usr/bin/env python
#
# Description: 
# 
import numpy as np
from math import acos, pi, cos, sin, sqrt
from math import acos, pi, cos, sin, atan2
import requests
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

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
    elems_dict['e'] = eccentricity_from_vectors(r, v, mu)[1]
    elems_dict['i'] = inclination(r, v)
    elems_dict['cap_ohm'] = ra_o_an(r, v)
    elems_dict['low_ohm'] = arg_of_periapse(r, v, mu)
    elems_dict['f'] = true_anom_f(r, v, mu)

    return elems_dict

def radial_vect_from_orbital_elems(a, e, inc, cap_ohm, low_ohm, f):

    i, j, k = np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])

    # work out radial distance using standard ellipse formula
    e = np.linalg.norm(e)
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

def create_json_obj(data:list):

    mu = float(1.32712440018E+11)

    for object in data:

        # print(object)

        new_dict_obj = {}
        
        # create time
        new_dict_obj['time'] = object[0]

        # create X, Y, Z
        x_index = object[1].find('X =')+3
        new_dict_obj['X']  = float(object[1][x_index:x_index+22])
        
        y_index = object[1].find('Y =')+3
        new_dict_obj['Y']  = float(object[1][y_index:y_index+22])
        
        z_index = object[1].find('Z =')+3
        new_dict_obj['Z']  = float(object[1][z_index:z_index+22])
        
        # create VX, VY, VZ
        x_index = object[2].find('VX')+3
        new_dict_obj['VX']  = float(object[2][x_index:x_index+22])
        
        y_index = object[2].find('VY')+3
        new_dict_obj['VY']  = float(object[2][y_index:y_index+22])
        
        z_index = object[2].find('VZ')+3
        new_dict_obj['VZ']  = float(object[2][z_index:z_index+22])

        # create orbital elements around the sun
        r = np.array([new_dict_obj['X'], new_dict_obj['Y'], new_dict_obj['Z']])
        v = np.array([new_dict_obj['VX'], new_dict_obj['VY'], new_dict_obj['VZ']])

        new_dict_obj['a'] = a_from_vectors(r, v, mu)
        new_dict_obj['e'] = eccentricity_from_vectors(r, v, mu)[1]
        new_dict_obj['i'] = inclination(r, v)
        new_dict_obj['cap_ohm'] = ra_o_an(r, v)
        new_dict_obj['low_ohm'] = arg_of_periapse(r, v, mu)
        new_dict_obj['f'] = true_anom_f(r, v, mu)

        # print(new_dict_obj)

    return new_dict_obj

def jpl_body_request(body_code, time='2000-Jan-01 00:00:00'):

    api_url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    api_url += "?format=json&EPHEM_TYPE=VECTORS&OBJ_DATA=NO"
    api_url += f"&COMMAND='{body_code}'"
    api_url += f"&CENTER='500@0'"
    api_url += f"&TLIST='{time}'"

    # Help query:

    try:
        response = requests.get(api_url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
                       
            # The response content will contain the data returned by the API
            data = response.json()  # Assuming the response is in JSON format

            with open('request_output.json', "w") as json_file:
                json.dump(data, json_file, indent=4)

            # Check if actually have emphemeris data returned:
            if "error" in data:
                print(f"Error in emphemerides request: {data['error']}")
                return 0
            
            else:
                # extract name of target body
                name_start = data['result'].find('Target body name:')
                target_name = data['result'][name_start+18:name_start+49].rstrip()

                elements_start = data['result'].find('$$SOE')
                elements_end = data['result'].find('$$EOE')

                data['result'] = data['result'][elements_start+6:elements_end]
                data['result'] = data['result'].split('\n')

                vector_set = [data['result']]
                body_dict = create_json_obj(vector_set)
                body_dict['Name'] = target_name

        else:
            print(f"Request failed with status code: {response.status_code}")
            return 0

    except requests.exceptions.RequestException as e:
        print(f"Request exception: {e}")

    return body_dict

def plot_vectors(bodies:list):

    #  Create a figure and axis
    fig, ax = plt.subplots()
    
    # Add the Sun
    target = patches.Circle([0, 0],  
                            2e7, 
                            fill=True, 
                            color=(0.8, 0.8, 0),
                            label='Sun')

    ax.add_patch(target)
    

    for data in bodies:
  
        x = float(data['X'])
        y = float(data['Y'])

        r_vect = np.array([x, y])
        r_mag = np.linalg.norm(r_vect)
        
        # Plot target orbit
        f = np.linspace(0, 360, 360)
        r_vectors_check = np.array([radial_vect_from_orbital_elems(data['a'], 
                                                                  data['e'], 
                                                                  data['i'], 
                                                                  data['cap_ohm'], 
                                                                  data['low_ohm'], 
                                                                  pos) for pos in f])
        
        x, y = [pos[0] for pos in r_vectors_check], [pos[1] for pos in r_vectors_check]

        target = patches.Circle(r_vect, 
                                1e7, 
                                fill=True, 
                                color=(random.random(), random.random(), random.random()),
                                label=f"{data['name']}, {data['Time']}")

        ax.add_patch(target)
        ax.plot(x, y, color='k', linestyle='--', linewidth=0.5)

    # Set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')

    # Set the axis limits to show the entire ellipse
    plt.autoscale()
    ax.legend()

    # Optional: Add labels and title
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plotting Solar Orbits at time: ')

    # Show the plot
    plt.grid()
    plt.show()

    return 0


def main():

    mu = float(1.32712440018E+11)
    r0 = np.array([160077584.6301427, -51942023.9865322, -10038863.59733499])
    v0 = np.array([2.582958403490801, 32.05362233994889, 0.4121945888753515])

    # elems = orbital_elems_from_vectors(r, v, mu)

    # # print(elems)

    # r_vector_check = radial_vect_from_orbital_elems(elems['a'], elems['e'], elems['i'], elems['cap_ohm'], elems['low_ohm'], elems['f'])
    # # print(r_vector_check)

    # print(v_vis_viva(1.515e+8, 2.457e+8, mu))

    rf, vf = r_and_v_as_function_of_t(mu, r0, v0, 3600*780)
    print(rf, vf)
    
    elems_f = orbital_elems_from_vectors(rf, vf, mu)

    for key, value in elems_f.items():
        print(f"{key}: {value}")

    pass


if __name__ == '__main__':
    main()

    # test commit 