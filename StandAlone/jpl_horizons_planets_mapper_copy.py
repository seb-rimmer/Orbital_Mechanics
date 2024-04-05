#!/usr/bin/env python
#
# Description: 
#
# This file maps the position of major celestial bodies at the start and end of the Hera mission interplanetary trajectory.

import numpy as np
import requests
import json
import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

from math import sqrt , atan2, pi, cos, acos

import functions

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
    e_vector = eccentricity_from_vectors(r_vector, v_vector, mu)[0]

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

    e_vector = eccentricity_from_vectors(r_vector, v_vector, mu)[0]
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

    F_value = 1 - (a / np.linalg.norm(r0_vector) * (1 - math.cos(E_t - E0)))
    G_value = dt - ((a ** 3 / mu) ** 0.5) * ((E_t - E0) - math.sin(E_t - E0))

    r_t1 = F_value * np.array(r0_vector) + G_value * np.array(v0_vector)

    F_dot = f_dot_value_for_v_t1(a, mu, r0_vector, r_t1, E0, E_t)
    G_dot = g_dot_value_for_v_t1(a, mu, r_t1, dt, E0, E_t)

    v_t1 = F_dot * r0_vector + G_dot * v0_vector

    return r_t1, v_t1

def create_json_obj(name,
                    start_time,
                    end_time,
                    x,y,z,
                    vx, vy, vz,
                    ):

    mu = float(1.32712440018E+11)

    new_dict_obj = {}

    if start_time is str and end_time is str:
        datetime_object_start = datetime.datetime.strptime(object[0], '%Y-%b-%d %H:%M:%S.%f')
        datetime_object_end = datetime.datetime.strptime(object[1], '%Y-%b-%d %H:%M:%S.%f')
    else:
        datetime_object_start = start_time
        datetime_object_end = end_time

    # Add name
    new_dict_obj['name'] = name

    # Add time to object
    new_dict_obj['start_time'] = datetime_object_start
    new_dict_obj['end_time'] = datetime_object_end
    
    # create X, Y, Z
    new_dict_obj['X']  = x
    new_dict_obj['Y']  = y
    new_dict_obj['Z']  = z
    
    # create VX, VY, VZ
    new_dict_obj['VX']  = vx
    new_dict_obj['VY']  = vy
    new_dict_obj['VZ']  = vz

    # create orbital elements around the sun
    r = np.array([x, y, z])
    v = np.array([vx, vy, vz])

    new_dict_obj['a'] = functions.a_from_vectors(r, v, mu)
    new_dict_obj['e'] = np.linalg.norm(functions.eccentricity_from_vectors(r, v, mu))
    new_dict_obj['i'] = functions.inclination(r, v)
    new_dict_obj['cap_ohm'] = functions.ra_o_an(r, v)
    new_dict_obj['low_ohm'] = functions.arg_of_periapse(r, v, mu)
    new_dict_obj['f'] = functions.true_anom_f(r, v, mu)

    # print(new_dict_obj)

    return new_dict_obj


def plot_vectors(bodies:list, canonical):

    if canonical:
        divider = 149597871
        radius = 149597871/10
    else:
        divider = 1
        radius = 0.05

    #  Create a figure and axis
    fig, ax = plt.subplots()
    
    center = (0, 0)

    for data in bodies:

        # Plot things as an ellipse, so need to calculate semi-minor axis:
        e = data['e']
        a = data['a'] / divider
        b = sqrt( a**2 - (e*a)**2 )

        x = float(data['X']) / divider
        y = float(data['Y']) / divider

        r_vect = np.array([x, y])
        r_mag = np.linalg.norm(r_vect)
        
        # Plot target orbit
        f = np.linspace(0, 360, 360)
        r_vectors_check = np.array([functions.radial_vect_from_orbital_elems(data['a'], 
                                                                  data['e'], 
                                                                  data['i'], 
                                                                  data['cap_ohm'], 
                                                                  data['low_ohm'], 
                                                                  pos) for pos in f])
        
        x, y = [pos[0]/divider for pos in r_vectors_check], [pos[1]/divider for pos in r_vectors_check]

        target = patches.Circle(r_vect, 
                                1e7, 
                                fill=True, 
                                color=data['color'],
                                label=f"{data['name']} ,  {data['Time']}")

        ax.add_patch(target)
        ax.plot(x, y, color='k', linestyle='--', linewidth=0.5)

    # Set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')

    # Set the axis limits to show the entire ellipse
    # lim = r_mag
    # ax.set_xlim(-lim, lim)
    # ax.set_ylim(-lim, lim)
    ax.legend()

    # Optional: Add labels and title
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plotting Solar Orbits')

    # Show the plot
    plt.grid()
    plt.show()
    plt.savefig('bodies_at_depart_and_arrive.png')

    return 0

def jpl_horizons_request(body, time):

    api_url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    api_url += "?format=json&EPHEM_TYPE=VECTORS&OBJ_DATA=NO"
    api_url += f"&COMMAND='{body}'"
    api_url += f"&CENTER='*@sun'"
    api_url += f"&TLIST='{time}'"

    response = requests.get(api_url)

    if response.status_code == 200:
                    # The response content will contain the data returned by the API
                    return response.json()  # Assuming the response is in JSON format
    else:
        return 0


def main():

    # Define the time span:
    # Get the current time
    current_time = datetime.datetime.now()

    # Baseline launch date of 12th Oct 2024
    launch_date = datetime.datetime(2024, 10, 12)
    arrival_date = datetime.datetime(2026, 10, 16)

    # Format the current time as a string
    # now = current_time.strftime("%Y-%m-%d %H:%M:%S")
    
    # bodies
    # Earth - 3
    # Mars - 4
    # Didymos/Dimorhpos Barycentre - 20065803

    bodies = [3, 20065803]
    colors_dict = {3:'green', 4:'orange', 20065803:'darkgrey'}
    mu_sun = float(1.32712440018E+11)

    body_positions = []
    time = launch_date
    time_label = " at launch"
    
    for i in range(0, 2):

        for body in bodies:
            
            # Help query:

            try:
                
                data = jpl_horizons_request(body, time)

                # Check if the request was successful (status code 200)
                if data != 0:

                    try:
                        if data['error']:
                            print(f"Error in request for : {data['error']} ")
                            if "No ephemeris for target" in data['error']:
                                print(f"Need to do our own propagation for target ... ")
                                print(f"Executing propagator function ... ") 

                                # Send request for start time
                                didymos = body_positions[2]
                                time = arrival_date
                                r0_vector = np.array([didymos['X'], didymos['Y'], didymos['Z']])
                                v0_vector = np.array([didymos['VX'], didymos['VY'], didymos['VZ']])                         
                                
                                # propagate
                                dt_datetime = time - launch_date
                                dt_seconds = dt_datetime.days * 3600

                                r1_vector, v1_vector = r_and_v_as_function_of_t(mu_sun, r0_vector, v0_vector, dt_seconds)

                                print(r1_vector)

                                body_dict = create_json_obj('Didymos (arrival)',
                                                            launch_date,
                                                            arrival_date,
                                                            r1_vector[0],
                                                            r1_vector[1],
                                                            r1_vector[2],
                                                            v1_vector[0],
                                                            v1_vector[1],
                                                            v1_vector[2]
                                                            )

                                body_dict['Time'] = time_label
                                body_dict['color'] = colors_dict[body]
                                body_positions += [body_dict]

                    except KeyError:                     
                        
                        # write to json file
                        with open('request_output.json', "w") as json_file:
                            json.dump(data, json_file)
                        
                        # write text request to file
                        with open('request_output.txt', "w") as text_file:
                            text_file.write(data['result'])

                        # extract name of target body
                        name_start = data['result'].find('Target body name:')
                        target_name = data['result'][name_start+18:name_start+49].rstrip()
                        
                        # extract time
                        t_index_start = data['result'].find('Start time')
                        t_index_end = data['result'].find('Stop  time')

                        # Extract the time string
                        time_string_start =data['result'][t_index_start+23:t_index_start+48].strip()  # Extracts everything after the ':'
                        time_string_end =data['result'][t_index_end+23:t_index_end+48].strip()  # Extracts everything after the ':'

                        # Convert the time string to a datetime object


                        # new_dict_obj['start_time'] = t_stamp
                        # '2460595.500000000 = A.D. 2024-Oct-12 00:00:00.0000 TDB '
                        
                        # extract elements
                        elements_start = data['result'].find('$$SOE')
                        elements_end = data['result'].find('$$EOE')

                        data['result'] = data['result'][elements_start+6:elements_end]
                        data['result'] = f'{time_string_end}\n' + data['result']
                        data['result'] = f'{time_string_start}\n' + data['result']
                        data['result'] = data['result'].split('\n')

                        vector_set = data['result']

                        x_index = vector_set[3].find('X =')+3                       
                        y_index = vector_set[3].find('Y =')+3                        
                        z_index = vector_set[3].find('Z =')+3
                        vx_index = vector_set[4].find('VX')+3
                        vy_index = vector_set[4].find('VY')+3                        
                        vz_index = vector_set[4].find('VZ')+3

                        body_dict = create_json_obj(target_name,
                                                    vector_set[0],
                                                    vector_set[1],
                                                    float(vector_set[3][x_index:x_index+22]),
                                                    float(vector_set[3][y_index:y_index+22]),
                                                    float(vector_set[3][z_index:z_index+22]),
                                                    float(vector_set[4][vx_index:vx_index+22]),
                                                    float(vector_set[4][vy_index:vy_index+22]),
                                                    float(vector_set[4][vz_index:vz_index+22])
                                                    )
                        
                        body_dict['Time'] = time_label
                        body_dict['color'] = colors_dict[body]

                        body_positions += [body_dict]
                        print(f"Body: {body_dict['name']}")

                        # with open('output/bodies_output.json', 'a') as output:
                        #     json.dump(body_dict, output, indent=4)

                        # Load existing JSON data
                        with open('output/bodies_output.json', 'r') as f:
                            existing_data = json.load(f)

                        # Append new JSON object to existing data
                        existing_data.append(body_dict)

                        # Write the updated data back to the file
                        with open('output/bodies_output.json', 'w') as f:
                            json.dump(existing_data, f, indent=4)


                else:
                    print(f"Request failed with status code: {response.status_code}")

            except requests.exceptions.RequestException as e:
                print(f"Request exception: {e}")

            # Vector checks
            # r_mag = sqrt(body_dict['X']**2 + body_dict['Y']**2)
            # r_ellipse = (body_dict['a'] * (1 - body_dict['e']**2) ) / \
            #             ( 1 + body_dict['e']*cos(body_dict['f'] * pi/180))
            
            # theta = atan2(body_dict['Y'], body_dict['X']) * 180/pi
            

            # print(f"{body_dict['X']}")
            # print(f"{body_dict['Y']}")
            # print(f"{body_dict['Z']}")
            # print(f"{body_dict['VX']}")
            # print(f"{body_dict['VY']}")
            # print(f"{body_dict['VZ']}")

            # print(f'r_mag from vectors = {r_mag}')
            # print(f'r from ellipse eq  = {r_ellipse}')
            
            # print(f'theta = {theta}')
            # print(f"arg periapse = {body_dict['low_ohm']}")

        
        time = arrival_date
        time_label = " at arrival"

    plot_vectors(body_positions, canonical=False)

    return 0


if __name__ == '__main__':
    main()
    