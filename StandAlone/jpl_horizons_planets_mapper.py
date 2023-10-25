#!/usr/bin/env python
#
# Description: 
# 

import numpy as np
import requests
import json
import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

from math import sqrt 

import functions

def create_json_obj(data:list):

    mu = float(1.32712440018E+11)

    for object in data:

        print(object)

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

        new_dict_obj['a'] = functions.a_from_vectors(r, v, mu)
        new_dict_obj['e'] = np.linalg.norm(functions.eccentricity_from_vectors(r, v, mu))
        new_dict_obj['i'] = functions.inclination(r, v)
        new_dict_obj['cap_ohm'] = functions.ra_o_an(r, v)
        new_dict_obj['low_ohm'] = functions.arg_of_periapse(r, v, mu)
        new_dict_obj['f'] = functions.true_anom_f(r, v, mu)

        print(new_dict_obj)

    return new_dict_obj


def plot_vectors(bodies:list, canonical):

    if canonical:
        divider = 149597871
    else:
        divider = 1

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
        target_orbit = patches.Circle(center, r_mag, fill=False, color='k', linestyle='--')
        
        ellipse_centre = (-a*e, 0)
        ellipse = patches.Ellipse(ellipse_centre, 2*a, 2*b, angle=data['cap_ohm'], fill=False, color='b', linestyle='--')

        target = patches.Circle(r_vect, 
                                0.2, 
                                fill=True, 
                                color=(random.random(), random.random(), random.random()),
                                label=data['Name'])

        ax.add_patch(ellipse)
        # ax.add_patch(target_orbit)
        ax.add_patch(target)
        ax.plot(ellipse_centre, marker='*')

    # Set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')

    # Set the axis limits to show the entire ellipse
    lim = r_mag
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.legend()

    # Optional: Add labels and title
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plotting Solar Orbits')

    # Show the plot
    plt.grid()
    plt.show()

    return 0


def main():

    # Define the time span:
    # Get the current time
    current_time = datetime.datetime.now()

    # Format the current time as a string
    now = current_time.strftime("%Y-%m-%d %H:%M:%S")
    
    # bodies = [i for i in range(3, 5)]
    bodies = [3, 20065803 ]
    body_positions = []

    for body in bodies:
        
        api_url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
        api_url += "?format=json&EPHEM_TYPE=VECTORS&OBJ_DATA=NO"
        api_url += f"&COMMAND='{body}'"
        api_url += f"&CENTER='*@sun'"
        api_url += f"&TLIST='{now}'"

        # Help query:

        try:
            response = requests.get(api_url)

            # Check if the request was successful (status code 200)
            if response.status_code == 200:
                # The response content will contain the data returned by the API
                data = response.json()  # Assuming the response is in JSON format
        
                # print(f'{}')
                # for line in data['result']:
                print(data['result'])
                print("------------ END OF REQUEST --------------\n")

                with open('request_output.json', "w") as json_file:
                    json.dump(data, json_file)
                
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

                body_positions += [body_dict]

            else:
                print(f"Request failed with status code: {response.status_code}")

        except requests.exceptions.RequestException as e:
            print(f"Request exception: {e}")

        
    plot_vectors(body_positions, canonical=True)


if __name__ == '__main__':
    main()
    