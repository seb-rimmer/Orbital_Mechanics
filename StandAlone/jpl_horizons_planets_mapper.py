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

def create_json_obj(data:list):

    for object in data:

        print(object)

        new_dict_obj = {}
        
        # create time
        new_dict_obj['time'] = object[0]

        # create X, Y, Z
        x_index = object[1].find('X =')+3
        new_dict_obj['X']  = object[1][x_index:x_index+22]
        
        y_index = object[1].find('Y =')+3
        new_dict_obj['Y']  = object[1][y_index:y_index+22]
        
        z_index = object[1].find('Z =')+3
        new_dict_obj['Z']  = object[1][z_index:z_index+22]
        
        # create VX, VY, VZ
        x_index = object[1].find('VX =')+5
        new_dict_obj['VX']  = object[1][x_index:x_index+22]
        
        y_index = object[1].find('VY =')+5
        new_dict_obj['VY']  = object[1][y_index:y_index+22]
        
        z_index = object[1].find('VZ =')+5
        new_dict_obj['VZ']  = object[1][z_index:z_index+22]

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

        x = float(data['X']) / divider
        y = float(data['Y']) / divider

        r_vect = np.array([x, y])
        r_mag = np.linalg.norm(r_vect)

        # Plot target orbit
        target_orbit = patches.Circle(center, r_mag, fill=False, color='k', linestyle='--')
        target = patches.Circle(r_vect, 
                                0.3, 
                                fill=True, 
                                color=(random.random(), random.random(), random.random()),
                                label=data['Name'])

        ax.add_patch(target_orbit)
        ax.add_patch(target)

    # Set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')

    # Set the axis limits to show the entire ellipse
    lim = r_mag * 1.2
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
    
    bodies = [i for i in range(1, 6)]
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
    