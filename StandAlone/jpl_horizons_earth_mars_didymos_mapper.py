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

from math import sqrt , atan2, pi, sin, cos, tan

import functions

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
        e = np.linalg.norm(data['e'])
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
                                color=(random.random(), random.random(), random.random()),
                                label=data['Name'])

        ax.add_patch(target)
        ax.plot(x, y, color='k', linestyle='--', linewidth=0.5)

    # Set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')

    # Set the axis limits to show the entire ellipse
    lim = r_mag
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.legend()

    # Add earth periphelion
    # r = bodies[0]['a'] * ( 1 - bodies[0]['e'])
    # theta = bodies[0]['low_ohm'] + bodies[0]['cap_ohm']
    # ax.plot([0, r * cos(theta * pi/180)], [0, r*sin(theta * pi/180)])

    # Add didymos aphelion
    # r = bodies[2]['a'] * ( 1 + bodies[2]['e'])
    # theta = bodies[2]['low_ohm'] + bodies[2]['cap_ohm'] + 180
    # ax.plot([0, r * cos(theta * pi/180)], [0, r*sin(theta * pi/180)])

    # Optional: Add labels and title
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plotting Solar Orbits at time: ')

    # Show the plot
    plt.grid()
    plt.show()

    return 0


def main():


    now = '2019-Oct-26 00:00:00'

    bodies = [3, 4, 20065803 ]
    body_positions = []

    for body in bodies:
        
        body_dict = functions.jpl_body_request(body, now)

        body_positions += [body_dict]

        # Vector checks
        r_mag = sqrt(body_dict['X']**2 + body_dict['Y']**2)
        r_ellipse = (body_dict['a'] * (1 - body_dict['e']**2) ) / \
                    ( 1 + body_dict['e']*cos(body_dict['f'] * pi/180))
        
        theta = atan2(body_dict['Y'], body_dict['X']) * 180/pi
        
        # print(f"Body: {body_dict['Name']}")
        # print(f"{body_dict['X']}")
        # print(f"{body_dict['Y']}")
        # print(f"{body_dict['Z']}")
        # print(f"{body_dict['VX']}")
        # print(f"{body_dict['VY']}")
        # print(f"{body_dict['VZ']}")

    plot_vectors(body_positions, canonical=False)


if __name__ == '__main__':
    main()
    