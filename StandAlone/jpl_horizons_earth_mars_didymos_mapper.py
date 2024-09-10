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
        r_vectors_check = np.array([functions.radial_vect_from_orbital_elems(data['a'], 
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
                                label=data['Name'])

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

    launch =  '2024-Oct-12 00:00:00'
    arrival = '2026-Dec-01 00:00:00'

    bodies = [3, 4, 20065803]
    colors = ['blue', 'orange', 'grey']
    body_positions = []

    for body in bodies:
        
        # Add body at launch
        body_dict = functions.jpl_body_request(body, launch)
        if body_dict != 0:
            body_positions += [body_dict]

        # Add body at arrival
        # body_dict = functions.jpl_body_request(body, arrival)
        # if body_dict != 0:
        #     body_positions += [body_dict]

    plot_vectors(body_positions)


if __name__ == '__main__':
    main()
    