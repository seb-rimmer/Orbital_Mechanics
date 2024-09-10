"""

File: lamberts_between_actual_bodies.py
Author: Seb
Date: 2024-09-10
Description: This script will calcualte the dV between two bodies in the solar system using a lamberts transfer.

            Inputs - Starting and Target bodies, date of departure, date of arrival
            Ouputs -  dV required for transfer, information to plot the results

"""
import json
import functions
from math import sqrt, acos, pi
import numpy as np

def main():

    # Accept inputs for start, target bodies and time
    launch =  '2024-Oct-12 00:00:00'
    arrival = '2026-Dec-01 00:00:00'

    # quote JPL horizons for bodies info and positions
    # Path to the JSON file
    file_path = r"C:\Users\Sebastian.Rimmer\OneDrive - ESA\Documents\Learning\Orbital_Mechanics\output\bodies_output_static.json"

    # Load the JSON file into a dictionary
    with open(file_path, 'r') as f:
        data_dict = json.load(f)
        earth_launch = data_dict[0]
        mars_launch = data_dict[1]
        didymos_launch = data_dict[2]
        earth_arrival = data_dict[3]
        mars_arrival = data_dict[4]
        didymos_arrival = data_dict[5]

    # Plot where bodies are at start and arrival
    functions.plot_vectors([earth_launch, didymos_launch, earth_arrival, didymos_arrival])

    # extract r1, r2, separation anlge between two bodies. Make clear this is an approximation in the x-y plane
    xe, ye, ze = earth_launch['X'], earth_launch['Y'], earth_launch['Z']
    xd, yd, zd = didymos_arrival['X'], didymos_arrival['Y'], didymos_arrival['Z']

    r_earth = np.array([xe, ye, ze])
    r_didymos = np.array([xd, yd, zd])

    r1 = np.linalg.norm(r_earth)
    r2 = np.linalg.norm(r_didymos)

    theta_cos = np.dot(r_earth, r_didymos) / (r1*r2)
    theta = acos(theta_cos) * 180/pi
    
    print(f"R1 (Sun-Earth departure distance) : {r1:3.2f}\n"
          f"R2 (Sun-Didymos arrival distance) : {r2:3.2f}\n"
          f"theta (Earth-Didymos Separation angle) : {theta:3.2f}\n")
 
    # feed info to lamberts solution solver, prograde solutions only

    # accept data returned, translate into diagrams

    return 0

if __name__ == '__main__':
    main()