import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from matplotlib.widgets import Slider, Button, RadioButtons

import cv2

import json
import datetime
import jpl_horizons_planets_mapper_copy

def update(slider, scat_e, x_e, y_e, scat_d, x_d, y_d,scat_m, x_m, y_m):

    frame = int(slider)

    # for each frame, update the data stored on each artist.
    x_t_e = x_e[frame-1]
    y_t_e = y_e[frame-1]
    x_t_d = x_d[frame-1]
    y_t_d = y_d[frame-1]
    x_t_m = x_m[frame-1]
    y_t_m = y_m[frame-1]
    
    # update the scatter plot:
    data_e = np.stack([x_t_e, y_t_e]).T
    data_d = np.stack([x_t_d, y_t_d]).T
    data_m = np.stack([x_t_m, y_t_m]).T

    scat_e.set_offsets(data_e)
    scat_d.set_offsets(data_d)
    scat_m.set_offsets(data_m)

    return (scat_e, scat_d)

def main():

    ## Parse data from json file
    # ----------------------------

    # Path to the JSON file
    file_path = r"C:\Users\Sebastian.Rimmer\OneDrive - ESA\Documents\Learning\Orbital_Mechanics\output\bodies_output.json"

    # Load the JSON file into a dictionary
    with open(file_path, 'r') as f:
        data_dict = json.load(f)
        earth = data_dict[0]
        mars = data_dict[1]
        didymos = data_dict[2]

    ## Create data for plot
    # ----------------------------
    
    mu_sun = float(1.32712440018E+11)   
    start = datetime.datetime.strptime(earth['start_time'], '%Y-%b-%d %H:%M:%S.%f')
    end = datetime.datetime(2026, 12, 16)

    dt_days = (end - start).days
    r0_vector_earth,  v0_vector_earth = np.array([earth['X'], earth['Y'], earth['Z']]), np.array([earth['VX'], earth['VY'], earth['VZ']])
    r0_vector_didymos,  v0_vector_didymos = np.array([didymos['X'], didymos['Y'], didymos['Z']]), np.array([didymos['VX'], didymos['VY'], didymos['VZ']])
    r0_vector_mars,  v0_vector_mars = np.array([mars['X'], mars['Y'], mars['Z']]), np.array([mars['VX'], mars['VY'], mars['VZ']])
    
    x_e, y_e = [], []
    x_d, y_d = [], []
    x_m, y_m = [], []

    for i in range(dt_days):
        
        r1_vector_e, v1_vector_e = jpl_horizons_planets_mapper_copy.r_and_v_as_function_of_t(mu_sun, r0_vector_earth, v0_vector_earth, 24*3600)
        r1_vector_d, v1_vector_d = jpl_horizons_planets_mapper_copy.r_and_v_as_function_of_t(mu_sun, r0_vector_didymos, v0_vector_didymos, 24*3600)
        r1_vector_m, v1_vector_m = jpl_horizons_planets_mapper_copy.r_and_v_as_function_of_t(mu_sun, r0_vector_mars, v0_vector_mars, 24*3600)

        x_e.append(r1_vector_e[0])
        y_e.append(r1_vector_e[1])
        x_d.append(r1_vector_d[0])
        y_d.append(r1_vector_d[1])
        x_m.append(r1_vector_m[0])
        y_m.append(r1_vector_m[1])

        r0_vector_earth, v0_vector_earth = r1_vector_e, v1_vector_e
        r0_vector_didymos, v0_vector_didymos = r1_vector_d, v1_vector_d
        r0_vector_mars, v0_vector_mars = r1_vector_m, v1_vector_m
    
    x_e, y_e, = np.array(x_e), np.array(y_e)
    x_d, y_d, = np.array(x_d), np.array(y_d)
    x_m, y_m, = np.array(x_m), np.array(y_m)

    ## Plot
    # ----------------------------
    
    # Define the figure
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.25)

    plt.axis('equal')
    plt.grid(True)

    # Plotting orbit tracks, Sun, initial positions
    # ----------------------------
    ax.plot(x_e, y_e, '--', color='k', linewidth=0.5)
    ax.plot(x_d, y_d, '--', color='k', linewidth=0.5)
    ax.plot(x_m, y_m, '--', color='k', linewidth=0.5)
    circle = patches.Circle((2.5, 20), radius=149597871/10, edgecolor='orange', facecolor='orange')
    ax.add_patch(circle)
    scat_e = ax.scatter(x_e[0], y_e[0], color='blue')
    scat_d = ax.scatter(x_d[0], y_d[0], color='gray')
    scat_m = ax.scatter(x_m[0], y_m[0], color='orange')

    # Slider bullshits
    # ----------------------------
    ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sDays = Slider(ax_slider, 'Days', 0, 735, valinit=735, valstep=1) # this needs to be initialised at 735 to fit in full range
    sDays.on_changed(lambda val: update(val, scat_e, x_e, y_e, scat_d, x_d, y_d, scat_m, x_m, y_m))

    plt.show()

    return 0

if __name__ == '__main__':
    main()
    