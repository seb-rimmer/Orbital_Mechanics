import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from matplotlib.widgets import Slider, Button, RadioButtons

import cv2

import json
import datetime
import functions

def update(slider, scat_e, x_e, y_e, scat_d, x_d, y_d,scat_m, x_m, y_m, text_annotation):


    date_0 = datetime.datetime(2019, 10, 26)
    selected_date = date_0 + datetime.timedelta(days=int(slider))
    text_annotation.set_text(selected_date.strftime("%Y-%m-%d"))

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
    file_path = r"C:\Users\Sebastian.Rimmer\OneDrive - ESA\Documents\Learning\Orbital_Mechanics\output\bodies_output_static.json"

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
    end = datetime.datetime(2028, 12, 16)

    dt_days = (end - start).days + 1
    r0_vector_earth,  v0_vector_earth = np.array([earth['X'], earth['Y'], earth['Z']]), np.array([earth['VX'], earth['VY'], earth['VZ']])
    r0_vector_didymos,  v0_vector_didymos = np.array([didymos['X'], didymos['Y'], didymos['Z']]), np.array([didymos['VX'], didymos['VY'], didymos['VZ']])
    r0_vector_mars,  v0_vector_mars = np.array([mars['X'], mars['Y'], mars['Z']]), np.array([mars['VX'], mars['VY'], mars['VZ']])
    
    x_e, y_e = [], []
    x_d, y_d = [], []
    x_m, y_m = [], []
    theta_j2000_e, theta_j2000_d = np.array([]), np.array([])
    for i in range(dt_days):
        
        r1_vector_e, v1_vector_e = functions.r_and_v_as_function_of_t(mu_sun, r0_vector_earth, v0_vector_earth, 24*3600)
        r1_vector_d, v1_vector_d = functions.r_and_v_as_function_of_t(mu_sun, r0_vector_didymos, v0_vector_didymos, 24*3600)
        r1_vector_m, v1_vector_m = functions.r_and_v_as_function_of_t(mu_sun, r0_vector_mars, v0_vector_mars, 24*3600)

        x_e.append(r1_vector_e[0])
        y_e.append(r1_vector_e[1])
        x_d.append(r1_vector_d[0])
        y_d.append(r1_vector_d[1])
        x_m.append(r1_vector_m[0])
        y_m.append(r1_vector_m[1])

        r0_vector_earth, v0_vector_earth = r1_vector_e, v1_vector_e
        r0_vector_didymos, v0_vector_didymos = r1_vector_d, v1_vector_d
        r0_vector_mars, v0_vector_mars = r1_vector_m, v1_vector_m

        # true anomalies
        f_e, f_d = functions.true_anom_f(r1_vector_e, v1_vector_e, mu_sun), functions.true_anom_f(r1_vector_d, v1_vector_d, mu_sun)
        w_e, w_d = functions.arg_of_periapse(r1_vector_e, v1_vector_e, mu_sun), functions.arg_of_periapse(r1_vector_d, v1_vector_d, mu_sun)
        Ohm_e, Ohm_d = functions.ra_o_an(r1_vector_e, v1_vector_e), functions.ra_o_an(r1_vector_d, v1_vector_d)

        # theta in ecliptic j2000 frame, for reference
        sum_e = f_e + w_e + Ohm_e -720
        if sum_e < 0:
            sum_e += 360

        sum_d = f_d + w_d + Ohm_d
        if sum_d > 360:
            sum_d -= 360
        if sum_d > 360:
            sum_d -= 360

        theta_j2000_e = np.append(theta_j2000_e, sum_e)
        theta_j2000_d = np.append(theta_j2000_d, sum_d)

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
    ax_slider = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    sDays = Slider(ax_slider, 'Date', 0, dt_days, valinit=dt_days, valstep=1)

    # Add a text annotation to display the selected date
    text_annotation = ax_slider.text(1.05, 0.5, '', transform=ax_slider.transAxes, va='center')
    text_annotation.set_text((start + datetime.timedelta(days=int(sDays.val))).strftime("%Y-%m-%d"))

    sDays.on_changed(lambda val: update(val, scat_e, x_e, y_e, scat_d, x_d, y_d, scat_m, x_m, y_m, text_annotation))
    sDays.valtext.set_visible(False)

    # Calculating the distance for transfer opps
    # ----------------------------
    # Comparing circular orbit with time period of elliptical orbit. For visual sake, lets just plot sun distances for now
    plt.figure()

    separation_theta = np.array([])
    for i in range(dt_days):
        
        xi_e, yi_e = x_e[i], y_e[i] 
        xi_d, yi_d = x_d[i], y_d[i]
   
        r_earth_theta = np.sqrt(xi_e**2 + yi_e**2)
        r_didymos_theta = np.sqrt(xi_d**2 + yi_d**2)
        r_diff = np.sqrt(abs(xi_e-xi_d)**2 + abs(yi_e-yi_d)**2)
        theta_diff = (r_earth_theta**2 + r_didymos_theta**2 -  r_diff**2) / ( 2 * r_earth_theta * r_didymos_theta)
        theta_diff = np.arccos(theta_diff) * 180/np.pi

        separation_theta = np.append(separation_theta, theta_diff)

    plt.plot(pd.date_range(start=start, end=end, freq='D'), separation_theta)

    # for i in range(1, 10):
    #     plt.plot([365*i, 365*i], [0, 180], 'red')
    plt.grid(True)
    plt.ylim([0, 180])
    
   # Set x-axis label format
    plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m'))
    plt.gca().xaxis.set_major_locator(plt.matplotlib.dates.MonthLocator(interval=3))

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Add labels and title
    plt.xlabel('Date')
    plt.ylabel('Value')
    plt.title('Earth-Didymos separation angle')


    plt.show()

    # To be udapted to be representative of something similar too ... 
    # https://stackoverflow.com/questions/44985966/managing-dynamic-plotting-in-matplotlib-animation-module/44989063#44989063
    
    return 0

if __name__ == '__main__':
    main()
    