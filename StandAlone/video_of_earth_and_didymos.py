import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import cv2

import json
import datetime
import jpl_horizons_planets_mapper_copy

def update(frame, x, y, fig, ax, output_video):
    ax.clear()
    ax.plot(x[:frame], y[:frame])
    ax.set_xlim(-149597871*2, 149597871*2)
    ax.set_ylim(-149597871*2, 149597871*2)
    ax.set_title('Frame {}'.format(frame))
    # plt.tight_layout()
    ax.set_aspect('equal')
    fig.canvas.draw()

    data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    data = cv2.cvtColor(data, cv2.COLOR_RGB2BGR)
    output_video.write(data)

def main():

    ## Parse data from json file
    # ----------------------------

    # Path to the JSON file
    file_path = r"C:\Users\Sebastian.Rimmer\OneDrive - ESA\Documents\Learning\Orbital_Mechanics\output\bodies_output.json"

    # Load the JSON file into a dictionary
    with open(file_path, 'r') as f:
        data_dict = json.load(f)
        earth = data_dict[0]
        didymos = data_dict[1]

    # Print the dictionary
    # print(json.dumps(earth))
    # print(json.dumps(didymos))
    
    ## Create data for video plot
    # ----------------------------
    
    mu_sun = float(1.32712440018E+11)   
    start = datetime.datetime.strptime(earth['start_time'], '%Y-%b-%d %H:%M:%S.%f')
    end = datetime.datetime(2026, 12, 16)

    dt_days = (end - start).days
    r0_vector_earth,  v0_vector_earth = np.array([earth['X'], earth['Y'], earth['Z']]), np.array([earth['VX'], earth['VY'], earth['VZ']])
    r0_vector_didymos,  v0_vector_didymos = np.array([didymos['X'], didymos['Y'], didymos['Z']]), np.array([didymos['VX'], didymos['VY'], didymos['VZ']])
    
    x_e, y_e = [], []
    x_d, y_d = [], []

    for i in range(dt_days):
        
        r1_vector_e, v1_vector_e = jpl_horizons_planets_mapper_copy.r_and_v_as_function_of_t(mu_sun, r0_vector_earth, v0_vector_earth, 24*3600)
        r1_vector_d, v1_vector_d = jpl_horizons_planets_mapper_copy.r_and_v_as_function_of_t(mu_sun, r0_vector_didymos, v0_vector_didymos, 24*3600)

        x_e.append(r1_vector_e[0])
        y_e.append(r1_vector_e[1])
        x_d.append(r1_vector_d[0])
        y_d.append(r1_vector_d[1])

        r0_vector_earth, v0_vector_earth = r1_vector_e, v1_vector_e
        r0_vector_didymos, v0_vector_didymos = r1_vector_d, v1_vector_d
    
    x_e, y_e, = np.array(x_e), np.array(y_e)
    x_d, y_d, = np.array(x_d), np.array(y_d)

    ## Video plot
    # ----------------------------
    
    # Define the figure
    fig, ax = plt.subplots()
    ax.plot(x_e, y_e, '--', color='k', linewidth=0.5)
    ax.plot(x_d, y_d, '--', color='k', linewidth=0.5)
    circle = patches.Circle((2.5, 20), radius=149597871/10, edgecolor='orange', facecolor='orange')
    ax.add_patch(circle)
    scat_e = ax.scatter(x_e[0], y_e[0], color='blue')
    scat_d = ax.scatter(x_d[0], y_d[0], color='gray')
    
    start_day = 0
    title = ax.text(2e8, 2e8, f'time = {start_day} days')
    
    plt.axis('equal')
    plt.grid(True)
    # plt.show()

    ani = animation.FuncAnimation(fig=fig, func=update, fargs=(scat_e, x_e, y_e, scat_d, x_d, y_d, title), frames=dt_days, interval=10, repeat=False)
    plt.show()

    # ani.save(filename="tmp/pillow_example.apng", writer="pillow")

    return 0

def update(frame, scat_e, x_e, y_e, scat_d, x_d, y_d, title):

    # for each frame, update the data stored on each artist.
    x_t_e = x_e[frame]
    y_t_e = y_e[frame]
    x_t_d = x_d[frame]
    y_t_d = y_d[frame]
    
    # update the scatter plot:
    data_e = np.stack([x_t_e, y_t_e]).T
    data_d = np.stack([x_t_d, y_t_d]).T

    scat_e.set_offsets(data_e)
    scat_d.set_offsets(data_d)

    # update title
    title.set_text(f'time = {frame} days')

    return (scat_e, scat_d)

if __name__ == '__main__':
    main()
    