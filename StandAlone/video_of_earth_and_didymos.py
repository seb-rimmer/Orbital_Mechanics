import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
    end = datetime.datetime(2026, 10, 16)

    dt_days = (end - start).days
    r0_vector = np.array([earth['X'], earth['Y'], earth['Z']])
    v0_vector = np.array([earth['VX'], earth['VY'], earth['VZ']])  
    
    x, y = [], []
    for i in range(dt_days):
        
        r1_vector, v1_vector = jpl_horizons_planets_mapper_copy.r_and_v_as_function_of_t(mu_sun, r0_vector, v0_vector, 24*3600)
        
        # print(r1_vector)
        x.append(r1_vector[0])
        y.append(r1_vector[1])

        r0_vector, v0_vector = r1_vector, v1_vector

    ## Video plot
    # ----------------------------
    
    # Define the figure
    fig, ax = plt.subplots()
    
    # Initialize the video writer - saving as mp4 file broken for now ... 
    output_video = cv2.VideoWriter('xy_data_video.mp4', cv2.VideoWriter_fourcc(*'MJPEG'), 30, (640, 480), False)

    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=len(x), fargs=(x, y, fig, ax, output_video), interval=50)

    # Close the video writer
    output_video.release()

    plt.show()

    return 0


if __name__ == '__main__':
    main()
    