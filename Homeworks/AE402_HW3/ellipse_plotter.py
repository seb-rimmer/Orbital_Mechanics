import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import cos, sin, pi, sqrt, ceil, atan
import numpy as np

def plot_transfer(r2, sma, theta):

    #  Create a figure and axis
    fig, ax = plt.subplots()

    # Calcualte the ellipse parameters (center, width, height, and angle)
    r1 = 1
    r2 = np.linalg.norm(r2)

    # this must definitely be correct
    e = (sma - 1) / sma

    x_j = r2 * cos(theta) + sma*e
    y_j = r2 * sin(theta)

    b = sqrt((y_j)**2 * 1/(1 - (x_j**2 / sma**2)))
    
    print(f'a = {sma}')
    print(f'b = {b}')
    print(f'e = {e}')

    x_traj = np.linspace(r2 * cos(theta), 1, 100)
    y_traj = np.array([b*sqrt(1 - ( (x + sma*e) / sma)**2) for x in x_traj])

    center = (0, 0)
    earth_pos = (1, 0)
    target_pos = (cos(theta)*r2, sin(theta)*r2)
    ellipse_centre = (-sma*e, 0)

    # Plot Earth orbit
    earth_orbit = patches.Circle(center, r1, fill=False, color='k', linestyle='--')
    earth = patches.Circle(earth_pos, 0.1, fill=True, color='g')

    # Plot target orbit
    target_orbit = patches.Circle(center, r2, fill=False, color='k', linestyle='--')
    target = patches.Circle(target_pos, 0.1, fill=True, color='orange')

    ellipse = patches.Ellipse(ellipse_centre, 2*sma, 2*b, fill=False, color='b', linestyle='--')
    ax.plot(x_traj, y_traj, color='blue')

    # centre and focci positioning
    ax.plot(ellipse_centre[0], ellipse_centre[1],  marker='*')
    ax.plot(ellipse_centre[0]+sma*e, ellipse_centre[1],  marker='*', color='red')

    # Plot vectors
    ax.plot([0, cos(theta)*r2], [0, sin(theta)*r2], label='My Line', color='red', linestyle='-', marker='')
    ax.plot([0, 1], [0, 0], label='My Line', color='red', linestyle='-', marker='')

    # Add the ellipse to the axis
    ax.add_patch(ellipse)
    ax.add_patch(earth_orbit)
    # ax.add_patch(earth)
    ax.add_patch(target_orbit)
    # ax.add_patch(target)

    # Set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')

    # Set the axis limits to show the entire ellipse
    lim = (max(r1, r2)*1.2)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    # Optional: Add labels and title
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plotting an Ellipse')

    # Show the plot
    plt.grid()
    plt.show()

def main():
    # plot_transfer(5.2, 5.32, 147)

    return 0

if __name__ == '__main__':
    main()