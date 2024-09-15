import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import cos, sin, pi, sqrt, ceil, atan2
import numpy as np
from scipy.optimize import fsolve

def focci_solver(p, *data):
    p1, r1, p2, r2, a = data
    x, y = p

    return ((x - p1[0])**2 + (y - p1[1])**2 - (2*a - r1)**2, (x - p2[0])**2 + (y - p2[1])**2 - (2*a - r2)**2)

def plot_transfer(r2, sma, theta):

    #  Create a figure and axis
    fig, ax = plt.subplots()

    # Calcualte the ellipse parameters (center, width, height, and angle)
    r1 = 1
    r2 = np.linalg.norm(r2)

    # Calculaing eccentricity
    e = (sma - 1) / sma

    # Calculating semi-minor axis b using a point on the transfer ellipse (target body at r2)
    x_j = r2 * cos(theta) + sma*e
    y_j = r2 * sin(theta)
    b = sqrt((y_j)**2 * 1/(1 - (x_j**2 / sma**2)))
    
    # These ellipse parameters characterise the 2 ellipses corresponding to 4 transfer arcs.
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

    # ellipse = patches.Ellipse(ellipse_centre, 2*sma, 2*b, fill=False, color='b', linestyle='--')
    # ax.plot(x_traj, y_traj, color='blue')

    # Solve ellipses equations for the solutions of the 2 vacant focci
    p1 = (1, 0)
    p2 = (r2*cos(theta), r2*sin(theta))
    data = (p1, 1, p2, r2, sma)
    f1 = fsolve(focci_solver, (1, 1), args=data)
    f2 = fsolve(focci_solver, (-1, -1), args=data)

    # Calculating points for both possible transfer ellipses. 1 parameters, 2 points, 3 rotate
    e1 = np.linalg.norm(f1)/(2*sma)
    tilt_ellipse1 = atan2(f1[1], f1[0])
    b1 = sma * sqrt(1 - e1**2)

    e2 = np.linalg.norm(f2)/(2*sma)
    tilt_ellipse2 = atan2(f2[1], f2[0])
    b2 = sma * sqrt(1 - e2**2)

    # Make both angles positive
    if tilt_ellipse1 < 0:
        tilt_ellipse1 += 2*pi   
    if tilt_ellipse2 < 0:
        tilt_ellipse2 += 2*pi

    # Transfer ellipse 1 points
    actual_center1 = np.array(f1/2)
    angle1 = tilt_ellipse1
    theta1 = np.linspace(0, 2*pi, 10000)  # angle from 0 to 2*p
    x1, y1 = np.array([sma * cos(theta) for theta in theta1]), np.array([b1 * sin(theta) for theta in theta1])

    # Calculating rotation and translation matrix to rotate points on ellipse by tilt angle (between focii)
    R = np.array([[cos(angle1), -sin(angle1)], [sin(angle1), cos(angle1)]])
    t = np.array([actual_center1[0], actual_center1[1]])
    x1, y1 = np.matmul(R, np.array([x1, y1])) + t[:, np.newaxis]
    
    # Transfer ellipse 2 points
    actual_center2 = np.array(f2/2)
    angle2 = tilt_ellipse2
    theta2 = np.linspace(0, 2*pi, 10000)  # angle from 0 to 2*p
    x2, y2 = np.array([sma * cos(theta) for theta in theta2]), np.array([b2 * sin(theta) for theta in theta2])
    
    # Calculating rotation and translation matrix to rotate points on ellipse by tilt angle (between focii)
    R = np.array([[cos(angle2), -sin(angle2)], [sin(angle2), cos(angle2)]])
    t = np.array([actual_center2[0], actual_center2[1]])
    x2, y2 = np.matmul(R, np.array([x2, y2])) +  + t[:, np.newaxis]

    # centre and focci positioning
    # ax.plot(ellipse_centre[0], ellipse_centre[1],  marker='*')
    ax.plot(f1[0],f1[1],  marker='*', color='red')
    ax.plot(f2[0],f2[1],  marker='*', color='red')

    # Plot vectors
    ax.plot([0, cos(theta)*r2], [0, sin(theta)*r2], label='My Line', color='red', linestyle='-', marker='')
    ax.plot([0, 1], [0, 0], label='My Line', color='red', linestyle='-', marker='')

    # Add the ellipse to the axis
    ax.plot(x1, y1, color='blue', linestyle='--', marker='')
    ax.plot(x2, y2, color='blue', linestyle='--', marker='')
    ax.add_patch(earth_orbit)
    # ax.add_patch(earth)
    ax.add_patch(target_orbit)
    # ax.add_patch(target)

    # Add actual trajectory (index index points of transfer ellipses)
    # Will never use the retrograde option, so select and choose between prograde options
    # Solution is a bit of a hack, probably a more elegant solution
    
    # euclidean distances for start with first ellipse
    argp1 = np.argmin(np.sqrt((x1 - 1)**2 + (y1 - 0)**2))
    argp2 = np.argmin(np.sqrt((x1 - cos(theta)*r2)**2 + (y1 - sin(theta)*r2)**2))
    traj_ellipse = [x1, y1]

    if y1[argp1] < 0:
        argp1 = np.argmin(np.sqrt((x2 - 1)**2 + (y2 - 0)**2))
        argp2 = np.argmin(np.sqrt((x2 - cos(theta)*r2)**2 + (y2 - sin(theta)*r2)**2))
        traj_ellipse = [x2, y2]

    start, stop = min(argp1, argp2), max(argp1, argp2)
    ax.plot(traj_ellipse[0][start:stop], traj_ellipse[1][start:stop], color='red', linestyle='-', marker='')

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