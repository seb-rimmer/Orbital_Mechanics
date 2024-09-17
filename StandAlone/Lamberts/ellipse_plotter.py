import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import cos, sin, pi, sqrt, ceil, atan2
import numpy as np
from scipy.optimize import fsolve

class transfer_ellipse:
        
        f1, f2 = [1, 0], None

        def __init__(self, theta, r2, sma, name):
            
            self.name = name
            self.r1 = np.array([1, 0])
            self.r2 = r2
            self.p1 = (1, 0)
            self.p2 = (r2*cos(theta), r2*sin(theta))
            self.sma = sma
            self.e = np.linalg.norm(self.f1)/(2*sma)
            self.b = self.sma * sqrt(1 - self.e**2)

            self.calculate_second_focii()
            self.calcualte_transformation_parameters()
            self.calculate_cart_coords()
    
        def calculate_second_focii(self):
            
            points = np.add(self.p1, self.p2)
            if self.name == 'b':
                points = -points

            self.f2 = fsolve(focci_solver, points, args=(self.p1, 1, self.p2, self.r2, self.sma))

        def calcualte_transformation_parameters(self):
            
            self.tilt_ellipse = atan2(self.f2[1], self.f2[0])
            if self.tilt_ellipse < 0:
                self.tilt_ellipse += 2*pi

            self.ellipse_centre = np.array(self.f2/2)
            self.R = np.array([[cos(self.tilt_ellipse), -sin(self.tilt_ellipse)], [sin(self.tilt_ellipse), cos(self.tilt_ellipse)]])
            self.translate = np.array([self.ellipse_centre[0], self.ellipse_centre[1]])

        def calculate_cart_coords(self):
            theta1 = np.linspace(0, 2*pi, 10000)  # angle from 0 to 2*p
            x1, y1 = np.array([self.sma * cos(theta) for theta in theta1]), np.array([self.b * sin(theta) for theta in theta1])
            self.x_cart, self.y_cart = np.matmul(self.R, np.array([x1, y1])) + self.translate[:, np.newaxis]


def focci_solver(p, *data):
    p1, r1, p2, r2, a = data
    x, y = p

    return ((x - p1[0])**2 + (y - p1[1])**2 - (2*a - r1)**2, (x - p2[0])**2 + (y - p2[1])**2 - (2*a - r2)**2)

def theta_positive_from_x(x, y):
    theta = atan2(y, x)
    if theta < 0:
        theta += 2 * pi
    return theta

def plot_transfer(r2, sma, theta, tf, tm):

    #  Create a figure and axis
    fig, ax = plt.subplots()

    # Calcualte the ellipse parameters (center, width, height, and angle)
    r1 = 1
    r2 = np.linalg.norm(r2)

    center = (0, 0)
    earth_pos = (1, 0)
    target_pos = (cos(theta)*r2, sin(theta)*r2)

    # Plot Earth orbit
    earth_orbit = patches.Circle(center, r1, fill=False, color='k', linestyle='--')
    earth = patches.Circle(earth_pos, 0.1, fill=True, color='g')

    # Plot target orbit
    target_orbit = patches.Circle(center, r2, fill=False, color='k', linestyle='--')
    target = patches.Circle(target_pos, 0.1, fill=True, color='orange')
    
    # Declare and define transfer ellipses 

    # Ellipse A - defined as with the second vacant focii in the > pi part of r1 and r2 vectors
    transfer_ellipse_a = transfer_ellipse(theta, r2, sma, 'a')

    # Ellipse B - defined as with the second vacant focii in the > pi part of r1 and r2 vectors
    transfer_ellipse_b = transfer_ellipse(theta, r2, sma, 'b')
    
    # Plot vectors
    ax.plot([0, cos(theta)*r2], [0, sin(theta)*r2], label='My Line', color='red', linestyle='-', marker='')
    ax.plot([0, 1], [0, 0], label='My Line', color='red', linestyle='-', marker='')

    # Add the ellipse to the axis
    ax.plot(transfer_ellipse_a.x_cart, transfer_ellipse_a.y_cart, color='blue', linestyle='--', marker='')
    ax.plot(transfer_ellipse_b.x_cart, transfer_ellipse_b.y_cart, color='blue', linestyle='--', marker='')
    ax.plot(transfer_ellipse_a.f2[0],transfer_ellipse_a.f2[1],  marker='*', color='red')
    ax.plot(transfer_ellipse_b.f2[0],transfer_ellipse_b.f2[1],  marker='*', color='red')

    ax.add_patch(earth_orbit)
    # ax.add_patch(earth)
    ax.add_patch(target_orbit)
    # ax.add_patch(target)

    #TODO: Add how to choose the transfer arcs to plot in red, corresponding to tf

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