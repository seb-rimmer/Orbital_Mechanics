import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import cos, sin, pi, sqrt, ceil, atan2, acos
import numpy as np
from scipy.optimize import fsolve

class transfer_ellipse:
        
        def __init__(self, theta, r2, sma, f_star):
            
            self.r1 = np.array([1, 0])
            self.r2 = r2
            self.theta = theta
            self.p1 = (1, 0)
            self.p2 = (r2*cos(self.theta), r2*sin(self.theta))
            self.sma = sma
            self.f1 = np.array([0, 0])
            self.f2 = f_star
            self.e = np.linalg.norm(self.f2 - self.f1)/(2*sma)
            self.b = self.sma * sqrt(1 - self.e**2)

            self.calcualte_transformation_parameters()
            self.calculate_cart_coords()

        def calcualte_transformation_parameters(self):
            
            self.tilt_ellipse = atan2(self.f2[1], self.f2[0])

            if self.tilt_ellipse < 0 and self.tilt_ellipse > -pi:
                self.tilt_ellipse += pi 

            self.ellipse_centre = np.array(self.f2/2)
            self.R = np.array([
                                [cos(self.tilt_ellipse), -sin(self.tilt_ellipse)], 
                                [sin(self.tilt_ellipse),  cos(self.tilt_ellipse)]
                              ])
            self.translate = np.array([self.ellipse_centre[0], self.ellipse_centre[1]])

        def calculate_cart_coords(self):
            theta1 = np.linspace(0, 2*pi, 10000)  # angle from 0 to 2*p
            x1, y1 = np.array([self.sma * cos(theta) for theta in theta1]), np.array([self.b * sin(theta) for theta in theta1])
            self.x_cart_trans, self.y_cart_trans = np.matmul(self.R, np.array([x1, y1])) + self.translate[:, np.newaxis]
            self.x_cart, self.y_cart = x1, y1

        def transfer_arc_vals(self):
            
            a, b, e = self.sma, self.b, self.e
            h, k = self.ellipse_centre

            vect1 = np.array([1, 0])
            vect2 = np.array(-self.f2) / np.linalg.norm(self.f2)
            start_f = acos(np.dot(vect1, vect2))

            # Determine true anomaly range
            
            f = np.linspace(start_f, start_f+self.theta, 100)
            
            r = (a * (1 - e**2)) / (1 + e*np.cos(f))

            x_ellipse_plane, y_ellipse_plane = np.cos(f)*r, np.sin(f)*r

            self.R = np.array([
                                [cos(self.tilt_ellipse-pi), -sin(self.tilt_ellipse-pi)], 
                                [sin(self.tilt_ellipse-pi),  cos(self.tilt_ellipse-pi)]
                              ])
            x_cart_trans, y_cart_trans = np.matmul(self.R, np.array([x_ellipse_plane, y_ellipse_plane]))

            return x_cart_trans, y_cart_trans
            # return x_ellipse_plane, y_ellipse_plane
        
def focci_solver(p, *data):
    p1, r1, p2, r2, a = data
    x, y = p

    return ((x - p1[0])**2 + (y - p1[1])**2 - (2*a - r1)**2, (x - p2[0])**2 + (y - p2[1])**2 - (2*a - r2)**2)

def theta_positive_from_x(x, y):
    theta = atan2(y, x)
    if theta < 0:
        theta += 2 * pi
    return theta

def find_circle_intersection(x1, y1, r1, x2, y2, r2):
    
    # Distance between the centers
    d = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    l = (d**2 - r2**2 + r1**2) / (2*d)
    h = sqrt(r1**2 - l**2)
    
    f1_int1 = np.array([[l], [+h]])
    f2_int1 = np.array([[l], [-h]])

    R_theta = atan2((y2- y1), (x2 - x1))
    R= np.array([ [np.cos(R_theta), -np.sin(R_theta)],
                  [np.sin(R_theta), np.cos(R_theta)]
                ])
    T = np.array([1, 0])

    f1 = np.matmul(R, f1_int1) + T[:, np.newaxis]
    f2 = np.matmul(R, f2_int1) + T[:, np.newaxis]

    return f1.flatten(), f2.flatten()

def plot_transfer(r2, sma, theta, tf, tm):

    # Calcualte the ellipse parameters (center, width, height, and angle)
    r1 = 1
    r2 = np.linalg.norm(r2)

    #  Create a figure and axis
    plt.ion()
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    plt.grid()
    lim = max(1, r2)*1.2
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    plt.draw()

    center = (0, 0)
    earth_pos = (1, 0)
    target_pos = (cos(theta)*r2, sin(theta)*r2)

    # Plot Earth orbit
    earth_orbit = patches.Circle(center, r1, fill=False, color='k', linestyle='--')

    # Plot target orbit
    target_orbit = patches.Circle(center, r2, fill=False, color='k', linestyle='--')
    
    ax.add_patch(earth_orbit)
    ax.add_patch(target_orbit)
    plt.draw()

    # Plot bodies position vectors
    ax.plot([0, cos(theta)*r2], [0, sin(theta)*r2], label='My Line', color='red', linestyle='-', marker='')
    ax.plot([0, 1], [0, 0], label='My Line', color='red', linestyle='-', marker='')
    # 
    # Plot ellipse property 2*a -r circles
    # ax.add_patch(patches.Circle(earth_pos, 2*sma-r1, fill=False, color='gray'))
    # ax.add_patch(patches.Circle(target_pos, 2*sma-r2, fill=False, color='gray'))

    # Calculate two vacant focii points based on intersection of above circles
    f1, f2 = find_circle_intersection(earth_pos[0], earth_pos[1], 2*sma-r1, target_pos[0], target_pos[1], 2*sma-r2)
    ax.plot(f1[0], f1[1], '*', color='r')
    # ax.plot(f2[0], f2[1], '*', color='m')

    # Declare and define transfer ellipses 

    # Ellipse A - defined as with the second vacant focii in the > pi part of r1 and r2 vectors
    transfer_ellipse_a = transfer_ellipse(theta, r2, sma, f1)

    ax.plot(transfer_ellipse_a.x_cart_trans, transfer_ellipse_a.y_cart_trans, color='b', linestyle='--', label='Transfer ellipse A')
    plt.draw()

    # Ellipse B - defined as with the second vacant focii in the > pi part of r1 and r2 vectors
    transfer_ellipse_b = transfer_ellipse(theta, r2, sma, f2)

    # ax.plot(transfer_ellipse_b.x_cart_trans, transfer_ellipse_b.y_cart_trans, color='m', linestyle='--', label='Transfer ellipse B')
    plt.draw()

    #TODO: Add how to choose the transfer arcs to plot in red, corresponding to tf
    # The two retrograde arcs correspond to:
    #       - Long arc of Ellipse A (where second focii is in <pi region between r1 and r2)
    #       - Short arc of Ellipse B (focii in >pi region between r1 r2 vectors)
    
    if tf > tm: # long arcs

        # Always going to choose prograde option 
        if transfer_ellipse_a.long_arc_vals[0][0] > 0:
            trans_ellipse = transfer_ellipse_a
        else:
            trans_ellipse = transfer_ellipse_b
    
    else:    # short arcs

        # Always going to choose prograde option
        angle = atan2(transfer_ellipse_a.f2[1], transfer_ellipse_a.f2[0])
        if angle < 0:
            angle =  2*pi + angle
        if angle > theta:
            trans_ellipse = transfer_ellipse_a
        else:
            trans_ellipse = transfer_ellipse_b

    x_transfer, y_transfer = trans_ellipse.transfer_arc_vals()

    # Plot long arc of ellipse a
    ax.plot(x_transfer, y_transfer, 'r')

    # Optional: Add labels and title
    ax.legend()
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plotting an Ellipse')

    # Show the plot
    plt.ioff()
    plt.show()

def main():
    # plot_transfer(5.2, 5.32, 147)

    return 0

if __name__ == '__main__':
    main()