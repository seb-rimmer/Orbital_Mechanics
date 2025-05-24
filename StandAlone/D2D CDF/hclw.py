## Script that contains HCLW dynamics
import functions
import math
from math import sin, cos
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

def triad(ax):
    
    # Define the origin and direction of each axis
    origin = np.array([0, 0, 0])  # Origin at (0, 0, 0)

    # Define the axis vectors (unit vectors in X, Y, Z directions)
    x_axis = np.array([1, 0, 0])  # X-axis
    y_axis = np.array([0, 1, 0])  # Y-axis
    z_axis = np.array([0, 0, 1])  # Z-axis

    # Plot the arrows for the X, Y, and Z axes
    ax.quiver(*origin, *x_axis, color='r', length=100, arrow_length_ratio=0.1)
    ax.quiver(*origin, *y_axis, color='g', length=100, arrow_length_ratio=0.1)
    ax.quiver(*origin, *z_axis, color='b', length=100, arrow_length_ratio=0.1)

    return 0

def starting_conditions(r_pco, n, phi):

    x0 = r_pco*cos(phi)
    y0 = -2*r_pco*sin(phi)
    z0 = 2*r_pco*cos(phi)
    x0dot = -r_pco*n*sin(phi)
    y0dot = -2*r_pco*n*cos(phi)
    z0dot = -2*r_pco*n*sin(phi)

    return [x0, y0, z0, x0dot, y0dot, z0dot]

mu_earth = 3.986*10**14
r_earth  = 6378*1e3

# Assumed altitude for Mosaic formation from slides
altitude = 8000*1e3
sma = r_earth+altitude
e = 0
n = math.sqrt(mu_earth/(sma**3)) 

P_s = functions.orbital_period(sma, mu_earth)

# Define HCLW System A matrix
A = np.array([[0,      0,     0,    1,   0,   0], 
                [0,      0,     0,    0,   1,   0], 
                [0,      0,     0,    0,   0,   1],
                [3*n**2, 0,     0,    0, 2*n,   0], 
                [0,      0,     0, -2*n,   0,   0], 
                [0,      0, -n**2,    0,   0,   0]])

# Propagate A matrix
dt = 10
num_steps = 500
A_d = expm(A*dt)

# Initialize state history
r_pco = 100
# phase = math.pi/4
phase = 0
[x0, y0, z0, x0dot, y0dot, z0dot] = starting_conditions(r_pco, n, phase)

n = A.shape[0]
x_hist = np.zeros((num_steps + 1, n))
x_hist[0] = [x0, y0, z0, x0dot, y0dot, z0dot]  # Set initial state

# Propagate state through time steps
for k in range(num_steps):
    x_hist[k+1] = A_d @ x_hist[k]

# Create figure and first axis
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')

# Plot first dataset
x, y, z = x_hist[:, 0], x_hist[:, 1], x_hist[:, 2]

ax1.plot(x, y, z, 'b-o', label="Dataset 1")
ax1.plot(0, 0, 'r-o', label="Dataset 2")
ax1.set_xlabel("Z-axis")
ax1.set_ylabel("Y-axis")
ax1.set_zlabel("X-axis")

# plot triad of HCLW frame
triad(ax1)

ax1.axis('equal')

# States plot
fig2, ax2 = plt.subplots(6, 1, figsize=(6, 10))  # 3 rows, 1 column
t = range(0, num_steps+1, 1)
ax2[0].plot(t, x_hist[:, 0], color='r')
ax2[0].set_title('X state (radial)')
ax2[0].grid(True)
ax2[1].plot(t, x_hist[:, 1], color='r')
ax2[1].set_title('Y state (along-track)')
ax2[1].grid(True)
ax2[2].plot(t, x_hist[:, 2], color='r')
ax2[2].set_title('Z state (cross-track)')
ax2[2].grid(True)
ax2[3].plot(t, x_hist[:, 3], color='r')
ax2[3].set_title('X dot state (radial)')
ax2[3].grid(True)
ax2[4].plot(t, x_hist[:, 4], color='r')
ax2[4].set_title('Y dot state (along-track)')
ax2[4].grid(True)
ax2[5].plot(t, x_hist[:, 5], color='r')
ax2[5].set_title('Z dot state (cross-track)')
ax2[5].grid(True)
ax2[5].set_xlabel("Time")


# Plot nadir-pointing
fig3, ax3 = plt.subplots()  # 3 rows, 1 column
ax3.plot(x_hist[:, 1], x_hist[:, 2], color='r')
ax3.set_title('Y-Z nadir plane state (radial)')
ax3.set_xlabel("Y-axis (along-track)")
ax3.set_ylabel("Z-axis (cross-track)")
ax3.grid(True)
ax3.set_ylim([-200, 200])
ax3.set_xlim([-200, 200])

fig4, ax4 = plt.subplots()  # 3 rows, 1 column
ax4.plot(x_hist[:, 2], x_hist[:, 0], color='r')
ax4.set_title('X-Z nadir plane state (radial)')
ax4.set_ylabel("X-axis (radial)")
ax4.set_xlabel("Z-axis (cross-track)")
ax4.grid(True)
ax4.set_ylim([-200, 200])
ax4.set_xlim([-200, 200])
ax4.axis('equal')


plt.show()