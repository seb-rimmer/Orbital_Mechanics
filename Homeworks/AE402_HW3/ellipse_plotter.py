import matplotlib.pyplot as plt
import matplotlib.patches as patches
from math import cos, sin, pi, sqrt

#  Create a figure and axis
fig, ax = plt.subplots()

# Define the ellipse parameters (center, width, height, and angle)
theta_d = 147
theta_r = theta_d * pi/180
re = 1
rj = 5.2
sma = 5.32

a, b, c = 1, rj*cos(theta_r)/sma, (rj - sma)/sma
dis = b**2 - 4*a*c
e1, e2 = (-b + sqrt(dis)) / (2*a), (-b - sqrt(dis)) / (2*a)
e = max(e1, e2)

center = (0, 0)
earth_pos = (1, 0)
jupiter_pos = (cos(theta_r)*rj, sin(theta_r)*rj)
ellipse_centre = (-sma*e+1, 0)

# Plot Earth orbit
earth_orbit = patches.Circle(center, re, fill=False, color='k', linestyle='--')
earth = patches.Circle(earth_pos, 0.2, fill=True, color='g')

# Plot Jupiter orbit
jupiter_orbit = patches.Circle(center, rj, fill=False, color='k', linestyle='--')
jupiter = patches.Circle(jupiter_pos, 0.5, fill=True, color='orange')

# Create an Ellipse patch
ellipse = patches.Ellipse(ellipse_centre, 2*sma, sma * sqrt(1 - e**2), fill=False, color='b')


# Plot vectors
ax.plot([0, cos(theta_r)*rj], [0, sin(theta_r)*rj], label='My Line', color='red', linestyle='-', marker='')
ax.plot([0, 1], [0, 0], label='My Line', color='red', linestyle='-', marker='')

# Add the ellipse to the axis
ax.add_patch(ellipse)
ax.add_patch(earth_orbit)
ax.add_patch(earth)
ax.add_patch(jupiter_orbit)
ax.add_patch(jupiter)

# Set the aspect ratio of the plot to be equal
ax.set_aspect('equal')

# Set the axis limits to show the entire ellipse
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)

# Optional: Add labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Plotting an Ellipse')

# Show the plot
plt.grid()
plt.show()