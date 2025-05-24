import numpy as np

a, b = 5.328, 3.106
theta = -3.077
h, k = -4.32, -0.2764

x = np.linspace(-4.36, 1, 100)

xtest = np.array([1, 0, -4.36])
x= xtest

# Constants for equation to give y coordinates of rotated and translated ellipse
A = (np.sin(theta)**2)/(a**2) + (np.cos(theta)**2)/(b**2)
B = (2*(x-h)*np.cos(theta)*np.sin(theta))/(a**2) - (2*(x-h)*np.cos(theta)*np.sin(theta))/(b**2)
C = (((x-h)*np.cos(theta))**2)/(a**2) + (((x-h)*np.sin(theta))**2)/(b**2)

# determinant check > 0
det = B**2 - 4*A*(C-1)

y = k + (-B + np.sqrt(det))/2*A

print(y)