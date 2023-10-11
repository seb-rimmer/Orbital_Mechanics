# Lambert's Problem - Transfer a SC from one point to another in a specified time 
# 
# Step 1 - Compute minimum transfer time
# Step 2 - Compute the principal values of alhpa and beta
# Step 3 - Compute the minimum semi-major axis for the transfer ellipse
# Step 4 - Compute the time corresponding to the minimum semimajor axis transfer ellipse
# Step 5 - Determine values of alpha and beta
# Step 6 - Numerically solve lambert's time equation for semimajor axis
# Step 7 - Compute terminal velocity vectors v1 and v2
# Step 8 - Compute total dV for transfer
# Step 9 - Plot trajectories, start points and end points

import numpy as np

from math import sqrt, cos, sin, tan, acos, asin, pi
from scipy.optimize import fsolve 

def angle_between_vectors(r1, r2):
    
    return 0

def lamberts_time_equation(a, alpha, beta, mu, tf):
    
    return sqrt(a**3 / mu) * (alpha - beta - sin(alpha) +sin(beta)) - tf

# def lamberts_time_equation(a, s, c, mu, tf):

#     s  = 1.8381;        % Pre-calculated space-triangle semi-perimeter
#     c  = 1.6007;        % Pre-calculated chord length
#     tf = (410/365.25) * 2 * pi;     % ToF of 410 days

#     # Updated value for alpha because our t_f is greater than t_m 
#     alpha = 2*pi - 2*asin(sqrt(s/(2*a)));
#     beta  = - 2*asin(sqrt((s-c)/(2*a)));

#     return tf-(a^(3/2))*(alpha-beta-sin(alpha)+sin(beta));

def main():

    # System constants
    TU = (3.137 * 10**7) / (2*pi)
    au = 149_597_871  # km for 1 AU for canonical units
    DU = 1 * au
    mu = 1

    # Starting vector
    # r0 = np.array([0, 0, 0,])
    r0 = 1.0
    u_1 = np.linalg.norm(r0)
    v_dep = 0

    # Final vector
    # rf = np.array([0, 0, 0,])
    rf = 5.2
    u_2 = np.linalg.norm(rf)
    v_arr = 0

    # Time of flight
    tf_days = 524
    tf = (tf_days * 24 * 3600) / TU

    # Parameters (no out of plane transfers)
    # r1 = np.linalg.norm(r0[0:1])
    # r2 = np.linalg.norm(rf[0:1])
    # theta = angle_between_vectors(r1, r2)
    r1 = r0
    r2 = rf    
    theta = 147.0 * pi/180

    # chord
    c = sqrt(r1**2 + r2**2 - 2*r1*r2*cos(theta))
    u_c = (rf - r0) / c

    s = (r1 + r2 + c) / 2

    # Step 1 - Compute minimum transfer time
    t_p = sqrt(2) / (3 * sqrt(mu)) * (s**1.5 - np.sign(sin(theta))*(s - c)**1.5)
    
    # Step 3 - Compute the minimum semi-major axis for the transfer ellipse
    a_m = s/2

    # Step 2 - Compute the principal values of alpha and beta from minimum sma
    alpha_m = 2 * asin(sqrt(s / (2*a_m)))
    beta_m  = 2 * asin(sqrt((s-c) / (2*a_m)))

    # Step 4 - Compute the time corresponding to the minimum semimajor axis transfer ellipse
    t_m = sqrt((s**3) / 8) * (pi - beta_m + sin(beta_m))

    # Step 5 - Determine values of alpha and beta
    if theta <= pi:
        beta = beta_m
    else:
        beta = -beta_m

    if tf <= t_m:
        alpha = alpha_m
    else:
        alpha = 2*pi - alpha_m

    # Step 6 - Numerically solve lambert's time equation for semimajor axis
    #
    # return sqrt(a**3 / mu) * (alpha - beta - sin(alpha) +sin(beta)) - tf
    
    solution = fsolve(lamberts_time_equation, x0=0, args=(alpha, beta, mu, tf))
    
    if solution.success:
        result = solution.x
        print(f'Solution: x = {result[0]}')
    else:
        print('No solution found.')

    # Step 7 - Compute terminal velocity vectors v1 and v2
    
    # Work out A and B constants for velocity vectors
    A = sqrt(1/(4*a)) * 1/tan(alpha * 0.5)
    B = sqrt(1/(4*a)) * 1/tan(beta * 0.5)

    v1 = (B+A) * u_c + (B-A) * u_1      # departure velocity vector
    v2 = (B+A) * u_c - (B-A) * u_2      # arrival velocity vector

    # Step 8 - Compute total dV for transfer
    dv1 = abs(v1 - v_dep)
    dv2 = abs(v_arr - v2)

    # Step 9 - Plot trajectories, start points and end points
    

if __name__ == '__main__':
    main()
