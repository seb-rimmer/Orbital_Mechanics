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
# 
# problem_1_lamberts_refactored.py
# 

import numpy as np
import matplotlib

from math import sqrt, cos, sin, tan, acos, asin, pi
from scipy.optimize import fsolve 

# my imports
from ellipse_plotter import plot_transfer

def angle_between_vectors(r1, r2):
    
    return 0

def lamberts_time_equation_solver_lower_branch(vars, tf, s, c):

    a, alpha, beta = [x for x in vars]
    
    a2 = 2*asin(sqrt(s/(2*a)))
    b2 = 2*asin(sqrt((s-c)/(2*a)))

    return [alpha - a2,
            beta - b2,
            tf - (a**(3/2))*(a2-b2-(sin(a2)-sin(b2)))]

def lamberts_time_equation_solver_upper_branch(vars, tf, s, c):

    a, alpha, beta = [x for x in vars]
    
    a2 = 2*pi - 2*asin(sqrt(s/(2*a)))
    b2 = 2*asin(sqrt((s-c)/(2*a)))

    return [alpha - a2,
            beta - b2,
            tf - (a**(3/2))*(a2-b2-(sin(a2)-sin(b2)))]

def lamberts_time_equation_tf_TU(a, alpha, beta):

    return (a**(3/2))*(alpha-beta-sin(alpha)+sin(beta))

def main():

    # System constants
    TU = (3.137 * 10**7) / (2*pi)
    AU = 149_597_871  # km for 1 AU for canonical units
    DU = 1 * AU
    mu = 1

    # Some examples to test
    # example var can be j, v, or m

    example ='j'
    
    if example == 'j':
        # Jupiter
        theta_d = 147
        r0 = 1
        rf = 5.2
        tf_days = 524
    
    elif example == 'v':
        # Venus
        theta_d = 135
        r0 = 1
        rf = 0.723
        tf_days = 337
    
    elif example == 'm':
        # Mars
        theta_d = 75
        r0 = 1
        rf = 1.524
        tf_days = 115

    # Angle between vectors
    theta = theta_d * pi/180

    # Starting vector
    r0 = np.array([r0, 0])
    v_dep = np.array([0, 1])

    # Final vector
    rf = np.array([-rf*cos(pi - theta), rf*sin(pi-theta)])
    u_2 = cos(theta)
    v_arr = np.sqrt(1/np.linalg.norm(rf)) * np.array([-sin(theta), cos(theta)])

    # Unit vectors
    u_1 = r0 / np.linalg.norm(r0)
    u_2 = rf / np.linalg.norm(rf)

    # Time of flight
    tf = (tf_days/365.25) * 2 * pi

    r1 = np.linalg.norm(r0)
    r2 = np.linalg.norm(rf)

    # chord
    c = sqrt(r1**2 + r2**2 - 2*r1*r2*cos(theta))
    u_c = (rf - r0) / c

    s = (r1 + r2 + c) / 2

    # Step 1 - Compute minimum transfer time
    t_p = sqrt(2) / (3 * sqrt(mu)) * (s**1.5 - np.sign(sin(theta))*(s - c)**1.5)
    
    # Step 3 - Compute the minimum semi-major axis for the transfer ellipse
    a_m = s/2

    # Step 2 - Compute the principal values of alpha and beta from minimum sma
    alpha0 = 2 * asin( sqrt(  (s)   / (2*a_m) ) ) 
    beta0  = 2 * asin( sqrt(  (s-c) / (2*a_m) ) )

    # Step 4 - Compute the time corresponding to the minimum semimajor axis transfer ellipse
    tm = sqrt((s**3) / 8) * (pi - beta0 + sin(beta0))

    # Step 5 - Determine initial values of alpha and beta
    if theta <= pi:
        beta = beta0
    else:
        beta = -beta0

    if tf <= tm:
        alpha = alpha0
    else:
        alpha = 2*pi - alpha0

    # Step 6 - Numerically solve lambert's time equation for semimajor axis
    #
    # lamberts_time_equation(a, s, c, tf)
    # initial guess of a needs to be a bit further away from a_m so we don't hit a 
    # singularity with asin() calcualtions
    x0 = [a_m*1.2, alpha0, beta0]

    if tf > tm:
        solutions = fsolve(lamberts_time_equation_solver_upper_branch, x0, args=(tf, s, c))
    else:
        solutions = fsolve(lamberts_time_equation_solver_lower_branch, x0, args=(tf, s, c))
    
    a, alpha, beta = [x for x in solutions]
    
    print(f'Solution: x = {solutions}')

    # Step 7 - Compute terminal velocity vectors v1 and v2
    
    # Work out A and B constants for velocity vectors
    A = sqrt(1/(4*a)) * 1/tan(alpha/2)
    B = sqrt(1/(4*a)) * 1/tan(beta/2)

    v1 = (B+A) * u_c + (B-A) * u_1      # departure velocity vector
    v2 = (B+A) * u_c - (B-A) * u_2      # arrival velocity vector

    # Step 8 - Compute total dV for transfer
    dv1 = np.linalg.norm(v1 - v_dep)
    dv2 = np.linalg.norm(v_arr - v2)
    
    print(f"Departure dV = {dv1} DU/TU  /  dV = {dv1 * DU/TU} km/s")
    print(f"Arrival dV   = {dv2} DU/TU  /  dV = {dv2 * DU/TU} km/s")
    print(f"Total dV     = {dv1 + dv2} DU/TU  /  dV = {(dv1+dv2) * DU/TU} km/s")

    # Step 9 - Plot trajectories, start points and end points
    plot_transfer(rf, a, theta)

if __name__ == '__main__':
    main()
