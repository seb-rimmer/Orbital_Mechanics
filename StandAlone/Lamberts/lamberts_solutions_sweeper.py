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
import matplotlib.pyplot as plt
import argparse


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

def lamberts_equation(a, s, c):

    alph = 2*pi - 2*asin(sqrt(s/(2*a)))
    beta = 2*asin(sqrt((s-c)/(2*a)))
    tof = a**(3/2) * (alph - beta - (sin(alph)-sin(beta)))

    return tof

def newton_lambert_solver(a0, s, c, tof, branch, biggerThetaTransfer):

    # Use Newton's method to return a value for variable a to within precision
    # specified by delta, for the function tof = f(a)
    alph = 2*asin(sqrt(s/(2*a0)))
    beta = 2*asin(sqrt((s-c)/(2*a0)))
    zeta = alph - beta

    if branch:
        alph = alph
    else:
        alph = 2*pi - alph

    if biggerThetaTransfer:
        beta = -beta

    # g(a0)
    g_a0 = (3*zeta - sin(alph) + sin(beta)) * (sin(zeta) + sin(alph) - sin(beta) ) - 8*(1-cos(zeta))

    # f'(a0)
    f_d_a0 = (0.5 * sqrt(a0)) * (sin(zeta) + sin(alph) - sin(beta))**-1 * g_a0

    # f(a0)
    f_a0 = a0**(3/2) * (alph - beta - (sin(alph)-sin(beta))) - tof

    a1 = a0 - f_a0/f_d_a0

    return a1, alph, beta

def main():

    parser = argparse.ArgumentParser(description="A simple Python script that shows lamberts problem trajectories")
    
    # Add command line arguments
    parser.add_argument('--example', 
                        type = str, 
                        choices = ['v', 'm', 'j'],
                        help = "Which Lambert's example you want to use, out of Mars (m), Jupiter (j) or Venus (v) transfer")
    
    # parser.add_argument('--output', type=str, help="Output file")
    # parser.add_argument('--verbose', action='store_true', help="Enable verbose mode")
    
    args = parser.parse_args()
    
    # Access the argument values
    example = args.example
    # output_file = args.output
    # verbose = args.verbose

    # System constants
    TU = (3.137 * 10**7) / (2*pi)
    AU = 149_597_871  # km for 1 AU for canonical units
    DU = 1 * AU
    mu = 1

    # Some examples to test
    # example var can be j, v, or m
    example = 'm'    
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
        tf_days = 400
    
    elif example == 'm':
        # Mars
        theta_d = 175
        r0 = 1
        rf = 1.524
        tf_days = 200

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

    # Time of flight converted to TU
    tf = (tf_days/365.25) * 2 * pi

    r1 = np.linalg.norm(r0)
    r2 = np.linalg.norm(rf)

    # chord
    c = sqrt(r1**2 + r2**2 - 2*r1*r2*cos(theta))
    u_c = (rf - r0) / c

    # spacetrianlge semiperimeter
    s = (r1 + r2 + c) / 2

    # Step 3 - Compute the minimum semi-major axis for the transfer ellipse
    a_m = s/2

    # Step 1 - Compute minimum transfer time (days)
    tp = 1/3 * sqrt(2/mu) * (s**1.5 - np.sign(sin(theta))*(s - c)**1.5)
    
    if tf_days > tp:
        # elliptic transfer orbit exists.
        print(f'Elliptic transfer ellipse possible; min possible transfer time is {tp:3.2f} TU, chosen transfer time is: {tf_days * 2*pi/365.25:3.2f} TU.')
    else:
        # orbit transfer must be parabolic - we don't do that here
        print('Desired transfer time is less than minimum possible transfer time, transfer not possible.')

    # Step 2 - Compute the principal values of alpha and beta from minimum sma
    alpha0 = 2 * asin( sqrt(  (s)   / (2*a_m) ) ) 
    beta0  = 2 * asin( sqrt(  (s-c) / (2*a_m) ) )

    # Step 4 - Compute the transfer time corresponding to the minimum semimajor axis transfer ellipse
    tm = sqrt(a_m**3 / mu) * (alpha0 - beta0 - (sin(alpha0) - sin(beta0)))

    print(f"Minimum energy ellipse sma = {a_m}, corresponding tof: {tm}")

    # Do sweep and plot for semimajor axis vs tof
    # a_m = a_m * 1.1
    a_sweep = np.linspace(a_m, 1.5, 100)
    a_sweep = np.concatenate((np.linspace(a_m, a_m*1.1, 100), np.linspace(a_m*1.1, a_m*1.5, 10)))

    t_fx_top = np.array([])
    t_fx_bottom = np.array([])

    # Sweep bottom half
    for a in a_sweep:
        alpha0 = 2 * asin( sqrt(  (s)   / (2*a) ) ) 
        beta0  = 2 * asin( sqrt(  (s-c) / (2*a) ) )
        tm_i = sqrt(a**3 / mu) * (alpha0 - beta0 - (sin(alpha0) - sin(beta0)))
        t_fx_bottom = np.append(t_fx_bottom, tm_i)

    # Sweep top half
    for a in a_sweep:
        alpha0 = 2*pi -  2 * asin( sqrt(  (s)   / (2*a) ) ) 
        beta0  = 2 * asin( sqrt(  (s-c) / (2*a) ) )
        tm_i = sqrt(a**3 / mu) * (alpha0 - beta0 - (sin(alpha0) - sin(beta0)))
        t_fx_top = np.append(t_fx_top, tm_i)

    # convert time to years
    # t_fx_top = np.divide(t_fx_top, 2*pi)
    # t_fx_bottom = np.divide(t_fx_bottom, 2*pi)

    #  Create a figure and axis
    fig, ax = plt.subplots()
    plt.plot(a_sweep, t_fx_bottom)
    plt.plot(a_sweep, t_fx_top)    

    # plot additional data for tm, am etc
    plt.plot([a_m, a_m], [0, tm], 'r--', label='Minimum energy sma')
    plt.plot([0, a_m], [tm, tm], 'r--', label='ToF corresponding to minimum energy sma')

    plt.plot([0, 2], [tf, tf], 'r--', label='Chosen ToF')

    ax.set_xlim(a_m*0.9,a_m* 1.5)
    ax.set_ylim(0, 12)

    # Optional: Add labels and title
    plt.xlabel('Semi-major axis (AU)')
    plt.ylabel('Time of flight (TU)')
    plt.title('Plotting semi-major axis sweep')
    plt.legend()

    # Show the plot
    plt.grid()
    # plt.show()

    # Step 5 - Determine initial values of alpha and beta
    biggerThetaTransfer = (theta <= pi)
    if not biggerThetaTransfer:
        beta = beta0
    else:
        beta = -beta0

    # which branch; upper is false, lowwer is true
    branch = (tf <= tm)
    if branch:
        alpha = alpha0
    else:
        alpha = 2*pi - alpha0

    # Step 6 - Numerically solve lambert's time equation for semimajor axis
    #
    # legacy fsolve replaced with newton's method
    # x0 = [a_m*1.2, alpha0, beta0]
    # if tf > tm:
    #     solutions = fsolve(lamberts_time_equation_solver_upper_branch, x0, args=(tf, s, c))
    # else:
    #     solutions = fsolve(lamberts_time_equation_solver_lower_branch, x0, args=(tf, s, c))
    a0 = a_m*1.00001
    a1 = 10*a0
    diff = abs(a1-a0)
    print(a0, s/(a0*2), s, c)
    while diff > 1e-6:
        
        # a0 = a1
        a1, alpha, beta = newton_lambert_solver(a0, s, c, tf, branch, biggerThetaTransfer)
        print(a1, s/(a1*2), s, c, alpha, beta)

        new_tf = lamberts_equation(a1, s, c)
        
        # print(f"a: {a1}, TU: {new_tof}, diff = {abs(a1-a0)}, diff to fsolve result: {abs(a1-1.2321040153545382)}")
        
        diff = abs(a1-a0)
        
        a0=a1

    a = a1
    # print(f'Solution: x = {[a1, alpha, beta]}')

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

    # calculate eccentricity for transfer orbit (from Eq 5.40 in prussing)
    p = sin((alpha+beta)/2)**2 * (4*a *(s - r1)*(s - r2)) / (c**2)
    e = sqrt(1 - p/a) 
    plot_transfer(rf, a, theta)

if __name__ == '__main__':
    main()