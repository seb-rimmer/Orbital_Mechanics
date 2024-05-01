from math import sqrt, cos, sin, tan, acos, asin, pi
import numpy as np

def lamberts_equation(a, s, c):

    alph = 2*asin(sqrt(s/(2*a)))
    beta = 2*asin(sqrt((s-c)/(2*a)))
    tof = a**(3/2) * (alph - beta - (sin(alph)-sin(beta)))

    return tof

def newton_lambert_solver(a0, s, c, tof):

    # Use Newton's method to return a value for variable a to within precision
    # specified by delta, for the function tof = f(a)

    alph = 2*asin(sqrt(s/(2*a0)))
    beta = 2*asin(sqrt((s-c)/(2*a0)))
    zeta = alph - beta

    # g(a0)
    g_a0 = (3*zeta - sin(alph) + sin(beta)) * (sin(zeta) + sin(alph) - sin(beta) ) - 8*(1-cos(zeta))

    # f'(a0)
    f_d_a0 = (0.5 * sqrt(a0)) * (sin(zeta) + sin(alph) - sin(beta))**-1 * g_a0

    # f(a0)
    f_a0 = a0**(3/2) * (alph - beta - (sin(alph)-sin(beta))) - tof

    a1 = a0 - f_a0/f_d_a0

    return a1, alph, beta

def main():

    tof_TU  = 1.9782787414802256
    s = 2.0578793172534886
    c = 1.5917586345069774
    a0 = s/2 * 1.00001
    a1 = 10*a0    
    # check lamberts eq

    # print(f'Check result for tof is 1.9798 TU - result : {lamberts_equation(a_target, s, c)}')
    # print(f'Check result for tof is 1.5495986122986443 TU - result : {lamberts_equation(1.9671577306446713, s, c)}')
    # print(f'Check result for tof is 1.9798 TU - result : {lamberts_equation(a_target, s, c)}')

    diff = abs(a1-a0)
    while diff > 1e-6:
        # a0 = a1
        a1, alp, beta = newton_lambert_solver(a0, s, c, tof_TU)
        new_tof = lamberts_equation(a1, s, c)

        print(f"a: {a1}, TU: {new_tof}, diff = {abs(a1-a0)}, diff to fsolve result: {abs(a1-1.2321040153545382)}")
        diff = abs(a1-a0)
        
        a0=a1
    

    return 0


if __name__ == '__main__':
    main()