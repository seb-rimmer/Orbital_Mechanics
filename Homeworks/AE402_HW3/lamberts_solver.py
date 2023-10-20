from math import sin, cos, asin, acos, pi, sqrt
from scipy.optimize import fsolve 


def lamberts_time_equation_solver_lower_branch(vars, tf, s, c):

    a, alpha, beta = [x for x in vars]
    
    a2 = 2*asin(sqrt(s/(2*a)))
    b2 = 2*asin(sqrt((s-c)/(2*a)))

    return [alpha - a2,
            beta - b2,
            tf - (a**(3/2))*(a2-b2-(sin(a2)-sin(b2)))]

def main():

    tf = (115/365) * 2*pi
    s, c = 2.0578793172534886, 1.5917586345069774

    am     = 1.0289396586267443 * 1.2
    alpha0 = 2 * asin( sqrt(  (s)   / (2*am) ) )  # 3.141592653589793
    beta0  = 2 * asin( sqrt(  (s-c) / (2*am) ) )  # 0.9920327545964835

    solution = fsolve(lamberts_time_equation_solver_lower_branch,
                      [am, alpha0, beta0], 
                      args=(tf, s, c))
    
    print(solution)

if __name__ == '__main__':
    main()
