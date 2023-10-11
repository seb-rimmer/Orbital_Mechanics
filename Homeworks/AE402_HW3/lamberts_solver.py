from math import sin, cos, asin, acos, pi, sqrt
from scipy.optimize import fsolve 


def lamberts_time_equation(a):

    s  = 1.8381
    c  = 1.6007
    tf = 9.06

    # % Updated value for alpha because our t_f is greater than t_m

    alpha = 2*pi - 2*asin(sqrt(s/(2*a)))
    beta  = - 2*asin(sqrt((s-c)/(2*a)))

    return tf - (a**(3/2))*(alpha-beta-sin(alpha)+sin(beta));

def main():

    solution = fsolve(lamberts_time_equation, x0=1.1)
    print(solution)
    print("check : {}")

    # if solution.success:
    #     result = solution.x
    #     print(f'Solution: x = {result[0]}')
    # else:
    #     print('No solution found.')

if __name__ == '__main__':
    main()
