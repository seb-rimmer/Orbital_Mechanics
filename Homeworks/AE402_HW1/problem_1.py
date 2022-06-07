"""
At time 𝑡0 in units for which 𝜇 = 1 (so-called canonical units),
the following data is given for a two-body problem:
𝒓 = [0 2 0] and 𝒗 = [−1 1 0]/√2

    a) Calculate the vectors 𝒉 and 𝒆 and verify that 𝒉 ∙ 𝒆 = 0.
    b) Write the polar equation 𝑟(𝑓) for this conic orbit.
    c) What is the value of 𝑓 at 𝑡0 ?
    d) Determine the speed of the spacecraft when 𝑟 = 32.
    e) Determine the value of 𝑓 when 𝑟 = 32.
    f) For the Earth 𝜇 = 398600 km 3 /s^2 . Determine the value in
       kilometers of the length unit (LU) for which 𝜇 = 1 LU 3 /s 2 .
"""

import numpy as np
import math as m


def main():

    r = [0, 2, 0]
    v = [-1, 1, 0]

    r_arr = np.array(r)
    v_arr = np.array(v)  / 2**0.5

    # h vector is cross product of r and v (r dot)
    h_arr = np.cross(r_arr, v_arr)

    # e vector is cross product of v and h minus normalised r vector
    e_arr = np.cross(v_arr, h_arr) - r_arr/np.linalg.norm(r_arr)

    print(f"H vector is: {h_arr}\n")
    print(f"e vector is: {e_arr}\n")

    print(f"Cross product of e {e_arr} and h {h_arr} is: {np.dot(h_arr, e_arr)}")

    print(f"")

    return 0


if __name__ == '__main__':
    main()