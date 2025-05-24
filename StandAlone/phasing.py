### Script to perform phasing manoeuvre between two orbits

# Assume both sc are starting in same orbit

# Set time constraint for manoeuvre

import matplotlib.pyplot as plt
from math import tan, atan2, sin, sqrt, pi
from functions import orbital_period, convert_seconds
import numpy as np

# TODO: Write out assumptions

def phasing_man(N):
    
    # Standard grav parameters, m^3 / s^2
    mu_sun   = 1.327*10**20
    mu_earth = 3.986*10**14
    mu_moon  = 4.904*10**12

    r_moon = 1.737*10**6

    # In-plane orbital parameters of target
    # sma_target = r_moon + 200_000
    sma_target = 10200 * 1e3
    e_target = 0.33
    mu_target = mu_earth

    # Phase angle of chaser
    phase_chaser = 90 * pi/180

    # Orbital period
    P_target = orbital_period(sma_target, mu_target)
    print(f'Orbital period of target orbit: {P_target/60:.2f} mins')

    # Time separation
    ecc_anom = 2 * atan2((sqrt(1 - e_target)/sqrt(1 + e_target)), 1/tan(phase_chaser/2) )
    # print(f'{ecc_anom * 180/pi:.2f}')
    
    dt = sqrt(sma_target**3 / mu_target) * (ecc_anom - e_target*sin(ecc_anom))
    # print(f'{dt} seconds or {dt/60:.2f} minutes')
    
    # Phasing orbit parameters - assume catch up in N orbits of target orbit
    P_phasing = P_target - dt/N
    # print(f'P phasing is {P_phasing/60:.2f} mins')
    sma_phasing = (mu_target * (P_phasing/(2*pi))**2)**(1/3)
    print(sma_phasing/1e3)

    # V target at periapse
    r_peri_target = sma_target - sma_target*e_target
    h_target = sqrt(sma_target * mu_target * (1 - e_target**2))
    v_target_peri = h_target / r_peri_target

    # V phasing at periapse
    e_phasing = (sma_phasing - r_peri_target) / sma_phasing
    r_peri_phasing = sma_phasing - sma_phasing*e_phasing
    h_phasing = sqrt(sma_phasing * mu_target * (1 - e_phasing**2))
    v_target_phasing = h_phasing / r_peri_phasing

    # print(v_target_peri/1e3, v_target_phasing/1e3)

    # Dv for phasing manoeuvre
    dv = 2*abs(v_target_phasing - v_target_peri)
    print(f'dV for phasing is {dv /1e3:.4f} km/s')
    d, h, m, s = convert_seconds(N*P_phasing)
    print(f'Total phasing manoeuvre time: {d} days, {h} hs, {m} min, {s} s (or {(N*P_phasing/P_target):.3f} target orbits)')
    
    return [dv, N*P_phasing/P_target]

def main():

    dv_t = []
    N = 10
    for n in range(1, N):
        dv_t.append(phasing_man(n))
    dv_t = np.array(dv_t)

    # Create figure and first axis
    fig, ax1 = plt.subplots()

    # Plot first dataset
    ax1.plot(range(1, N), dv_t[:, 0], 'b-o', label="Dataset 1")
    ax1.set_xlabel("X-axis")
    ax1.set_ylabel("Y1-axis", color='blue')
    ax1.tick_params(axis='y', colors='blue')

    # Create second y-axis (twin)
    ax2 = ax1.twinx()
    ax2.plot(range(1, N), dv_t[:, 1], 'r-s', label="Dataset 2")
    ax2.set_ylabel("Y2-axis", color='red')
    ax2.tick_params(axis='y', colors='red')

    # Show the plot
    plt.grid(True)
    plt.show()
    return 0

if __name__ == '__main__':
    main()
