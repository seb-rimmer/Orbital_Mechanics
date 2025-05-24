#!/usr/bin/env python
#
# Description: Calcualte dV for Hohmann transfers from Earth to other places in the solar system
#              using the Vis-Viva formula. To be used as rough values.
# 
from math import pi, sqrt

def orbital_period(a, mu):

    return 2*pi * sqrt(a**3 / mu)

def hohmann(r1, r2, mu):

    # SMA of Hohmann transfer is from simple geometry
    aH = (r1 + r2)/2

    # Transfer time is half the period of the transfer ellipse
    tH = pi * sqrt(aH**3 / mu)

    eH = (r2-r1) / (r1+r2) 

    # dVA - What we need minus what we have
    v1 = sqrt(mu/r1)
    vH_depart = sqrt(mu * (2/r1 - 1/aH))
    dVA = vH_depart - v1

    # dVB - What we need minus what we have
    v2 = sqrt(mu/r2)
    vH_arrive = sqrt(mu * (2/r2 - 1/aH))
    dVB = v2 - vH_arrive

    dvTotal = dVA + dVB

    return dvTotal, tH

def main():

    # System constants

    # Standard grav parameters, m^3 / s^2
    mu_sun   = 1.327*10**20
    mu_earth = 3.986*10**14
    mu_moon  = 4.904*10**12

    AU = 149_597_871  # km for 1 AU for canonical units
    TU = sqrt(AU**3 / mu_sun)

    r1 = 1.737*10**6 + 200_000
    r2 = r1 + 10_000
    mu = mu_moon
    result = hohmann(r1, r2, mu)
    

    # orbital period 
    orbital_period_s_r1 = orbital_period(r1, mu_moon)
    orbital_period_s_r2 = orbital_period(r2, mu_moon)

    print(f'Orbital period at r1: {orbital_period_s_r1//3600:.0f} hrs {(orbital_period_s_r1 - (orbital_period_s_r1//3600)*3600) /60:.1f} mins')
    print(f'Orbital period at r1: {orbital_period_s_r2//3600:.0f} hrs {(orbital_period_s_r2 - (orbital_period_s_r2//3600)*3600) /60:.1f} mins')

    dv_ms = result[0]
    dT_s = result[1]
    dT_hrs = dT_s / 3600

    print(f"Delta-v required to get from orbit at {r1/1000} km to orbit at {r2/1000}km on a Hohmann transfer is: {dv_ms:.4f} m/s ")
    print(f"Transfer will take: {dT_hrs:.4f} hours")


    # if mu == 1:
    
    #     dv_canon = result[0]
    #     dv_kms = dv_canon * AU / TU

    #     dT_canon = result[1]
    #     dT_days = dT_canon * TU / (24 * 3600)

    #     print(f"Delta-v required to get from Earth to orbit at {r2} AU on a Hohmann transfer is: \n {dv_canon:.4f} AU/TU \n {dv_kms:.4f} km/s ")
    #     print(f"Transfer will take: {dT_canon:.4f} TU , or {dT_days:.4f} days, {dT_days/365:.2f} years.")

    # else:
    #     dv_kms = result[0]

    #     dT_s = result[1]
    #     dT_hrs = dT_s / 36006

    #     print(f"Delta-v required to get from orbit at {r1} km to orbit at {r2}km on a Hohmann transfer is: \n {dv_kms:.4f} km/s ")
    #     print(f"Transfer will take: {dT_hrs:.4f} hours")
    

if __name__ == '__main__':
    main()