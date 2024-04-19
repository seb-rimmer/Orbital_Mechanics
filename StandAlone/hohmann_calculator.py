#!/usr/bin/env python
#
# Description: Calcualte dV for Hohmann transfers from Earth to other places in the solar system
#              using the Vis-Viva formula. To be used as rough values.
# 
from math import pi, sqrt

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

    mu_sun = 1.327e11
    AU = 149_597_871  # km for 1 AU for canonical units
    TU = sqrt(AU**3 / mu_sun)

    r1 = 1
    r2 = 1.64
    mu = 1
    result = hohmann(r1, r2, mu)

    if mu == 1:
    
        dv_canon = result[0]
        dv_kms = dv_canon * AU / TU

        dT_canon = result[1]
        dT_days = dT_canon * TU / (24 * 3600)

        print(f"Delta-v required to get from Earth to orbit at {r2} AU on a Hohmann transfer is: \n {dv_canon:.4f} AU/TU \n {dv_kms:.4f} km/s ")
        print(f"Transfer will take: {dT_canon:.4f} TU , or {dT_days:.4f} days, {dT_days/365:.2f} years.")

    else:
        dv_kms = result[0]

        dT_s = result[1]
        dT_hrs = dT_s / 3600

        print(f"Delta-v required to get from orbit at {r1} km to orbit at {r2}km on a Hohmann transfer is: \n {dv_kms:.4f} km/s ")
        print(f"Transfer will take: {dT_hrs:.4f} hours")
    

if __name__ == '__main__':
    main()