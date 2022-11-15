"""
    r_v_from_tle.py created by Seb at 19:31 07/12/2022

    -----------------------------------------------------------------------------
    This file takes a Phoenix TLE from the Celestrak web address below and converts
    it to radius and velocity vectors using the sgp4 python library.

    - Option to manual input a TLE, or scrape the most recent TLE from NORAD
    - additional option to convert directly to orbital elements, using the problem_1 python script

    Faraday Phoenix TLE at: https://celestrak.org/NORAD/elements/gp.php?CATNR=48924
    Python SGP4 documentation at: https://pypi.org/project/sgp4/

"""

from sgp4.api import Satrec
from sgp4.api import jday


def main():

    # FARADAY_PHOENIX TLE, manual
    s = '1 48924U 21059AX  22193.48066053  .00003018  00000+0  16384-3 0  9999'
    t = '2 48924  97.5591 323.3340 0018265  15.0996  49.7083 15.15174404 5774'

    # create satellite object from the TLE
    satellite = Satrec.twoline2rv(s, t)

    # set the julian date, year/month/day/hour/minute/second format
    jd, fr = jday(2022, 7, 12, 20, 2, 0)

    # return error, radial vector (km) and velocity vector (km/s) using sgp4 model
    e, r, v = satellite.sgp4(jd, fr)

    # convert to orbital elements
    print("e : " + str(e) + "\n" + "r : " + str(r) + "\n" + "v : " + str(v))


if __name__ == '__main__':
    main()