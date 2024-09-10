#!/usr/bin/env python
#
# Description: 
# 

import numpy as np
import requests
import json
import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

from math import sqrt , atan2, pi, sin, cos, tan

import functions

def main():

    launch =  '2024-Oct-12 00:00:00'
    arrival = '2026-Dec-01 00:00:00'

    bodies = [3, 4, 20065803]
    colors = ['blue', 'orange', 'grey']
    body_positions = []

    for body in bodies:
        
        # Add body at launch
        body_dict = functions.jpl_body_request(body, launch)
        if body_dict != 0:
            body_positions += [body_dict]

        # Add body at arrival
        body_dict = functions.jpl_body_request(body, arrival)
        if body_dict != 0:
            body_positions += [body_dict]

    functions.plot_vectors(body_positions)
    



if __name__ == '__main__':
    main()
    