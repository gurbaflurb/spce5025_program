from keplarianElements import KeplerianElements
from groundSite import GroundSite
import keHelperFunctions

import math
import datetime
import dateutil
import logging
from pprint import pprint

import tabulate
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

def main():

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    vectors_file = 'vectors.yaml'
    vector_data = keHelperFunctions.read_in_yaml(vectors_file)

    ke1 = KeplerianElements(vector_data['vectors'][f'vector1']['x_pos'],
                            vector_data['vectors']['vector1']['y_pos'],
                            vector_data['vectors']['vector1']['z_pos'],
                            vector_data['vectors']['vector1']['x_velocity'],
                            vector_data['vectors']['vector1']['y_velocity'],
                            vector_data['vectors']['vector1']['z_velocity'])


    # Epoch provided: 2012-10-04/050044.523
    epoch = datetime.datetime(2012, 10, 4, 0, 0, 0) # Double check this, the epoch provided isn't the same as what I have in here.

    # Problem 1
    print('----- Problem 1 -----')
    logger.info('Keplarian Elements')
    ke1.print_ke()

    # Compute the following:
    # Compute the Delta-V required to change the inclination by 5 degrees
    # Apply impulsive burn to vector and compute post-burn Keplarian elements to verify results
    # Delta-V required to change SMA by 26km
    # Verify by applying impulsive Delta-V in velocity direction and compute new Keplerian elements
    # Delta-SMA resulting from Delta-V=1m/s
    # Phase rate between given orbit and one with period of 728.0 minutes
    # Time required for phase angle to change by 30 degrees


if __name__ == '__main__':
    main()
