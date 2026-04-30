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

    ke1 = KeplerianElements(vector_data['vectors']['vector1']['x_pos'],
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

    print()
    # Compute the following:
    # Compute the Delta-V required to change the inclination by 5 degrees
    p1_delta_v = keHelperFunctions.estimate_plane_change_maneuver(ke1, 5)
    logger.info(f'Delta-V required to change inclination by 5 degrees: {p1_delta_v}')
    print()

    # Apply impulsive burn to vector and compute post-burn Keplarian elements to verify results
    logger.info('Computing verification of Delta-V...')
    p1_ke2, p1_delta_v_check = keHelperFunctions.estimate_plane_change_burn_ke(ke1, 5)
    logger.info(f'Delta-V check for changing inclination by 5 degrees: {p1_delta_v_check}')
    logger.info(f'Delta-V check for changing inclination by 5 degrees (Magnitude): {np.linalg.norm(p1_delta_v_check)}')
    logger.info(f'New Keplarian Elements post 5 degree burn')
    p1_ke2.print_ke()
    print()
    logger.info(f'SMA Diff: {ke1.semi_major_axis - p1_ke2.semi_major_axis}')
    logger.info(f'Eccentricity Diff: {ke1.eccentricity - p1_ke2.eccentricity}')
    logger.info(f'Inclination Diff: {math.degrees(ke1.inclination - p1_ke2.inclination)}')
    logger.info(f'RAAN Diff: {math.degrees(ke1.raan - p1_ke2.raan)}')
    logger.info(f'AoP Diff: {math.degrees(ke1.aop - p1_ke2.aop)}')
    logger.info(f'True Anomaly Diff: {math.degrees(ke1.nu - p1_ke2.nu)}')
    logger.info(f'TP Diff: {ke1.tp - p1_ke2.tp}')
    logger.info(f'Apogee Radii Diff: {ke1.apogee_radii - p1_ke2.apogee_radii}')
    logger.info(f'Perigee Radii Diff: {ke1.perigee_radii - p1_ke2.perigee_radii}')

    


    # Delta-V required to change SMA by 26km
    # Verify by applying impulsive Delta-V in velocity direction and compute new Keplerian elements


    # Delta-SMA resulting from Delta-V=1m/s


    # Phase rate between given orbit and one with period of 728.0 minutes


    # Time required for phase angle to change by 30 degrees


if __name__ == '__main__':
    main()
