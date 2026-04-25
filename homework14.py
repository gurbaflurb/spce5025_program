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


    # Epoch provided: 06/01/2012 00:00:00
    epoch = datetime.datetime(2024, 6, 1, 0, 0, 0)


    # Problem 1
    logger.info('----- Problem 1 -----')
    logger.info('Keplarian Elements')
    ke1.print_ke()

    # Using Time of Flight approach, estimate the times of the next five
    # 
    # 1. Node Crossings (Use J2000 vectors)
    # 2. Perigee Passages
    # 3. Eclipse centers (orbit midnight)
    #   - For each eclipse, estimate the duration of the eclipse
    #   - You may need to compute a new sun vector for each estimate
    # Verify results against the tabulated data
    #
    # Output to a CSV in the following format:
    #
    # Sun Percentage values
    # Pct = 1      - Full Illumination
    # 0 < Pct < 1  - Partially illuminated
    # Pct = 0      - Fully Eclipsed
    #
    # DATE | X | Y | Z | Xecef | Yecef | Zecef | A (SMA) | E (Eccentricity) | I (Inclination) | RAAN | ARGP | TA | MA | SVE | Orbit Beta Angle | SunPct | Eclipse (Denoted by *)
    #

    print()
    print(f'Next Node Passes')
    node_crossings = keHelperFunctions.estimate_node_crossing_times(ke1, ke1.aop, 5)

    for i in node_crossings:
        pass_time = epoch + datetime.timedelta(seconds=i)
        logger.info(f'{pass_time.date()} {pass_time.time()}')

    exit(0)

    logger.info('Next Perigee Passes')
    

    time_to_apogee = keHelperFunctions.compute_time_to_apogee(ke1)
    time_to_perigee = keHelperFunctions.compute_time_to_perigee(ke1)
    logger.info(f'Time to Apogee: {time_to_apogee/60} minutes')
    logger.info(f'Time to Perigee: {time_to_perigee/60} minutes')

    ttp = ke1.determine_time_to_angle(2*math.pi)
    tta = ke1.determine_time_to_angle(math.pi)
    logger.info(f'Time to 2pi: {ttp}')
    logger.info(f'Time to pi: {tta}')
    

    p1_jd = keHelperFunctions.convert_date_to_jd(epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, epoch.second)

    sun_vector = keHelperFunctions.determine_sun_vector_lf(p1_jd)

    noon, midnight = keHelperFunctions.compute_noon_and_midnight(sun_vector, ke1.raan, ke1.inclination, ke1.aop)

    logger.info(f'Noon Vector: {noon}')
    logger.info(f'Midnight Vector: {midnight}')


if __name__ == '__main__':
    main()
