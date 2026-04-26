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
    epoch = datetime.datetime(2012, 6, 1, 0, 0, 0)

    # Problem 1
    print('----- Problem 1 -----')
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

    print('\n')
    print(f'----- Part 1. Next Ascending Node Passes -----')
    node_crossings = keHelperFunctions.estimate_node_crossing_times(ke1, ke1.aop, 5)

    j = 1
    for i in node_crossings:
        pass_time = epoch + datetime.timedelta(seconds=i)
        logger.info(f'{j} - {pass_time.date()} {pass_time.time()}')
        j += 1

    print('\n')
    print('----- Part 2. Next Perigee Passes -----')

    time_to_perigee = keHelperFunctions.compute_time_to_perigee(ke1)
    logger.info(f'Time to Perigee: {time_to_perigee} seconds')

    logger.info('Predicted Perigee crossing times')
    pass_time = epoch + datetime.timedelta(seconds=time_to_perigee)
    logger.info(f'1 - {pass_time}')

    for i in range(1, 5):
        pass_time = epoch  + datetime.timedelta(seconds=(time_to_perigee + ke1.tp*i))
        logger.info(f'{i+1} - {pass_time}')
    
    print('\n')

    print('----- Part 3. Next Eclipse Midnight Times -----')

    p1_jd = keHelperFunctions.convert_date_to_jd(epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, epoch.second)

    sun_vector = keHelperFunctions.determine_sun_vector_lf(p1_jd)

    noon, midnight = keHelperFunctions.compute_noon_and_midnight(sun_vector, ke1.raan, ke1.inclination, ke1.aop)

    noon_nu, midnight_nu = keHelperFunctions.compute_noon_and_midnight_true_anomaly(noon)
    # midnight_nu = math.atan2(midnight[1], midnight[0])

    logger.info(f'Midnight Nu: {math.degrees(midnight_nu)}')
    logger.info(f'Noon Nu: {math.degrees(noon_nu)}')

    # Time of flight to midnight
    tof_to_midnight = keHelperFunctions.compute_time_to_midnight(ke1, midnight_nu)
    logger.info(f'TOF to Midnight: {tof_to_midnight} seconds')
    logger.info(f'TOF to Midnight: {tof_to_midnight/60} minutes')


    beta_angle = keHelperFunctions.compute_keplarian_beta_angle(sun_vector, ke1)
    logger.info(f'Beta Angle: {math.degrees(beta_angle)}')

    apparent_angular_radius = keHelperFunctions.compute_angular_radius(np.linalg.norm(ke1.r_vector)/1000)
    logger.info(f'App Rad: {math.degrees(apparent_angular_radius)}')

    epsilon = keHelperFunctions.estimate_eclipse_duration(beta_angle, apparent_angular_radius)
    logger.info(f'Epsilon: {math.degrees(epsilon)}')

    eclipse_duration = keHelperFunctions.compute_eclipse_duration(ke1, epsilon)
    logger.info(f'ECL dur: {eclipse_duration} seconds')
    logger.info(f'ECL dur: {eclipse_duration/60} minutes')

    eclipse_midnight_time = epoch +datetime.timedelta(seconds=tof_to_midnight)
    eclipse_midnight_start_time = eclipse_midnight_time - datetime.timedelta(seconds=eclipse_duration/2)
    eclipse_midnight_end_time = eclipse_midnight_time + datetime.timedelta(seconds=eclipse_duration/2)

    logger.info(f'Eclipse Start Time: {eclipse_midnight_start_time}')
    logger.info(f'Eclipse Midnight Time: {eclipse_midnight_time}')
    logger.info(f'Eclipse End Time: {eclipse_midnight_end_time}')

    print('\n')

    logger.info(f'Estimating next four eclipses...')

    eclipse_data = []
    eclipse_data_headers = ['Eclipse Number', 'Start Time', 'Midnight Time', 'End Time', 'Beta Angle', 'Duration (min)']

    eclipse_midnights = []
    eclipse_midnights.append([eclipse_midnight_time, keHelperFunctions.convert_date_to_jd(eclipse_midnight_time.year, eclipse_midnight_time.month, eclipse_midnight_time.day, eclipse_midnight_time.hour, eclipse_midnight_time.minute, eclipse_midnight_time.second)])
    
    eclipse_data.append([1, eclipse_midnight_start_time, eclipse_midnight_time, eclipse_midnight_end_time, math.degrees(beta_angle), eclipse_duration/60])
    
    for i in range(1,5):
        new_date = eclipse_midnight_time + datetime.timedelta(seconds=ke1.tp*(i))
        eclipse_midnights.append([new_date, keHelperFunctions.convert_date_to_jd(new_date.year, new_date.month, new_date.day, new_date.hour, new_date.minute, new_date.second)])


    for i in range(1,5):
        cur_sun_vector = keHelperFunctions.determine_sun_vector_lf(eclipse_midnights[i][1])
        current_eclipse_midnight = eclipse_midnights[i][0]
        cur_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(cur_sun_vector, ke1)
        cur_epsilon = keHelperFunctions.estimate_eclipse_duration(cur_beta_angle, apparent_angular_radius)
        cur_eclipse_duration = keHelperFunctions.compute_eclipse_duration(ke1, cur_epsilon)

        cur_eclipse_midnight_start_time = current_eclipse_midnight - datetime.timedelta(seconds=cur_eclipse_duration/2)
        cur_eclipse_midnight_end_time = current_eclipse_midnight + datetime.timedelta(seconds=cur_eclipse_duration/2)

        eclipse_data.append([i+1, cur_eclipse_midnight_start_time, current_eclipse_midnight, cur_eclipse_midnight_end_time, math.degrees(cur_beta_angle), cur_eclipse_duration/60])

    print(tabulate.tabulate(eclipse_data, eclipse_data_headers))


if __name__ == '__main__':
    main()
