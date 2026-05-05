import math
import datetime
import dateutil
import logging
from pprint import pprint

import tabulate
import matplotlib.pyplot as plt
import numpy as np

from keplarianElements import KeplerianElements
from groundSite import GroundSite
import keHelperFunctions


logger = logging.getLogger(__name__)

def main():

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    vectors_file = 'vectors.yaml'
    vector_data = keHelperFunctions.read_in_yaml(vectors_file)

    target = KeplerianElements(vector_data['vectors']['vector1']['x_pos'],
                               vector_data['vectors']['vector1']['y_pos'],
                               vector_data['vectors']['vector1']['z_pos'],
                               vector_data['vectors']['vector1']['x_velocity'],
                               vector_data['vectors']['vector1']['y_velocity'],
                               vector_data['vectors']['vector1']['z_velocity'])
    
    chaser = KeplerianElements(vector_data['vectors']['vector2']['x_pos'],
                               vector_data['vectors']['vector2']['y_pos'],
                               vector_data['vectors']['vector2']['z_pos'],
                               vector_data['vectors']['vector2']['x_velocity'],
                               vector_data['vectors']['vector2']['y_velocity'],
                               vector_data['vectors']['vector2']['z_velocity'])
    

    # Epoch provided: 2012-05-01 00:00:00
    epoch = datetime.datetime(2012, 5, 1, 0, 0, 0)
    logger.info(f'Provided Epoch: {epoch}')

    dgsa = GroundSite(vector_data['ground_site']['DGSA']['geodetic_lat'],
                      vector_data['ground_site']['DGSA']['lon'],
                      vector_data['ground_site']['DGSA']['height'])

    vtsa = GroundSite(vector_data['ground_site']['VTSA']['geodetic_lat'],
                      vector_data['ground_site']['VTSA']['lon'],
                      vector_data['ground_site']['VTSA']['height'])



    # Problem 1
    print('----- Problem 1 -----')
    logger.info('Target Keplarian Elements')
    target.print_ke()
    print()
    logger.info('Chaser Keplarian Elements')
    chaser.print_ke()
    print()

    # Problem 2
    print('----- Problem 2 -----')
    logger.info(f'Chaser Orbital Period: {chaser.tp} seconds')
    logger.info(f'Target Orbital Period: {target.tp} seconds')
    print()

    # Problem 3
    print('----- Problem 3 -----')
    p3_sma_diff = target.semi_major_axis - chaser.semi_major_axis
    logger.info(f'SMA Difference between target and chaser: {p3_sma_diff} meters')
    print()

    # Problem 4
    print('----- Problem 4 -----')
    # Delta-V required to change SMA by problem 3 difference
    print('----- Problem 3 -----')
    logger.info('Computing Delta-V to change SMA by problem 3 SMA difference')
    p4_delta_v = keHelperFunctions.estimate_in_plane_burn(chaser.tp, p3_sma_diff)
    logger.info(f'Delta-V to adjust chaser orbital period to target orbital period: {p4_delta_v} m/s')
    print()

    # Problem 5
    print('----- Problem 5 -----')
    p5_phase_angle = dgsa.determine_angle_between_two_sv(target.r_vector, chaser.r_vector)
    logger.info(f'Phase between Target and Chaser: {math.degrees(p5_phase_angle)}')
    
    # chaser_cur_pos = chaser.r_vector
    # chaser_cur_vel = chaser.r_dot_vector

    # target_cur_pos = target.r_vector
    # target_cur_vel = target.r_dot_vector

    # for i in range(1,1001):
    #     chaser_new_pos, chaser_new_vel = keHelperFunctions.keplarian_rk4(chaser_cur_pos, chaser_cur_vel, 60, chaser.mu)
    #     target_new_pos, target_new_vel = keHelperFunctions.keplarian_rk4(target_cur_pos, target_cur_vel, 60, target.mu)

    #     new_phase_angle = dgsa.determine_angle_between_two_sv(target_new_pos, chaser_new_pos)

    #     chaser_cur_pos = chaser_new_pos
    #     chaser_cur_vel = chaser_new_vel
    #     target_cur_pos = target_new_pos
    #     target_cur_vel = target_new_vel

    #     print(f'New Phase Angle: {math.degrees(new_phase_angle)} deg at {epoch + datetime.timedelta(seconds=i*60)}')



    print()


    # Problem 6
    print('----- Problem 6 -----')
    p6_phase_rate_angle = keHelperFunctions.compute_phase_rate(target.tp, chaser.tp)
    logger.info(f'Current phase angle at epoch: {p6_phase_rate_angle} radians/second')
    logger.info(f'Current phase angle at epoch: {math.degrees(p6_phase_rate_angle)} deg/second')
    logger.info(f'Current phase angle at epoch: {math.degrees(p6_phase_rate_angle)*3600} deg/hour')
    logger.info(f'Current phase angle at epoch: {(math.degrees(p6_phase_rate_angle)*3600)*24} deg/day')
    print()


    # Problem 7
    print('----- Problem 7 -----')


    # Problem 8
    print('----- Problem 8 -----')


    # Problem 9
    print('----- Problem 9 -----')


    # Problem 10
    print('----- Problem 10 -----')


    # Problem 11
    print('----- Problem 11 -----')
    

    # Problem 12
    print('----- Problem 12 -----')

    # Problem 13
    print('----- Problem 13 -----')

    # Problem 14
    print('----- Problem 14 -----')

    # Problem 15
    print('----- Problem 15 -----')

    # Problem 16
    print('----- Problem 16 -----')
    # End epoch provided: 2012-05-02 00:00:00
    p16_end_epoch = datetime.datetime(2012, 5, 2, 0, 0, 0)

    # Problem 17
    print('----- Problem 17 -----')

    # Problem 18
    print('----- Problem 18 -----')

    # Problem 19
    print('----- Problem 19 -----')

    # Problem 20
    print('----- Problem 20 -----')


    # EXTRA CREDIT SECTION IF I HAVE TIME!







if __name__ == '__main__':
    main()
