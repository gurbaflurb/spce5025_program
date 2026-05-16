import argparse
import math
import datetime
import dateutil
import logging
from pprint import pprint

import tabulate
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

from keplarianElements import KeplerianElements
from groundSite import GroundSite
import keHelperFunctions



def main(args):

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    else:
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
    logger.info('Computing Delta-V to change SMA by problem 3 SMA difference')
    p4_delta_v = keHelperFunctions.estimate_in_plane_burn(chaser.tp, p3_sma_diff)
    logger.info(f'Delta-V to adjust chaser orbital period to target orbital period: {p4_delta_v} m/s')
    print()

    # Problem 5
    print('----- Problem 5 -----')
    p5_phase_angle = dgsa.determine_angle_between_two_sv(target.r_vector, chaser.r_vector)
    logger.info(f'Phase between Target and Chaser: {math.degrees(p5_phase_angle)}')
    
    print()


    # Problem 6
    print('----- Problem 6 -----')
    p6_phase_rate_angle = keHelperFunctions.compute_phase_rate(target.tp, chaser.tp)
    logger.info(f'Current phase angle at epoch: {p6_phase_rate_angle} radians/second')
    logger.info(f'Current phase angle at epoch: {math.degrees(p6_phase_rate_angle)} deg/second')
    logger.info(f'Current phase angle at epoch: {math.degrees(p6_phase_rate_angle)*3600} deg/hour')
    logger.info(f'Current phase angle at epoch: {(math.degrees(p6_phase_rate_angle)*3600)*24} deg/day')
    print()


    logger.info('Calculating Hohmann Burns')

    burn1, burn2 = keHelperFunctions.get_hohmann_transfer_burns(chaser.mu, chaser.semi_major_axis, target.semi_major_axis)

    logger.info(f'Burn 1 to Intermediate Delta-v: {burn1}')
    logger.info(f'Burn 2 to Intermediate Delta-v: {burn2}')

    print()

    logger.info('Applying burn1 delta-v to chaser!')

    test_intermediate_chaser_ke = keHelperFunctions.apply_delta_v_ke(chaser, burn1)
    logger.info('----- Chaser Intermediate Orbit -----')
    test_intermediate_chaser_ke.print_ke()
    print()

    logger.info('Computing Intermediate Phase Rate Angle')
    intermediate_phase_rate_angle = keHelperFunctions.compute_phase_rate(target.tp, test_intermediate_chaser_ke.tp)

    logger.info('Computing distance covered over half the intermediate orbit')
    half_orbit_iteration = np.abs(test_intermediate_chaser_ke.tp/2 * math.degrees(intermediate_phase_rate_angle))
    logger.info(f'Degrees changed over half the intermediate orbit: {half_orbit_iteration} deg')

    half_orbit_angular_change = math.degrees(p5_phase_angle) - half_orbit_iteration
    logger.info(f'Phase Angle Angular distance we need to travel: {half_orbit_angular_change}')

    angular_time_to_first_burn = np.abs(half_orbit_angular_change/math.degrees(p6_phase_rate_angle))
    logger.info(f'Seconds to Angular first burn: {angular_time_to_first_burn}')

    first_burn_epoch = epoch + datetime.timedelta(seconds=angular_time_to_first_burn)
    logger.info(f'Time to first burn: {first_burn_epoch}')


    ###### For problem 7 and on, it makes more sense to simulate the entire scenario, save the data, and answer each question using that data ######
    
    # Set initial values
    target_data_headers = ['Step', 'Seconds Since Epoch', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD']


    # Propagate the target for two days (So we get ample data). 60 seconds per step. 172800 seconds
    target_two_day_orbit_data = []
    target_two_day_orbit_data.append([0, 0, target.r_vector[0], target.r_vector[1], target.r_vector[2], target.r_dot_vector[0], target.r_dot_vector[1], target.r_dot_vector[2]])

    target_cur_pos = target.r_vector
    target_cur_vel = target.r_dot_vector

    logger.info('Propagating Target for 86400 in 1 second intervals')
    for i in range(1, 86401):
        target_new_pos, target_new_vel = keHelperFunctions.keplarian_rk4(target_cur_pos, target_cur_vel, 1, target.mu)

        target_cur_pos = target_new_pos
        target_cur_vel = target_new_vel

        data_line = [i, i, target_new_pos[0], target_new_pos[1], target_new_pos[2], target_new_vel[0], target_new_vel[1], target_new_vel[2]]

        target_two_day_orbit_data.append(data_line)

    logger.info('Target Orbit Propagated 86400 seconds!')

    # Set a counter to count up the seconds since the epoch
    seconds_since_epoch = 0



    logger.info(f'Propagating Chaser Orbit out to first burn')

    chaser_initial_orbit_data = []
    chaser_initial_orbit_data.append([0, 0, chaser.r_vector[0], chaser.r_vector[1], chaser.r_vector[2], chaser.r_dot_vector[0], chaser.r_dot_vector[1], chaser.r_dot_vector[2]])

    chaser_cur_pos = chaser.r_vector
    chaser_cur_vel = chaser.r_dot_vector

    for i in range(1, int(angular_time_to_first_burn)+1):
        chaser_new_pos, chaser_new_vel = keHelperFunctions.keplarian_rk4(chaser_cur_pos, chaser_cur_vel, 1, chaser.mu)

        chaser_cur_pos = chaser_new_pos
        chaser_cur_vel = chaser_new_vel

        seconds_since_epoch += 1

        data_line = [i, i, chaser_new_pos[0], chaser_new_pos[1], chaser_new_pos[2], chaser_new_vel[0], chaser_new_vel[1], chaser_new_vel[2]]

        chaser_initial_orbit_data.append(data_line)


    logger.info('Chaser Orbit Propagated to first burn!')


    logger.info('Calculating Hohmann Burns')

    first_burn_target_data = target_two_day_orbit_data[int(seconds_since_epoch)]
    first_burn_target_ke = KeplerianElements(first_burn_target_data[2], first_burn_target_data[3], first_burn_target_data[4], first_burn_target_data[5], first_burn_target_data[6], first_burn_target_data[7])

    first_burn_chaser_data = chaser_initial_orbit_data[int(seconds_since_epoch)]
    first_burn_chaser_ke = KeplerianElements(first_burn_chaser_data[2], first_burn_chaser_data[3], first_burn_chaser_data[4], first_burn_chaser_data[5], first_burn_chaser_data[6], first_burn_chaser_data[7])

    logger.info('Applying burn1 delta-v to chaser!')

    intermediate_chaser_ke = keHelperFunctions.apply_delta_v_ke(first_burn_chaser_ke, burn1)
    logger.info('----- Chaser Intermediate Orbit -----')
    intermediate_chaser_ke.print_ke()
    print()
    logger.info('----- Target Orbit At Intermediate Burn -----')
    first_burn_target_ke.print_ke()
    print()

    intermediate_phase_angle = dgsa.determine_angle_between_two_sv(first_burn_target_ke.r_vector, intermediate_chaser_ke.r_vector)
    logger.info(f'Phase Angle between Intermediate Chaser and Target: {math.degrees(intermediate_phase_angle)}')
    print()

    intermediate_phase_rate_angle = keHelperFunctions.compute_phase_rate(first_burn_target_ke.tp, intermediate_chaser_ke.tp)
    logger.info(f'Current phase angle at intermediate: {intermediate_phase_rate_angle} radians/second')
    logger.info(f'Current phase angle at intermediate: {math.degrees(intermediate_phase_rate_angle)} deg/second')
    logger.info(f'Current phase angle at intermediate: {math.degrees(intermediate_phase_rate_angle)*3600} deg/hour')
    logger.info(f'Current phase angle at intermediate: {(math.degrees(intermediate_phase_rate_angle)*3600)*24} deg/day')
    print()

    seconds_to_burn2 = np.abs(math.degrees(intermediate_phase_angle)/math.degrees(intermediate_phase_rate_angle))

    logger.info(f'Seconds to Second Burn: {seconds_to_burn2}')

    second_burn_epoch = first_burn_epoch + datetime.timedelta(seconds=seconds_to_burn2)

    logger.info(f'Second Burn Epoch: {second_burn_epoch}')

    logger.info('Propagating Intermediate Orbit Out to Second Burn')

    chaser_intermediate_orbit_data = []
    # chaser_intermediate_orbit_data.append([0, 0, intermediate_chaser_ke.r_vector[0], intermediate_chaser_ke.r_vector[1], intermediate_chaser_ke.r_vector[2], intermediate_chaser_ke.r_dot_vector[0], intermediate_chaser_ke.r_dot_vector[1], intermediate_chaser_ke.r_dot_vector[2]])

    intermediate_chaser_cur_pos = intermediate_chaser_ke.r_vector
    intermediate_chaser_cur_vel = intermediate_chaser_ke.r_dot_vector

    iterations = np.abs((second_burn_epoch - first_burn_epoch).total_seconds())

    logger.info(f'Iterations until Second Burn Point: {iterations}')

    for i in range(1, int(iterations)+1):
        intermediate_chaser_new_pos, intermediate_chaser_new_vel = keHelperFunctions.keplarian_rk4(intermediate_chaser_cur_pos, intermediate_chaser_cur_vel, 1, intermediate_chaser_ke.mu)

        intermediate_chaser_cur_pos = intermediate_chaser_new_pos
        intermediate_chaser_cur_vel = intermediate_chaser_new_vel

        seconds_since_epoch += 1

        data_line = [i, seconds_since_epoch, intermediate_chaser_new_pos[0], intermediate_chaser_new_pos[1], intermediate_chaser_new_pos[2], intermediate_chaser_new_vel[0], intermediate_chaser_new_vel[1], intermediate_chaser_new_vel[2]]

        chaser_intermediate_orbit_data.append(data_line)

        

    logger.info('Intermediate Orbit Propagated to Second Burn!')
    print()


    second_burn_target_data = target_two_day_orbit_data[seconds_since_epoch]
    second_burn_target_ke = KeplerianElements(second_burn_target_data[2], second_burn_target_data[3], second_burn_target_data[4], second_burn_target_data[5], second_burn_target_data[6], second_burn_target_data[7])
    print(f'Target SMA: {second_burn_target_ke.semi_major_axis}')

    second_burn_chaser_data = chaser_intermediate_orbit_data[-1]
    second_burn_chaser_ke = KeplerianElements(second_burn_chaser_data[2], second_burn_chaser_data[3], second_burn_chaser_data[4], second_burn_chaser_data[5], second_burn_chaser_data[6], second_burn_chaser_data[7])
    print(f'Chaser SMA: {second_burn_chaser_ke.semi_major_axis}')



    logger.info('Applying new Delta-V to Chaser')

    final_chaser_ke = keHelperFunctions.apply_delta_v_ke(second_burn_chaser_ke, burn2)
    logger.info('----- Final Chaser Orbit -----')
    final_chaser_ke.print_ke()
    print()
    logger.info('----- Target Orbit At Second Burn -----')
    second_burn_target_ke.print_ke()

    new_phase_angle_deg = dgsa.determine_angle_between_two_sv(final_chaser_ke.r_vector, second_burn_target_ke.r_vector)
    print(f'Phase Angle at second Burn: {math.degrees(new_phase_angle_deg)}')

    print()

    iterations_to_next_day = 86401 - seconds_since_epoch

    logger.info(f'Propagating Final Orbit Out {iterations_to_next_day} seconds')

    chaser_final_orbit_data = []
    # chaser_final_orbit_data.append([0, 0, final_chaser_ke.r_vector[0], final_chaser_ke.r_vector[1], final_chaser_ke.r_vector[2], final_chaser_ke.r_dot_vector[0], final_chaser_ke.r_dot_vector[1], final_chaser_ke.r_dot_vector[2]])

    final_chaser_cur_pos = final_chaser_ke.r_vector
    final_chaser_cur_vel = final_chaser_ke.r_dot_vector

    for i in range(1, iterations_to_next_day):
        final_chaser_new_pos, final_chaser_new_vel = keHelperFunctions.keplarian_rk4(final_chaser_cur_pos, final_chaser_cur_vel, 1, final_chaser_ke.mu)

        final_chaser_cur_pos = final_chaser_new_pos
        final_chaser_cur_vel = final_chaser_new_vel

        seconds_since_epoch += 1

        data_line = [i, seconds_since_epoch, final_chaser_cur_pos[0], final_chaser_cur_pos[1], final_chaser_cur_pos[2], final_chaser_cur_vel[0], final_chaser_cur_vel[1], final_chaser_cur_vel[2]]

        chaser_final_orbit_data.append(data_line)

    logger.info(f'Final Orbit Propagated out {iterations_to_next_day} seconds!')
    print()

    logger.info(f'Calculating Distance between Chaser and Target')

    distance = np.linalg.norm(second_burn_target_ke.r_vector - final_chaser_ke.r_vector)

    logger.info(f'Distance from Target to Chaser: {distance} meters')



    print('----- Problem 17 -----')

    all_chaser_data = chaser_initial_orbit_data 

    for line in chaser_intermediate_orbit_data:
        all_chaser_data.append(line)

    for line in chaser_final_orbit_data:
        all_chaser_data.append(line)

    first_burn_data = all_chaser_data[int((first_burn_epoch - epoch).total_seconds())]
    
    first_burn_ke = KeplerianElements(first_burn_data[2], first_burn_data[3], first_burn_data[4], first_burn_data[5], first_burn_data[6], first_burn_data[7])

    first_burn_ecef_pos = keHelperFunctions.convert_eci_ecef(first_burn_ke.r_vector, int((first_burn_epoch - epoch).total_seconds()))

    relative_pos_dgsa =  dgsa.compute_relative_pos(first_burn_ecef_pos)
    relative_pos_vtsa =  vtsa.compute_relative_pos(first_burn_ecef_pos)

    dgsa_topo = keHelperFunctions.convert_ecef_topocentric(math.radians(dgsa.lat), math.radians(dgsa.lon), relative_pos_dgsa)
    vtsa_topo = keHelperFunctions.convert_ecef_topocentric(math.radians(vtsa.lat), math.radians(vtsa.lon), relative_pos_vtsa)

    dgsa_az, dgsa_el = keHelperFunctions.compute_azimuth_elevation(dgsa_topo)
    vtsa_az, vtsa_el = keHelperFunctions.compute_azimuth_elevation(vtsa_topo)

    print('First Burn Check!')
    print(f'DGSA AZ: {math.degrees(dgsa_az)}')
    print(f'DGSA EL: {math.degrees(dgsa_el)}')
    print(f'VTSA AZ: {math.degrees(vtsa_az)}')
    print(f'VTSA EL: {math.degrees(vtsa_el)}')
    print()

    print('----- Problem 18 -----')

    second_burn_data = all_chaser_data[int((second_burn_epoch - epoch).total_seconds())]
    
    second_burn_ke = KeplerianElements(second_burn_data[2], second_burn_data[3], second_burn_data[4], second_burn_data[5], second_burn_data[6], second_burn_data[7])

    second_burn_ecef_pos = keHelperFunctions.convert_eci_ecef(second_burn_ke.r_vector, (second_burn_epoch - epoch).total_seconds())

    relative_pos_dgsa =  dgsa.compute_relative_pos(second_burn_ecef_pos)
    relative_pos_vtsa =  vtsa.compute_relative_pos(second_burn_ecef_pos)

    dgsa_topo = keHelperFunctions.convert_ecef_topocentric(math.radians(dgsa.lat), math.radians(dgsa.lon), relative_pos_dgsa)
    vtsa_topo = keHelperFunctions.convert_ecef_topocentric(math.radians(vtsa.lat), math.radians(vtsa.lon), relative_pos_vtsa)

    dgsa_az, dgsa_el = keHelperFunctions.compute_azimuth_elevation(dgsa_topo)
    vtsa_az, vtsa_el = keHelperFunctions.compute_azimuth_elevation(vtsa_topo)

    print('Second Burn Check!')
    print(f'DGSA AZ: {math.degrees(dgsa_az)}')
    print(f'DGSA EL: {math.degrees(dgsa_el)}')
    print(f'VTSA AZ: {math.degrees(vtsa_az)}')
    print(f'VTSA EL: {math.degrees(vtsa_el)}')
    print()

    
    print(f'----- Problem 19 -----')
    initial_jd = keHelperFunctions.convert_date_to_jd(epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, epoch.second)
    initial_sun_vector = keHelperFunctions.determine_sun_vector_lf(initial_jd)
    
    initial_target_noon, initial_target_midnight = keHelperFunctions.compute_noon_and_midnight(initial_sun_vector, target.raan, target.inclination, target.aop)
    target_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(initial_sun_vector, target)

    target_apparent_angular_radius = keHelperFunctions.compute_angular_radius(np.linalg.norm(target.r_vector)/1000)

    target_epsilon = keHelperFunctions.estimate_eclipse_duration(target_beta_angle, target_apparent_angular_radius)

    eclipse_duration = keHelperFunctions.compute_eclipse_duration(target, target_epsilon)
    logger.info(f'Target ECL dur: {eclipse_duration} seconds')
    logger.info(f'Target ECL dur: {eclipse_duration/60} minutes')
    print()

    initial_chaser_noon, initial_chaser_midnight = keHelperFunctions.compute_noon_and_midnight(initial_sun_vector, chaser.raan, chaser.inclination, chaser.aop)
    chaser_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(initial_sun_vector, chaser)

    chaser_apparent_angular_radius = keHelperFunctions.compute_angular_radius(np.linalg.norm(chaser.r_vector)/1000)

    chaser_epsilon = keHelperFunctions.estimate_eclipse_duration(chaser_beta_angle, chaser_apparent_angular_radius)

    eclipse_duration = keHelperFunctions.compute_eclipse_duration(chaser, chaser_epsilon)
    logger.info(f'Initial Chaser ECL dur: {eclipse_duration} seconds')
    logger.info(f'Initial Chaser ECL dur: {eclipse_duration/60} minutes')
    print()



    print(f'----- Problem 20 -----')
    logger.info('Burn 1 Eclipse Calculations')
    burn1_jd = keHelperFunctions.convert_date_to_jd(first_burn_epoch.year, first_burn_epoch.month, first_burn_epoch.day, first_burn_epoch.hour, first_burn_epoch.minute, first_burn_epoch.second)

    burn1_sun_vector = keHelperFunctions.determine_sun_vector_lf(burn1_jd)

    burn1_noon, burn1_midnight = keHelperFunctions.compute_noon_and_midnight(burn1_sun_vector, first_burn_chaser_ke.raan, first_burn_chaser_ke.inclination, first_burn_chaser_ke.aop)

    burn1_noon_nu, burn1_midnight_nu = keHelperFunctions.compute_noon_and_midnight_true_anomaly(burn1_noon)

    burn1_tof_to_midnight = keHelperFunctions.compute_time_to_midnight(first_burn_chaser_ke, burn1_midnight_nu)

    burn1_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(burn1_sun_vector, first_burn_chaser_ke)

    burn1_apparent_angular_radius = keHelperFunctions.compute_angular_radius(np.linalg.norm(first_burn_chaser_ke.r_vector)/1000)

    burn1_epsilon = keHelperFunctions.estimate_eclipse_duration(burn1_beta_angle, burn1_apparent_angular_radius)

    burn1_eclipse_duration = keHelperFunctions.compute_eclipse_duration(first_burn_chaser_ke, burn1_epsilon)

    burn1_eclipse_midnight_time = epoch +datetime.timedelta(seconds=burn1_tof_to_midnight)
    burn1_eclipse_midnight_start_time = burn1_eclipse_midnight_time - datetime.timedelta(seconds=burn1_eclipse_duration/2)
    burn1_eclipse_midnight_end_time = burn1_eclipse_midnight_time + datetime.timedelta(seconds=burn1_eclipse_duration/2)

    logger.info(f'Estimating next 15 eclipses...')

    eclipse_data = []
    eclipse_data_headers = ['Eclipse Number', 'Start Time', 'Midnight Time', 'End Time', 'Beta Angle', 'Duration (min)']

    burn1_eclipse_midnights = []
    burn1_eclipse_midnights.append([burn1_eclipse_midnight_time, keHelperFunctions.convert_date_to_jd(burn1_eclipse_midnight_time.year, burn1_eclipse_midnight_time.month, burn1_eclipse_midnight_time.day, burn1_eclipse_midnight_time.hour, burn1_eclipse_midnight_time.minute, burn1_eclipse_midnight_time.second)])
    
    eclipse_data.append([1, burn1_eclipse_midnight_start_time, burn1_eclipse_midnight_time, burn1_eclipse_midnight_end_time, math.degrees(burn1_beta_angle), burn1_eclipse_duration/60])
    
    for i in range(1,15):
        new_date = burn1_eclipse_midnight_time + datetime.timedelta(seconds=chaser.tp*(i))
        burn1_eclipse_midnights.append([new_date, keHelperFunctions.convert_date_to_jd(new_date.year, new_date.month, new_date.day, new_date.hour, new_date.minute, new_date.second)])


    for i in range(1,15):
        cur_sun_vector = keHelperFunctions.determine_sun_vector_lf(burn1_eclipse_midnights[i][1])
        current_eclipse_midnight = burn1_eclipse_midnights[i][0]
        cur_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(cur_sun_vector, chaser)
        cur_epsilon = keHelperFunctions.estimate_eclipse_duration(cur_beta_angle, burn1_apparent_angular_radius)
        cur_eclipse_duration = keHelperFunctions.compute_eclipse_duration(chaser, cur_epsilon)

        cur_eclipse_midnight_start_time = current_eclipse_midnight - datetime.timedelta(seconds=cur_eclipse_duration/2)
        cur_eclipse_midnight_end_time = current_eclipse_midnight + datetime.timedelta(seconds=cur_eclipse_duration/2)

        eclipse_data.append([i+1, cur_eclipse_midnight_start_time, current_eclipse_midnight, cur_eclipse_midnight_end_time, math.degrees(cur_beta_angle), cur_eclipse_duration/60])

    print(tabulate.tabulate(eclipse_data, eclipse_data_headers))

    logger.info(f'Time of First Burn: {first_burn_epoch}')

    print()

    logger.info('Burn 2 Eclipse Calculations')
    burn2_jd = keHelperFunctions.convert_date_to_jd(second_burn_epoch.year, second_burn_epoch.month, second_burn_epoch.day, second_burn_epoch.hour, second_burn_epoch.minute, second_burn_epoch.second)

    burn2_sun_vector = keHelperFunctions.determine_sun_vector_lf(burn2_jd)

    burn2_noon, burn2_midnight = keHelperFunctions.compute_noon_and_midnight(burn2_sun_vector, second_burn_chaser_ke.raan, second_burn_chaser_ke.inclination, second_burn_chaser_ke.aop)

    burn2_noon_nu, burn2_midnight_nu = keHelperFunctions.compute_noon_and_midnight_true_anomaly(burn2_noon)

    burn2_tof_to_midnight = keHelperFunctions.compute_time_to_midnight(second_burn_chaser_ke, burn2_midnight_nu)

    burn2_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(burn2_sun_vector, second_burn_chaser_ke)

    burn2_apparent_angular_radius = keHelperFunctions.compute_angular_radius(np.linalg.norm(second_burn_chaser_ke.r_vector)/1000)

    burn2_epsilon = keHelperFunctions.estimate_eclipse_duration(burn2_beta_angle, burn2_apparent_angular_radius)

    burn2_eclipse_duration = keHelperFunctions.compute_eclipse_duration(second_burn_chaser_ke, burn2_epsilon)

    burn2_eclipse_midnight_time = epoch +datetime.timedelta(seconds=burn2_tof_to_midnight)
    burn2_eclipse_midnight_start_time = burn2_eclipse_midnight_time - datetime.timedelta(seconds=burn2_eclipse_duration/2)
    burn2_eclipse_midnight_end_time = burn2_eclipse_midnight_time + datetime.timedelta(seconds=burn2_eclipse_duration/2)

    logger.info(f'Estimating next 15 eclipses...')

    eclipse_data = []
    eclipse_data_headers = ['Eclipse Number', 'Start Time', 'Midnight Time', 'End Time', 'Beta Angle', 'Duration (min)']

    burn2_eclipse_midnights = []
    burn2_eclipse_midnights.append([burn2_eclipse_midnight_time, keHelperFunctions.convert_date_to_jd(burn2_eclipse_midnight_time.year, burn2_eclipse_midnight_time.month, burn2_eclipse_midnight_time.day, burn2_eclipse_midnight_time.hour, burn2_eclipse_midnight_time.minute, burn2_eclipse_midnight_time.second)])
    
    eclipse_data.append([1, burn2_eclipse_midnight_start_time, burn2_eclipse_midnight_time, burn2_eclipse_midnight_end_time, math.degrees(burn2_beta_angle), burn2_eclipse_duration/60])
    
    for i in range(1,15):
        new_date = burn2_eclipse_midnight_time + datetime.timedelta(seconds=second_burn_chaser_ke.tp*(i))
        burn1_eclipse_midnights.append([new_date, keHelperFunctions.convert_date_to_jd(new_date.year, new_date.month, new_date.day, new_date.hour, new_date.minute, new_date.second)])


    for i in range(1,15):
        cur_sun_vector = keHelperFunctions.determine_sun_vector_lf(burn1_eclipse_midnights[i][1])
        current_eclipse_midnight = burn1_eclipse_midnights[i][0]
        cur_beta_angle = keHelperFunctions.compute_keplarian_beta_angle(cur_sun_vector, chaser)
        cur_epsilon = keHelperFunctions.estimate_eclipse_duration(cur_beta_angle, burn1_apparent_angular_radius)
        cur_eclipse_duration = keHelperFunctions.compute_eclipse_duration(chaser, cur_epsilon)

        cur_eclipse_midnight_start_time = current_eclipse_midnight - datetime.timedelta(seconds=cur_eclipse_duration/2)
        cur_eclipse_midnight_end_time = current_eclipse_midnight + datetime.timedelta(seconds=cur_eclipse_duration/2)

        eclipse_data.append([i+1, cur_eclipse_midnight_start_time, current_eclipse_midnight, cur_eclipse_midnight_end_time, math.degrees(cur_beta_angle), cur_eclipse_duration/60])

    print(tabulate.tabulate(eclipse_data, eclipse_data_headers))

    logger.info(f'Time of Second Burn: {second_burn_epoch}')
    print()





    # EXTRA CREDIT SECTION IF I HAVE TIME!
    # Once the second burn has occured, plot the angular separation over time. I'm thinking over the period of a day
    # Thinking about it, if all goes well we should see a slight oscillation over a day
    print('----- Extra Credit -----')

    seconds_since_epoch_data = []
    phase_angle_data = []

    logger.info(f'Generating Phase Angles from {first_burn_epoch - datetime.timedelta(seconds=1000)} to {epoch + datetime.timedelta(seconds=86400)}')
    # for i in range(70000, 86401):
    for i in range(int((first_burn_epoch - datetime.timedelta(seconds=1000) - epoch).total_seconds()), 86401):
        cur_chaser_data = all_chaser_data[i]
        cur_target_data =  target_two_day_orbit_data[i]

        if (cur_chaser_data[1] != cur_target_data[1]):
            print('Seconds since epoch dont match!')
            print(f'Chaser Data: {cur_chaser_data[1]}')
            print(f'Target Data: {cur_target_data[1]}')
            
            exit(0)

        # logger.info(f'Seconds since Epoch: {i}')

        chaser_pos = [cur_chaser_data[2], cur_chaser_data[3], cur_chaser_data[4]]
        target_pos = [cur_target_data[2], cur_target_data[3], cur_target_data[4]]

        # logger.info(f'Chaser POS: {chaser_pos}')
        # logger.info(f'Target POS: {target_pos}')

        current_phase_angle = math.degrees(dgsa.determine_angle_between_two_sv(target_pos, chaser_pos))

        # if i % 20 == 0:
        #     logger.info(f'Current Phase Angle at {i}: {current_phase_angle}')

        phase_angle_data.append(current_phase_angle)
        seconds_since_epoch_data.append(i)

    logger.info('Plotting data to chart')
    keHelperFunctions.graph_sma_difference(f'Phase Angle Over Two Burns from\n{epoch + datetime.timedelta(seconds=70000)} to {epoch + datetime.timedelta(seconds=86400)}', 'Seconds Since Epoch', 'Phase Angle', seconds_since_epoch_data, phase_angle_data, 'final_exam_extracredit.png')

    logger.info('Writing Data to CSV')

    curated_data = []
    for line in all_chaser_data:
        if line[1] % 60 == 0:
            curated_data.append(line)

    keHelperFunctions.print_data_to_csv('problem_16_orbit_propagation_data.csv', target_data_headers, curated_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Turn on Debug verbosity. DEFAULT=False')
    
    main(parser.parse_args())
