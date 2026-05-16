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
    
    # # Basic proagator to follow the chaser and target over 60060 seconds. Built to check my work a bit
    # # 
    # chaser_cur_pos = chaser.r_vector
    # chaser_cur_vel = chaser.r_dot_vector

    # target_cur_pos = target.r_vector
    # target_cur_vel = target.r_dot_vector

    # end_time = 864000

    # logger.info('Propagating Chaser and Target a day at .1 second intervals to find burn 2')
    # for i in range(1, end_time+1):
    #     chaser_new_pos, chaser_new_vel = keHelperFunctions.keplarian_rk4(chaser_cur_pos, chaser_cur_vel, .1, chaser.mu)
    #     target_new_pos, target_new_vel = keHelperFunctions.keplarian_rk4(target_cur_pos, target_cur_vel, .1, target.mu)

    #     new_phase_angle = dgsa.determine_angle_between_two_sv(target_new_pos, chaser_new_pos)

    #     chaser_cur_pos = chaser_new_pos
    #     chaser_cur_vel = chaser_new_vel
    #     target_cur_pos = target_new_pos
    #     target_cur_vel = target_new_vel

        # temp_ke = KeplerianElements(chaser_new_pos[0], chaser_new_pos[1], chaser_new_pos[2], chaser_new_vel[0], chaser_new_vel[1], chaser_new_vel[2])

        # check_nu = math.degrees(temp_ke.nu)

        # print(f'New Phase Angle: {math.degrees(new_phase_angle)} deg at {epoch + datetime.timedelta(seconds=i*60)}')

        # if new_phase_angle < 0.0001:
        #     print(f'Phase Angle: {new_phase_angle}')
        #     print(f'Occured at time: {i}')
        #     break

        # if check_nu > 197.63 and check_nu < 197.70:
        #     logger.info(f'Chaser at time {i} .1 seconds since epoch')
        #     temp_ke.print_ke()
        #     break


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

















    # logger.info('Calculating seconds to intercept')
    # seconds_to_rendezvous = np.abs(math.degrees(p5_phase_angle)/math.degrees(p6_phase_rate_angle))
    # logger.info(f'Seconds to Rendevous: {seconds_to_rendezvous}')
    # seconds_to_burn1 = seconds_to_rendezvous - int(chaser.tp/2)
    # logger.info(f'Seconds to Burn 1: {seconds_to_burn1}')

    ###### For problem 7 and on, it makes more sense to simulate the entire scenario, save the data, and answer each question using that data ######
    
    # Set initial values
    target_data_headers = ['Step', 'Seconds Since Epoch', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD']


    # Propagate the target for two days (So we get ample data). 60 seconds per step. 172800 seconds
    target_two_day_orbit_data = []
    target_two_day_orbit_data.append([0, 0, target.r_vector[0], target.r_vector[1], target.r_vector[2], target.r_dot_vector[0], target.r_dot_vector[1], target.r_dot_vector[2]])

    target_cur_pos = target.r_vector
    target_cur_vel = target.r_dot_vector

    logger.info('Propagating Target for 86400 in 1 second intervals')
    for i in range(1, 86400):
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






    # This was just some code to check perigee. I am def doing something wrong here, but it does spit out a pretty close number. That being said I know i'm off because I did the same for apogee and found a similar discrepency

    # logger.info('Checking Intermediate Orbit Perigee')

    # intermediate_chaser_cur_pos = intermediate_chaser_ke.r_vector
    # intermediate_chaser_cur_vel = intermediate_chaser_ke.r_dot_vector

    # for i in range(1, 864000):
    #     intermediate_chaser_new_pos, intermediate_chaser_new_vel = keHelperFunctions.keplarian_rk4(intermediate_chaser_cur_pos, intermediate_chaser_cur_vel, .1, intermediate_chaser_ke.mu)

    #     intermediate_chaser_cur_pos = intermediate_chaser_new_pos
    #     intermediate_chaser_cur_vel = intermediate_chaser_new_vel

    #     temp_ke = KeplerianElements(intermediate_chaser_cur_pos[0], intermediate_chaser_cur_pos[1], intermediate_chaser_cur_pos[2], intermediate_chaser_cur_vel[0], intermediate_chaser_cur_vel[1], intermediate_chaser_cur_vel[2])

    #     temp_ta = math.degrees(temp_ke.nu)

    #     # if temp_ta < 0.2:
    #     #     temp_ke.print_ke()
    #     #     logger.info(f'Pergee of Intermediate: {np.linalg.norm(intermediate_chaser_cur_pos)}')
    #     #     break
    #     if temp_ta < 180.02 and temp_ta > 179.9999:
    #         temp_ke.print_ke()
    #         logger.info(f'Apogee of Intermediate: {np.linalg.norm(intermediate_chaser_cur_pos)}')
    #         break

    # exit(0)













    logger.info('Propagating Intermediate Orbit Out to Second Burn')

    chaser_intermediate_orbit_data = []
    chaser_intermediate_orbit_data.append([0, 0, intermediate_chaser_ke.r_vector[0], intermediate_chaser_ke.r_vector[1], intermediate_chaser_ke.r_vector[2], intermediate_chaser_ke.r_dot_vector[0], intermediate_chaser_ke.r_dot_vector[1], intermediate_chaser_ke.r_dot_vector[2]])

    intermediate_chaser_cur_pos = intermediate_chaser_ke.r_vector
    intermediate_chaser_cur_vel = intermediate_chaser_ke.r_dot_vector

    # logger.info('Computing distance covered over half the intermediate orbit')
    # half_orbit_iteration = np.abs(intermediate_chaser_ke.tp/2 * math.degrees(intermediate_phase_rate_angle))
    # logger.info(f'Degrees changed over half the intermediate orbit: {half_orbit_iteration} deg')
    # half_orbit_angular_change = math.degrees(p5_phase_angle) - half_orbit_iteration
    # logger.info(f'Phase Angle Angular distance we need to travel: {half_orbit_angular_change}')
    # angular_time_to_first_burn = np.abs(half_orbit_angular_change/math.degrees(p6_phase_rate_angle))
    # logger.info(f'Seconds to Angular first burn: {angular_time_to_first_burn}')


    # iterations = int(seconds_to_burn2 - chaser.tp/2) + 3 # Old method where I was kinda correct
    # iterations = int(seconds_to_burn2)
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
    chaser_final_orbit_data.append([0, 0, final_chaser_ke.r_vector[0], final_chaser_ke.r_vector[1], final_chaser_ke.r_vector[2], final_chaser_ke.r_dot_vector[0], final_chaser_ke.r_dot_vector[1], final_chaser_ke.r_dot_vector[2]])

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

    # print('Last line for chaser data')
    # print(all_chaser_data[-1])

    # for line in all_chaser_data:
        # print(line)

    # EXTRA CREDIT SECTION IF I HAVE TIME!
    # Once the second burn has occured, plot the angular separation over time. I'm thinking over the period of a day
    # Thinking about it, if all goes well we should see a slight oscillation over a day
    print('----- Extra Credit -----')




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Turn on Debug verbosity. DEFAULT=False')
    
    main(parser.parse_args())
