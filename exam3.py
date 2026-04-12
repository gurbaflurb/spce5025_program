from keplarianElements import KeplerianElements
from groundSite import GroundSite
import keHelperFunctions

import math
import datetime
import logging

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

    sv_1_altd = vector_data['vectors']['vector1']['alt']
    sv_1_latd = vector_data['vectors']['vector1']['latd']
    sv_1_lon = vector_data['vectors']['vector1']['lon']

    ke2 = KeplerianElements(vector_data['vectors'][f'vector2']['x_pos'],
                               vector_data['vectors']['vector2']['y_pos'],
                               vector_data['vectors']['vector2']['z_pos'],
                               vector_data['vectors']['vector2']['x_velocity'],
                               vector_data['vectors']['vector2']['y_velocity'],
                               vector_data['vectors']['vector2']['z_velocity'])

    sv_2_altd = vector_data['vectors']['vector2']['alt']
    sv_2_latd = vector_data['vectors']['vector2']['latd']
    sv_2_lon = vector_data['vectors']['vector2']['lon']

    pogo_gs = GroundSite(vector_data['ground_site']['POGO']['geodetic_lat'],
                         vector_data['ground_site']['POGO']['lon'],
                         vector_data['ground_site']['POGO']['height'],
                         vector_data['ground_site']['POGO']['X'],
                         vector_data['ground_site']['POGO']['Y'],
                         vector_data['ground_site']['POGO']['Z'],
                         vector_data['ground_site']['POGO']['radius'])

    # Epoch provided: 04/20/2012 00:00:00
    epoch = datetime.datetime(2012, 4, 20, 0, 0, 0)

    # We're given a couple starting notes and instructions before we get started
    #
    # a. Note that the input reference frame for the propagations is True of Date
    #    (TOD), and the output frame will also be TOD.
    #
    # b. The equations used to create a ground site vector from latitude, longitude,
    #    and altitude, return the vector in the Earth Centered, Earth Fixed (ECEF)
    #    reference frame.
    #
    # c. Azimuth and elevation from a ground site to a satellite are computed
    #    in the Topocentric reference frame
    #
    # d. The TOD→ECEF transformation is a Z-rotation with an angular magnitude
    #    corresponding to the Greenwich Hour Angle (GHA).
    #
    # e. For the final two problems it will probably be necessary to interpolate
    #    between time points of the propagation to find times where the
    #    satellite(s) cross a particular angular threshold. You can use linear
    #    interpolation to do this. (Note that you will get more accurate interpolation
    #    results if your trajectory results are available at a smaller step size.)
    #


    # Problem 1
    logger.info('----- Problem 1 -----')

    logger.info('Initializing Problem 1 Variables')
    p1_step_size = 200 # Seconds
    p1_number_of_steps = 217 # 12 hours, which is 216 steps, but I'm off by one, so for my implementation its 217

    # Using central body gravity only, propogate each SV defined above.
    # Output the results in the following table format:
    # Time from Epoch | X | Y | Z | XD | YD | ZD
    # I'm gonna add step to my implementation for tracking
    p1_headers = ['Step', 'Seconds from Epoch', 'X', 'Y', 'Z', 'XD', 'YD', 'ZD']

    p1_ke1 = ke1
    p1_ke2 = ke2

    p1_mu = ke1.mu

    p1_cur_date = epoch
    ke1_cur_pos = ke1.r_vector
    ke1_cur_vel = ke1.r_dot_vector
    ke2_cur_pos = ke2.r_vector
    ke2_cur_vel = ke2.r_dot_vector

    p1_sv1_pos_and_vel_vectors = []
    p1_sv2_pos_and_vel_vectors = []


    # Initial Positions vectors
    p1_sv1_pos_and_vel_vectors.append([0, 0, ke1_cur_pos[0], ke1_cur_pos[1],
                                       ke1_cur_pos[2], ke1_cur_vel[0],
                                       ke1_cur_vel[1], ke1_cur_vel[2]])
    p1_sv2_pos_and_vel_vectors.append([0, 0, ke2_cur_pos[0], ke2_cur_pos[1],
                                       ke2_cur_pos[2], ke2_cur_vel[0],
                                       ke2_cur_vel[1], ke2_cur_vel[2]])

    logger.info('Simulating SV positions using central body only...')
    for i in range(1, p1_number_of_steps):
        # Initialize some data
        current_step = i
        seconds_since_epoch = i * p1_step_size
        p1_cur_date = p1_cur_date + datetime.timedelta(seconds=p1_step_size)
        # print(f'{p1_cur_date.day}-{p1_cur_date.month}-{p1_cur_date.year} {p1_cur_date.hour}:{p1_cur_date.minute}:{p1_cur_date.second}')

        # Do the central body for ke1 first
        # new_pos, new_vel = keHelperFunctions.keplarian_rk4_perturbations(ke1_cur_pos, ke1_cur_vel, )  (ke1_cur_pos, ke1_cur_vel, p1_step_size, p1_mu)
        new_pos, new_vel = keHelperFunctions.keplarian_rk4(ke1_cur_pos, ke1_cur_vel, p1_step_size, p1_mu)
        ke1_cur_pos = new_pos
        ke1_cur_vel = new_vel

        p1_sv1_pos_and_vel_vectors.append([current_step, seconds_since_epoch, ke1_cur_pos[0], ke1_cur_pos[1], ke1_cur_pos[2], ke1_cur_vel[0], ke1_cur_vel[1], ke1_cur_vel[2]])

        # Now do the central body for ke2
        new_pos, new_vel = keHelperFunctions.keplarian_rk4(ke2_cur_pos, ke2_cur_vel, p1_step_size, p1_mu)
        ke2_cur_pos = new_pos
        ke2_cur_vel = new_vel

        p1_sv2_pos_and_vel_vectors.append([current_step, seconds_since_epoch, ke2_cur_pos[0], ke2_cur_pos[1], ke2_cur_pos[2], ke2_cur_vel[0], ke2_cur_vel[1], ke2_cur_vel[2]])
        if i == 64:
            logger.info(f'Step: {i}')
            logger.info(f'Seconds Since Epoch: {seconds_since_epoch}')
            logger.info(f'Time: {p1_cur_date}')
            logger.info(f'X1: {ke1_cur_pos[0]}')
            logger.info(f'Y1: {ke1_cur_pos[1]}')
            logger.info(f'Z1: {ke1_cur_pos[2]}')
            logger.info(f'XD1: {ke1_cur_vel[0]}')
            logger.info(f'YD1: {ke1_cur_vel[1]}')
            logger.info(f'ZD1: {ke1_cur_vel[2]}')
            logger.info(f'X2: {ke2_cur_pos[0]}')
            logger.info(f'Y2: {ke2_cur_pos[1]}')
            logger.info(f'Z2: {ke2_cur_pos[2]}')
            logger.info(f'XD2: {ke2_cur_vel[0]}')
            logger.info(f'YD2: {ke2_cur_vel[1]}')
            logger.info(f'ZD2: {ke2_cur_vel[2]}')


    logger.info('Writing SV positions and velocities to CSV')
    keHelperFunctions.print_data_to_csv('exam3_p1_sv1.csv', p1_headers, p1_sv1_pos_and_vel_vectors)
    keHelperFunctions.print_data_to_csv('exam3_p1_sv2.csv', p1_headers, p1_sv2_pos_and_vel_vectors)
    print()

    # Problem 2
    logger.info('----- Problem 2 -----')

    # All of this is for SV 1, not necessarily SV 2. But I think in problem 3, Having this info for SV2 may be helpful

    # Report the results in the following table
    # Event         | Time | Elevation (Deg) | Azimuth (Deg)
    # Rise          |      |                 |
    # Max Elevation |      |                 |
    # Set           |      |                 |
    #
    # All calculations up to this point are in ECI/TOD reference frame

    p2_headers = ['Step', 'Seconds from Epoch', 'Time (UTC)', 'Azimuth (Degrees)', 'Elevation (Degrees)', 'Event']
    sv1_pass_data = []
    sv2_pass_data = []
    p2_sv1_table = []

    logger.info('Finding all passes for SV1...')
    for vector in p1_sv1_pos_and_vel_vectors:
        cur_pos = [vector[2], vector[3], vector[4]]
        cur_time = epoch + datetime.timedelta(seconds=vector[1])
        cur_jd = keHelperFunctions.convert_date_to_jd(cur_time.year, cur_time.month, cur_time.day, cur_time.hour, cur_time.minute, cur_time.second)

        cur_ecef_pos = keHelperFunctions.convert_eci_ecef_gha(cur_pos, cur_jd)

        cur_relative_pos = pogo_gs.compute_relative_pos(cur_ecef_pos)

        topo_relative_position = keHelperFunctions.convert_ecef_topocentric(math.radians(pogo_gs.lat), math.radians(pogo_gs.lon), cur_relative_pos)

        cur_az, cur_el = pogo_gs.determine_az_el_to_sv(topo_relative_position)

        degree_az = math.degrees(cur_az)
        degree_el = math.degrees(cur_el)

        if degree_el > 0 and degree_el < 180:
            pass_data = [vector[0], vector[1], epoch + datetime.timedelta(seconds=vector[1]), degree_az, degree_el]
            sv1_pass_data.append(pass_data)


        if vector[1] == 12800:
            logger.info('STEP 64 AZIMUTH AND ELEVATION')
            logger.info(f'AZ for SV1 at step 64: {degree_az} deg')
            logger.info(f'EL for SV1 at step 64: {degree_el} deg')
    
    logger.info('Finding all passes for SV2...')
    for vector in p1_sv2_pos_and_vel_vectors:
        cur_pos = [vector[2], vector[3], vector[4]]
        cur_time = epoch + datetime.timedelta(seconds=vector[1])
        cur_jd = keHelperFunctions.convert_date_to_jd(cur_time.year, cur_time.month, cur_time.day, cur_time.hour, cur_time.minute, cur_time.second)

        cur_ecef_pos = keHelperFunctions.convert_eci_ecef_gha(cur_pos, cur_jd)

        cur_relative_pos = pogo_gs.compute_relative_pos(cur_ecef_pos)

        topo_relative_position = keHelperFunctions.convert_ecef_topocentric(math.radians(pogo_gs.lat), math.radians(pogo_gs.lon), cur_relative_pos)

        cur_az, cur_el = pogo_gs.determine_az_el_to_sv(topo_relative_position)

        degree_az = math.degrees(cur_az)
        degree_el = math.degrees(cur_el)

        if degree_el > 0 and degree_el < 180:
            pass_data = [vector[0], vector[1], epoch + datetime.timedelta(seconds=vector[1]), degree_az, degree_el]
            sv2_pass_data.append(pass_data)
    
    logger.info('Writing SV2 Azimuth and Elevation Data to CSV')
    keHelperFunctions.print_data_to_csv('exam3_p2_sv2_az_and_el.csv', p2_headers, sv2_pass_data)

    # This part only works since there is only one pass of SV1
    start_of_pass_event = sv1_pass_data[0]
    start_of_pass_event.append('Rise')
    p2_sv1_table.append(start_of_pass_event)
    highest_el = 0
    highest_data = []

    
    logger.info('Building Problem 2 Table')
    for measurement in sv1_pass_data:
        if measurement[4] > highest_el:
            highest_el = measurement[4]
            highest_data = measurement
    
    highest_data.append('Max Elevation')
    p2_sv1_table.append(highest_data)
    end_of_pass = sv1_pass_data[-1]
    end_of_pass.append('Set')
    p2_sv1_table.append(end_of_pass)

    print(tabulate.tabulate(p2_sv1_table, p2_headers), end='\n\n')


    # Problem 3
    logger.info('----- Problem 3 -----')

    # Part c
    # Compute the ECEF position vector for the POGO ground site, r_s_vector
    logger.info('Computing ECEF coordinates for POGO')
    pogo_ecef = pogo_gs.compute_ecef_coords()
    logger.info(pogo_ecef)
    logger.info('MAGNITUDE OF ECEF VECTOR')
    logger.info(np.linalg.norm(pogo_ecef))

    # Part f
    # For each time point that Satellite 1 is in view of the POGO ground site, compute the
    # angular separation theta between satellites
    logger.info('Computing Angles between SV1 and SV2 when SV1 IS IN VIEW')
    sv1_start_pass_step = sv1_pass_data[0][0]
    sv1_end_pass_step = sv1_pass_data[-1][0]

    p3_headers = ['Step', 'Seconds since Epoch', 'Angular Difference (deg)']

    p3_separation_angles = []

    for step in range(sv1_start_pass_step, sv1_end_pass_step+1):
        p3_sv1_eci_pos = [p1_sv1_pos_and_vel_vectors[step][2], p1_sv1_pos_and_vel_vectors[step][3], p1_sv1_pos_and_vel_vectors[step][4]]
        p3_sv2_eci_pos = [p1_sv2_pos_and_vel_vectors[step][2], p1_sv2_pos_and_vel_vectors[step][3], p1_sv2_pos_and_vel_vectors[step][4]]

        p3_cur_time = epoch + datetime.timedelta(seconds=p1_sv1_pos_and_vel_vectors[step][1])
        p3_cur_jd = keHelperFunctions.convert_date_to_jd(p3_cur_time.year, p3_cur_time.month, p3_cur_time.day, p3_cur_time.hour, p3_cur_time.minute, p3_cur_time.second)

        p3_sv1_ecef_pos = keHelperFunctions.convert_eci_ecef_gha(p3_sv1_eci_pos, p3_cur_jd)
        p3_sv2_ecef_pos = keHelperFunctions.convert_eci_ecef_gha(p3_sv2_eci_pos, p3_cur_jd)

        p3_sv1_relative_pos = pogo_gs.compute_relative_pos(p3_sv1_ecef_pos)
        p3_sv2_relative_pos = pogo_gs.compute_relative_pos(p3_sv2_ecef_pos)

        p3_sv1_topo_relative_position = keHelperFunctions.convert_ecef_topocentric(math.radians(pogo_gs.lat), math.radians(pogo_gs.lon), p3_sv1_relative_pos)
        p3_sv2_topo_relative_position = keHelperFunctions.convert_ecef_topocentric(math.radians(pogo_gs.lat), math.radians(pogo_gs.lon), p3_sv2_relative_pos)

        angle_between_sv1_sv2 = pogo_gs.determine_angle_between_two_sv(p3_sv1_topo_relative_position, p3_sv2_topo_relative_position)

        separation_angle_data = [step, p1_sv1_pos_and_vel_vectors[step][1], math.degrees(angle_between_sv1_sv2)]

        p3_separation_angles.append(separation_angle_data)

        if step == 64:
            logger.info(f'STEP 64 ANGLE BETWEEN SV: {math.degrees(angle_between_sv1_sv2)}')

    keHelperFunctions.print_data_to_csv('exam3_p3_angular_separation.csv', p3_headers, p3_separation_angles)

    # Using the orbit predictions from problem 1 for SV 1 and 2, complete the below

    # Part h
    # Estimate the magnitude of the minimum separation angle between the two satellites
    logger.info('Magnitude of the minimum separation angle between SV1 and SV2')

    # Hardcoded to 97 since I know its the loweest from the CSV
    p3_sv1_eci_pos = [p1_sv1_pos_and_vel_vectors[97][2], p1_sv1_pos_and_vel_vectors[97][3], p1_sv1_pos_and_vel_vectors[97][4]]
    p3_sv2_eci_pos = [p1_sv2_pos_and_vel_vectors[97][2], p1_sv2_pos_and_vel_vectors[97][3], p1_sv2_pos_and_vel_vectors[97][4]]

    p3_cur_time = epoch + datetime.timedelta(seconds=p1_sv1_pos_and_vel_vectors[97][1])
    p3_cur_jd = keHelperFunctions.convert_date_to_jd(p3_cur_time.year, p3_cur_time.month, p3_cur_time.day, p3_cur_time.hour, p3_cur_time.minute, p3_cur_time.second)

    p3_sv1_ecef_pos = keHelperFunctions.convert_eci_ecef_gha(p3_sv1_eci_pos, p3_cur_jd)
    p3_sv2_ecef_pos = keHelperFunctions.convert_eci_ecef_gha(p3_sv2_eci_pos, p3_cur_jd)

    p3_sv1_relative_pos = pogo_gs.compute_relative_pos(p3_sv1_ecef_pos)
    p3_sv2_relative_pos = pogo_gs.compute_relative_pos(p3_sv2_ecef_pos)

    p3_sv1_topo_relative_position = keHelperFunctions.convert_ecef_topocentric(math.radians(pogo_gs.lat), math.radians(pogo_gs.lon), p3_sv1_relative_pos)
    p3_sv2_topo_relative_position = keHelperFunctions.convert_ecef_topocentric(math.radians(pogo_gs.lat), math.radians(pogo_gs.lon), p3_sv2_relative_pos)

    logger.info(f'SV1 MINIMUM RFI TOPO COORDINATES: {p3_sv1_topo_relative_position}')
    logger.info(f'SV2 MINIMUM RFI TOPO COORDINATES: {p3_sv2_topo_relative_position}')

    r3_vector = np.add(p3_sv1_topo_relative_position, p3_sv2_topo_relative_position)

    r3_norm = np.linalg.norm(r3_vector)
    logger.info(f'MAGNITUDE AT LOWEST ANGLE AT STEP 97: {r3_norm}')

    logger.info(f'RFI START TIME {epoch + datetime.timedelta(seconds=19100)}')
    logger.info(f'RFI END TIME {epoch + datetime.timedelta(seconds=19900)}')

    logger.info(f'RFI MINIMUM TIME {epoch + datetime.timedelta(seconds=19400)}')




if __name__ == '__main__':
    main()
