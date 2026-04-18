from keplarianElements import KeplerianElements
from groundSite import GroundSite
import keHelperFunctions

import math
import datetime
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

    chief = KeplerianElements(vector_data['vectors'][f'vector1']['x_pos'],
                              vector_data['vectors']['vector1']['y_pos'],
                              vector_data['vectors']['vector1']['z_pos'],
                              vector_data['vectors']['vector1']['x_velocity'],
                              vector_data['vectors']['vector1']['y_velocity'],
                              vector_data['vectors']['vector1']['z_velocity'])

    deputy = KeplerianElements(vector_data['vectors'][f'vector2']['x_pos'],
                               vector_data['vectors']['vector2']['y_pos'],
                               vector_data['vectors']['vector2']['z_pos'],
                               vector_data['vectors']['vector2']['x_velocity'],
                               vector_data['vectors']['vector2']['y_velocity'],
                               vector_data['vectors']['vector2']['z_velocity'])


    # Epoch provided: 04/09/2024 00:00:00
    epoch = datetime.datetime(2024, 4, 9, 0, 0, 0)


    # Problem 1
    logger.info('----- Problem 1 -----')
    logger.info('Chief Keplarian Elements')
    chief.print_ke()
    print()

    logger.info('Deputy Keplarian Elements')
    deputy.print_ke()
    print()

    p1_r_vector_diff = chief.r_vector - deputy.r_vector
    p1_ke_diff = keHelperFunctions.compute_ke_diff(chief, deputy)
    logger.info('KE DIFF')
    logger.info(f'ECI X, Y, Z: {p1_ke_diff['eci_pos']}')
    logger.info(f'ECI XD, YD, ZD: {p1_ke_diff['eci_vel']}')
    logger.info(f'SMA: {p1_ke_diff['sma']}')
    logger.info(f'Eccentricity: {p1_ke_diff['eccentricity']}')
    logger.info(f'Inclination: {p1_ke_diff['inclination']}')
    logger.info(f'raan: {p1_ke_diff['raan']}')
    logger.info(f'AoP: {p1_ke_diff['aop']}')
    logger.info(f'TA: {p1_ke_diff['nu']}')
    logger.info(f'MA: {p1_ke_diff['mean_anomaly']}')
    logger.info(f'ARGLAT {p1_ke_diff['arglat']}')
    print()

    # Compute UVW differences using transport theorem (Slide 39)
    chief_uvw = keHelperFunctions.convert_coordinates_to_uvw(chief.r_vector, chief.r_dot_vector)
    chief_uvw_vel = keHelperFunctions.convert_vel_to_uvw(chief.r_vector, chief.r_dot_vector)

    deputy_uvw = keHelperFunctions.convert_coordinates_to_uvw(deputy.r_vector, deputy.r_dot_vector)
    deputy_uvw_vel = keHelperFunctions.convert_vel_to_uvw(deputy.r_vector, deputy.r_dot_vector)

    p1_rss_pos = keHelperFunctions.compute_rss_diff(chief.r_vector, deputy.r_vector)
    p1_rss_vel = keHelperFunctions.compute_rss_diff(chief.r_dot_vector, deputy.r_dot_vector)

    logger.info(f'RSS POS DIFF: {p1_rss_pos}')
    logger.info(f"RSS VEL DIFF: {p1_rss_vel}")
    print()
    logger.info(f'Chief UVW Coordinates: {chief_uvw}')
    logger.info(f'Chief UVW Velocity: {chief_uvw_vel}')
    print()
    logger.info(f'Deputy UVW Coordinates: {deputy_uvw}')
    logger.info(f'Deputy UVW Velocity: {deputy_uvw_vel}')
    print()

    # Compute Mean Motion for the chief satellite
    n = chief.mean_motion
    logger.info(f'Chief Mean Motion: {n}')

    # Using that n, evaluate the HCW matrix for t=1 second
    p1_transformation_matrix = keHelperFunctions.compute_hcw_matrix(chief_uvw, chief_uvw_vel, n, 1)
    logger.info('HCW Transformation Matrix')
    pprint(p1_transformation_matrix)
    hcw_pos, hcw_vel = keHelperFunctions.compute_hcw(chief_uvw, chief_uvw_vel, n, 1)
    print()

    logger.info(f'CHIEF HCW POS: {hcw_pos}')
    logger.info(f'CHIEF HCW VEL: {hcw_vel}')


    # Propagate Chief and Deputy for 6000 seconds at 60 second step size
    # Compute and plot UVW position difference at each step
    p1_chief_cur_pos = chief.r_vector
    p1_chief_cur_vel = chief.r_dot_vector

    p1_deputy_cur_pos = deputy.r_vector
    p1_deupty_cur_vel = deputy.r_dot_vector

    uvw_relative_pos_x_data = []
    uvw_relative_pos_y_data = []
    uvw_relative_pos_csv_data = []


    # Calulate the first UVW Diff outside of loop
    d_eci = p1_deputy_cur_pos - p1_chief_cur_pos
    # d_dot_eci = p1_deupty_cur_vel - p1_chief_cur_vel
    # angular_velocity = (np.cross(d_eci, d_dot_eci))/math.pow(np.linalg.norm(d_eci), 2)
    cur_uvw_transformation_matrix = keHelperFunctions.get_uvw_transformation_matrix(p1_chief_cur_pos, p1_chief_cur_vel)
    d_uvw = np.dot(cur_uvw_transformation_matrix, d_eci)
    # d_dot_uvw = np.dot(cur_uvw_transformation_matrix, np.subtract(d_dot_eci, np.cross(angular_velocity, d_eci)))

    uvw_relative_pos_csv_data.append([d_uvw[1], d_uvw[0]])
    uvw_relative_pos_x_data.append(d_uvw[1])
    uvw_relative_pos_y_data.append(d_uvw[0])


    for i in range(1, 100):
        chief_new_pos, chief_new_vel = keHelperFunctions.keplarian_rk4(p1_chief_cur_pos, p1_chief_cur_vel, 60, chief.mu)
        deputy_new_pos, deputy_new_vel = keHelperFunctions.keplarian_rk4(p1_deputy_cur_pos, p1_deupty_cur_vel, 60, deputy.mu)

        # Slide 46
        d_eci = deputy_new_pos - chief_new_pos
        # d_dot_eci = deputy_new_vel - chief_new_vel

        # angular_velocity = (np.cross(d_eci, d_dot_eci))/math.pow(np.linalg.norm(d_eci), 2)

        cur_uvw_transformation_matrix = keHelperFunctions.get_uvw_transformation_matrix(chief_new_pos, chief_new_vel)    

        d_uvw = np.dot(cur_uvw_transformation_matrix, d_eci)

        # d_dot_uvw = np.dot(cur_uvw_transformation_matrix, np.subtract(d_dot_eci, np.cross(angular_velocity, d_eci)))

        uvw_relative_pos_csv_data.append([d_uvw[1], d_uvw[0]])
        uvw_relative_pos_x_data.append(d_uvw[1])
        uvw_relative_pos_y_data.append(d_uvw[0])

        p1_chief_cur_pos = chief_new_pos
        p1_chief_cur_vel = chief_new_vel

        p1_deputy_cur_pos = deputy_new_pos
        p1_deupty_cur_vel = deputy_new_vel


    logger.info('Graphing UVW data')
    plt.title('Deputy to Chief Relative Position')
    plt.xlabel('RK4 Prop')
    plt.ylabel('dVp')
    plt.plot(uvw_relative_pos_x_data, uvw_relative_pos_y_data)
    plt.grid()
    plt.savefig('uvw_diff.png')
    plt.clf()

    logger.info('Writing UVW data to CSV')
    keHelperFunctions.print_data_to_csv('problem1_uvw_rk4_difference.csv', ['dVp', 'RK4 Prop'], uvw_relative_pos_csv_data)


    # Problem 2
    logger.info('----- Problem 2 -----')

    csv_data = keHelperFunctions.read_in_csv('Homework 13 Problem 2.csv')

    # 



if __name__ == '__main__':
    main()
