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

    # Problem 3
    print('----- Problem 3 -----')

    # Problem 4
    print('----- Problem 4 -----')

    # Problem 5
    print('----- Problem 5 -----')

    # Problem 6
    print('----- Problem 6 -----')

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
