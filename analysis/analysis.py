#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import sys
from util import *

if __name__ == '__main__':
    # Get arguments
    PATIENT_ID = sys.argv[1]
    try:
        dilate_radius = int(sys.argv[2])
    except IndexError:
        dilate_radius = 0
    except TypeError:
        print 'Please enter an integer for dilation/erosion radius.'
        raise

    # Write to CSV file
    write_resected_electrodes(PATIENT_ID,dilate_radius)
