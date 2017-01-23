#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import sys
from util import *

if __name__ == '__main__':
    # Get arguments
    try:
        patient_id = sys.argv[1]
        dilate_radius = int(sys.argv[2])
    except IndexError:
        dilate_radius = 0
    except TypeError:
        print 'Please enter an appropriate PATIENT_ID and an integer for dilation/erosion radius as inputs.'
        raise

    # Generate list of resected electrodes and write to CSV file
    write_resected_electrodes(patient_id, dilate_radius)

    # Load resected electrodes
    # Load electrodes to ignore
    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
    # For each clip, broken up in to epochs of length epoch_length
        # Compute multiband connectivity
        # Compute util.region_control == cres(t)
        # Save with appropriate name
    # Load
