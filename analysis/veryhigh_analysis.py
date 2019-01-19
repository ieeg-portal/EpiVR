#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import sys
import glob
import json
import time

from util import *
from veryhigh_util_connectivity import *
from veryhigh_util_virtual_resection import *

if __name__ == '__main__':
    # Get arguments
    try:
        patient_id = sys.argv[1]
    except TypeError:
        print 'Please enter an appropriate PATIENT_ID.'
        raise
    try:
        epoch_length = int(sys.argv[2]) # seconds
    except IndexError:
        epoch_length = 1
    try:
        run_null = int(sys.argv[3])
    except IndexError:
        run_null = 0

    try:
        force_flag = int(sys.argv[4])
    except IndexError:
        force_flag = 0

    try:
        starting_null_id = int(sys.argv[5])-1
    except IndexError:
        starting_null_id = 0

    # Compute multi band connectivity and store adjacency matricies
    # blah
    # pass
    compute_multiband_connectivity(patient_id, epoch_length)

    # Compute node-level virtual resection
    nodal_virtual_resection(patient_id, data=data)

    # Keep dilating and eroding resection zone to capture percentage of network nodes
    for dilate_radius in [0, -5, 5, 10, -10, 15, -15, 20, -20]:
        print dilate_radius
        # Compute virtual resection to get c_{res}(t)
        unique_idx = virtual_resection(patient_id, dilate_radius, data)
