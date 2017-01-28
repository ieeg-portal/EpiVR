#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import sys

from util import *
from util_connectivity import *
from util_virtual_resection import *

if __name__ == '__main__':
    # Get arguments
    try:
        patient_id = sys.argv[1]
        epoch_length = int(sys.argv[2]) # seconds
    except IndexError:
        epoch_length = 1
    except TypeError:
        print 'Please enter an appropriate PATIENT_ID.'
        raise

    # Compute multi band connectivity and store adjacency matricies
    compute_multiband_connectivity(patient_id, epoch_length)

    for dilate_radius in np.arange(-8,8):
        # Compute virtual resection to get c_{res}(t)
        unique_idx = virtual_resection(patient_id, dilate_radius, data)

        # Compute null models for virtual resection to get c_{null}(t)
        for unique_id, event_type, event_id in unique_idx:
            null_virtual_resection(patient_id, unique_id, event_type, event_id, dilate_radius, data)

        # Generate figures for result
        for unique_id, event_type, event_id in unique_idx:
            plot_experiment(patient_id,unique_id)
