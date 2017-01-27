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
        dilate_radius = int(sys.argv[2]) # voxels
        epoch_length = int(sys.argv[3]) # seconds
    except IndexError:
        dilate_radius = 0
        epoch_length = 1
    except TypeError:
        print 'Please enter an appropriate PATIENT_ID and integers for dilation/erosion radius and epoch length as inputs.'
        raise

    # Compute multi band connectivity and store adjacency matricies
    # computed_file_idx = compute_multiband_connectivity(patient_id, epoch_length)

    # Compute virtual resection to get c_{res}(t)
    # unique_idx = virtual_resection(patient_id, dilate_radius, data)

    # Compute null models for virtual resection to get c_{null}(t)
    # for unique_id, event_type, event_id in unique_idx:
    #     null_virtual_resection(patient_id, unique_id, event_type, event_id, dilate_radius, data)

    # Generate figures for result
    unique_id = 'b50a1916-d49f-4b05-b879-d22ad266f474'
    plot_experiment(patient_id,unique_id)
