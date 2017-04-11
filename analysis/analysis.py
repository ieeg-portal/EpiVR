#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

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
        starting_null_id = int(sys.argv[4])-1
    except IndexError:
        starting_null_id = 0

    # Check if connectivity already exists
    count_conn = 0
    for k,v in data["PATIENTS"][patient_id]['Events']['Ictal'].items():
        try:
            if(v["STATUS"] == 'ALL_DROPOUT'):
                continue
            else:
                count_conn += 1
        except Exception:
            count_conn += 1

    # if(not(os.path.isfile('%s/%s/aim3/%s.Ictal.3.multiband.npz'%(os.path.expanduser(data['COMP_DIR']),patient_id,patient_id)))):
    if(len(glob.glob('%s/%s/aim3/*multiband*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) != count_conn):
        # Compute multi band connectivity and store adjacency matricies
        blah
        pass
        # compute_multiband_connectivity(patient_id, epoch_length)

    # Keep dilating and eroding resection zone to capture percentage of network nodes
    # for dilate_radius in [0, -5, 5, 10, -10, 15, -15, 20, -20]:
    for dilate_radius in [0]:
        # Check if already exists
        unique_idx = []
        for fn in glob.glob(os.path.expanduser('%s/%s/aim3/*pipedef*'%(data['COMP_DIR'],patient_id))):
            pipedef = json.load(open(fn,'r'))
            if(pipedef['fconn'] == 'multiband+broadband' and pipedef['dilate_radius'] == dilate_radius):
                uid = fn.split('.')[-3]
                event_type = fn.split('.')[1]
                event_id = fn.split('.')[2]
                unique_idx.append((uid,event_type,event_id))
        if(len(unique_idx) == 0):
            # Compute virtual resection to get c_{res}(t)
            unique_idx = virtual_resection(patient_id, dilate_radius, data)

        # Compute null models for virtual resection to get c_{null}(t)
        if(run_null):
            for unique_id, event_type, event_id in unique_idx:
                print 'Running null for patient %s, event id %s, unique_id %s'%(patient_id, event_id, unique_id)
                start = time.time()
                null_virtual_resection(patient_id, unique_id, event_type, event_id, dilate_radius, data, starting_null_id=starting_null_id)
                print 'Time Took for 1 event for null: '
                print time.time() - start

        # Generate figures for result
        # unique_idx = [('9b266cab-c6cd-4b6d-88aa-098f1f8d9a26','Ictal','1'),('9b266cab-c6cd-4b6d-88aa-098f1f8d9a26','Ictal','2'),('9b266cab-c6cd-4b6d-88aa-098f1f8d9a26','Ictal','3')]
        for unique_id, event_type, event_id in unique_idx:
            # plot_experiment(patient_id,unique_id)
            print unique_id, event_type, event_id
            pass
