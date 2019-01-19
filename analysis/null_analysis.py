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
        starting_null_id = int(sys.argv[3])-1
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

    # NULL MODEL
    # Check if all NULL node-level virtual resections have been computed
    if(len(glob.glob('%s/%s/aim3/*nodenull*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) != count_conn):
        for event_type, events in data['PATIENTS'][patient_id]['Events'].items():
            for event_id in events.keys():
                try:
                    if(events[event_id]['STATUS'] == 'ALL_DROPOUT'):
                        continue # unusable clip
                except KeyError:
                    pass
                # Compute node-level null models of virtual resection
                null_nodal_virtual_resection(patient_id, event_type, event_id, data=data)

    # Keep dilating and eroding resection zone to capture percentage of network nodes
    for dilate_radius in [0, -5, 5, 10, -10, 15, -15, 20, -20]:
        ## NULL MODEL
        # Compute null models for virtual resection to get c_{null}(t)
        if(run_null):
            for unique_id, event_type, event_id in unique_idx:
                print 'Running null for patient %s, event id %s, unique_id %s'%(patient_id, event_id, unique_id)
                start = time.time()
                null_virtual_resection(patient_id, unique_id, event_type, event_id, dilate_radius, data, starting_null_id=starting_null_id)
                print 'Time Took for 1 event for null: '
                print time.time() - start

