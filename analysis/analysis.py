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
        force_flag = int(sys.argv[3])
    except IndexError:
        force_flag = 0

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

    # Check if all connectivity adjacency matrices have been computed
    if(len(glob.glob('%s/%s/aim3/*multiband*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) != count_conn) or force_flag:
        # Compute multi band connectivity and store adjacency matricies
        # blah
        # pass
        compute_multiband_connectivity(patient_id, epoch_length)

    # Check if all node-level virtual resections have been computed
    if(len(glob.glob('%s/%s/aim3/*noderes*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) != count_conn) or force_flag:
        # Compute node-level virtual resection
        nodal_virtual_resection(patient_id, data=data)

    # Keep dilating and eroding resection zone to capture percentage of network nodes
    # for dilate_radius in [0, -5, 5, 10, -10, 15, -15, 20, -20]:
    for dilate_radius in [0]:
        # Check if already exists
        unique_idx = []
        for fn in glob.glob(os.path.expanduser('%s/%s/aim3/*pipedef*'%(data['COMP_DIR'],patient_id))):
            pipedef = json.load(open(fn,'r'))
            try:
                if(pipedef['fconn'] == 'multiband+broadband' and pipedef['dilate_radius'] == dilate_radius):
                    uid = fn.split('.')[-3]
                    event_type = fn.split('.')[1]
                    event_id = fn.split('.')[2]
                    unique_idx.append((uid,event_type,event_id))
            except KeyError:
                continue
        if(len(unique_idx) == 0) or force_flag or len(glob.glob('%s/%s/aim3/*cres*1000*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) == 0:
        # if dilate_radius < 0:
            print dilate_radius
            # Compute virtual resection to get c_{res}(t)
            unique_idx = virtual_resection(patient_id, dilate_radius, data)

    # Compute control centrality of SOZ
    # Check if already exists
    if not len(glob.glob(os.path.expanduser('%s/%s/aim3/*sozres*'%(data['COMP_DIR'],patient_id)))) or force_flag:
        # Compute virtual resection to get c_{res}(t)
        unique_idx = soz_virtual_resection(patient_id, data)


