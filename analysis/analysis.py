#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import sys
from util import *

if __name__ == '__main__':
    try:
        PATIENT_ID = sys.argv[1]
        dilate_radius = 1
        if dilate_radius > 0:
            suffix = 'dilate_{0:d}'.format(dilate_radius)
        elif dilate_radius < 0:
            suffix = 'erode_{0:d}'.format(dilate_radius)
        else:
            suffix = '0'
        open('../data/%s_resected_electrodes_%s.csv' % (PATIENT_ID, suffix), 'w').write(
            '\n'.join(map(lambda x: ','.join(x),
                          get_resected_electrodes(PATIENT_ID, dilate_radius))))
    except Exception:
        raise
