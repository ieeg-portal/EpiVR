#!/usr/bin/python
'''
Analysis main script to build Virtual Resection application pipeline.
'''

import os
import sys

from util import *
import json


if __name__ == '__main__':
    try:
        PATIENT_ID = sys.argv[1]
        get_resected_electrodes(PATIENT_ID)
    except Exception:
        raise
