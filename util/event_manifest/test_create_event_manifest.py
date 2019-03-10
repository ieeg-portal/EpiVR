""" Unit test for create_event_manifest.py
"""

import os
import csv
import create_event_manifest
import pytest

if not os.path.isdir('tmp'):
    os.makedirs('tmp')
    os.makedirs('tmp/HUP116')
    os.makedirs('tmp/sample_pt')
create_event_manifest.create_event_manifest(
    r'/gdrive/public/USERS/lkini/3T_Subjects/HUP116/img',
    'tmp')
csvfile = open("tmp/HUP116/event_manifest.csv", "r")
reader = csv.reader(csvfile)
for row in reader:
    print(row)

with pytest.raises(OSError):
    create_event_manifest.create_event_manifest(
        r'/blah',
        'tmp')

create_event_manifest.create_event_manifest(
    r'tests/data/sample_pt/img',
    'tmp')

# Clean up test assets generated
os.system('rm -r tmp/')
