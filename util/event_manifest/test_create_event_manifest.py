""" Unit test for create_event_manifest.py
"""

import os
import csv
import create_event_manifest

if not os.path.isdir('tmp'):
    os.makedirs('tmp')
create_event_manifest.create_event_manifest(
    r'/gdrive/public/USERS/lkini/3T_Subjects/HUP116/img',
    'tmp')
csvfile = open("tmp/event_manifest.csv", "r")
reader = csv.reader(csvfile)
for row in reader:
    print(row)

# Clean up test assets generated
os.system('rm -r tmp/')
