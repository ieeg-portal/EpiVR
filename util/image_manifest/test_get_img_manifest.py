""" Unit test for get_img_manifest.py
"""

import os
import csv
import get_img_manifest

if not os.path.isdir('tmp'):
    os.makedirs('tmp')
get_img_manifest.get_img_manifest(
    r'/gdrive/public/USERS/lkini/3T_Subjects/HUP116/img',
    'tmp')
csvfile = open("image_manifest.csv", "r")
reader = csv.reader(csvfile)
for row in reader:
    print(row)

# Clean up test assets generated
os.system('rm tmp/image_manifest.csv')
