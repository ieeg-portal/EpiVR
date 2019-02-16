""" Unit test for get_img_manifest.py
        Description:
            This function takes in a file path containing all images for a
            given subject and outputs a csv file listing the date, contrast,
            and filepath. The images should be organized into directories
            where the directory name is the date the image was collected.
            There must be a subdirectory 'nii' containing all the nifti
            files. Note that this function only searches for .nii and
            .nii.gz files.

        Args:
            input_path (str): The path location containing the Nifti files
            based on /gdrive/public/USERS/lkini/3T_Subjects/HUP*/img/

        Returns:
            csv_file (object): date, contrast, full_path_location

"""

import os
import csv
import get_img_manifest

get_img_manifest.get_img_manifest(
    r'/gdrive/public/USERS/lkini/3T_Subjects/HUP116/img')
csvfile = open("image_manifest.csv", "r")
reader = csv.reader(csvfile)
for row in reader:
    print(row)

# Clean up test assets generated
os.system('rm image_manifest.csv')
