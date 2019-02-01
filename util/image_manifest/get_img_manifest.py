"""get_img_manifest.py
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
import sys
import csv


def get_img_manifest(input_path):
    # open a csv file for writing the data out
    with open('image_manifest.csv', mode='w') as img_mani:
        img_mani = csv.writer(img_mani, delimiter=',', quotechar='"')

        for date in os.listdir(input_path):  # loop through dir
            if os.path.isdir(os.path.join(input_path, date)):  # only if dir
                datepath = os.path.join(input_path, date)  # get full path
                datepath = os.path.join(datepath, 'nii')  # add nii subdir
                for file in os.listdir(datepath):  # get all images present
                    filepath = os.path.join(datepath, file)
                    # only get nifti files
                    if file.endswith('.nii.gz') or file.endswith('.nii'):
                        filepath_components = file.split('_')  # use _ as delimiter
                        modality_type = filepath_components[1:-1]  # remove ends
                        modality_type = '_'.join(modality_type)  # rejoin w/o ends
                        if modality_type.startswith('FLAIR') or modality_type.startswith('T1') or modality_type.startswith('T2') or modality_type.startswith('DWI'):
                            modality = 'MRI'
                            modality_type = modality_type.split('_')[0]
                        elif modality_type.startswith('PET'):
                            modality = 'PET'
                            modality_type = 'FDG' # assume all lohiths imaging data has only FDG PET
                        elif modality_type.startswith('CT'):
                            modality = 'CT'
                            modality_type = 'TO BE MANUALLY FILLED OUT AT A LATER TIME'
                        img_mani.writerow([date, modality, modality_type, filepath])
    return img_mani


# Gather our code in a main() function
def main():
    get_img_manifest(sys.argv[1])


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
