"""create_event_manifest.py
        Description:
            This function takes in a file path containing all images for a
            given subject and outputs a csv file listing the date, event,
            and event description. This represents major timeline of events
            such as date of implant, explant, resection, etc.

        Args:
            input_path (str): The path location containing the csv files
                based on /gdrive/public/USERS/lkini/3T_Subjects/HUP*/img/
            output_path (str): The path location in which to output the event manifest file
                based on /gdrive/public/DATA/Human_Data/VirtualResection

        Returns:
            csv_file (object): date, event, event_description

"""
import os
import sys
import csv


def create_event_manifest(input_path, output_path='/gdrive/public/DATA/Human_Data/VirtualResection/'):
    if not os.path.exists(output_path) or not os.path.exists(input_path):
        raise OSError('Make sure input path %s and output path %s exist' % (input_path, output_path))

    # Get patient ID
    patient_id = input_path.split('/')[-2]

    # open a csv file for writing the data out
    with open(os.path.join(output_path, patient_id, 'event_manifest.csv'), mode='w') as img_mani:
        img_mani = csv.writer(img_mani, delimiter=',', quotechar='"')

        img_mani.writerow(["Date", "Event", "Event Description"])

        if os.path.isfile(os.path.join(input_path, '%s_img_t2.csv' % patient_id)):
            lines = open(os.path.join(input_path, '%s_img_t2.csv' % patient_id), 'r').readlines()
            assert len(lines) == 2
            header = lines[0].replace('\n', '').split(',')
            dates = lines[1].replace('\n', '').split(',')
            event_data = {}
            for ii, event in enumerate(header):
                event_data[event] = dates[ii]
            for event_key, event, event_description in [
                ('PATIENT_PREOP_T1_DATE', 'Reference', 'Reference Pre-Op T1 date'),
                ('PATIENT_POSTIMPLANT_CT_DATE', 'Implant Date', 'Date of implant surgery based on when CT was taken'),
                ('PATIENT_RESECTION_DATE', 'Resection Date', '')]:
                try:
                    event_date = event_data[event]
                except KeyError:
                    continue
                img_mani.writerow([event_date, event, event_description])
    return img_mani


# Gather our code in a main() function
def main():
    try:
        create_event_manifest(sys.argv[1], sys.argv[2])
    except IndexError:
        create_event_manifest(sys.argv[1])


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
