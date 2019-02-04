"""DATAjson_parse.py

This function is designed to parse the DATA.json file and generate separate
annotations.json and clinical.json files for each patient. Currently the
code makes a directory for each patient and outputs the files there. You may
want to update to change the output methodology.

Example:
    python DATAjson_parse.py
    NOTE: must be in directory containing DATA.json

Outputs:
    current_dir
        patient_dir (patient_dir is named based on patient name in DATA.json)
            annotations.json
            clinical.json

Todo:
    * Update output for borel usage

Alterations:
    TCA 2/1/19 - initialized
"""
import json
import os

# load in the data file
with open('DATA.json') as f:
    data = json.load(f)

# loop through each pateint and generate annotations.json & clinical.json
for patient in data['PATIENTS']:
    try:  # try/except allows code to run dispite not processing CHOP data

        ann = dict()  # setup annotation variable and populate
        ann['patient'] = patient
        ann['Ictal'] = {}
        ann['Interictal'] = {}
        current_patient = data['PATIENTS'][patient]['Events']['Ictal']

        # convert dict_keys into an array of ints (sequential events)
        my_events = []
        for event in sorted(current_patient.viewkeys()):
            my_events.append(int(event))
        my_events = sorted(my_events)

        # loop through all events listed
        for event in my_events:
            current_event = current_patient[str(event)]
            subfields = {
                    'FILE': current_event['FILE'],
                    'SeizureType': current_event['SeizureType'],
                    'EMU_Report_Event_Number': current_event['EMU_Report_Event\
                        _Number'],
                    'SeizureEEC': current_event['SeizureEEC'],
                    'SeizureUEO': current_event['SeizureUEO'],
                    'SeizureEnd': current_event['SeizureEnd'],
                    'SEIZURE_ONSET_ELECTRODES': current_event['SEIZURE_ONSET_E\
                        LECTRODES'],
                    'IGNORE_ELECTRODES': data['PATIENTS'][patient]['IGNORE_ELE\
                        CTRODES'],
                    'SeizurePhenotype': current_event['SeizurePhenotype'],
                }
            if event == 1000:
                ann['Interictal']['1'] = subfields
            else:
                ann['Ictal'][str(event)] = subfields

        # set up clinical variable
        clin = {
                'PATIENT': patient,
                'Outcome': data['PATIENTS'][patient]['Outcome'],
                'Sex': data['PATIENTS'][patient]['Sex'],
                'AgeOnset': data['PATIENTS'][patient]['AgeOnset'],
                'AgeSurgery': data['PATIENTS'][patient]['AgeSurgery'],
                'SeizureOnset': data['PATIENTS'][patient]['SeizureOnset'],
                'Lesion Status': data['PATIENTS'][patient]['Lesion Status'],
                'Pathology': data['PATIENTS'][patient]['Pathology'],
                'Resection Type': data['PATIENTS'][patient]['Resection Type'],
                'Implant Type': 'ECoG'
            }

        # Write out .json files
        os.mkdir(patient)
        outname = os.path.join(patient, 'annotations.json')
        with open(outname, 'w') as outfile:
            json.dump(ann, outfile, indent=4)
        outname = os.path.join(patient, 'clinical.json')
        with open(outname, 'w') as outfile:
            json.dump(clin, outfile, indent=4)

    # This code does not process data from CHOP
    except:
        print('Did not process ' + patient)
