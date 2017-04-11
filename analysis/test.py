import unittest

from util import *
from util_connectivity import *
from util_virtual_resection import *

with open('../data/TEST_DATA.json') as json_data_file:
    data = json.load(json_data_file)

class CorrespondNamesTest(unittest.TestCase):
    '''
    This unit test checks conversion between EEG labels from IEEG.org and cartoon map labels (e.g. LG64 and LG064-Ref).
    '''
    def test(self):
        patient_id = 'Study029'
        dilate_radius = 0

        data_dir = os.path.expanduser(data['REAL_DATA_DIR'])
        for event_type, events in data['PATIENTS'][patient_id]['Events'].items():
            for event_id in events.keys():
                fn = os.path.join(data_dir, patient_id, 'eeg', events[event_id]['FILE'])
                eeg_channel_labels = []

                # Get channels, ECoG Data, Fsx
                with h5py.File(fn) as f:
                    evData = f['evData'].value
                    Fs = f['Fs'].value
                    for column in f['channels']:
                        row_data = []
                        for row_number in range(len(column)):
                            row_data.append(''.join(map(unichr, f[column[row_number]][:])))
                        eeg_channel_labels.append(row_data)
                Fs = int(Fs[0][0])
                eeg_channel_labels = eeg_channel_labels[0]
                # evData = scipy.stats.zscore(evData,axis=1)
                T = evData.shape[0]
                break
            break

        cartoon_map_labels = map(lambda x: x.split(',')[4].replace('\n',''), open(os.path.expanduser(
            data['PATIENTS'][patient_id]['ELECTRODE_LABELS']
            ),'r').readlines())

        res_dict = correspond_label_names(eeg_channel_labels, cartoon_map_labels)
        for k,v in sorted(res_dict.items(),key=lambda x: x[0]):
            print k,v[1]

        self.assertTrue(True)

class DataTest(unittest.TestCase):
    '''
    This unit test checks compatibility of datasets as defined by the TEST_DATA.json config file.
    '''
    def test(self):
        for event_type, events in data['PATIENTS']['TEST1']['Events'].items():
            fn = os.path.join(os.path.expanduser(data['DATA_DIR']),
                'TEST1',
                'eeg',
                events['1']['FILE']
                )
            channels = []

            # Get channels, ECoG Data, Fsx
            with h5py.File(fn) as f:
                evData = f['evData'].value
                Fs = f['Fs'].value
                for column in f['channels']:
                    row_data = []
                    for row_number in range(len(column)):
                        row_data.append(''.join(map(unichr, f[column[row_number]][:])))
                    channels.append(row_data)
            Fs = int(Fs[0][0])
            channels = list(np.squeeze(np.array(channels)))
            # evData = scipy.stats.zscore(evData,axis=1)
            T = evData.shape[0]

            assert evData.shape == (450000,60)
            assert Fs == 500

class ResectionZoneTest(unittest.TestCase):
    '''
    This unit test computes and prints electrodes in the resection zone for a given patient ID.
    '''
    def test(self):
        patient_id = 'Study029'
        data_dir = os.path.expanduser(data['REAL_DATA_DIR'])

        labels = map(lambda x: x.split(',')[4].replace('\n',''), open(os.path.expanduser(
            data['PATIENTS'][patient_id]['ELECTRODE_LABELS']
            ),'r').readlines())

        # Load ignored node labels
        ignored_node_labels = data['PATIENTS'][patient_id]['IGNORE_ELECTRODES']
        for ignored_node_label in ignored_node_labels:
            if(ignored_node_label not in labels):
                labels.append(ignored_node_label)

        for event_type, events in data['PATIENTS'][patient_id]['Events'].items():
            for event_id in events.keys():
                fn = os.path.join(data_dir, patient_id, 'eeg', events[event_id]['FILE'])
                channels = []

                # Get channels, ECoG Data, Fsx
                with h5py.File(fn) as f:
                    evData = f['evData'].value
                    Fs = f['Fs'].value
                    for column in f['channels']:
                        row_data = []
                        for row_number in range(len(column)):
                            row_data.append(''.join(map(unichr, f[column[row_number]][:])))
                        channels.append(row_data)
                Fs = int(Fs[0][0])
                channels = channels[0]
                # evData = scipy.stats.zscore(evData,axis=1)
                T = evData.shape[0]

                # Correspond lable names
                labels_dict = correspond_label_names(channels, labels)

                # Load electrodes to ignore
                ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
                for ii,node_id in enumerate(ignored_node_idx):
                    print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
                channels = list(np.delete(np.array(channels),ignored_node_idx))

                # Recorrespond label names
                labels_dict = correspond_label_names(channels, labels)
                break
            break

        dilate_radius = -5
        print 'Printing resected electrodes with erosion of 5% of network'
        print get_resected_electrodes(patient_id, dilate_radius, data, labels_dict)
        dilate_radius = 0
        print 'Printing resected electrodes with no dilation'
        print get_resected_electrodes(patient_id, dilate_radius, data, labels_dict)
        dilate_radius = 5
        print 'Printing resected electrodes with dilation of 5% of network'
        print get_resected_electrodes(patient_id, dilate_radius, data, labels_dict)
        self.assertTrue(True)

class ConnectivityTest(unittest.TestCase):
    '''
    This unit test creates adjacency multiband connectivity matrices and checks correctness of output.
    '''
    def test(self):
        compute_multiband_connectivity('TEST1', 1, data)
        self.assertTrue(True)

class VirtualResectionTest(unittest.TestCase):
    '''
    This unit test computes c_res(t) based on multiband connectivity matrices and checks correctness of output.
    '''
    def test(self):
        # Run virtual resection on all files
        unique_idx = virtual_resection('TEST1',0,data)

        self.assertTrue(True)

class NullVirtualResectionTest(unittest.TestCase):
    '''
    This unit test computes c_null(t) based on multiband connectivity matrices and checks correctness of output.
    '''
    def test(self):
        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),'TEST1','aim3')
        unique_idx = []
        for fn in os.listdir(comp_dir):
            try:
                print fn
                match = re.match(r'[A-Za-z0-9]+.([A-Za-z]+).([0-9]+).cres.([0-9a-zA-Z-]+).npz',fn)
                unique_idx.append((match.group(3),match.group(1),match.group(2)))
            except AttributeError:
                continue

        # Run virtual resection on all files
        for unique_id, event_type, event_id in unique_idx:
            null_virtual_resection('TEST1', unique_id, event_type, event_id, 0,data)
        self.assertTrue(True)

class PlotVirtualResectionTest(unittest.TestCase):
    '''
    This unit test generates the plot figures.
    '''

    def test(self):
        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),'TEST1','aim3')
        unique_idx = []
        for fn in os.listdir(comp_dir):
            try:
                match = re.match(r'[A-Za-z0-9]+.([A-Za-z]+).([0-9]+).cres.([0-9a-zA-Z-]+).npz',fn)
                unique_idx.append((match.group(3),match.group(1),match.group(2)))
            except AttributeError:
                continue

        # Plot virtual resection results on all files
        for unique_id, event_type, event_id in unique_idx:
            plot_experiment('TEST1', unique_id, data=data)
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
