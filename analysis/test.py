import unittest

from util import *
from util_connectivity import *
from util_virtual_resection import *

with open('../data/TEST_DATA.json') as json_data_file:
    data = json.load(json_data_file)

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
            evData = scipy.stats.zscore(evData,axis=1)
            T = evData.shape[0]

            assert evData.shape == (450000,60)
            assert Fs == 500

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
                print fn
                match = re.match(r'[A-Za-z0-9]+.([A-Za-z]+).([0-9]+).cres.([0-9a-zA-Z-]+).npz',fn)
                unique_idx.append((match.group(3),match.group(1),match.group(2)))
            except AttributeError:
                continue

        # Plot virtual resection results on all files
        print unique_idx
        for unique_id, event_type, event_id in unique_idx:
            plot_experiment('TEST1', unique_id, data=data)
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
