#!/usr/bin/python
'''
Utility module to compute broadband and multiband connectivity using Echobase.
'''

from util import *

np.random.seed(sum(map(ord, "aesthetics")))

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

warnings.filterwarnings('ignore')

## FUNCTIONAL CONNECTIVITY METHODS
def _helper_compute_multiband_connectivity(job):
    """
    Helper function for compute_multiband_connectivity.
    Parameters
    ----------
        jobs: tuple,
            Contains job details for each ictal block to process.

    Returns
    -------
        unique_idx: list
            Saves all adjacency matrices in different bands as npz files in comp_dir with unique ids generated using UUID4. These IDs along with each event type (e.g. 'ICTAL') and event id (e.g. '1') are output as a list of tuples.
    """
    # Get job details
    epoch,epoch_length,Fs,evData = job

    # Get clip
    data_clip = evData[epoch*epoch_length*Fs:(epoch+1)*epoch_length*Fs,:]

    # Compute multiband connectivity
    adj_alphatheta, adj_beta, adj_lowgamma, adj_highgamma = multiband_conn(data_clip, Fs)

    # Compute broadband connectivity
    adj_broadband_CC = broadband_conn(data_clip, Fs)

    return (epoch,adj_alphatheta, adj_beta, adj_lowgamma, adj_highgamma, adj_broadband_CC)

def compute_multiband_connectivity(patient_id, epoch_length=1, data=data):
    """
    Function for computing multiband connectivity adjacency matrices
    Parameters
    ----------
        patient_id: str
            Patient ID in DATA.json

        epoch_length: int
            Duration of epoch in seconds used to compute the adjacency matrix. The total number of adjacency matrices saved will be block length (sec) /epoch_length where block length = 900 seconds (15 minutes) of recording for e.g.

        data: dict
            Dictionary of data loaded from DATA.json (or TEST_DATA.json during unit tests). Default is DATA.json.
    Returns
    -------
        None
            Saves all adjacency matrices in different bands as npz files in comp_dir.
    """

    n_proc = 40
    pool = Pool(n_proc)

    # Generate list of cartoon map labels
    labels = map(lambda x: x.split(',')[4].replace('\n',''), open(os.path.expanduser(
        data['PATIENTS'][patient_id]['ELECTRODE_LABELS']
        ),'r').readlines())

    # Get path
    comp_dir = os.path.expanduser(data['COMP_DIR'])
    data_dir = os.path.expanduser(data['DATA_DIR'])

    # Load ignored node labels
    ignored_node_labels = data['PATIENTS'][patient_id]['IGNORE_ELECTRODES']
    for ignored_node_label in ignored_node_labels:
        if(ignored_node_label not in labels):
            labels.append(ignored_node_label)

    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
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
            evData = scipy.stats.zscore(evData,axis=1)
            T = evData.shape[0]

            # Correspond lable names
            labels = correspond_label_names(channels, labels)

            # Load electrodes to ignore
            ignored_node_idx = map(lambda x: labels[x][0],ignored_node_labels)
            for ii,node_id in enumerate(ignored_node_idx):
                print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
            evData = np.delete(evData, ignored_node_idx, axis=1)

            # For each clip, broken up in to epochs of length epoch_length
            epochs = int(T/(epoch_length*Fs))
            all_adj_alphatheta = np.zeros((evData.shape[1],evData.shape[1],epochs))
            all_adj_beta = np.zeros((evData.shape[1],evData.shape[1],epochs))
            all_adj_lowgamma = np.zeros((evData.shape[1],evData.shape[1],epochs))
            all_adj_highgamma = np.zeros((evData.shape[1],evData.shape[1],epochs))
            all_adj_broadband_CC = np.zeros((evData.shape[1],evData.shape[1],epochs))

            # Create parallel jobs for each block
            jobs = []
            for epoch in range(epochs):
                jobs.append((epoch, epoch_length, Fs, evData))

            return_list = pool.map(_helper_compute_multiband_connectivity, jobs)

            # Process results from Pool
            return_list = sorted(return_list, key=lambda x: x[0])
            for epoch in range(epochs):
                assert return_list[epoch][0] == epoch
                all_adj_alphatheta[:,:,epoch] = return_list[epoch][1]
                all_adj_beta[:,:,epoch] = return_list[epoch][2]
                all_adj_lowgamma[:,:,epoch] = return_list[epoch][3]
                all_adj_highgamma[:,:,epoch] = return_list[epoch][4]
                all_adj_broadband_CC[:,:,epoch] = return_list[epoch][5]

            # Save with appropriate name
            print 'Writing adjacency matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
            adj_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id))
            np.savez(open(adj_fn,'w'), all_adj_alphatheta=all_adj_alphatheta, all_adj_beta=all_adj_beta, all_adj_lowgamma=all_adj_lowgamma, all_adj_highgamma=all_adj_highgamma, epoch_length=epoch_length, all_adj_broadband_CC=all_adj_broadband_CC)
