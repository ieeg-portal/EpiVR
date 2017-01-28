#!/usr/bin/python
'''
Utility module to compute virtual resection.
'''

from util import *

np.random.seed(sum(map(ord, "aesthetics")))

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

warnings.filterwarnings('ignore')

## VIRTUAL RESECTION METHODS
def get_resected_electrodes(patient_id, dilate_radius=0, data=data):
    """
    Utility function to compute and output labels of all resected electrodes.
    Parameters
    ----------
        patient_id: str
            Patient ID

        dilate_radius: int
            Amount of dilation/erosion to apply to resection image.

    Returns
    -------
        resected_electrodes: list of tuples
            Each tuple will contain electrode IDs and their label text
    """

    # Initialize output
    resected_electrodes = []

    # Get resection image
    try:
        resection_nii = nib.load(os.path.expanduser(data['PATIENTS'][patient_id]['RESECTION_IMAGE']))
        resection_data = resection_nii.get_data()
        # Threshold and binarize the resection mask
        resection_data[resection_data >= 0.8] = 1
        resection_data[resection_data < 0.8] = 0
        resection_affine = resection_nii.get_affine()
    except KeyError:
        assert isinstance(patient_id, str)
        print 'The PATIENT ID %s is not configured in the DATA JSON config file. ' \
              'Please use a PATIENT ID that is included in the study.' % patient_id
        raise

    # Load electrode labels and generate electrode image
    electrodes = {}
    ele = np.zeros(resection_data.shape)
    radius = 2  # Radius of electrodes
    try:
        # Generate electrode image
        lines = open(os.path.expanduser(data['PATIENTS'][patient_id]['ELECTRODE_LABELS']), 'r').readlines()
        for line in lines:
            electrode_id = int(line.split(',')[0])
            X = int(float(line.split(',')[1]))
            Y = int(float(line.split(',')[2]))
            Z = int(float(line.split(',')[3]))
            electrode_label = line.split(',')[4]
            ele[max(X - radius, 0):min(X + radius, ele.shape[0]), max(Y - radius, 0):min(Y + radius, ele.shape[1]),
            max(Z - radius, 0):min(Z + radius, ele.shape[2])] = electrode_id

            # Keep filling in electrodes dictionary
            electrodes[electrode_id] = electrode_label.replace('\n','')
    except ValueError:
        print 'The electrode labels is incorrectly formatted. ' \
              'Make sure X, Y, Z are the second, third and forth columns ' \
              'and make sure all coordinates are floating point values.'
        raise
    except KeyError:
        print 'Make sure the data config file has the patient %s and their ELECTRODE_LABELS file path is set.' % patient_id
        raise

    # Erode/Dilate the resection image appropriately
    if dilate_radius > 0:
        resection_data = morphology.binary_dilation(resection_data, iterations=np.abs(dilate_radius))
    elif dilate_radius < 0:
        resection_data = morphology.binary_erosion(resection_data, iterations=np.abs(dilate_radius))

    # Apply the resection image mask and determine which electrodes are still present
    ele_masked = np.multiply(ele, resection_data)
    resected_electrode_ids = np.unique(ele_masked[ele_masked > 0])
    for resected_electrode_id in resected_electrode_ids:
        resected_electrodes.append((str(int(resected_electrode_id)), electrodes[resected_electrode_id]))

    return resected_electrodes

def write_resected_electrodes(patient_id, dilate_radius=0, data=data):
    """
    Utility function to write resected electrode labels to CSV.
    Parameters
    ----------
        patient_id: str
            Patient ID

        dilate_radius: int
            Amount of dilation/erosion to apply to resection image.

    Returns
    -------
        path to the resected electrodes csv file
            The resected electrodes csv file will contain the electrode IDs and their labels.
    """

    # Create suffix for output filename
    if dilate_radius > 0:
        suffix = 'dilate_{}'.format(dilate_radius)
    elif dilate_radius < 0:
        suffix = 'erode_{}'.format(dilate_radius)
    else:
        suffix = '0'

    # Check if already created
    comp_dir = os.path.expanduser(data["COMP_DIR"])
    if os.path.isfile(os.path.join(comp_dir, patient_id, 'aim1/%s_resected_electrodes_%s.csv' % (patient_id, suffix))):
        return os.path.join(comp_dir, patient_id, 'aim1/%s_resected_electrodes_%s.csv' % (patient_id, suffix))

    # Write to CSV file
    open(os.path.join(comp_dir, patient_id, 'aim1/%s_resected_electrodes_%s.csv' % (patient_id, suffix)), 'w').write(
        '\n'.join(map(lambda x: ','.join(x),
                      get_resected_electrodes(patient_id, dilate_radius, data))))

    return os.path.join(comp_dir, patient_id, 'aim1/%s_resected_electrodes_%s.csv' % (patient_id, suffix))

def region_control(adj, node_list):
    """
    Function for computing control centrality of a subregion, aka resection zone (change in synchronizability)
    Parameters
    ----------
        adj: ndarray, shape (N, N)
            Undirected, symmetric adjacency matrix with N nodes

        node_list: list
            List of node indices to simultaneously remove from the network
    Returns
    -------
        control_vec: ndarray, shape (N,)
            Control centrality based on delta_sync for each of N variates
    """

    # Standard param checks
    errors.check_type(adj, np.ndarray)
    errors.check_dims(adj, 2)
    if not (adj == adj.T).all():
        raise Exception('Adjacency matrix is not undirected and symmetric.')

    # Get data attributes
    n_node = adj.shape[0]

    # Get the original synchronizability
    base_sync = synchronizability(adj)

    adj_lesion = lesion.node_lesion(adj, node_list)
    lesion_sync = synchronizability(adj_lesion)

    return (lesion_sync-base_sync) / base_sync

def _null_region_control(job):
    """
    Function for computing control centrality of a subregion, aka resection zone (change in synchronizability)
    Parameters
    ----------
        jobs: tuple
            Job to run parallely for null computation purposes.

    Returns
    -------
        control_vec: ndarray, shape (N,)
            Control centrality based on delta_sync for each of N variates
    """

    epochs, all_adj_alphatheta, all_adj_beta, all_adj_lowgamma, all_adj_highgamma, all_adj_broadband_CC, resected_node_idx = job

    # Perform resection of network
    control_centrality_alphatheta = np.zeros((epochs,))
    control_centrality_beta = np.zeros((epochs,))
    control_centrality_lowgamma = np.zeros((epochs,))
    control_centrality_highgamma = np.zeros((epochs,))
    control_centrality_broadband_CC = np.zeros((epochs,))

    for epoch in range(epochs):
        control_centrality_alphatheta[epoch] = region_control(all_adj_alphatheta[:,:,epoch],resected_node_idx)
        control_centrality_beta[epoch] = region_control(all_adj_beta[:,:,epoch],resected_node_idx)
        control_centrality_lowgamma[epoch] = region_control(all_adj_lowgamma[:,:,epoch],resected_node_idx)
        control_centrality_highgamma[epoch] = region_control(all_adj_highgamma[:,:,epoch],resected_node_idx)
        control_centrality_broadband_CC[epoch] = region_control(all_adj_broadband_CC[:,:,epoch],resected_node_idx)

    return (control_centrality_alphatheta, control_centrality_beta, control_centrality_lowgamma, control_centrality_highgamma, control_centrality_broadband_CC, resected_node_idx)

def virtual_resection(patient_id, dilate_radius=0, data=data):
    """
    Function for computing c_resection(t).
    Parameters
    ----------
        patient_id: str
            Patient ID in DATA.json

        unique_id: str
            Unique UUID code for npz files to load adjacency matrices

        dilate_radius: int
            Radius of dilation (>0) or erosion (<2) of the resection estimate.

        data: dict
            Dictionary of data loaded from DATA.json (or TEST_DATA.json during unit tests). Default is DATA.json.
    Returns
    -------
        None
            Saves all adjacency matrices in different bands as npz files in comp_dir.
    """

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

    # Generate list of resected electrodes and write to CSV file
    resected_electrodes_fn = write_resected_electrodes(patient_id, dilate_radius, data)

    # Load resected electrodes
    resected_nodes = map(lambda x: int(x.split(',')[0]), open(resected_electrodes_fn,'r').readlines())
    resected_node_labels = map(lambda x: x.split(',')[1].replace('\n',''), open(resected_electrodes_fn,'r').readlines())

    # Create output UUID codx
    unique_idx = []

    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
    for event_type, events in data['PATIENTS'][patient_id]['Events'].items():
        unique_id = str(uuid.uuid4())
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
            labels_dict = correspond_label_names(channels, labels)

            # Load electrodes to ignore
            ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
            for ii,node_id in enumerate(ignored_node_idx):
                print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
            channels = list(np.delete(np.array(channels),ignored_node_idx))

            # Recorrespond label names
            labels_dict = correspond_label_names(channels, labels)

            # Map the resected electrodes to channels
            clean_resected_node_labels = []
            for resected_node_label in resected_node_labels:
                if resected_node_label in ignored_node_labels:
                    continue
                else:
                    clean_resected_node_labels.append(resected_node_label)
            resected_node_idx = map(lambda x: labels_dict[x][0], clean_resected_node_labels)
            for ii,node_id in enumerate(resected_node_idx):
                print 'Virtually resecting node label: %s because label %s is in the resection zone'%(channels[node_id],resected_node_labels[ii])

            # For each clip, load up adjacency matrices
            adj_file = np.load(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id)))

            epoch_length = int(adj_file['epoch_length'])
            all_adj_alphatheta = adj_file['all_adj_alphatheta']
            all_adj_beta = adj_file['all_adj_beta']
            all_adj_lowgamma = adj_file['all_adj_lowgamma']
            all_adj_highgamma = adj_file['all_adj_highgamma']
            all_adj_broadband_CC = adj_file['all_adj_broadband_CC']
            epochs = int(T/(epoch_length*Fs))

            assert all_adj_alphatheta.shape[2] == epochs
            assert all_adj_beta.shape[2] == epochs
            assert all_adj_lowgamma.shape[2] == epochs
            assert all_adj_highgamma.shape[2] == epochs
            assert all_adj_broadband_CC.shape[2] == epochs

            # Perform resection of network
            control_centrality_alphatheta = np.zeros((epochs,))
            control_centrality_beta = np.zeros((epochs,))
            control_centrality_lowgamma = np.zeros((epochs,))
            control_centrality_highgamma = np.zeros((epochs,))
            control_centrality_broadband_CC = np.zeros((epochs,))

            for epoch in range(epochs):
                control_centrality_alphatheta[epoch] = region_control(all_adj_alphatheta[:,:,epoch],resected_node_idx)
                control_centrality_beta[epoch] = region_control(all_adj_beta[:,:,epoch],resected_node_idx)
                control_centrality_lowgamma[epoch] = region_control(all_adj_lowgamma[:,:,epoch],resected_node_idx)
                control_centrality_highgamma[epoch] = region_control(all_adj_highgamma[:,:,epoch],resected_node_idx)
                control_centrality_broadband_CC[epoch] = region_control(all_adj_broadband_CC[:,:,epoch],resected_node_idx)

            # Save with appropriate name
            print 'Writing c_res(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
            cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.npz'%(patient_id,event_type,event_id,unique_id))
            np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma, control_centrality_broadband_CC=control_centrality_broadband_CC)
            pipeline_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.pipedef.json'%(patient_id,event_type,event_id,unique_id))
            pipedef = {'fconn':'multiband+broadband', 'epoch_length':epoch_length, 'dilate_radius':dilate_radius}
            with open(pipeline_fn, 'w') as fp:
                json.dump(pipedef, fp)
            unique_idx.append((unique_id,event_type,event_id))
    return unique_idx

def null_virtual_resection(patient_id, unique_id, event_type, event_id, dilate_radius=0, data=data):
    """
    Function for computing c_null(t).
    Parameters
    ----------
        patient_id: str
            Patient ID in DATA.json

        unique_id: str
            Unique UUID code for npz files to load adjacency matrices

        dilate_radius: int
            Radius of dilation (>0) or erosion (<2) of the resection estimate.

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

    # Generate list of resected electrodes and write to CSV file
    resected_electrodes_fn = write_resected_electrodes(patient_id, dilate_radius, data)

    # Load resected electrodes
    resected_nodes = map(lambda x: int(x.split(',')[0]), open(resected_electrodes_fn,'r').readlines())
    resected_node_labels = map(lambda x: x.split(',')[1].replace('\n',''), open(resected_electrodes_fn,'r').readlines())

    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
    fn = os.path.join(data_dir, patient_id, 'eeg', data['PATIENTS'][patient_id]['Events'][event_type][event_id]['FILE'])
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
    labels_dict = correspond_label_names(channels, labels)

    # Load electrodes to ignore
    ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
    for ii,node_id in enumerate(ignored_node_idx):
        print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
    channels = list(np.delete(np.array(channels),ignored_node_idx))

    # Recorrespond label names
    labels_dict = correspond_label_names(channels, labels)

    # Map the resected electrodes to channels
    clean_resected_node_labels = []
    for resected_node_label in resected_node_labels:
        if resected_node_label in ignored_node_labels:
            continue
        else:
            clean_resected_node_labels.append(resected_node_label)
    resected_node_idx = map(lambda x: labels_dict[x][0], clean_resected_node_labels)
    for ii,node_id in enumerate(resected_node_idx):
        print 'Virtually resecting node label: %s because label %s is in the resection zone'%(channels[node_id],resected_node_labels[ii])

    # For each clip, load up adjacency matrices
    adj_file = np.load(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id)))


    epoch_length = int(adj_file['epoch_length'])
    all_adj_alphatheta = adj_file['all_adj_alphatheta']
    all_adj_beta = adj_file['all_adj_beta']
    all_adj_lowgamma = adj_file['all_adj_lowgamma']
    all_adj_highgamma = adj_file['all_adj_highgamma']
    all_adj_broadband_CC = adj_file['all_adj_broadband_CC']
    epochs = int(T/(epoch_length*Fs))

    assert all_adj_alphatheta.shape[2] == epochs
    assert all_adj_beta.shape[2] == epochs
    assert all_adj_lowgamma.shape[2] == epochs
    assert all_adj_highgamma.shape[2] == epochs
    assert all_adj_broadband_CC.shape[2] == epochs

    # Create parallel jobs for computation
    jobs = []
    for perm_iter in range(1000):
        permuted_resected_node_idx = list(np.random.choice(np.arange(len(channels)), len(resected_node_idx)))
        jobs.append((epochs, all_adj_alphatheta, all_adj_beta, all_adj_lowgamma, all_adj_highgamma, all_adj_broadband_CC, permuted_resected_node_idx))

    return_list = pool.map(_null_region_control,jobs)

    # Save with appropriate name
    for ii,result in enumerate(return_list):
        print 'Writing c_null(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
        control_centrality_alphatheta, control_centrality_beta, control_centrality_lowgamma, control_centrality_highgamma, control_centrality_broadband_CC, permuted_resected_node_idx = result
        cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cnull.%i.%s.npz'%(patient_id,event_type,event_id,ii+1,unique_id))
        np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma, control_centrality_broadband_CC
            =control_centrality_broadband_CC, permuted_resected_node_idx=permuted_resected_node_idx)
