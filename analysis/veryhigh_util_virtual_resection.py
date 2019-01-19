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
def get_resected_electrodes(patient_id, dilate_radius=0, data=data, labels_dict=None):
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
        resection_data[resection_data > 0.2] = 1
        resection_data[resection_data <= 0.2] = 0
        resection_affine = resection_nii.get_affine()
    except KeyError:
        assert isinstance(patient_id, str)
        print 'The PATIENT ID %s is not configured in the DATA JSON config file. ' \
              'Please use a PATIENT ID that is included in the study.' % patient_id
        raise

    # Load electrode labels and generate electrode image
    electrodes = {}
    ele = {}
    # ele = np.zeros(resection_data.shape)
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
            ele[electrode_id] = [X,Y,Z]
            # ele[max(X - radius, 0):min(X + radius, ele.shape[0]), max(Y - radius, 0):min(Y + radius, ele.shape[1]),
            # max(Z - radius, 0):min(Z + radius, ele.shape[2])] = electrode_id

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

    # Compute coordinates of resection zone
    hull = np.array(np.where(resection_data)).T

    # Determine which points are in the resection zone
    resected_electrode_ids = np.array(ele.keys())[in_hull(np.array(ele.values()), hull)]

    # Determine IDs that are included in the network
    network_node_labels = labels_dict.keys()
    resected_node_labels = map(lambda x: electrodes[x], resected_electrode_ids)

    clean_resected_node_labels = []
    for resected_node_label in resected_node_labels:
        if resected_node_label not in network_node_labels:
            continue
        else:
            clean_resected_node_labels.append(resected_node_label)

    resected_node_ids = map(lambda x: labels_dict[x][0], clean_resected_node_labels)

    network_resection_fraction = len(resected_node_ids)/np.float(len(network_node_labels))
    resection_fraction = network_resection_fraction

    # Erode/Dilate the resection image appropriately
    dilate_radius = dilate_radius/100.0
    if dilate_radius > 0:
        while resection_fraction < network_resection_fraction + dilate_radius:
            resection_data = morphology.binary_dilation(resection_data, iterations=1)
            # Compute coordinates of resection zone
            hull = np.array(np.where(resection_data)).T

            # Determine which points are in the resection zone
            resected_electrode_ids = np.array(ele.keys())[in_hull(np.array(ele.values()), hull)]
            # Determine IDs that are included in the network
            resected_node_labels = map(lambda x: electrodes[x], resected_electrode_ids)

            clean_resected_node_labels = []
            for resected_node_label in resected_node_labels:
                if resected_node_label not in network_node_labels:
                    continue
                else:
                    clean_resected_node_labels.append(resected_node_label)

            resected_node_ids = map(lambda x: labels_dict[x][0], clean_resected_node_labels)

            resection_fraction = len(resected_node_ids)/np.float(len(network_node_labels))
            print resection_fraction, resected_node_ids, clean_resected_node_labels

    elif dilate_radius < 0:
        while resection_fraction > max(network_resection_fraction + dilate_radius,0):
            resection_data = morphology.binary_erosion(resection_data, iterations=1)
            # Compute coordinates of resection zone
            hull = np.array(np.where(resection_data)).T

            # Determine which points are in the resection zone
            resected_electrode_ids = np.array(ele.keys())[in_hull(np.array(ele.values()), hull)]
            # Determine IDs that are included in the network
            resected_node_labels = map(lambda x: electrodes[x], resected_electrode_ids)

            clean_resected_node_labels = []
            for resected_node_label in resected_node_labels:
                if resected_node_label not in network_node_labels:
                    continue
                else:
                    clean_resected_node_labels.append(resected_node_label)

            resected_node_ids = map(lambda x: labels_dict[x][0], clean_resected_node_labels)

            resection_fraction = len(resected_node_ids)/np.float(len(network_node_labels))
            print resection_fraction, resected_node_ids, clean_resected_node_labels

    resected_electrodes = zip(map(lambda x: str(x),resected_node_ids),clean_resected_node_labels)

    return resected_electrodes

def write_resected_electrodes(patient_id, dilate_radius=0, data=data, labels_dict=None):
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
                      get_resected_electrodes(patient_id, dilate_radius, data, labels_dict))))

    return os.path.join(comp_dir, patient_id, 'aim1/%s_resected_electrodes_%s.csv' % (patient_id, suffix))

def base_synchronizability(adj):
    """
    Function for computing base synchronizability of the entire network
    Parameters
    ----------
        adj: ndarray, shape (N, N)
            Undirected, symmetric adjacency matrix with N nodes

    Returns
    -------
        base_sync: float
            Base synchronizability of the network
    """

    # Standard param checks
    errors.check_type(adj, np.ndarray)
    errors.check_dims(adj, 2)
    if not (adj == adj.T).all():
        raise Exception('Adjacency matrix is not undirected and symmetric.')
    if(np.isnan(adj).any()):
        return np.nan

    # Get data attributes
    n_node = adj.shape[0]

    # Get the original synchronizability
    base_sync = synchronizability(adj)

    return base_sync

def region_control(adj, node_list, base_sync=None):
    """
    Function for computing control centrality of a subregion, aka resection zone (change in synchronizability)
    Parameters
    ----------
        adj: ndarray, shape (N, N)
            Undirected, symmetric adjacency matrix with N nodes

        node_list: list
            List of node indices to simultaneously remove from the network

        base_sync: ndarray
            Base synchronization
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
    if(np.isnan(adj).any()):
        return np.nan

    # Get data attributes
    n_node = adj.shape[0]

    # Get the original synchronizability
    if(base_sync == None):
        base_sync = synchronizability(adj)

    adj_lesion = lesion.node_lesion(adj, node_list)
    lesion_sync = synchronizability(adj_lesion)

    return (lesion_sync-base_sync) / base_sync

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

    # Create output UUID codx
    unique_idx = []

    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
    for event_type, events in data['PATIENTS'][patient_id]['Events'].items():

        # Find the unique id
        for fn in os.listdir(os.path.join(comp_dir,patient_id,'aim3')):
            if 'pipedef.json' in fn:
                with open(os.path.join(comp_dir,patient_id,'aim3',fn)) as pipedef_file:
                    pipedef = json.load(pipedef_file)
                    if pipedef['dilate_radius'] == dilate_radius:
                        unique_id = fn.split('.')[-3]
                        assert unique_id == fn.split('.')[4]
                        break

        for event_id in events.keys():
            cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.npz'%(patient_id,event_type,event_id,unique_id))
            tmp = np.load(cres_fn)
            control_centrality_alphatheta = tmp['control_centrality_alphatheta']
            control_centrality_beta = tmp['control_centrality_beta']
            control_centrality_lowgamma = tmp['control_centrality_lowgamma']
            control_centrality_highgamma = tmp['control_centrality_highgamma']

            control_centrality_broadband_CC = tmp['control_centrality_broadband_CC']
            non_control_centrality_alphatheta = tmp['non_control_centrality_alphatheta']
            non_control_centrality_beta = tmp['non_control_centrality_beta']
            non_control_centrality_lowgamma = tmp['non_control_centrality_lowgamma']
            non_control_centrality_highgamma = tmp['non_control_centrality_highgamma']

            non_control_centrality_broadband_CC = tmp['non_control_centrality_broadband_CC']
            base_sync_alphatheta = tmp['base_sync_alphatheta']
            base_sync_beta = tmp['base_sync_beta']
            base_sync_lowgamma = tmp['base_sync_lowgamma']
            base_sync_highgamma = tmp['base_sync_highgamma']

            base_sync_broadband_CC = tmp['base_sync_broadband_CC']


            try:
                if(events[event_id]['STATUS'] == 'ALL_DROPOUT'):
                        continue # unusable clip
            except KeyError:
                pass

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

            # Correspond label names
            labels_dict = correspond_label_names(channels, labels)

            # Load electrodes to ignore
            ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
            for ii,node_id in enumerate(ignored_node_idx):
                print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
            channels = list(np.delete(np.array(channels),ignored_node_idx))

            # Recorrespond label names
            labels_dict = correspond_label_names(channels, labels)

            # Generate list of resected electrodes and write to CSV file
            try:
                if(dilate_radius == 0):
                    resected_node_labels = data['PATIENTS'][patient_id]['RESECTED_ELECTRODES']
                elif(dilate_radius > 0 or dilate_radius < 0):
                    resected_node_labels = data['PATIENTS'][patient_id]['RESECTED_ELECTRODES']
                    for fringe_node_label in data['PATIENTS'][patient_id]['RESECTED_FRINGE_ELECTRODES']:
                        resected_node_labels.append(fringe_node_label)
                else:
                    return []
            except Exception:
                resected_electrodes_fn = write_resected_electrodes(patient_id, dilate_radius, data, labels_dict)

                # Load resected electrodes
                try:
                    resected_nodes = map(lambda x: int(x.split(',')[0]), open(resected_electrodes_fn,'r').readlines())
                    resected_node_labels = map(lambda x: x.split(',')[1].replace('\n',''), open(resected_electrodes_fn,'r').readlines())
                except IndexError:
                    print 'ERROR! Resected electrodes %s does not have any electrodes. Skipping'%(resected_electrodes_fn)
                    return []

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

            # Map the NON-resected electrodes to channels
            all_node_idx = map(lambda x: labels_dict[x][0], labels_dict.keys())
            non_resected_node_idx = []
            for idx in all_node_idx:
                if(idx in resected_node_idx):
                    continue
                else:
                    non_resected_node_idx.append(idx)
            # non_resected_node_idx = np.array(non_resected_node_idx)

            # For each clip, load up adjacency matrices
            adj_file = np.load(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id)))

            epoch_length = int(adj_file['epoch_length'])
            all_adj_veryhigh = adj_file['all_adj_veryhigh']
            epochs = int(T/(epoch_length*Fs))

            assert all_adj_veryhigh.shape[2] == epochs

            # Perform resection of network
            control_centrality_veryhigh = np.zeros((epochs,))

            for epoch in range(epochs):
                if(np.isnan(all_adj_veryhigh[:,:,epoch]).any() or len(resected_node_idx) >= len(channels)-1):
                    control_centrality_veryhigh[epoch] = np.nan
                else:
                    control_centrality_veryhigh[epoch] = region_control(all_adj_veryhigh[:,:,epoch],resected_node_idx)

            # Perform resection of non-ROZ network
            non_control_centrality_veryhigh = np.zeros((epochs,))

            for epoch in range(epochs):
                if(np.isnan(all_adj_veryhigh[:,:,epoch]).any() or len(non_resected_node_idx) >= len(channels)-1):
                    non_control_centrality_veryhigh[epoch] = np.nan
                else:
                    non_control_centrality_veryhigh[epoch] = region_control(all_adj_veryhigh[:,:,epoch],non_resected_node_idx)

            # Compute base synchronizability of network
            base_sync_veryhigh = np.zeros((epochs,))

            for epoch in range(epochs):
                if(np.isnan(all_adj_veryhigh[:,:,epoch]).any()):
                    base_sync_veryhigh[epoch] = np.nan
                else:
                    base_sync_veryhigh[epoch] = base_synchronizability(all_adj_veryhigh[:,:,epoch])

            # Save with appropriate name
            print 'Writing c_res(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
            cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.npz'%(patient_id,event_type,event_id,unique_id))
            np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma, control_centrality_veryhigh=control_centrality_veryhigh, control_centrality_broadband_CC=control_centrality_broadband_CC, non_control_centrality_alphatheta=non_control_centrality_alphatheta, non_control_centrality_beta=non_control_centrality_beta, non_control_centrality_lowgamma=non_control_centrality_lowgamma, non_control_centrality_highgamma=non_control_centrality_highgamma, non_control_centrality_veryhigh=non_control_centrality_veryhigh, non_control_centrality_broadband_CC=non_control_centrality_broadband_CC, base_sync_alphatheta=base_sync_alphatheta, base_sync_beta=base_sync_beta, base_sync_lowgamma=base_sync_lowgamma, base_sync_highgamma=base_sync_highgamma, base_sync_veryhigh=base_sync_veryhigh, base_sync_broadband_CC=base_sync_broadband_CC,labels=labels_dict)
            pipeline_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.pipedef.json'%(patient_id,event_type,event_id,unique_id))
            timestamp = datetime.datetime.now()
            pipedef = {'fconn':'multiband+broadband', 'epoch_length':epoch_length, 'dilate_radius':dilate_radius, 'time_stamp':timestamp.strftime('%Y-%m-%d %H:%M:%S')}
            with open(pipeline_fn, 'w') as fp:
                json.dump(pipedef, fp)
            unique_idx.append((unique_id,event_type,event_id))
    return unique_idx


def nodal_virtual_resection(patient_id, data=data):

    """
    Function for computing c_resection(t).
    Parameters
    ----------
        patient_id: str
            Patient ID in DATA.json

        data: dict
            Dictionary of data loaded from DATA.json (or TEST_DATA.json during unit tests). Default is DATA.json.
    Returns
    -------
        None
            Saves all control centrality matrices in different bands as npz files in comp_dir.
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


    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
    for event_type, events in data['PATIENTS'][patient_id]['Events'].items():
        for event_id in events.keys():
            try:
                if(events[event_id]['STATUS'] == 'ALL_DROPOUT'):
                    continue # unusable clip
            except KeyError:
                pass

            if os.path.exists(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.noderes.npz'%(patient_id,event_type,event_id))):
                tmp = np.load(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.noderes.npz'%(patient_id,event_type,event_id)))
                control_centrality_alphatheta = tmp['control_centrality_alphatheta']
                control_centrality_beta = tmp['control_centrality_beta']
                control_centrality_lowgamma = tmp['control_centrality_lowgamma']
                control_centrality_highgamma = tmp['control_centrality_highgamma']

                control_centrality_broadband_CC = tmp['control_centrality_broadband_CC']
                base_sync_alphatheta = tmp['base_sync_alphatheta']
                base_sync_beta = tmp['base_sync_beta']
                base_sync_lowgamma = tmp['base_sync_lowgamma']
                base_sync_highgamma = tmp['base_sync_highgamma']

                base_sync_broadband_CC = tmp['base_sync_broadband_CC']

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
                T = evData.shape[0]

                # Correspond label names
                labels_dict = correspond_label_names(channels, labels)

                # Load electrodes to ignore
                ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
                for ii,node_id in enumerate(ignored_node_idx):
                    print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
                channels = list(np.delete(np.array(channels),ignored_node_idx))

                # Recorrespond label names
                labels_dict = correspond_label_names(channels, labels)

                # For each clip, load up adjacency matrices
                adj_file = np.load(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id)))

                epoch_length = int(adj_file['epoch_length'])
                all_adj_alphatheta = adj_file['all_adj_alphatheta']
                all_adj_beta = adj_file['all_adj_beta']
                all_adj_lowgamma = adj_file['all_adj_lowgamma']
                all_adj_highgamma = adj_file['all_adj_highgamma']
                all_adj_veryhigh = adj_file['all_adj_veryhigh']
                all_adj_broadband_CC = adj_file['all_adj_broadband_CC']
                epochs = int(T/(epoch_length*Fs))
                num_nodes = len(channels)

                assert all_adj_alphatheta.shape[0] == num_nodes
                assert all_adj_alphatheta.shape[1] == num_nodes
                assert all_adj_beta.shape[0] == num_nodes
                assert all_adj_beta.shape[1] == num_nodes
                assert all_adj_lowgamma.shape[0] == num_nodes
                assert all_adj_lowgamma.shape[1] == num_nodes
                assert all_adj_highgamma.shape[0] == num_nodes
                assert all_adj_highgamma.shape[1] == num_nodes
                assert all_adj_veryhigh.shape[0] == num_nodes
                assert all_adj_veryhigh.shape[1] == num_nodes
                assert all_adj_broadband_CC.shape[0] == num_nodes
                assert all_adj_broadband_CC.shape[1] == num_nodes
                assert all_adj_alphatheta.shape[2] == epochs
                assert all_adj_beta.shape[2] == epochs
                assert all_adj_lowgamma.shape[2] == epochs
                assert all_adj_highgamma.shape[2] == epochs
                assert all_adj_veryhigh.shape[2] == epochs
                assert all_adj_broadband_CC.shape[2] == epochs

                # Perform resection of network
                control_centrality_veryhigh = np.zeros((num_nodes,epochs,))

                for epoch in range(epochs):
                    if(np.isnan(all_adj_veryhigh[:,:,epoch]).any()):
                        control_centrality_veryhigh[:,epoch] = np.nan
                    else:
                        control_centrality_veryhigh[:,epoch] = node_control(all_adj_veryhigh[:,:,epoch])

                # Compute base synchronizability of network
                base_sync_veryhigh = np.zeros((epochs,))

                for epoch in range(epochs):
                    if(np.isnan(all_adj_veryhigh[:,:,epoch]).any()):
                        base_sync_veryhigh[epoch] = np.nan
                    else:
                        base_sync_veryhigh[epoch] = base_synchronizability(all_adj_veryhigh[:,:,epoch])

                # Save with appropriate name
                print 'Writing nodal c_res(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
                cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.noderes.npz'%(patient_id,event_type,event_id))
                np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma, control_centrality_veryhigh=control_centrality_veryhigh, control_centrality_broadband_CC=control_centrality_broadband_CC, base_sync_alphatheta=base_sync_alphatheta, base_sync_beta=base_sync_beta, base_sync_lowgamma=base_sync_lowgamma, base_sync_highgamma=base_sync_highgamma, base_sync_veryhigh=base_sync_veryhigh, base_sync_broadband_CC=base_sync_broadband_CC,labels=labels_dict)


            else:
                # blah

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
                T = evData.shape[0]

                # Correspond label names
                labels_dict = correspond_label_names(channels, labels)

                # Load electrodes to ignore
                ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
                for ii,node_id in enumerate(ignored_node_idx):
                    print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
                channels = list(np.delete(np.array(channels),ignored_node_idx))

                # Recorrespond label names
                labels_dict = correspond_label_names(channels, labels)

                # For each clip, load up adjacency matrices
                adj_file = np.load(os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id)))

                epoch_length = int(adj_file['epoch_length'])
                all_adj_alphatheta = adj_file['all_adj_alphatheta']
                all_adj_beta = adj_file['all_adj_beta']
                all_adj_lowgamma = adj_file['all_adj_lowgamma']
                all_adj_highgamma = adj_file['all_adj_highgamma']
                all_adj_veryhigh = adj_file['all_adj_veryhigh']
                all_adj_broadband_CC = adj_file['all_adj_broadband_CC']
                epochs = int(T/(epoch_length*Fs))
                num_nodes = len(channels)

                assert all_adj_alphatheta.shape[0] == num_nodes
                assert all_adj_alphatheta.shape[1] == num_nodes
                assert all_adj_beta.shape[0] == num_nodes
                assert all_adj_beta.shape[1] == num_nodes
                assert all_adj_lowgamma.shape[0] == num_nodes
                assert all_adj_lowgamma.shape[1] == num_nodes
                assert all_adj_highgamma.shape[0] == num_nodes
                assert all_adj_highgamma.shape[1] == num_nodes
                assert all_adj_veryhigh.shape[0] == num_nodes
                assert all_adj_veryhigh.shape[1] == num_nodes
                assert all_adj_broadband_CC.shape[0] == num_nodes
                assert all_adj_broadband_CC.shape[1] == num_nodes
                assert all_adj_alphatheta.shape[2] == epochs
                assert all_adj_beta.shape[2] == epochs
                assert all_adj_lowgamma.shape[2] == epochs
                assert all_adj_highgamma.shape[2] == epochs
                assert all_adj_veryhigh.shape[2] == epochs
                assert all_adj_broadband_CC.shape[2] == epochs

                # Perform resection of network
                control_centrality_alphatheta = np.zeros((num_nodes,epochs,))
                control_centrality_beta = np.zeros((num_nodes,epochs,))
                control_centrality_lowgamma = np.zeros((num_nodes,epochs,))
                control_centrality_highgamma = np.zeros((num_nodes,epochs,))
                control_centrality_veryhigh = np.zeros((num_nodes,epochs,))
                control_centrality_broadband_CC = np.zeros((num_nodes,epochs,))

                for epoch in range(epochs):
                    if(np.isnan(all_adj_alphatheta[:,:,epoch]).any()):
                        control_centrality_alphatheta[:,epoch] = np.nan
                    else:
                        control_centrality_alphatheta[:,epoch] = node_control(all_adj_alphatheta[:,:,epoch])
                    if(np.isnan(all_adj_beta[:,:,epoch]).any()):
                        control_centrality_beta[:,epoch] = np.nan
                    else:
                        control_centrality_beta[:,epoch] = node_control(all_adj_beta[:,:,epoch])
                    if(np.isnan(all_adj_lowgamma[:,:,epoch]).any()):
                        control_centrality_lowgamma[:,epoch] = np.nan
                    else:
                        control_centrality_lowgamma[:,epoch] = node_control(all_adj_lowgamma[:,:,epoch])
                    if(np.isnan(all_adj_highgamma[:,:,epoch]).any()):
                        control_centrality_highgamma[:,epoch] = np.nan
                    else:
                        control_centrality_highgamma[:,epoch] = node_control(all_adj_highgamma[:,:,epoch])
                    if(np.isnan(all_adj_veryhigh[:,:,epoch]).any()):
                        control_centrality_veryhigh[:,epoch] = np.nan
                    else:
                        control_centrality_veryhigh[:,epoch] = node_control(all_adj_veryhigh[:,:,epoch])
                    if(np.isnan(all_adj_broadband_CC[:,:,epoch]).any()):
                        control_centrality_broadband_CC[:,epoch] = np.nan
                    else:
                        control_centrality_broadband_CC[:,epoch] = node_control(all_adj_broadband_CC[:,:,epoch])

                # Compute base synchronizability of network
                base_sync_alphatheta = np.zeros((epochs,))
                base_sync_beta = np.zeros((epochs,))
                base_sync_lowgamma = np.zeros((epochs,))
                base_sync_highgamma = np.zeros((epochs,))
                base_sync_veryhigh = np.zeros((epochs,))
                base_sync_broadband_CC = np.zeros((epochs,))

                for epoch in range(epochs):
                    if(np.isnan(all_adj_alphatheta[:,:,epoch]).any()):
                        base_sync_alphatheta[epoch] = np.nan
                    else:
                        base_sync_alphatheta[epoch] = base_synchronizability(all_adj_alphatheta[:,:,epoch])
                    if(np.isnan(all_adj_beta[:,:,epoch]).any()):
                        base_sync_beta[epoch] = np.nan
                    else:
                        base_sync_beta[epoch] = base_synchronizability(all_adj_beta[:,:,epoch])
                    if(np.isnan(all_adj_lowgamma[:,:,epoch]).any()):
                        base_sync_lowgamma[epoch] = np.nan
                    else:
                        base_sync_lowgamma[epoch] = base_synchronizability(all_adj_lowgamma[:,:,epoch])
                    if(np.isnan(all_adj_highgamma[:,:,epoch]).any()):
                        base_sync_highgamma[epoch] = np.nan
                    else:
                        base_sync_highgamma[epoch] = base_synchronizability(all_adj_highgamma[:,:,epoch])
                    if(np.isnan(all_adj_veryhigh[:,:,epoch]).any()):
                        base_sync_veryhigh[epoch] = np.nan
                    else:
                        base_sync_veryhigh[epoch] = base_synchronizability(all_adj_veryhigh[:,:,epoch])
                    if(np.isnan(all_adj_broadband_CC[:,:,epoch]).any()):
                        base_sync_broadband_CC[epoch] = np.nan
                    else:
                        base_sync_broadband_CC[epoch] = base_synchronizability(all_adj_broadband_CC[:,:,epoch])

                # Save with appropriate name
                print 'Writing nodal c_res(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
                cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.noderes.npz'%(patient_id,event_type,event_id))
                np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma, control_centrality_veryhigh=control_centrality_veryhigh, control_centrality_broadband_CC=control_centrality_broadband_CC, base_sync_alphatheta=base_sync_alphatheta, base_sync_beta=base_sync_beta, base_sync_lowgamma=base_sync_lowgamma, base_sync_highgamma=base_sync_highgamma, base_sync_veryhigh=base_sync_veryhigh, base_sync_broadband_CC=base_sync_broadband_CC,labels=labels_dict)

