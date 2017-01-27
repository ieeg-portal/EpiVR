#!/usr/bin/python
'''
Utility module run Echobase operations on any given
'''

import os
import sys
import json
import h5py
import re
import uuid

from multiprocessing import Pool
import warnings

import nibabel as nib
import matplotlib.pyplot as plt

import numpy as np
import numpy as np
import pandas as pd
import scipy as sp
import scipy.stats
from scipy.ndimage import morphology
from scipy.signal import convolve

import seaborn as sns

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from Echobase.Common import errors
from Echobase.Network.Metrics.globaltopo import synchronizability
from Echobase.Network.Transforms import lesion
from Echobase.Pipelines.ecog_network import *

np.random.seed(sum(map(ord, "aesthetics")))

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

warnings.filterwarnings('ignore')

## Utility Functions
def correspond_label_names(eeg_channel_labels, cartoon_map_labels):
    '''
    Function for corresponding channel labels from ECoG (usually set as a montage) and electrode labels on clinical cartoon maps.
    Parameters
    ----------
        eeg_channel_labels: list of str
            Channel labels as found on IEEG.org.

        cartoon_map_labels: list of str
            Electrode labels from cartoon map. Does not have to be complete set of labels.

    Returns
    -------
        dict
            Dictionary contains as key the cartoon map label and its corresponding EEG channel label along with index in eeg_channel_labels
    '''
    labels = {}
    for ii,channel_label in enumerate(eeg_channel_labels):
        for cartoon_label in cartoon_map_labels:
            for p in re.findall(r'([A-Za-z]+)[ _-]*([0-9]{1,2})',channel_label):
                if p[0] == re.match(r'([A-Za-z]+)([0-9]+)',cartoon_label).group(1) and str(int(p[1]))==re.match(r'([A-Za-z]+)([0-9]+)',cartoon_label).group(2):
                    labels[cartoon_label] = (ii,channel_label)
    return labels


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h


def plot_experiment(patient_id, unique_id, data=data):
    '''
    Utility function to plot the results of an experiment given unique_id.
    Parameters
    ----------
        unique_id: str,
            Patient ID

    Returns
    -------
        TBD: TBD
            TBD
    '''
    files = {'cres':[],'cnull':[],'pipedef':[]}

    # Find all the appropriate files
    comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')
    for fn in os.listdir(comp_dir):
        if(unique_id in fn):
            if('cres' in fn and '.npz' in fn):
                files['cres'].append(os.path.join(comp_dir,fn))
            elif('cres' in fn and '.json' in fn):
                files['pipedef'].append(os.path.join(comp_dir,fn))
            elif('null' in fn):
                files['cnull'].append(os.path.join(comp_dir,fn))



    # Load all cres data for each block
    cres_data = np.zeros((
        len(files['cres']),
        np.load(files['cres'][0])['control_centrality_highgamma'].shape[0]
        ))
    xticklabels = []
    for fn_iter,fn in enumerate(sorted(files['cres'])):
        cres_data[fn_iter,:] = np.load(fn)['control_centrality_highgamma']

    # Additionally smooth cres_data
    cres_smooth_data = np.zeros(cres_data.shape)
    for clip_iter in range(cres_data.shape[0]):
        smooth_window = 5
        kernel = 1/np.float(smooth_window)*np.ones((smooth_window,1))
        kernel = kernel.flatten()
        cres_smooth_data[clip_iter,:] = convolve(cres_data[clip_iter,:],kernel,'same')

    # Load all cnull data across all blocks
    cnull_data = np.array(())
    for fn_iter, fn in enumerate(files['cnull']):
        cnull_data = np.hstack((cnull_data,np.load(fn)['control_centrality_highgamma']))

    cnull_mean, ci5, ci95 = mean_confidence_interval(data=cnull_data)
    print cnull_mean, ci5, ci95

    # Convert all cres_data clips into dataframe
    cres_df = []
    for clip in range(1,4):
        for epoch in range(900):
            if epoch < 150:
                time_window = '-5 min. to -2.5 min'
            elif epoch < 300:
                time_window = '-2.5 min. to seizure onset'
            elif epoch < 450:
                time_window = 'seizure onset to 2.5 min'
            elif epoch < 600:
                time_window = '2.5 min. to 5 min'
            elif epoch < 750:
                time_window = '5 min. to 7.5 min'
            else:
                time_window = '7.5 min. to 10 min'
            cres_df.append(['Ictal Clip %i'%clip,time_window,epoch,cres_data[clip-1,epoch],patient_id])

    columns = ['Clip','TimeWindow','Epoch','Cres','PatientID']
    cres_df = pd.DataFrame(np.array(cres_df),columns=columns)
    cres_df.Clip = cres_df.Clip.astype("category")
    cres_df.PatientID = cres_df.PatientID.astype("category")
    cres_df.Cres = cres_df.Cres.astype('float64')
    cres_df.TimeWindow = cres_df.TimeWindow.astype("category", categories=['-5 min. to -2.5 min','-2.5 min. to seizure onset','seizure onset to 2.5 min','2.5 min. to 5 min','5 min. to 7.5 min','7.5 min. to 10 min'], ordered=True)

    # Convert all cres_smooth_data clips into dataframe
    cres_smooth_df = []
    for clip in range(1,4):
        for epoch in range(900):
            if epoch < 150:
                time_window = '-5 min. to -2.5 min'
            elif epoch < 300:
                time_window = '-2.5 min. to seizure onset'
            elif epoch < 450:
                time_window = 'seizure onset to 2.5 min'
            elif epoch < 600:
                time_window = '2.5 min. to 5 min'
            elif epoch < 750:
                time_window = '5 min. to 7.5 min'
            else:
                time_window = '7.5 min. to 10 min'
            cres_smooth_df.append(['Ictal Clip %i'%clip,time_window,epoch,cres_smooth_data[clip-1,epoch],patient_id])

    columns = ['Clip','TimeWindow','Epoch','Cres','PatientID']
    cres_smooth_df = pd.DataFrame(np.array(cres_smooth_df),columns=columns)
    cres_smooth_df.Clip = cres_smooth_df.Clip.astype("category")
    cres_smooth_df.PatientID = cres_df.PatientID.astype("category")
    cres_smooth_df.Cres = cres_smooth_df.Cres.astype('float64')
    cres_smooth_df.TimeWindow = cres_smooth_df.TimeWindow.astype("category", categories=['-5 min. to -2.5 min','-2.5 min. to seizure onset','seizure onset to 2.5 min','2.5 min. to 5 min','5 min. to 7.5 min','7.5 min. to 10 min'], ordered=True)


    # Get pipeline details
    pipedef = json.load(open(files['pipedef'][0],'r'))

    # Generate a figure
    fig = plt.figure(1,figsize=(18,4))
    epoch_length = int(pipedef['epoch_length'])
    plt.plot(np.arange(-5*60/epoch_length,10*60/epoch_length),cres_data.T)
    plt.plot([0, 0],[-3, 3], alpha=0.75)
    fig.axes[0].text(5, -0.75, 'Seizure onset')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Control Centrality')
    plt.title('%s Control Centrality across 3 ictal clips\n Null Confidence Interval: [%.3E %.3E]'%(patient_id, ci5,ci95))
    plt.ylim([-1,1])
    plt.savefig('%s/%s_%s_CC_highgamma_time_plot.png'%(comp_dir,patient_id,unique_id))

    # Define the color palette
    color_palette_sz = [(0.18431373554117539, 0.47266437376246734, 0.7116493828156415), (0.94071511310689593, 0.60991928098248505, 0.48127645896930321), (0.75617071810890646, 0.21038062695194693, 0.22352941947824814), (0.76617071810890646, 0.21038062695194693, 0.22352941947824814), (0.77617071810890646, 0.21038062695194693, 0.22352941947824814), (0.78617071810890646, 0.21038062695194693, 0.22352941947824814)]

    # Plot factor plot
    # # f, ax = plt.subplots()
    # ax = sns.factorplot(x="Cres", y="Clip", hue="TimeWindow", row="PatientID", data=cres_df, orient='h', size=2, aspect=3.5, palette=color_palette_sz, kind="violin", cut=0, bw=.2)
    # sns.despine(offset=10, trim=True)
    # ax.fig.set_alpha = 0.25
    # ax.fig.set_size_inches((8,18))
    # ax.fig.axes[0].fill_between([ci5, ci95],[-1, 3],[-1, 3],facecolor='gray',alpha=0.5)
    # plt.title('Patient %s Dilation Radius %i - High Gamma (Weakly synchronizing Engel Ib)'%(fn.split('.')[0].split('/')[-1],pipedef['dilate_radius']))
    # plt.show()

    ax = sns.factorplot(x="Cres", y="Clip", hue="TimeWindow", row="PatientID", data=cres_smooth_df, orient='h', size=2, aspect=3.5, palette=color_palette_sz, kind="violin", cut=0, bw=.2)
    sns.despine(offset=10, trim=True)
    ax.fig.set_alpha = 0.25
    ax.fig.set_size_inches((8,18))
    # ax.fig.axes[0].fill_between([ci5, ci95],[-1, 3],[-1, 3],facecolor='gray',alpha=0.5)

    plt.title('Patient %s Dilation Radius %i - High Gamma (Weakly synchronizing Engel Ib)'%(fn.split('.')[0].split('/')[-1],pipedef['dilate_radius']))
    plt.savefig('%s/%s_%s_CC_highgamma_distributions.png'%(comp_dir,patient_id,unique_id))
    for ictal_clip in range(cres_smooth_data.shape[0]):
        m,l, h = mean_confidence_interval(cres_smooth_data[ictal_clip,:])
        print m, l, h



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
            Saves all adjacency matrices in different bands as npz files in comp_dir with unique ids generated using UUID4. These IDs are output as a list.
    """

    epoch,epoch_length,Fs,evData = job

    # Get clip
    data_clip = evData[epoch*epoch_length*Fs:(epoch+1)*epoch_length*Fs,:]

    # Compute multiband connectivity
    adj_alphatheta, adj_beta, adj_lowgamma, adj_highgamma = multiband_conn(data_clip, Fs)

    return (epoch,adj_alphatheta, adj_beta, adj_lowgamma, adj_highgamma)

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

            # Create parallel jobs for each block
            jobs = []
            for epoch in range(epochs):
                jobs.append((epoch, epoch_length, Fs, evData))

            return_list = pool.map(_helper_compute_multiband_connectivity, jobs)

            return_list = sorted(return_list, key=lambda x: x[0])
            for epoch in range(epochs):
                assert return_list[epoch][0] == epoch
                all_adj_alphatheta[:,:,epoch] = return_list[epoch][1]
                all_adj_beta[:,:,epoch] = return_list[epoch][2]
                all_adj_lowgamma[:,:,epoch] = return_list[epoch][3]
                all_adj_highgamma[:,:,epoch] = return_list[epoch][4]

            # Save with appropriate name
            print 'Writing adjacency matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
            adj_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id))
            np.savez(open(adj_fn,'w'), all_adj_alphatheta=all_adj_alphatheta, all_adj_beta=all_adj_beta, all_adj_lowgamma=all_adj_lowgamma, all_adj_highgamma=all_adj_highgamma, epoch_length=epoch_length)

## VIRTUAL RESECTION METHODS
def get_resected_electrodes(patient_id, dilate_radius=0, data=data):
    '''
    Utility function to compute and output labels of all resected electrodes.
    :param patient_id: Patient ID (str)
    :param dilate_radius: Amount of dilation/erosion to apply to resection image.
    :return: list of tuples: (electrode IDs, electrode label)
    '''

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
    '''

    Utility function to write resected electrode labels to CSV.
    :param patient_id: Patient ID (str)
    :param dilate_radius: Amount of dilation/erosion to apply to resection image.
    :return: str Filename of resected electrodes CSV
    '''

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

    epochs, all_adj_alphatheta, all_adj_beta, all_adj_lowgamma, all_adj_highgamma, resected_node_idx = job

    # Perform resection of network
    control_centrality_alphatheta = np.zeros((epochs,))
    control_centrality_beta = np.zeros((epochs,))
    control_centrality_lowgamma = np.zeros((epochs,))
    control_centrality_highgamma = np.zeros((epochs,))

    for epoch in range(epochs):
        control_centrality_alphatheta[epoch] = region_control(all_adj_alphatheta[:,:,epoch],resected_node_idx)
        control_centrality_beta[epoch] = region_control(all_adj_beta[:,:,epoch],resected_node_idx)
        control_centrality_lowgamma[epoch] = region_control(all_adj_lowgamma[:,:,epoch],resected_node_idx)
        control_centrality_highgamma[epoch] = region_control(all_adj_highgamma[:,:,epoch],resected_node_idx)

    return (control_centrality_alphatheta, control_centrality_beta, control_centrality_lowgamma, control_centrality_highgamma, resected_node_idx)

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
            epochs = int(T/(epoch_length*Fs))

            assert all_adj_alphatheta.shape[2] == epochs
            assert all_adj_beta.shape[2] == epochs
            assert all_adj_lowgamma.shape[2] == epochs
            assert all_adj_highgamma.shape[2] == epochs

            # Perform resection of network
            control_centrality_alphatheta = np.zeros((epochs,))
            control_centrality_beta = np.zeros((epochs,))
            control_centrality_lowgamma = np.zeros((epochs,))
            control_centrality_highgamma = np.zeros((epochs,))

            for epoch in range(epochs):
                control_centrality_alphatheta[epoch] = region_control(all_adj_alphatheta[:,:,epoch],resected_node_idx)
                control_centrality_beta[epoch] = region_control(all_adj_beta[:,:,epoch],resected_node_idx)
                control_centrality_lowgamma[epoch] = region_control(all_adj_lowgamma[:,:,epoch],resected_node_idx)
                control_centrality_highgamma[epoch] = region_control(all_adj_highgamma[:,:,epoch],resected_node_idx)

            # Save with appropriate name
            print 'Writing c_res(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
            cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.npz'%(patient_id,event_type,event_id,unique_id))
            np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma)
            pipeline_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cres.%s.pipedef.json'%(patient_id,event_type,event_id,unique_id))
            pipedef = {'fconn':'multiband', 'epoch_length':epoch_length, 'dilate_radius':dilate_radius}
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
    epochs = int(T/(epoch_length*Fs))

    assert all_adj_alphatheta.shape[2] == epochs
    assert all_adj_beta.shape[2] == epochs
    assert all_adj_lowgamma.shape[2] == epochs
    assert all_adj_highgamma.shape[2] == epochs

    # Create parallel jobs for computation
    jobs = []
    for perm_iter in range(1000):
        permuted_resected_node_idx = list(np.random.choice(np.arange(len(channels)), len(resected_node_idx)))
        jobs.append((epochs, all_adj_alphatheta, all_adj_beta, all_adj_lowgamma, all_adj_highgamma, permuted_resected_node_idx))

    return_list = pool.map(_null_region_control,jobs)

    # Save with appropriate name
    for ii,result in enumerate(return_list):
        print 'Writing c_null(t) matrices for patient %s event %s %s'%(patient_id,event_type,event_id)
        control_centrality_alphatheta, control_centrality_beta, control_centrality_lowgamma, control_centrality_highgamma, permuted_resected_node_idx = result
        cres_fn = os.path.join(comp_dir,patient_id,'aim3','%s.%s.%s.cnull.%i.%s.npz'%(patient_id,event_type,event_id,ii+1,unique_id))
        np.savez(open(cres_fn,'w'), control_centrality_alphatheta=control_centrality_alphatheta, control_centrality_beta=control_centrality_beta, control_centrality_lowgamma=control_centrality_lowgamma, control_centrality_highgamma=control_centrality_highgamma, permuted_resected_node_idx=permuted_resected_node_idx)
