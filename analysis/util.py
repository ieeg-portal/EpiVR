#!/usr/bin/python
'''
Utility module run Echobase operations on ECoG network.
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
    """
    Function for computing confidence interval.
    Parameters
    ----------
        data: 1D array
            Data array to compute the CI.

        confidence: float
            Percent confidence

    Returns
    -------
        tuple
            The mean, lower bound and higher bound for the confidence interval.
    """
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
        patient_id: str,
            Patient ID

        unique_id: str,
            Unique UUID4 experiment ID.

    Returns
    -------
        None
            All plot figures are saved as png, svg in the COMP_DIR specified by data json file.
    '''

    # Load patient outcome
    try:
        outcome = np.float(data['PATIENTS'][patient_id]['Outcome'])
        if(outcome == 1.1):
            outcome = 'IA'
        elif(outcome == 1.2):
            outcome = 'IB'
        elif(outcome == 1.3):
            outcome = 'IC'
        elif(outcome == 1.4):
            outcome = 'ID'
        elif(outcome < 3):
            outcome = 'II'
        elif(outcome < 4):
            outcome = 'III'
        else:
            outcome = 'IV'
    except:
        outcome = 'NA'

    # Find all the appropriate files
    files = {'cres':[],'cnull':[],'pipedef':[]}
    comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')
    for fn in os.listdir(comp_dir):
        if(unique_id in fn):
            if('cres' in fn and '.npz' in fn):
                files['cres'].append(os.path.join(comp_dir,fn))
            elif('cres' in fn and '.json' in fn):
                files['pipedef'].append(os.path.join(comp_dir,fn))
            elif('null' in fn):
                files['cnull'].append(os.path.join(comp_dir,fn))

    # Generate figures for different types of connectivity measures
    for fconn in ['alphatheta','beta','lowgamma','highgamma', 'broadband_CC']:

        # Load all cres data for each block
        cres_data = np.zeros((
            len(files['cres']),
            np.load(files['cres'][0])['control_centrality_%s'%fconn].shape[0]
            ))
        xticklabels = []
        for fn_iter,fn in enumerate(sorted(files['cres'])):
            cres_data[fn_iter,:] = np.load(fn)['control_centrality_%s'%fconn]

        # Set details of dataset
        num_ictal_events = cres_data.shape[0]
        num_epochs = cres_data.shape[1]

        # Additionally smooth cres_data
        cres_smooth_data = np.zeros(cres_data.shape)
        for clip_iter in range(cres_data.shape[0]):
            smooth_window = 5
            kernel = 1/np.float(smooth_window)*np.ones((smooth_window,1))
            kernel = kernel.flatten()
            cres_smooth_data[clip_iter,:] = convolve(cres_data[clip_iter,:],kernel,'same')

        # Load all cnull data across all blocks
        cnull_data = np.zeros((num_ictal_events,num_epochs*1000))
        for event_id in range(1,num_ictal_events+1):
            tmp_data = np.array(())
            for fn_iter, fn in enumerate(files['cnull']):
                prefix = '%s.%s.%s'%(patient_id,'ICTAL',event_id)
                if(prefix.lower() in fn.lower()):
                    tmp_data = np.hstack((tmp_data,np.load(fn)['control_centrality_%s'%fconn]))
            cnull_data[event_id-1,:] = tmp_data

        # Compute 95% confidence interval
        cnull_mean, ci5, ci95 = mean_confidence_interval(data=cnull_data.flatten())
        print cnull_mean, ci5, ci95

        # Convert all cres_data clips into dataframe
        cres_df = []
        for clip in range(1,num_ictal_events+1):
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

        for clip in range(1,num_ictal_events+1):
            time_window = 'Null'
            for epoch in range(cnull_data.shape[1]):
                cres_df.append(['Ictal Clip %i'%clip,time_window,epoch,cnull_data[clip-1, epoch],patient_id])

        columns = ['Clip','TimeWindow','Epoch','Cres','PatientID']
        cres_df = pd.DataFrame(np.array(cres_df),columns=columns)
        cres_df.Clip = cres_df.Clip.astype("category")
        cres_df.PatientID = cres_df.PatientID.astype("category")
        cres_df.Cres = cres_df.Cres.astype('float64')
        cres_df.TimeWindow = cres_df.TimeWindow.astype("category", categories=['-5 min. to -2.5 min','-2.5 min. to seizure onset','seizure onset to 2.5 min','2.5 min. to 5 min','5 min. to 7.5 min','7.5 min. to 10 min','Null'], ordered=True)

        # Convert all cres_smooth_data clips into dataframe
        cres_smooth_df = []
        for clip in range(1,num_ictal_events+1):
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

        for clip in range(1,num_ictal_events+1):
            time_window = 'Null'
            for epoch in range(cnull_data.shape[1]):
                cres_smooth_df.append(['Ictal Clip %i'%clip,time_window,epoch,cnull_data[clip-1, epoch],patient_id])

        columns = ['Clip','TimeWindow','Epoch','Cres','PatientID']
        cres_smooth_df = pd.DataFrame(np.array(cres_smooth_df),columns=columns)
        cres_smooth_df.Clip = cres_smooth_df.Clip.astype("category")
        cres_smooth_df.PatientID = cres_df.PatientID.astype("category")
        cres_smooth_df.Cres = cres_smooth_df.Cres.astype('float64')
        cres_smooth_df.TimeWindow = cres_smooth_df.TimeWindow.astype("category", categories=['-5 min. to -2.5 min','-2.5 min. to seizure onset','seizure onset to 2.5 min','2.5 min. to 5 min','5 min. to 7.5 min','7.5 min. to 10 min','Null'], ordered=True)

        # Perform K-S 2 sample test on each window to null
        f = open('%s/%s_%s_CC_%s_distributions_stats.csv'%(comp_dir, patient_id, unique_id, fconn),'w')
        stats_txt = ''
        for event_id in range(1,num_ictal_events+1):
            for time_window in ['-5 min. to -2.5 min','-2.5 min. to seizure onset','seizure onset to 2.5 min','2.5 min. to 5 min','5 min. to 7.5 min','7.5 min. to 10 min']:
                D,p = scipy.stats.ks_2samp(cres_smooth_df[cres_smooth_df.TimeWindow == '-5 min. to -2.5 min'][cres_smooth_df.Clip=='Ictal Clip %i'%event_id].Cres,cnull_data[event_id-1,:])
                m,l,h = mean_confidence_interval(cres_smooth_df[cres_smooth_df.TimeWindow == time_window][cres_smooth_df.Clip == 'Ictal Clip %i'%event_id].Cres)
                if(p < 0.001):
                    stats_txt += 'Ictal Clip %i, %s, %0.6E, %0.6E, %0.6E, ***,%0.6E'%(event_id, time_window, m,l,h, p)
                elif(p < 0.005):
                    stats_txt += 'Ictal Clip %i, %s, %0.6E, %0.6E, %0.6E, **,%0.6E'%(event_id, time_window, m,l,h, p)
                elif(p < 0.05):
                    stats_txt += 'Ictal Clip %i, %s, %0.6E, %0.6E, %0.6E, *,%0.6E'%(event_id, time_window, m,l,h, p)
                else:
                    stats_txt += 'Ictal Clip %i, %s, %0.6E, %0.6E, %0.6E, NS,%0.6E'%(event_id, time_window, m,l,h, p)
                stats_txt += '\n'
        f.write(stats_txt)

        # Get pipeline details
        pipedef = json.load(open(files['pipedef'][0],'r'))

        # Generate a figure
        fig = plt.figure(1,figsize=(18,4))
        epoch_length = int(pipedef['epoch_length'])
        plt.plot(np.arange(-5*60/epoch_length,10*60/epoch_length),cres_smooth_data.T)
        plt.plot([0, 0],[-3, 3], alpha=0.1)
        fig.axes[0].text(5, -0.75, 'Seizure onset')
        plt.xlabel('Time (seconds)')
        plt.ylabel('Control Centrality')
        plt.title('%s Control Centrality across 3 ictal clips using %s\n Null Confidence Interval: [%.3E %.3E]'%(patient_id, fconn.upper(), ci5,ci95))
        plt.ylim([-1,1])
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_time_plot.png'%(comp_dir,patient_id,unique_id, fconn))
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_time_plot.svg'%(comp_dir,patient_id,unique_id, fconn))

        # Define the color palette
        color_palette_sz = [(0.18431373554117539, 0.47266437376246734, 0.7116493828156415), (0.94071511310689593, 0.60991928098248505, 0.48127645896930321), (0.75617071810890646, 0.21038062695194693, 0.22352941947824814), (0.76617071810890646, 0.21038062695194693, 0.22352941947824814), (0.77617071810890646, 0.21038062695194693, 0.22352941947824814), (0.78617071810890646, 0.21038062695194693, 0.22352941947824814), (0.125, 0.75, 0.125)]

        ax = sns.factorplot(x="Cres", y="Clip", hue="TimeWindow", row="PatientID", data=cres_smooth_df, orient='h', size=2, aspect=3.5, palette=color_palette_sz, kind="violin", cut=0, bw=.2)
        sns.despine(offset=10, trim=True)
        ax.fig.set_alpha = 0.25
        ax.fig.set_size_inches((8,18))
        # ax.fig.axes[0].fill_between([ci5, ci95],[-1, 3],[-1, 3],facecolor='gray',alpha=0.5)
        xlow = min(min(cres_smooth_data.flatten())*0.8,min(cres_smooth_data.flatten())*1.2)
        xhigh = max(max(cres_smooth_data.flatten())*0.8,max(cres_smooth_data.flatten())*1.2)

        plt.xlim([xlow,xhigh])

        # Characterize the plots approximately
        if(np.mean(cres_smooth_data.flatten()) < 0):
            if(np.mean(cres_smooth_data.flatten()) > -0.10):
                subnetwork_characterization = 'Weakly Synchronizing'
            else:
                subnetwork_characterization = 'Strongly Synchronizing'
        else:
            if(np.mean(cres_smooth_data.flatten()) > 0.10):
                subnetwork_characterization = 'Weakly Desynchronizing'
            else:
                subnetwork_characterization = 'Strongly Desynchronizing'

        plt.title('Patient %s Dilation Radius %i - %s (%s Engel %s)'%(patient_id,pipedef['dilate_radius'], fconn.upper(), subnetwork_characterization, outcome))
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_distributions.png'%(comp_dir,patient_id,unique_id, fconn))
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_distributions.svg'%(comp_dir,patient_id,unique_id, fconn))


