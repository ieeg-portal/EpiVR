import warnings

import os
import pandas as pd
import sys
import glob
import json
import time
import h5py, re

import nibabel as nib
import numpy as np
import pandas as pd
import scipy as sp
import scipy.stats
from scipy.ndimage import morphology
from scipy.signal import convolve

# Import data science modules
from scipy import interp, ndimage

# Import graphing modules
import matplotlib
import matplotlib.pyplot as plt

# Import machine learning modules
from sklearn import *
from sklearn.svm import *
from sklearn.ensemble import *
from sklearn.metrics import *
from sklearn.linear_model import *
from sklearn.model_selection import *
from sklearn.feature_selection import *
from sklearn.preprocessing import *
from sklearn.tree import *
from sklearn.externals import joblib

sys.path.append(os.path.expanduser('~/gdrive/aim3/code'))
from Echobase.Common import errors
from Echobase.Network.Metrics.globaltopo import synchronizability
from Echobase.Network.Transforms import lesion
from Echobase.Pipelines.ecog_network import *
from Echobase.Statistics.FDA.fda import *
from analysis.util import *
from analysis.util_connectivity import *
from analysis.util_virtual_resection import *

with open(os.path.expanduser('~/gdrive/aim3/code/data/DATA.json')) as json_data_file:
    data = json.load(json_data_file)

def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx+ny-2
    return (np.nanmean(x) - np.nanmean(y)) / np.sqrt(((nx-1)*np.nanstd(x,ddof=1)**2+(ny-1)*np.nanstd(y,ddof=1)**2)/dof)

def gather_cres_results(dilate_radius, fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all results.

    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        fconn: str,
            Connectivity metric
    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with cres for given dilation radius in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue

        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')
        try:
            if data['PATIENTS'][patient_id]['Cohort'] != 'DEVELOPMENT' and data['PATIENTS'][patient_id]['Cohort'] != 'VALIDATION':
                continue
        except Exception:
            continue

        # Find pipedef file
        for fn in os.listdir(comp_dir):
            if('pipedef' in fn):
                # Open pipedef
                pipedef = json.load(open('%s/%s'%(comp_dir,fn),'r'))
                # determine if correction dilation
                try:
                    if(np.float(pipedef['dilate_radius']) == np.float(dilate_radius)):
                        unique_id = fn.split('.')[4]
                        results[patient_id] = {}
                        break
                except KeyError:
                    continue

        # Open all cres
        try:
            for fn in os.listdir(comp_dir):
                if('cres.%s'%(unique_id) in fn and 'pipedef' not in fn):
                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type and seizure_type != '??'):
                        continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['control_centrality_%s'%fconn]
        except:
            continue

    return results

def gather_sozres_results(fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all results.

    Parameters
    ----------
        fconn: str,
            Connectivity metric
    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with soz cres for given dilation radius in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue

        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')
        try:
            if data['PATIENTS'][patient_id]['Cohort'] != 'DEVELOPMENT' and data['PATIENTS'][patient_id]['Cohort'] != 'VALIDATION':
                continue
        except Exception:
            continue

        # Open all sozres
        try:
            for fn in os.listdir(comp_dir):
                if('sozres' in fn):
                    if patient_id not in results.keys():
                        results[patient_id] = {}

                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type and seizure_type != '??'):
                        continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['control_centrality_%s'%fconn]
        except:
            continue

    return results

def gather_cnonres_results(dilate_radius, fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all results.

    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        fconn: str,
            Connectivity metric
    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with cres for given dilation radius in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue

        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')
        try:
            if data['PATIENTS'][patient_id]['Cohort'] != 'DEVELOPMENT' and data['PATIENTS'][patient_id]['Cohort'] != 'VALIDATION':
                continue
        except Exception:
            continue

        # Find pipedef file
        for fn in os.listdir(comp_dir):
            if('pipedef' in fn):
                # Open pipedef
                pipedef = json.load(open('%s/%s'%(comp_dir,fn),'r'))
                # determine if correction dilation
                if(np.float(pipedef['dilate_radius']) == np.float(dilate_radius)):
                    unique_id = fn.split('.')[4]
                    results[patient_id] = {}
                    break

        # Open all cres
        try:
            for fn in os.listdir(comp_dir):
                if('cres.%s'%(unique_id) in fn and 'pipedef' not in fn):
                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type and seizure_type != '??'):
                        continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['non_control_centrality_%s'%fconn]
        except:
            continue

    return results

def gather_base_sync_results(dilate_radius, fconn = 'highgamma', skip_chop = False, skip_mayo = False):
    '''
    Utility function to output a dictionary of all results.

    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        fconn: str,
            Connectivity metric

        skip_chop: bool,
            Boolean to skip all CHOP patients

        skip_mayo: bool,
            Boolean to skip all Mayo Clinical Study patients
    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with cres for given dilation radius in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue
        if(skip_chop and 'CHOP' in patient_id):
            continue
        if(skip_mayo and 'Study' in patient_id):
            continue
        try:
            if data['PATIENTS'][patient_id]['Cohort'] != 'DEVELOPMENT' and data['PATIENTS'][patient_id]['Cohort'] != 'VALIDATION':
                continue
        except Exception:
            continue
        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')

        # Find pipedef file
        for fn in os.listdir(comp_dir):
            if('pipedef' in fn):
                # Open pipedef
                pipedef = json.load(open('%s/%s'%(comp_dir,fn),'r'))
                # determine if correction dilation
                if(np.float(pipedef['dilate_radius']) == np.float(dilate_radius)):
                    unique_id = fn.split('.')[4]
                    results[patient_id] = {}
                    break

        try:
            for fn in os.listdir(comp_dir):
                if('cres.%s'%(unique_id) in fn and 'pipedef' not in fn):
                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type and seizure_type != '??'):
                        continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['base_sync_%s'%fconn]
        except:
            continue

    return results

def get_norm(dictionary_values, width=-1):
    '''
    Utility function to output a dictionary of all results normalized to standard length.

    Parameters
    ----------
        dictionary_values: dict,
            Dictionary of patients, with each patient having a dictionary of events with clips.

        width: int,
            Total number of bins across entire set of epochs (pre-seizure and seizure).

    -------
        dictionary_values_norm: dict,
            Results that contain as key patient id, and has as value another dictionary with cres
            for given dilation radius in each clip normalized to a standard length.
    '''
    dictionary_values_norm = {}
    if width == -1:
        width = 1E100
        for pt, clips in dictionary_values.items():
            for event_id, clip in clips.items():
                if clip.shape[0] < width:
                    width = clip.shape[0]

    for pt, clips in dictionary_values.items():
        dictionary_values_norm[pt] = {}
        for event_id, clip in clips.items():
            xp = np.linspace(-1.0,1.0,clip.shape[0])
            dictionary_values_norm[pt][event_id] = np.interp(np.linspace(-1.0,1,width),xp,clip)
    return dictionary_values_norm


def gather_adj_results(fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all results.

    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        fconn: str,
            Connectivity metric
    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with cres for given dilation radius in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue

        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')
        try:
            if data['PATIENTS'][patient_id]['Cohort'] != 'DEVELOPMENT' and data['PATIENTS'][patient_id]['Cohort'] != 'VALIDATION':
                continue
        except Exception:
            continue

        # Open all adj
        try:
            for fn in os.listdir(comp_dir):
                if('multiband' in fn):
                    if(patient_id not in results.keys()):
                        results[patient_id] = {}
                    print fn
                    clip_id = fn.split('.')[2]

                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type and seizure_type != '??'):
                        continue
                    results[patient_id][fn.split('.')[2]] = np.load('%s/%s'%(comp_dir,fn))['all_adj_%s'%fconn]
        except:
            print fn+'!!!!'
            continue

    return results

def get_outcome(outcome):
    """
    Function for determing poor and favorable surgical outcome.
    Parameters
    ----------
        outcome: str
            Surgical outcome as either Engel or ILAE.

    Returns
    -------
        str
            Returns either good or bad.
    """
    switcher = {
        '1': 'Good',
        '1.1': 'Good',
        '1.2': 'Good',
        '1.3': 'Good',
        '1.4': 'Good',
        '1A': 'Good',
        '1B': 'Good',
        '1C': 'Good',
        '1D': 'Good',
        'IA': 'Good',
        'IB': 'Good',
        'IC': 'Good',
        'ID': 'Good',
        '2': 'Poor',
        '2.1': 'Poor',
        '2.2': 'Poor',
        '2.3': 'Poor',
        '2.4': 'Poor',
        '2a': 'Poor',
        '2b': 'Poor',
        '2c': 'Poor',
        '2d': 'Poor',
        '3': 'Poor',
        '4': 'Poor',
        'II': 'Poor',
        'III': 'Poor',
        'IV': 'Poor',
        'ILAE1': 'Good',
        'ILAE2': 'Good',
        'ILAE3': 'Poor',
        'ILAE4': 'Poor',
        'ILAE5': 'Poor'
    }

    return switcher.get(outcome, "Good")


def get_resected_node_dx(patient_id,dilate_radius=0):
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
    event_type, events = data['PATIENTS'][patient_id]['Events'].items()[0]
    event_id = events.keys()[0]
    try:
        if(events[event_id]['STATUS'] == 'ALL_DROPOUT'):
            event_id = events.keys()[1]
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

    # Correspond lable names
    labels_dict = correspond_label_names(channels, labels)

    # Load electrodes to ignore
    ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
    for ii,node_id in enumerate(ignored_node_idx):
        pass
        # print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
    channels = list(np.delete(np.array(channels),ignored_node_idx))

    # Recorrespond label names
    labels_dict = correspond_label_names(channels, labels)

    # Generate list of resected electrodes and write to CSV file
    try:
        if(dilate_radius == 0):
            resected_node_labels = data['PATIENTS'][patient_id]['RESECTED_ELECTRODES']
        elif(dilate_radius > 0):
            resected_node_labels = data['PATIENTS'][patient_id]['RESECTED_ELECTRODES']
            for fringe_node_label in data['PATIENTS'][patient_id]['RESECTED_FRINGE_ELECTRODES']:
                resected_node_labels.append(fringe_node_label)
        else:
            blah
    except Exception:
        resected_electrodes_fn = write_resected_electrodes(patient_id, dilate_radius, data, labels_dict)

        # Load resected electrodes
        try:
            resected_nodes = map(lambda x: int(x.split(',')[0]), open(resected_electrodes_fn,'r').readlines())
            resected_node_labels = map(lambda x: x.split(',')[1].replace('\n',''), open(resected_electrodes_fn,'r').readlines())
        except IndexError:
            print 'ERROR! Resected electrodes %s does not have any electrodes. Skipping'%(resected_electrodes_fn)
            blah

    # Map the resected electrodes to channels
    clean_resected_node_labels = []
    for resected_node_label in resected_node_labels:
        if resected_node_label in ignored_node_labels:
            continue
        else:
            clean_resected_node_labels.append(resected_node_label)
    resected_node_idx = map(lambda x: labels_dict[x][0], clean_resected_node_labels)
    for ii,node_id in enumerate(resected_node_idx):
        pass
        # print 'Virtually resecting node label: %s because label %s is in the resection zone'%(channels[node_id],resected_node_labels[ii])
    return resected_node_idx, channels

