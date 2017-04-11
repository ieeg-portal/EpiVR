        # Additionally smooth cres_data
        cres_smooth_data = np.zeros(cres_data.shape)
        for clip_iter in range(cres_data.shape[0]):
            smooth_window = 5
            kernel = 1/np.float(smooth_window)*np.ones((smooth_window,1))
            kernel = kernel.flatten()
            cres_smooth_data[clip_iter,:] = convolve(cres_data[clip_iter,:],kernel,'same')



from sklearn import svm
from sklearn import metrics
from sklearn import cross_validation
import random
import matplotlib.mlab as mlab
import pandas as pd

# Create overall feature matrix and labels
feat_matrix = np.vstack((g_z[:,900:1800],b_z[:,900:1800]))
plt.figure(figsize=(10,8))
plt.imshow(feat_matrix,extent=[0,1,0,1])
title = 'Resection Zone $s_t$ for individual clips'
plt.title(title)
plt.clim([-2,2])
plt.xticks([])
plt.yticks([])
plt.show()




import os
import numpy as np

fn1 = 'HUP068_T1_19991230_electrode_labels.csv'
fn2 = '~/gdrive/aim1/results/HUP068/aim1/HUP068_T1_19991230_electrode_labels.csv'

fn1 = os.path.expanduser(fn1)
fn2 = os.path.expanduser(fn2)

lines1 = open(fn1,'r').readlines()
lines2 = open(fn2,'r').readlines()

ele1 = {}
for line in lines1:
    x = line.split(',')[1]
    y = line.split(',')[2]
    z = line.split(',')[3]
    label = line.split(',')[4]
    if(label == ''):
        continue
    ele1[label] = [float(x),float(y),float(z)]

ele2 = {}
for line in lines2:
    x = line.split(',')[1]
    y = line.split(',')[2]
    z = line.split(',')[3]
    label = line.split(',')[4]
    if(label == ''):
        continue
    ele2[label] = [float(x),float(y),float(z)]


for label in ele1.keys():
    a1 = np.array(ele1[label])
    a2 = np.array(ele2[label])
    print label, np.sqrt((a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 + (a1[2]-a2[2])**2)










import os
import glob
import nibabel as nib, numpy as np

data_dir = os.path.expanduser('~/gdrive/3T_Subjects')
comp_dir = os.path.expanduser('~/gdrive/aim1/results')
pt = 'HUP068'
fn1 = 'HUP068_T1_19991230_electrode_labels.csv'
fn2 = os.path.expanduser('~/gdrive/aim1/results/HUP068/aim1/HUP068_T1_19991230_electrode_labels.csv')


parameters = {}
lines = open('%s/%s/img/%s_img_t2.csv'%(data_dir,pt,pt)).readlines()
headers = lines[0].split(',')
values = lines[1].split(',')
for ii,header in enumerate(headers):
    parameters[header.replace('\n','')] = values[ii].replace('\n','')
t1_date = parameters['PATIENT_PREOP_T1_DATE']

# Get real files
if(not(os.path.isfile('%s/%s/img/%s/nii/%s_T1_%s.nii.gz'%(data_dir,pt,t1_date,pt,t1_date)))):
    if(not(os.path.isfile('%s/%s/img/%s/nii/%s_T1_%s_1.nii.gz'%(data_dir,pt,t1_date,pt,t1_date)))):
        if(not(os.path.isfile('%s/%s/img/%s/nii/%s_T1_POST_%s.nii.gz'%(data_dir,pt,t1_date,pt,t1_date)))):
            raise BaseException
        else:
            t1 = '%s/%s/img/%s/nii/%s_T1_POST_%s.nii.gz'%(data_dir,pt,t1_date,pt,t1_date)
    else:
        t1 = '%s/%s/img/%s/nii/%s_T1_%s_1.nii.gz'%(data_dir,pt,t1_date,pt,t1_date)
else:
    t1 = '%s/%s/img/%s/nii/%s_T1_%s.nii.gz'%(data_dir,pt,t1_date,pt,t1_date)

print pt, t1

# Make ele1
lines = open(fn1,'r').readlines()

t1_nii = nib.load(t1)
t1_data = t1_nii.get_data()
t1_affine = t1_nii.get_affine()

ele = np.zeros(t1_data.shape)
for line in lines:
    if(line.split(',')[1] == ''):
        continue
    eid = int(float(line.split(',')[0]))
    x = int(float(line.split(',')[1]))
    y = int(float(line.split(',')[2]))
    z = int(float(line.split(',')[3]))
    radius = 2
    ele[x-radius:x+radius,y-radius:y+radius,z-radius:z+radius] = eid

nib.save(nib.Nifti1Image(ele,t1_affine),'test1.nii.gz')

# Make ele1
lines = open(fn2,'r').readlines()

t1_nii = nib.load(t1)
t1_data = t1_nii.get_data()
t1_affine = t1_nii.get_affine()

ele = np.zeros(t1_data.shape)
for line in lines:
    eid = int(float(line.split(',')[0]))
    x = int(float(line.split(',')[1]))
    y = int(float(line.split(',')[2]))
    z = int(float(line.split(',')[3]))
    radius = 2
    ele[x-radius:x+radius,y-radius:y+radius,z-radius:z+radius] = eid

nib.save(nib.Nifti1Image(ele,t1_affine),'test2.nii.gz')














param['time_band'] = 5.
param['n_taper'] = 9
param['AlphaTheta_Band'] = [5., 15.]
param['Beta_Band'] = [15., 25.]
param['LowGamma_Band'] = [30., 40.]
param['HighGamma_Band'] = [95., 105.]

from mtspec import mt_coherence, mtspec


for k1 in range(0,evData.shape[0],100):
    k2 = k1+100

    n_samp, n_chan = evData[k1:k2,:].shape
    for n1 in range(n_chan):
        for n2 in range(n_chan):
            if(n1 == n2):
                continue
            print k1,n1,n2
            out = mt_coherence(1.0/Fs,evData[k1:k2, n1],evData[k1:k2, n2],5,9, int(n_samp/2.), 0.95, iadapt=1, cohe=True, freq=True)
















### DEBUGGING CONNECTIVITY STUFF

import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

patient_id = 'CHOP08'
event_id = '1'


dilate_radius = 0
epoch_length = 1

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

event_type, events = data['PATIENTS'][patient_id]['Events'].items()[0]
print 'Computing multiband connectivity for clip %s'%event_id
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

# Get smaller window if ictal clip has dropout issues
try:
    if(events[event_id]['STATUS'] == 'ALL_DROPOUT'):
        pass
    elif(events[event_id]['STATUS'] == 'DROPOUT'):
        start_epoch = 240
        end_epoch = 360-1
    elif(events[event_id]['STATUS'] == 'LATE_DROPOUT'):
        start_epoch = 0
        end_epoch = 600-1
    else:
        start_epoch = 0
        end_epoch = epochs-1
except Exception:
    start_epoch = 0
    end_epoch = epochs-1

# Create parallel jobs for each block
jobs = []
for epoch in range(start_epoch,end_epoch+1):
    jobs.append((epoch, epoch_length, Fs, evData))
job = jobs[756]

# Get job details
epoch,epoch_length,Fs,evData = job

# Get clip
data_clip = evData[epoch*epoch_length*Fs:(epoch+1)*epoch_length*Fs,:]
N = data_clip.shape[1]

data = data_clip
fs = Fs

# Parameter set
param = {}
param['time_band'] = 5.
param['n_taper'] = 9
param['AlphaTheta_Band'] = [5., 15.]
param['Beta_Band'] = [15., 25.]
param['LowGamma_Band'] = [30., 40.]
param['HighGamma_Band'] = [95., 105.]

# Build pipeline
data_hat = reref.common_avg_ref(data)
data = data_hat
time_band = param['time_band']
n_taper = param['n_taper']
cf = param['AlphaTheta_Band']

# Get data attributes
n_samp, n_chan = data.shape
triu_ix, triu_iy = np.triu_indices(n_chan, k=1)

# Initialize adjacency matrix
adj = np.zeros((n_chan, n_chan))




from __future__ import division
import numpy as np
from mtspec import mt_coherence, mtspec
from scipy.signal import coherence

n1 = 2
n2 = 3

print n1,n2
out = mt_coherence(1.0/fs,
                   data[:, n1],
                   data[:, n2],
                   time_band,
                   n_taper,
                   int(n_samp/2.), 0.95,
                   iadapt=1,
                   cohe=True, freq=True)












############ COMPARING WITH OLD RESULTS

import os
import sys
import glob
import json

import numpy as np
import pandas as pd
import scipy.io as io
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as plt_cm
from matplotlib import rcParams
import scipy.stats as stats

path_CoreData = '/gdrive/public/USERS/lkini/resectionProject/CORE.ieeg_resection.multiinst'
path_RsrchData = '/gdrive/public/USERS/lkini/resectionProject/RSRCH.vresect_validation'
path_InpData = path_RsrchData + '/e01-FuncNetw'
path_ExpData = path_RsrchData + '/e02-VirtualResect'

for path in [path_CoreData, path_RsrchData, path_ExpData]:
    if not os.path.exists(path):
        print('Path: {}, does not exist'.format(path))
        os.makedirs(path)

        # Get clinical_metadata file
df_meta = h5py.File('{}/clinical_metadata.mat'.format(path_CoreData), 'r')
meta_subj = [''.join(unichr(c) for c in df_meta[r])
             for r in df_meta['subject']['ID'][:, 0]]

# Get subject ids for processed data
subj_ids = [path.split('/')[-1] for path in glob.glob('{}/Subjects/*'.format(path_CoreData))]

# Populate dict
subj_dict = {}
for subj_id in subj_ids:
    # Retrieve channel indices for resection zone
    meta_ix = np.flatnonzero(np.array(meta_subj) == subj_id)
    if len(meta_ix) == 0:
        continue

    meta_ref = df_meta['subject']['Channels'][meta_ix, 0]
    if len(df_meta[meta_ref].shape) < 2:
        continue
    chan_lbl = [''.join(unichr(c) for c in df_meta[rr])
                for rr in df_meta[meta_ref][0, :]]
    all_chan_ix = np.arange(len(chan_lbl))

    meta_ref = df_meta['subject']['Channels_Sz'][meta_ix, 0]
    if len(df_meta[meta_ref].shape) < 2:
        soz_chan=[]
    else:
        soz_chan = [''.join(unichr(c) for c in df_meta[rr])
                for rr in df_meta[meta_ref][0, :]]

    meta_ref = df_meta['subject']['Channels_Resect'][meta_ix, 0]
    if len(df_meta[meta_ref].shape) < 2:
        continue
    resect_lbl = [''.join(unichr(c) for c in df_meta[rr])
                  for rr in df_meta[meta_ref][0, :]]
    resect_ix = np.array([chan_ix for chan_ix in all_chan_ix
                          if chan_lbl[chan_ix] in resect_lbl])

    meta_ref = df_meta['subject']['Outcome'][meta_ix, 0]
    outcome = df_meta[meta_ref][0,0]

    subj_dict[subj_id] = {'dyne': [],
                          'resect_ix': resect_ix,
                          'channel_ix': all_chan_ix,
                          'chan_lbl':chan_lbl,
                          'soz_chan':soz_chan,
                          'outcome': outcome}

    # Append dyne path for each subject event
    for dyne_path in glob.glob('{}/{}.*.*.dyne_log.csv'.format(path_InpData, subj_id)):
        subj_dict[subj_id]['dyne'].append({'log': dyne_path,
                                           'out': '{}.dyne_output.hdf'.format(dyne_path.split('.dyne_log')[0])})



proc_list = [dyne_ev
             for subj_id in subj_dict.keys()
             for dyne_ev in subj_dict[subj_id]['dyne']]

# Study the data
data = []
outcome = []
subjx = []
for proc_item in proc_list:
    try:
        subj_id, epoch_id, event_id, _, _ = proc_item['log'].split('/')[-1].split('.')
        file_out = '{}/{}.{}.{}.real_virtual_resect.npz'.format(path_ExpData, subj_id, epoch_id, event_id)

        adj = np.load(file_out)
        sync_out = adj['sync']
        vr_rz_out = adj['vr_rz']
        nvr_rz_out = adj['vr_nrz']
        print subj_id, epoch_id, event_id
        data.append(np.hstack((sync_out,vr_rz_out,nvr_rz_out)))
        subjx.append(subj_id+' ' + event_id)
        outcome.append(subj_dict[subj_id]['outcome'])
    except IOError:
        pass
data = np.array(data)
outcome = np.array(outcome)
subjx = np.array(subjx)


g = data[outcome<2,:]
b = data[outcome>=2,:]
print subjx[outcome<2]
print subjx[outcome>=2]
g_z = []
b_z = []
for k in g:
    g_z.append(np.hstack((stats.zscore(k[0:900],axis=0),stats.zscore(k[900:1800],axis=0),stats.zscore(k[1800:],axis=0))))
for k in b:
    b_z.append(np.hstack((stats.zscore(k[0:900],axis=0),stats.zscore(k[900:1800],axis=0),stats.zscore(k[1800:],axis=0))))
g_z = np.array(g_z)
b_z = np.array(b_z)

# Create overall feature matrix and labels
feat_matrix = np.vstack((g_z[:,900:1800],b_z[:,900:1800]))
# plt.figure(figsize=(10,8))
plt.imshow(feat_matrix,extent=[0,1,0,1])
title = 'Resection Zone $s_t$ for individual clips'
plt.title(title)
plt.clim([-2,2])
plt.xticks([])
plt.yticks([])
plt.show()

res = []
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP073/aim3/HUP073.Ictal.2.cres.16d31703-7226-42e0-a9d4-b78b805cd174.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP073/aim3/HUP073.Ictal.4.cres.16d31703-7226-42e0-a9d4-b78b805cd174.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP073/aim3/HUP073.Ictal.3.cres.16d31703-7226-42e0-a9d4-b78b805cd174.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP073/aim3/HUP073.Ictal.5.cres.16d31703-7226-42e0-a9d4-b78b805cd174.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP073/aim3/HUP073.Ictal.1.cres.16d31703-7226-42e0-a9d4-b78b805cd174.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP45/aim3/CHOP45.Ictal.2.cres.b9a7014c-a045-48d3-b7ff-29d03c0a380c.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP45/aim3/CHOP45.Ictal.1.cres.b9a7014c-a045-48d3-b7ff-29d03c0a380c.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP45/aim3/CHOP45.Ictal.4.cres.b9a7014c-a045-48d3-b7ff-29d03c0a380c.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP45/aim3/CHOP45.Ictal.3.cres.b9a7014c-a045-48d3-b7ff-29d03c0a380c.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP45/aim3/CHOP45.Ictal.5.cres.b9a7014c-a045-48d3-b7ff-29d03c0a380c.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP068/aim3/HUP068.Ictal.5.cres.4b42cec8-2dba-4f48-8af7-7aed985fb0b8.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP068/aim3/HUP068.Ictal.1.cres.4b42cec8-2dba-4f48-8af7-7aed985fb0b8.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP068/aim3/HUP068.Ictal.3.cres.4b42cec8-2dba-4f48-8af7-7aed985fb0b8.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP068/aim3/HUP068.Ictal.4.cres.4b42cec8-2dba-4f48-8af7-7aed985fb0b8.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP068/aim3/HUP068.Ictal.2.cres.4b42cec8-2dba-4f48-8af7-7aed985fb0b8.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP087/aim3/HUP087.Ictal.2.cres.4b5385bf-edd0-47cf-9c32-5807cc4b3a23.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP087/aim3/HUP087.Ictal.1.cres.4b5385bf-edd0-47cf-9c32-5807cc4b3a23.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP064/aim3/HUP064.Ictal.1.cres.ed1ac3bc-a8c9-4707-873a-7e4a955a655f.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP065/aim3/HUP065.Ictal.3.cres.f66c25d9-e974-40f9-8947-e58d3ac12fdc.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP065/aim3/HUP065.Ictal.1.cres.f66c25d9-e974-40f9-8947-e58d3ac12fdc.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP065/aim3/HUP065.Ictal.2.cres.f66c25d9-e974-40f9-8947-e58d3ac12fdc.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP078/aim3/HUP078.Ictal.2.cres.daba8715-5bc0-4a5a-87b3-c14879b18248.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP078/aim3/HUP078.Ictal.4.cres.daba8715-5bc0-4a5a-87b3-c14879b18248.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP078/aim3/HUP078.Ictal.1.cres.daba8715-5bc0-4a5a-87b3-c14879b18248.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP078/aim3/HUP078.Ictal.3.cres.daba8715-5bc0-4a5a-87b3-c14879b18248.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/HUP078/aim3/HUP078.Ictal.5.cres.daba8715-5bc0-4a5a-87b3-c14879b18248.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP42/aim3/CHOP42.Ictal.3.cres.ab7991b9-8f96-4b97-bfff-7bd96e0196aa.npz'))['control_centrality_broadband_CC'])
res.append(np.nan*np.ones((900,)))
res.append(np.nan*np.ones((900,)))
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP40/aim3/CHOP40.Ictal.2.cres.7ed66a4e-7ed4-44aa-b77a-e52f44fdb02f.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP40/aim3/CHOP40.Ictal.1.cres.7ed66a4e-7ed4-44aa-b77a-e52f44fdb02f.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP46/aim3/CHOP46.Ictal.1.cres.de89c35c-811d-48ee-aa2f-f3b555428ce7.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP46/aim3/CHOP46.Ictal.3.cres.de89c35c-811d-48ee-aa2f-f3b555428ce7.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP46/aim3/CHOP46.Ictal.2.cres.de89c35c-811d-48ee-aa2f-f3b555428ce7.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP38/aim3/CHOP38.Ictal.4.cres.49dd712c-f805-49fa-8db8-6904e9c4f1cb.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP38/aim3/CHOP38.Ictal.2.cres.49dd712c-f805-49fa-8db8-6904e9c4f1cb.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP14/aim3/CHOP14.Ictal.1.cres.dc77bd6d-adfe-45d9-9766-0b4431d8942c.npz'))['control_centrality_broadband_CC'])
res.append(np.nan*np.ones((900,)))
res.append(np.nan*np.ones((900,)))
res.append(np.nan*np.ones((900,)))
res.append(np.nan*np.ones((900,)))
res.append(np.nan*np.ones((900,)))
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP08/aim3/CHOP08.Ictal.2.cres.6a1f5022-7c92-494a-827c-0fdd82390f4a.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP08/aim3/CHOP08.Ictal.7.cres.6a1f5022-7c92-494a-827c-0fdd82390f4a.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP08/aim3/CHOP08.Ictal.1.cres.6a1f5022-7c92-494a-827c-0fdd82390f4a.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP08/aim3/CHOP08.Ictal.5.cres.6a1f5022-7c92-494a-827c-0fdd82390f4a.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.11.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.4.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.7.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.10.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.1.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.12.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.6.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.9.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP27/aim3/CHOP27.Ictal.8.cres.ea1bad6e-f680-4ccd-9a54-58b646970fd2.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP20/aim3/CHOP20.Ictal.1.cres.76f20b0e-b34a-4f10-aad4-c71b2c23a67e.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP20/aim3/CHOP20.Ictal.4.cres.76f20b0e-b34a-4f10-aad4-c71b2c23a67e.npz'))['control_centrality_broadband_CC'])
res.append(np.load(os.path.expanduser('~/gdrive/aim3/results/CHOP20/aim3/CHOP20.Ictal.3.cres.76f20b0e-b34a-4f10-aad4-c71b2c23a67e.npz'))['control_centrality_broadband_CC'])

res = np.array(res)
newvr_g_z = scipy.stats.zscore(res[outcome<2,:], axis=1)
newvr_b_z = scipy.stats.zscore(res[outcome>=2,:], axis=1)
newvr_feat_matrix = np.vstack((newvr_g_z,newvr_b_z))

plt.imshow(newvr_feat_matrix,extent=[0,1,0,1])
title = 'Resection Zone $s_t$ for individual clips'
plt.title(title)
plt.clim([-2,2])
plt.xticks([])
plt.yticks([])
plt.show()
















    outcome_type = 'Good'
    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    y_labels = []
    all_data = np.array(())
    for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
        outcome = data['PATIENTS'][patient_id]['Outcome']
        if(get_outcome(outcome) == outcome_type):
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                try:
                    all_data = np.hstack((all_data,np.reshape(clip_data,(900,1))))
                except:
                    all_data = np.reshape(clip_data,(900,1))
                y_labels.append('%s.clip.%s'%(patient_id,clip))

    # Apply z score normalization
    for k in range(all_data.shape[1]):
        tmp = all_data[:,k].squeeze()
        all_data[~np.isnan(tmp),k] = np.abs(scipy.stats.zscore(tmp[~np.isnan(tmp)]))

    # plt.figure()
    # plt.errorbar(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1),yerr=yerr, color=colors, capsize=5, capthick=2)
    plt.plot(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1), alpha=0.5)
    plt.hold(True)

    outcome_type = 'Poor'

    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    y_labels = []
    all_data = np.array(())
    for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
        outcome = data['PATIENTS'][patient_id]['Outcome']
        if(get_outcome(outcome) == outcome_type):
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                try:
                    all_data = np.hstack((all_data,np.reshape(clip_data,(900,1))))
                except:
                    all_data = np.reshape(clip_data,(900,1))
                y_labels.append('%s.clip.%s'%(patient_id,clip))

    # Apply z score normalization
    for k in range(all_data.shape[1]):
        tmp = all_data[:,k].squeeze()
        all_data[~np.isnan(tmp),k] = np.abs(scipy.stats.zscore(tmp[~np.isnan(tmp)]))

    # plt.figure()
    # plt.errorbar(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1),yerr=yerr, color=colors, capsize=5, capthick=2)
    plt.plot(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1), alpha=0.5)

    plt.show()







for fconn in ['broadband_CC','alphatheta','beta','lowgamma','highgamma']:
    for outcome_type in ['Good', 'Poor', 'Both']:
        for dilate_radius in [0, -5, 5]:
            plot_all_cres_time(dilate_radius, outcome_type, fconn)


for fconn in ['broadband_CC','alphatheta','beta','lowgamma','highgamma']:
    for dilate_radius in [0, -5, 5]:
        plot_all_cres_box(dilate_radius, fconn)


























############################
### COMPUTE cc ROZ and % of network resected

import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
    if(patient_id == 'TEST1'):
        continue
        # Initialize output
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
        unique_id = str(uuid.uuid4())
        for event_id in events.keys():
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

            except Exception:
                resected_electrodes_fn = write_resected_electrodes(patient_id, dilate_radius, data, labels_dict)

                # Load resected electrodes
                try:
                    resected_nodes = map(lambda x: int(x.split(',')[0]), open(resected_electrodes_fn,'r').readlines())
                    resected_node_labels = map(lambda x: x.split(',')[1].replace('\n',''), open(resected_electrodes_fn,'r').readlines())
                except IndexError:
                    print 'ERROR! Resected electrodes %s does not have any electrodes. Skipping'%(resected_electrodes_fn)

            # Map the resected electrodes to channels
            clean_resected_node_labels = []
            for resected_node_label in resected_node_labels:
                if resected_node_label in ignored_node_labels:
                    continue
                else:
                    clean_resected_node_labels.append(resected_node_label)
            resected_node_idx = map(lambda x: labels_dict[x][0], clean_resected_node_labels)
            break
        break
    network_resection_fraction = len(resected_node_idx)/np.float(len(labels))
    print patient_id + ' %0.4f'%network_resection_fraction







############################
### COMPUTE cc ROZ and % of network resected

import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

for patient_id in sorted(os.listdir(os.path.expanduser(data['COMP_DIR']))):
    if(patient_id == 'TEST1'):
        continue
        # Initialize output
    resected_electrodes = []

    # Get resection image
    try:
        resection_nii = nib.load(os.path.expanduser(data['PATIENTS'][patient_id]['RESECTION_IMAGE']))
        outcome = data['PATIENTS'][patient_id]['Outcome']
        resection_data = resection_nii.get_data()
        # Threshold and binarize the resection mask
        resection_data[resection_data > 0.2] = 1
        resection_data[resection_data <= 0.2] = 0
        resection_affine = resection_nii.get_affine()
        print patient_id + ' %0.4f cc'%(np.prod(resection_nii.header.get_zooms())*np.count_nonzero(resection_data[:])*0.001) + ' %s'%(get_outcome(outcome))
    except KeyError:
        print patient_id + ' No image'




############################
# Do t test on cc ROZ and outcome
good = [14.2879,54.4922,71.0793,14.1635,34.3072,36.5059,12.156,32.4287,14.1764,17.3397,40.8839,18.0851,28.2097,44.0665,20.4636,9.5827,19.3115]
poor = [53.1368,12.604,21.2765,12.6188,43.1842]
df = pd.DataFrame(np.vstack((zip(good,['good']*len(good)),zip(poor,['poor']*len(poor)))),columns=['cc','o'])
df.cc = df.cc.astype(np.float)

ax = sns.boxplot(x='o',y='cc',data=df, palette='Set2');
plt.xlabel('Surgical Outcome')
plt.ylabel('Volume Resected (cc.)')
plt.title('Comparison of resection volumes (cc.) between Surgical Outcomes')
plt.show()



















##############
# Figure out Chop and mayo VR issues
import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

patient_id = 'CHOP08'
event_id = '1'


dilate_radius = 0
epoch_length = 1

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
        blah
        pass # unusable clip
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
    print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
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
except Exception:
    resected_electrodes_fn = write_resected_electrodes(patient_id, dilate_radius, data, labels_dict)

    # Load resected electrodes
    try:
        resected_nodes = map(lambda x: int(x.split(',')[0]), open(resected_electrodes_fn,'r').readlines())
        resected_node_labels = map(lambda x: x.split(',')[1].replace('\n',''), open(resected_electrodes_fn,'r').readlines())
    except IndexError:
        print 'ERROR! Resected electrodes %s does not have any electrodes. Skipping'%(resected_electrodes_fn)

# Map the resected electrodes to channels
clean_resected_node_labels = []
for resected_node_label in resected_node_labels:
    if resected_node_label in ignored_node_labels:
        pass
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
    print epoch
    if(np.isnan(all_adj_alphatheta[:,:,epoch]).any()):
        control_centrality_alphatheta[epoch] = np.nan
    else:
        control_centrality_alphatheta[epoch] = region_control(all_adj_alphatheta[:,:,epoch],resected_node_idx)
    if(np.isnan(all_adj_beta[:,:,epoch]).any()):
        control_centrality_beta[epoch] = np.nan
    else:
        control_centrality_beta[epoch] = region_control(all_adj_beta[:,:,epoch],resected_node_idx)
    if(np.isnan(all_adj_lowgamma[:,:,epoch]).any()):
        control_centrality_lowgamma[epoch] = np.nan
    else:
        control_centrality_lowgamma[epoch] = region_control(all_adj_lowgamma[:,:,epoch],resected_node_idx)
    if(np.isnan(all_adj_highgamma[:,:,epoch]).any()):
        control_centrality_highgamma[epoch] = np.nan
    else:
        control_centrality_highgamma[epoch] = region_control(all_adj_highgamma[:,:,epoch],resected_node_idx)
    if(np.isnan(all_adj_broadband_CC[:,:,epoch]).any()):
        control_centrality_broadband_CC[epoch] = np.nan
    else:
        control_centrality_broadband_CC[epoch] = region_control(all_adj_broadband_CC[:,:,epoch],resected_node_idx)






############ p hacking
iter_id = 1
dilate_radius = 0
results = []

for skip_chop in [True, False]:
    for zscore in [True, False]:
        for feature in ['abs_cres','cres','cres_only_after_EEC','normalized_cres_ai_abs','normalized_cres_ai_no_abs']:
            for fconn in ['broadband_CC','alphatheta','beta','lowgamma','highgamma']:
                for window in [30,60,120,300]:
                    p1 = plot_all_cres_box(str(iter_id), dilate_radius, skip_chop = skip_chop, skip_mayo = True, skip_hup = False, zscore = zscore, feature = feature, fconn = fconn, window = window)
                    p2 = plot_all_patient_cres_box(str(iter_id), dilate_radius, skip_chop = skip_chop, skip_mayo = True, skip_hup = False, zscore = zscore, feature = feature, fconn = fconn, window = window)
                    results.append({'iter_id':iter_id,'skip_chop':skip_chop,'zscore':zscore,'feature':feature,'fconn':fconn,'p1':p1,'p2':p2,'window':window})
                    iter_id += 1
                    print iter_id

txt = ''
for res in results:
    txt += ','.join([str(res['iter_id']),str(res['skip_chop']),str(res['zscore']),res['feature'],res['fconn'],str(res['window']),str(res['p1']),str(res['p2'])]) + '\n'

open('res.csv','w').write(txt)







### GENERATE combined HUP 65/86 figure

import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *


patient_idx = ['HUP065','HUP086']
outcomex = ['Good','Poor']
fconn = 'broadband_CC'
window = 30


fig, ax = plt.subplots(5,1,sharex=True)
iter_id = 0
for patient_id_ii,patient_id in enumerate(patient_idx):
    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    times = np.arange(-300.0,600,window)*1.0/60

    for clip_id in all_cres[patient_id].keys():
        cres = all_cres[patient_id][clip_id]
        if('abs' in feature and 'no_abs' not in feature):
            cres = np.abs(cres)
        # plt.subplot(iter_id,1)
        box = ax[iter_id].boxplot(np.reshape(cres,(cres.shape[0]/window,window)).T);
        ax[iter_id].set_title('Patient %s (%s Surgical Outcome) - Clip %s'%(patient_id,outcomex[patient_id_ii],clip_id))
        ax[iter_id].set_ylim([-0.1,0.5])
        for patch in box['whiskers']:
            patch.set_c('black')

        for patch, time in zip(box['boxes'], times):
            if(time < 0):
                patch.set_c('black')
                patch.set_fillstyle('full')
            else:
                patch.set_c('red')
                patch.set_fillstyle('full')
        plt.xticks(map(lambda x: int(x+1), range(times.shape[0]))[::2],times[::2])
        # plt.yticks([])
        plt.xlabel('Time (min)')

        plt.hold(True)
        iter_id += 1
# plt.title('Patient %s (Poor Outcome) - Cres(t) in broadband cross-correlation connectivity \nover 30 second windows'%(patient_id))
plt.show()






### GENERATE SEM plot
import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *


for fconn in ['broadband_CC','alphatheta','beta','lowgamma','highgamma']:
    window = 5
    dilate_radius = 0
    zscore = False
    skip_chop = False

    # EVENT Event-based
    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    all_good_data = np.array(())
    all_poor_data = np.array(())
    for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
        if(skip_chop and 'CHOP' in patient_id):
            continue
        if(skip_mayo and 'Study' in patient_id):
            continue
        if(skip_hup and 'HUP' in patient_id):
            continue
        outcome = data['PATIENTS'][patient_id]['Outcome']
        if(get_outcome(outcome) == 'Good'):
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                # if(np.isnan(clip_data).any()):
                #     continue
                if(clip_data.shape[0] == 901):
                    clip_data = clip_data[:900]
                try:
                    all_good_data = np.hstack((all_good_data,np.reshape(clip_data,(900,1))))
                except:
                    all_good_data = np.reshape(clip_data,(900,1))
        if(get_outcome(outcome) == 'Poor'):
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                # if(np.isnan(clip_data).any()):
                #     continue
                try:
                    all_poor_data = np.hstack((all_poor_data,np.reshape(clip_data,(900,1))))
                except:
                    all_poor_data = np.reshape(clip_data,(900,1))

    # Apply z score normalization
    if(zscore):
        for k in range(all_good_data.shape[1]):
            tmp = all_good_data[:,k].squeeze()
            all_good_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
        for k in range(all_poor_data.shape[1]):
            tmp = all_poor_data[:,k].squeeze()
            all_poor_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])

    g = all_good_data.shape[1]
    b = all_poor_data.shape[1]
    p_clip = curve_test(np.hstack((all_good_data,all_poor_data)),np.arange(0,g),np.arange(g,g+b))



    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    all_good_data = np.array(())
    all_poor_data = np.array(())
    for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
        if(skip_chop and 'CHOP' in patient_id):
            continue
        if(skip_mayo and 'Study' in patient_id):
            continue
        if(skip_hup and 'HUP' in patient_id):
            continue
        outcome = data['PATIENTS'][patient_id]['Outcome']
        avg_data = np.array(())
        for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
            if(clip_data.shape[0] == 901):
                clip_data = clip_data[:900]
            try:
                avg_data = np.hstack((avg_data,np.reshape(clip_data,(900,1))))
            except Exception:
                avg_data = np.reshape(clip_data,(900,1))

        avg_data = np.nanmean(avg_data,axis=1)
        if(get_outcome(outcome) == 'Good'):
            try:
                all_good_data = np.hstack((all_good_data,np.reshape(avg_data,(900,1))))
            except Exception:
                all_good_data = np.reshape(avg_data,(900,1))
        if(get_outcome(outcome) == 'Poor'):
            try:
                all_poor_data = np.hstack((all_poor_data,np.reshape(avg_data,(900,1))))
            except Exception:
                all_poor_data = np.reshape(avg_data,(900,1))

    # Apply z score normalization
    if(zscore):
        for k in range(all_good_data.shape[1]):
            tmp = all_good_data[:,k].squeeze()
            all_good_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
        for k in range(all_poor_data.shape[1]):
            tmp = all_poor_data[:,k].squeeze()
            all_poor_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])

    g = all_good_data.shape[1]
    b = all_poor_data.shape[1]
    p_patient = curve_test(np.hstack((all_good_data,all_poor_data)),np.arange(0,g),np.arange(g,g+b))




    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    good_cres = []
    bad_cres = []

    for patient_id in all_cres.keys():
        if('CHOP' in patient_id and skip_chop):
            continue
        # Load patient outcome
        outcome = get_outcome(data['PATIENTS'][patient_id]['Outcome'])

        avg_data = np.array(())
        for clip_id in all_cres[patient_id].keys():
            cres = all_cres[patient_id][clip_id]
            if(cres.shape[0] == 901):
                cres = cres[:900]
            try:
                avg_data = np.hstack((avg_data,np.reshape(cres,(900,1))))
            except Exception:
                avg_data = np.reshape(cres,(900,1))
        avg_data = np.nanmean(avg_data,axis=1)

        if(outcome == 'Good'):
            good_cres.append(cres)
        else:
            bad_cres.append(cres)
    good_cres = np.array(good_cres)
    bad_cres = np.array(bad_cres)

    # Apply z score normalization
    if(zscore):
        for k in range(good_cres.shape[1]):
            tmp = good_cres[:,k].squeeze()
            good_cres[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
        for k in range(bad_cres.shape[1]):
            tmp = bad_cres[:,k].squeeze()
            bad_cres[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])

    if(feature == 'abs_cres'):
        good = []
        poor = []
        for w in np.arange(0,good_cres.shape[1]+1,window):
            good.append(np.nanmean(np.abs(good_cres[:,w:w+window]),axis=1))
            poor.append(np.nanmean(np.abs(bad_cres[:,w:w+window]),axis=1))
        good = np.array(good)
        poor = np.array(poor)
    elif(feature == 'cres'):
        good = []
        poor = []
        for w in np.arange(0,good_cres.shape[1]+1,window):
            good.append(np.nanmean(good_cres[:,w:w+window],axis=1))
            poor.append(np.nanmean(bad_cres[:,w:w+window],axis=1))
        good = np.array(good)
        poor = np.array(poor)
    else:
        raise

    times = np.arange(-300.0,600,window)*1.0/60

    cres = []
    for k in good:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(good.T,nan_policy='omit')
    error = error[:-1]
    plt.plot(times,cres,'b-')
    plt.fill_between(times,cres-error,cres+error,facecolor='blue',alpha=0.25)
    plt.hold(True)

    cres = []
    for k in poor:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(poor.T,nan_policy='omit')
    error = error[:-1]
    plt.plot(times,cres,'r-')
    plt.fill_between(times,cres-error,cres+error,facecolor='red',alpha=0.25)

    # plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('$c_res(t)$')
    plt.xlim([-5.0,10.0])
    plt.title('$c_{res}(t)$ in %s connectivity of all patients ($p$ = %0.3f)'%(fconn,p_patient))
    plt.legend(['Good Outcome','Poor Outcome'])
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    fig.savefig(os.path.expanduser('~/gdrive/aim3/fig/FDA_patient_%s.png'%(fconn)),dpi=100,bbox_inches='tight')
    plt.hold(False)



    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    all_good_data = np.array(())
    all_poor_data = np.array(())
    for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
        if(skip_chop and 'CHOP' in patient_id):
            continue
        if(skip_mayo and 'Study' in patient_id):
            continue
        if(skip_hup and 'HUP' in patient_id):
            continue
        outcome = data['PATIENTS'][patient_id]['Outcome']
        if(get_outcome(outcome) == 'Good'):
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                # if(np.isnan(clip_data).any()):
                #     continue
                if(clip_data.shape[0] == 901):
                    clip_data = clip_data[:900]
                try:
                    all_good_data = np.hstack((all_good_data,np.reshape(clip_data,(900,1))))
                except:
                    all_good_data = np.reshape(clip_data,(900,1))
        if(get_outcome(outcome) == 'Poor'):
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                # if(np.isnan(clip_data).any()):
                #     continue
                try:
                    all_poor_data = np.hstack((all_poor_data,np.reshape(clip_data,(900,1))))
                except:
                    all_poor_data = np.reshape(clip_data,(900,1))

    # Apply z score normalization
    if(zscore):
        for k in range(all_good_data.shape[1]):
            tmp = all_good_data[:,k].squeeze()
            all_good_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
        for k in range(all_poor_data.shape[1]):
            tmp = all_poor_data[:,k].squeeze()
            all_poor_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])


    if(feature == 'abs_cres'):
        good = []
        poor = []
        for w in np.arange(0,good_cres.shape[1]+1,window):
            good.append(np.nanmean(np.abs(all_good_data[w:w+window,:]),axis=0))
            poor.append(np.nanmean(np.abs(all_poor_data[w:w+window,:]),axis=0))
        good = np.array(good)
        poor = np.array(poor)
    elif(feature == 'cres'):
        good = []
        poor = []
        for w in np.arange(0,good_cres.shape[1]+1,window):
            good.append(np.nanmean(all_good_data[w:w+window,:],axis=0))
            poor.append(np.nanmean(all_poor_data[w:w+window,:],axis=0))
        good = np.array(good)
        poor = np.array(poor)
    else:
        raise

    times = np.arange(-300.0,600,window)*1.0/60

    cres = []
    for k in good:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(good.T,nan_policy='omit')
    error = error[:-1]
    plt.plot(times,cres,'b-')
    plt.fill_between(times,cres-error,cres+error,facecolor='blue',alpha=0.25)
    plt.hold(True)

    cres = []
    for k in poor:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(poor.T,nan_policy='omit')
    error = error[:-1]
    plt.plot(times,cres,'r-')
    plt.fill_between(times,cres-error,cres+error,facecolor='red',alpha=0.25)

    # plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('$c_{res}(t)$')
    plt.xlim([-5.0,10.0])
    plt.title('$c_{res}(t)$ in %s connectivity of all clips ($p$ = %0.3f)'%(fconn,p_clip))
    plt.legend(['Good Outcome','Poor Outcome'])
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    fig.savefig(os.path.expanduser('~/gdrive/aim3/fig/FDA_clip_%s.png'%(fconn)),dpi=100,bbox_inches='tight')
    plt.hold(False)




###########################################################################################################
for skip_chop in [True,False]:
    for fconn in ['broadband_CC','alphatheta','beta','lowgamma','highgamma']:
        # All cres
        all_cres = gather_results(dilate_radius, fconn)

        comp_dir = os.path.expanduser(data['COMP_DIR'])

        all_good_data = np.array(())
        all_poor_data = np.array(())
        for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
            if(skip_chop and 'CHOP' in patient_id):
                continue
            if(skip_mayo and 'Study' in patient_id):
                continue
            if(skip_hup and 'HUP' in patient_id):
                continue
            outcome = data['PATIENTS'][patient_id]['Outcome']
            avg_data = np.array(())
            for clip, clip_data in sorted(clips.items(), key=lambda x:x[0]):
                if(clip_data.shape[0] == 901):
                    clip_data = clip_data[:900]
                try:
                    avg_data = np.hstack((avg_data,np.reshape(clip_data,(900,1))))
                except Exception:
                    avg_data = np.reshape(clip_data,(900,1))

            avg_data = np.nanmean(avg_data,axis=1)
            if(get_outcome(outcome) == 'Good'):
                try:
                    all_good_data = np.hstack((all_good_data,np.reshape(avg_data,(900,1))))
                except Exception:
                    all_good_data = np.reshape(avg_data,(900,1))
            if(get_outcome(outcome) == 'Poor'):
                try:
                    all_poor_data = np.hstack((all_poor_data,np.reshape(avg_data,(900,1))))
                except Exception:
                    all_poor_data = np.reshape(avg_data,(900,1))

        # Apply z score normalization
        if(zscore):
            for k in range(all_good_data.shape[1]):
                tmp = all_good_data[:,k].squeeze()
                all_good_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
            for k in range(all_poor_data.shape[1]):
                tmp = all_poor_data[:,k].squeeze()
                all_poor_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
        if(feature == 'abs_cres'):
            good = np.nanmean(np.abs(all_good_data[300:300+window,:]),axis=0) - np.nanmean(np.abs(all_good_data[300-window:300,:]),axis=0)
            poor = np.nanmean(np.abs(all_poor_data[300:300+window,:]),axis=0) - np.nanmean(np.abs(all_poor_data[300-window:300,:]),axis=0)
        elif(feature == 'cres'):
            good = np.nanmean(all_good_data[300:300+window,:],axis=0) - np.nanmean(all_good_data[300-window:300,:],axis=0)
            poor = np.nanmean(all_poor_data[300:300+window,:],axis=0) - np.nanmean(all_poor_data[300-window:300,:],axis=0)
        elif(feature == 'cres_only_after_EEC'):
            good = np.nanmean(all_good_data[300:300+window,:],axis=0)
            poor = np.nanmean(all_poor_data[300:300+window,:],axis=0)
        elif(feature == 'normalized_cres_ai_abs'):
            top = np.nanmean(np.abs(all_good_data[300:300+window,:]),axis=0) - np.nanmean(np.abs(all_good_data[300-window:300,:]),axis=0)
            bottom = np.nanmean(np.abs(all_good_data[300:300+window,:]),axis=0) + np.nanmean(np.abs(all_good_data[300-window:300,:]),axis=0)
            good = top/bottom
            top = np.nanmean(np.abs(all_poor_data[300:300+window,:]),axis=0) - np.nanmean(np.abs(all_poor_data[300-window:300,:]),axis=0)
            bottom = np.nanmean(np.abs(all_poor_data[300:300+window,:]),axis=0) + np.nanmean(np.abs(all_poor_data[300-window:300,:]),axis=0)
            poor = top/bottom
        elif(feature == 'normalized_cres_ai_no_abs'):
            top = np.nanmean(all_good_data[300:300+window,:],axis=0) - np.nanmean(all_good_data[300-window:300,:],axis=0)
            bottom = np.nanmean(all_good_data[300:300+window,:],axis=0) + np.nanmean(all_good_data[300-window:300,:],axis=0)
            good = top/bottom
            top = np.nanmean(all_poor_data[300:300+window,:],axis=0) - np.nanmean(all_poor_data[300-window:300,:],axis=0)
            bottom = np.nanmean(all_poor_data[300:300+window,:],axis=0) + np.nanmean(all_poor_data[300-window:300,:],axis=0)
            poor = top/bottom
        else:
            good = all_good_data[62,:]
            poor = all_poor_data[62,:]

        good = good[~np.isnan(good)]
        poor = poor[~np.isnan(poor)]

        # plt.figure(1)
        plt.boxplot([good, poor], labels=['Good','Poor'])
        # for patch in box['whiskers']:
        #     patch.set_c('black')
        s,p = scipy.stats.ttest_ind(good,poor)
        # Permutation test
        all_data = np.hstack((good,poor))
        labels = np.hstack(([1]*len(good),[0]*len(poor)))
        perm_p = []
        for k in range(1000):
            perm = np.random.permutation(labels)
            rand_good = all_data[perm == 1]
            rand_poor = all_data[perm == 0]
            s,rand_p = scipy.stats.ttest_ind(rand_good, rand_poor)
            perm_p.append(rand_p)
        p = ((p >= np.array(perm_p)).sum() + 1.)/(1.0*len(perm_p))
        # return p
        plt.xlabel('Surgical Outcome')
        plt.xticks([1,2],['Good','Poor'])
        plt.ylabel('$\Delta c_{res}$')
        plt.title('mean c_res(1 min after) - mean c_res(1 min before) \nin %s connectivity \n p = %0.4f'%(fconn,p))
        plt.ylim([-0.1,0.1])
        fig = plt.gcf()
        fig.set_size_inches(18.5,10.5)
        plt.savefig(os.path.expanduser('~/gdrive/aim3/fig/box_patient_skip_chop_%s_%s.png'%(skip_chop,fconn)),dpi=100,bbox_inches='tight')
        # plt.hold(True)

        # plt.hold(False)
        plt.close()





############################################################################################################
### GENERATE SEM plot
import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *


for fconn in ['broadband_CC','alphatheta','beta','lowgamma','highgamma']:
    window = 5
    zscore = False
    skip_chop = True
    feature = 'cres'

    all_good = []
    all_poor = []
    for dilate_radius in [-5,0,5]:
        # All cres
        all_cres = gather_results(dilate_radius, fconn)

        good_cres = []
        bad_cres = []

        for patient_id in all_cres.keys():
            if('CHOP' in patient_id and skip_chop):
                continue
            # Load patient outcome
            outcome = get_outcome(data['PATIENTS'][patient_id]['Outcome'])

            avg_data = np.array(())
            for clip_id in all_cres[patient_id].keys():
                cres = all_cres[patient_id][clip_id]
                if(cres.shape[0] == 901):
                    cres = cres[:900]
                try:
                    avg_data = np.hstack((avg_data,np.reshape(cres,(900,1))))
                except Exception:
                    avg_data = np.reshape(cres,(900,1))
            avg_data = np.nanmean(avg_data,axis=1)

            if(outcome == 'Good'):
                good_cres.append(cres)
            else:
                bad_cres.append(cres)
        good_cres = np.array(good_cres)
        bad_cres = np.array(bad_cres)

        # Apply z score normalization
        if(zscore):
            for k in range(good_cres.shape[1]):
                tmp = good_cres[:,k].squeeze()
                good_cres[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
            for k in range(bad_cres.shape[1]):
                tmp = bad_cres[:,k].squeeze()
                bad_cres[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])

        if(feature == 'abs_cres'):
            good = []
            poor = []
            for w in np.arange(0,good_cres.shape[1]+1,window):
                good.append(np.nanmean(np.abs(good_cres[:,w:w+window]),axis=1))
                poor.append(np.nanmean(np.abs(bad_cres[:,w:w+window]),axis=1))
            good = np.array(good)
            poor = np.array(poor)
        elif(feature == 'cres'):
            good = []
            poor = []
            for w in np.arange(0,good_cres.shape[1]+1,window):
                good.append(np.nanmean(good_cres[:,w:w+window],axis=1))
                poor.append(np.nanmean(bad_cres[:,w:w+window],axis=1))
            good = np.array(good)
            poor = np.array(poor)
        else:
            raise

        times = np.arange(-300.0,600,window)*1.0/60

        cres = []
        for k in good:
            cres.append(np.nanmean(k[~np.isnan(k)]))
        cres = cres[:-1]
        all_good.append(cres)

        cres = []
        for k in poor:
            cres.append(np.nanmean(k[~np.isnan(k)]))
        cres = cres[:-1]
        all_poor.append(cres)

    all_good = np.array(all_good)
    all_poor = np.array(all_poor)
    error = scipy.stats.sem(all_good,nan_policy='omit')
    cres=  np.nanmean(all_good,axis=0)
    plt.plot(times,cres,'b-')
    plt.fill_between(times,cres-error,cres+error,facecolor='blue',alpha=0.25)
    plt.hold(True)

    error = scipy.stats.sem(all_poor,nan_policy='omit')
    cres=  np.nanmean(all_poor,axis=0)
    plt.plot(times,cres,'r-')
    plt.fill_between(times,cres-error,cres+error,facecolor='red',alpha=0.25)

    # plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('$c_{res}(t)$')
    plt.xlim([-5.0,10.0])
    plt.title('$c_{res}(t)$ in %s connectivity of all patients across -5%% to +5%% ROZ'%(fconn))
    plt.legend(['Good Outcome','Poor Outcome'])
    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    fig.savefig(os.path.expanduser('~/gdrive/aim3/fig/ROZ_delta_clip_%s.png'%(fconn)),dpi=100,bbox_inches='tight')
    plt.hold(False)































############################################################################################################
import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

comp_dir = data['COMP_DIR']
for patient_id in data['PATIENTS'].keys():
    try:
        tmp = data['PATIENTS'][patient_id]['MANUAL_RESECTED_ELECTRODES']
    except Exception:
        print patient_id
        try:
            fn = glob.glob(os.path.expanduser('%s/%s/aim1/%s_resected_electrodes_0.csv'%(comp_dir,patient_id,patient_id)))[0]
            data['PATIENTS'][patient_id]['MANUAL_RESECTED_ELECTRODES'] = map(lambda x: x.split(',')[1].replace('\n',''), open(fn,'r').readlines())
        except Exception:
            pass
        continue






############################
## Get csv for dc.js
import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

for patient_id in sorted(os.listdir(os.path.expanduser(data['COMP_DIR']))):
    if(patient_id == 'TEST1'):
        continue
        # Initialize output
    resected_electrodes = []

    # Get resection image
    try:
        outcome = data['PATIENTS'][patient_id]['Outcome']
        gender = data['PATIENTS'][patient_id]['Sex']
        age = data['PATIENTS'][patient_id]['AgeSurgery']

        print ','.join([patient_id,gender,age,get_outcome(outcome)])
    except KeyError:
        print patient_id + ' No image'




#####################
# stacked bar chart of demo
demo = pd.read_csv('~/gdrive/tmp/aim3/demographics.csv')


Series1 = pd.Series([demo[demo.study_site=='CHOP'][demo.outcome=='Good'].shape[0],demo[demo.study_site=='HUP'][demo.outcome=='Good'].shape[0],demo[demo.study_site=='Mayo'][demo.outcome=='Good'].shape[0]])

Series2 = pd.Series([demo[demo.study_site=='CHOP'][demo.outcome=='Poor'].shape[0],demo[demo.study_site=='HUP'][demo.outcome=='Poor'].shape[0],demo[demo.study_site=='Mayo'][demo.outcome=='Poor'].shape[0]])


import matplotlib as mpl
import seaborn as sns
tot = Series1+Series2
inst = pd.Series(['CHOP','HUP','Mayo'])

#Set general plot properties
sns.set_style("white")
sns.set_context({"figure.figsize": (24, 10)})

tmp = pd.DataFrame([inst,Series1,Series2,tot]).transpose()
tmp.columns = ['Inst','Good','Poor','Total']

#Plot 1 - background - "total" (top) series

sns.barplot(x = tmp.Inst, y = tmp.Total, palette = [(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), (1.0, 0.4980392156862745, 0.054901960784313725), (0.17254901960784313, 0.6274509803921569, 0.17254901960784313)])

#Plot 2 - overlay - "bottom" series
bottom_plot = sns.barplot(x = tmp.Inst, y = tmp.Good, palette = [(0.6823529411764706, 0.7803921568627451, 0.9098039215686274), (1.0, 0.7333333333333333, 0.47058823529411764), (0.596078431372549, 0.8745098039215686, 0.5411764705882353)])

#Optional code - Make plot look nicer
sns.despine(left=True)
bottom_plot.set_ylabel("Count")
bottom_plot.set_xlabel("Institution")

#Set fonts to consistent 16pt size
for item in ([bottom_plot.xaxis.label, bottom_plot.yaxis.label] +
             bottom_plot.get_xticklabels() + bottom_plot.get_yticklabels()):
    item.set_fontsize(16)

plt.show()





##########################33
##### CREATE NEW_DATA.json
import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
from util_plot import *

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

comp_dir = data['COMP_DIR']
# for patient_id in sorted(data['PATIENTS'].keys()):
for patient_id in ['HUP064','HUP065','HUP068','HUP086','HUP087','HUP088']:
    if('CHOP' in patient_id or 'Study' in patient_id):
        continue
    try:
        tmp = data['PATIENTS'][patient_id]['MANUAL_RESECTED_ELECTRODES']
    except Exception:
        print patient_id
        # Generate list of cartoon map labels
        labels = map(lambda x: x.split(',')[4].replace('\n',''), open(os.path.expanduser(
            data['PATIENTS'][patient_id]['ELECTRODE_LABELS']
            ),'r').readlines())

        fn = '%s-interictal-block-1.mat'%patient_id
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
        labels = correspond_label_names(channels, labels)

        fn = glob.glob(os.path.expanduser('%s/%s/aim1/%s_resected_electrodes_0.csv'%(comp_dir,patient_id,patient_id)))[0]
        manual_resected_electrodes = map(lambda x: x.split(',')[1].replace('\n',''), open(fn,'r').readlines())
        channels_manual_resected_electrodes = map(lambda x: labels[x][1],manual_resected_electrodes)
        manual_ignore_electrodes = data['PATIENTS'][patient_id]['IGNORE_ELECTRODES']
        channels_manual_ignore_electrodes = []
        for ele in manual_ignore_electrodes:
            if ele not in labels.keys():
                continue
            channels_manual_ignore_electrodes.append(labels[ele][1])
        data['PATIENTS'][patient_id]['MANUAL_RESECTED_ELECTRODES'] = channels_manual_resected_electrodes
        data['PATIENTS'][patient_id]['MANUAL_IGNORE_ELECTRODES'] = channels_manual_ignore_electrodes


import json
with open('../data/NEW_DATA.json', 'w') as fp:
    json.dump(data, fp)
