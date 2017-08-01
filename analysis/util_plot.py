'''
Utility module to plot results from the Virtual Resection project.
'''

from util import *
from Echobase.Statistics.FDA.fda import *
from Echobase.Common import errors
from Echobase.Sigproc import reref, prewhiten, filters
# from rsrch_vresect_validation.setup import *

np.random.seed(sum(map(ord, "aesthetics")))

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

warnings.filterwarnings('ignore')

def plot_eeg(fig_fn, patient_id, event_type, event_id, data=data, sep=0.75, lw=0.5):
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


    events = data['PATIENTS'][patient_id]['Events'][event_type]
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

    # Correspond lable names
    labels = correspond_label_names(channels, labels)

    # Load electrodes to ignore
    ignored_node_idx = map(lambda x: labels[x][0],ignored_node_labels)
    for ii,node_id in enumerate(ignored_node_idx):
        print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
    evData = np.delete(evData, ignored_node_idx, axis=1)

        # Parameter set
    param = {}
    param['Notch_60Hz'] = {'wpass': [58.0, 62.0],
                           'wstop': [59.0, 61.0],
                           'gpass': 0.1,
                           'gstop': 60.0}
    param['HPF_5Hz'] = {'wpass': [5.0],
                        'wstop': [4.0],
                        'gpass': 0.1,
                        'gstop': 60.0}
    param['LPF_115Hz'] = {'wpass': [115.0],
                          'wstop': [120.0],
                          'gpass': 0.1,
                          'gstop': 60.0}
    param['LPF_50Hz'] = {'wpass': [50.0],
                          'wstop': [55.0],
                          'gpass': 0.1,
                          'gstop': 60.0}
    param['XCorr'] = {'tau': 0.25}

    # Build pipeline
    data_hat = reref.common_avg_ref(evData)
    data_hat = prewhiten.ar_one(data_hat)
    data_hat = filters.elliptic(data_hat, Fs, **param['Notch_60Hz'])
    data_hat = filters.elliptic(data_hat, Fs, **param['HPF_5Hz'])
    if(Fs > 230):
        data_hat = filters.elliptic(data_hat, Fs, **param['LPF_115Hz'])
    else:
        data_hat = filters.elliptic(data_hat, Fs, **param['LPF_50Hz'])

    num_channels = data_hat.shape[1]

    ev_min = np.min(data_hat[:])
    ev_max = np.max(data_hat[:])
    ev_range = ev_max-ev_min

    ev_iter = 0
    plt.figure(1,dpi=1200)
    for channel in range(num_channels):
        plt.plot(np.linspace(-T/(2*Fs),T/(2*Fs),data_hat.shape[0]),data_hat[:,channel]+sep*ev_iter*ev_range, color='k', linewidth=lw)
        ev_iter += 1
        plt.hold(True)
    plt.grid(True)
    frame1 = plt.gcf()
    frame1.axes[0].get_yaxis().set_visible(False)
    plt.xlim([-T/(2*Fs),T/(2*Fs)])
    plt.ylim([ev_min,ev_max+sep*ev_iter*ev_range])

    plt.xlabel('Time (s)')
    plt.ylabel('EEG (uV)')

    # plt.show()
    plt.savefig(fig_fn,bbox_inches='tight')


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
    outcome = get_outcome(data['PATIENTS'][patient_id]['Outcome'])

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

        # Load all cnull data across all blocks
        cnull_data = np.zeros((num_ictal_events,num_epochs*1000))
        if(len(files['cnull'])  == 1000*num_ictal_events):
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
            for epoch in np.arange(240,360):
                if epoch >= 240 and epoch < 255:
                    time_window = '-1 min. to -45 sec.'
                elif epoch >= 255 and epoch < 270:
                    time_window = '-45 sec. to -30 sec.'
                elif epoch >= 270 and epoch < 285:
                    time_window = '-30 sec. to -15 sec.'
                elif epoch >= 285 and epoch < 300:
                    time_window = '-15 sec. to seizure onset'
                elif epoch >= 300 and epoch < 315:
                    time_window = 'seizure onset to +15 sec.'
                elif epoch >= 315 and epoch < 330:
                    time_window = '+15 sec. to +30 sec.'
                elif epoch >= 330 and epoch < 345:
                    time_window = '+30 sec. to +30 sec.'
                elif epoch >= 345 and epoch < 360:
                    time_window = '+45 sec. to +1 min.'

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
        cres_df.TimeWindow = cres_df.TimeWindow.astype("category", categories=['-1 min. to -45 sec.','-45 sec. to -30 sec.','-30 sec. to -15 sec.','-15 sec. to seizure onset','seizure onset to +15 sec.','+15 sec. to +30 sec.','+30 sec. to +30 sec.','+45 sec. to +1 min.','Null'], ordered=True)

        # Perform K-S 2 sample test on each window to null
        if(len(files['cnull']) == 1000*num_ictal_events):
            f = open('%s/%s_%s_CC_%s_distributions_stats.csv'%(comp_dir, patient_id, unique_id, fconn),'w')
            stats_txt = ''
            for event_id in range(1,num_ictal_events+1):
                for time_window in ['-1 min. to -45 sec.','-45 sec. to -30 sec.','-30 sec. to -15 sec.','-15 sec. to seizure onset','seizure onset to +15 sec.','+15 sec. to +30 sec.','+30 sec. to +30 sec.','+45 sec. to +1 min.']:
                    D,p = scipy.stats.ks_2samp(cres_df[cres_df.TimeWindow == time_window][cres_df.Clip=='Ictal Clip %i'%event_id].Cres,cnull_data[event_id-1,:])
                    m,l,h = mean_confidence_interval(cres_df[cres_df.TimeWindow == time_window][cres_df.Clip == 'Ictal Clip %i'%event_id].Cres)
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
        plt.plot(np.arange(-5*60/epoch_length,10*60/epoch_length),cres_data.T, alpha=0.65)
        plt.plot([0, 0],[-3, 3], alpha=0.1)
        fig.axes[0].text(5, -0.75, 'Seizure onset')
        plt.xlabel('Time (seconds)')
        plt.ylabel('Control Centrality')
        plt.title('%s Control Centrality across all ictal clips using %s\n Null Confidence Interval: [%.3E %.3E]'%(patient_id, fconn.upper(), ci5,ci95))
        plt.ylim([-1,1])
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_time_plot.png'%(comp_dir,patient_id,unique_id, fconn))
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_time_plot.svg'%(comp_dir,patient_id,unique_id, fconn))
        plt.clf()

        # Define the color palette
        color_palette_sz = [
        (0.18431373554117539, 0.47266437376246734, 0.7116493828156415),
        (0.18431373554117539, 0.47266437376246734, 0.7116493828156415),
        (0.18431373554117539, 0.47266437376246734, 0.7116493828156415),
        (0.94071511310689593, 0.60991928098248505, 0.48127645896930321),
        (0.75617071810890646, 0.21038062695194693, 0.22352941947824814),
        (0.76617071810890646, 0.21038062695194693, 0.22352941947824814),
        (0.77617071810890646, 0.21038062695194693, 0.22352941947824814),
        (0.78617071810890646, 0.21038062695194693, 0.22352941947824814),
        (0.125, 0.75, 0.125)]

        ax = sns.factorplot(x="Cres", y="Clip", hue="TimeWindow", row="PatientID", data=cres_df, orient='h', size=2, aspect=3.5, palette=color_palette_sz, kind="violin", cut=0, bw=.2)
        sns.despine(offset=10, trim=True)
        ax.fig.set_alpha = 0.25
        ax.fig.set_size_inches((8,18))
        # ax.fig.axes[0].fill_between([ci5, ci95],[-1, 3],[-1, 3],facecolor='gray',alpha=0.5)
        xlow = min(min(cres_data.flatten())*0.8,min(cres_data.flatten())*1.2)
        xhigh = max(max(cres_data.flatten())*0.8,max(cres_data.flatten())*1.2)

        plt.xlim([xlow,xhigh])

        # Characterize the plots approximately
        if(np.mean(cres_data.flatten()) < 0):
            if(np.mean(cres_data.flatten()) > -0.10):
                subnetwork_characterization = 'Weakly Synchronizing'
            else:
                subnetwork_characterization = 'Strongly Synchronizing'
        else:
            if(np.mean(cres_data.flatten()) > 0.10):
                subnetwork_characterization = 'Weakly Desynchronizing'
            else:
                subnetwork_characterization = 'Strongly Desynchronizing'

        plt.title('Patient %s Dilation Radius %i - %s (%s Outcome %s)'%(patient_id,pipedef['dilate_radius'], fconn.upper(), subnetwork_characterization, outcome))
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_distributions.png'%(comp_dir,patient_id,unique_id, fconn))
        plt.savefig('%s/%s_%s_ControlCentrality_res_%s_distributions.svg'%(comp_dir,patient_id,unique_id, fconn))

def plot_all_cres_heatmap(dilate_radius, outcome_type, fconn = 'highgamma'):
    '''
    Utility function to plot all c_res results given a dilation radius as a heatmap.
    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        outcome_type: str,
            Outcome, either "Good" or "Poor"
    Returns
    -------
        None
            All plot figures are saved as png, svg in the COMP_DIR specified by data json file.
    '''
    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    # Convert all cres_data clips into dataframe
    all_cres_df = []
    for patient_id in all_cres.keys():
        outcome = data['PATIENTS'][patient_id]['Outcome']
        if(get_outcome(outcome) == outcome_type):
            for clip in range(1,len(all_cres[patient_id].keys())+1):
                for epoch in np.arange(200,500):
                    if epoch >= 240 and epoch < 255:
                        time_window = '-1 min. to -45 sec.'
                    elif epoch >= 255 and epoch < 270:
                        time_window = '-45 sec. to -30 sec.'
                    elif epoch >= 270 and epoch < 285:
                        time_window = '-30 sec. to -15 sec.'
                    elif epoch >= 285 and epoch < 300:
                        time_window = '-15 sec. to seizure onset'
                    elif epoch >= 300 and epoch < 315:
                        time_window = 'seizure onset to +15 sec.'
                    elif epoch >= 315 and epoch < 330:
                        time_window = '+15 sec. to +30 sec.'
                    elif epoch >= 330 and epoch < 345:
                        time_window = '+30 sec. to +30 sec.'
                    elif epoch >= 345 and epoch < 360:
                        time_window = '+45 sec. to +1 min.'
                    else:
                        time_window = '-'

                    cres = all_cres[patient_id].values()[clip-1][epoch]
                    all_cres_df.append([patient_id+' Ictal Clip %i'%clip,time_window,epoch,cres])

    columns = ['PatientID_Clip','TimeWindow','Epoch','Cres']
    all_cres_df = pd.DataFrame(np.array(all_cres_df),columns=columns)
    # all_cres_df.Clip = all_cres_df.Clip.astype("category")
    all_cres_df.PatientID_Clip = all_cres_df.PatientID_Clip.astype("category")
    all_cres_df.Cres = all_cres_df.Cres.astype('float64')
    all_cres_df.TimeWindow = all_cres_df.TimeWindow.astype("category", categories=['-1 min. to -45 sec.','-45 sec. to -30 sec.','-30 sec. to -15 sec.','-15 sec. to seizure onset','seizure onset to +15 sec.','+15 sec. to +30 sec.','+30 sec. to +30 sec.','+45 sec. to +1 min.','-'], ordered=True)
    all_cres_rect = all_cres_df.pivot("PatientID_Clip","Epoch","Cres")
    # # #### PLOT HEATMAP
    # sns.heatmap(all_cres_rect, square=True, xticklabels=50)
    # # plt.yticks(rotation=0)
    # plt.show()

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

    print all_data.shape
    for k in range(all_data.shape[1]):
        tmp = all_data[:,k].squeeze()
        all_data[~np.isnan(tmp),k] = scipy.stats.zscore(tmp[~np.isnan(tmp)])
    fig, ax = plt.subplots()
    ax.imshow(all_data.T,extent=[0,1,0,0.1],cmap='jet',clim=[-2,2])
    # ax.set_yticklabels(y_labels)
    ax.grid(False)
    # plt.clim([-2,2])
    plt.xticks([])
    plt.yticks([])
    plt.title('Cres for patients with %s outcome'%outcome_type)
    plt.savefig('%s/../fig/ControlCentrality_res_%s_dilation_%s_Outcome_%s_heat_plot.png'%(comp_dir, str(dilate_radius), fconn, outcome_type),bbox_inches='tight')

def plot_all_cres_time(dilate_radius, outcome_type, fconn = 'highgamma'):
    '''
    Utility function to plot all c_res results given a dilation radius as a time plot with error bars.
    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        outcome_type: str,
            Outcome, either "Good" or "Poor"
    Returns
    -------
        None
            All plot figures are saved as png, svg in the COMP_DIR specified by data json file.
    '''
    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    if(outcome_type == 'Both'):
        # Plot for good outcome
        y_labels = []
        all_data = np.array(())
        for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
            outcome = data['PATIENTS'][patient_id]['Outcome']
            if(get_outcome(outcome) == 'Good'):
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

        plt.plot(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1), alpha=0.75)
        plt.hold(True)

        # Plot for poor outcome
        y_labels = []
        all_data = np.array(())
        for patient_id,clips in sorted(all_cres.items(), key=lambda x:x[0]):
            outcome = data['PATIENTS'][patient_id]['Outcome']
            if(get_outcome(outcome) == 'Poor'):
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

        plt.plot(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1), alpha=0.75)

    else:
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
        plt.plot(np.arange(-5,10,1/60.0),np.nanmean(all_data,axis=1))

    plt.title('Cres for patients with %s outcome'%outcome_type)
    plt.savefig('%s/../fig/ControlCentrality_res_%s_dilation_%s_Outcome_%s_time_plot.png'%(comp_dir, str(dilate_radius), fconn, outcome_type),bbox_inches='tight')
    plt.close()

def plot_all_cres_box(fn, dilate_radius, skip_chop = True, skip_mayo = True, skip_hup = False, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):
    '''
    Utility function to plot all c_res results given a dilation radius as a time plot with error bars.
    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        fconn: str
            Connectivity measure
    Returns
    -------
        None
            All plot figures are saved as png, svg in the COMP_DIR specified by data json file.
    '''
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

    good = good[~np.isnan(good)]
    poor = poor[~np.isnan(poor)]

    plt.boxplot([good, poor], labels=['Good','Poor'])
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
    return p
    # plt.title('mean |c_res(1 min after)| - |c_res(1 min before)| p = %0.4f'%p)
    # plt.show()
    # plt.savefig('%s/../fig/ControlCentrality_res_%s_dilation_%s_box_plot.png'%(comp_dir, str(dilate_radius), fconn),bbox_inches='tight')
    # plt.savefig(fn,bbox_inches='tight')
    # plt.close()

def plot_all_patient_cres_box(fn, dilate_radius, skip_chop = False, skip_mayo = False, skip_hup = False, zscore = False, feature = 'special', fconn = 'highgamma', window = 30):
    '''
    Utility function to plot all c_res results averaged within each patient given a dilation radius as a time plot with error bars.
    Parameters
    ----------
        dilate_radius: str,
            Dilation radius

        fconn: str
            Connectivity measure
    Returns
    -------
        None
            All plot figures are saved as png, svg in the COMP_DIR specified by data json file.
    '''
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

    plt.boxplot([good, poor], labels=['Good','Poor'])
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
    plt.ylabel('Network feature')
    plt.title('mean zscore[c_res(1 min after)] - mean zscore[c_res(1 min before)] p = %0.4f'%p)
    plt.show()
    # plt.savefig('%s/../fig/ControlCentrality_res_%s_dilation_%s_box_plot.png'%(comp_dir, str(dilate_radius), fconn),bbox_inches='tight')
    # plt.savefig(fn,bbox_inches='tight')
    # plt.close()

def gather_results(dilate_radius, fconn = 'highgamma'):
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
                    # if('CPS' not in seizure_type):
                    #     continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['control_centrality_%s'%fconn]
        except:
            continue

    return results

def gather_sync_results(dilate_radius, fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all results.

    Parameters
    ----------
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
                    # if('CPS' not in seizure_type):
                    #     continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['base_sync_%s'%fconn]
        except:
            continue

    return results


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

        # Open all adj
        try:
            for fn in os.listdir(comp_dir):
                if('multiband' in fn):
                    if(patient_id not in results.keys()):
                        results[patient_id] = {}
                    print fn
                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['EVENTS']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type):
                        continue
                    results[patient_id][fn.split('.')[2]] = np.load('%s/%s'%(comp_dir,fn))['all_adj_%s'%fconn]
        except:
            continue

    return results

def plot_clip_ticker_time(patient_id, clip_id, fn, dilate_radius, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):
    '''
    Function that plots ticker plot of cres
    '''

    # All cres
    all_cres = gather_results(dilate_radius, fconn)
    cres = all_cres[patient_id][clip_id]
    fig, ax = plt.subplots()
    times = np.arange(-300.0,600,window)*1.0/60
    box = ax.boxplot(np.reshape(cres,(cres.shape[0]/window,window)).T);
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')
    plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    fig = plt.gcf()
    fig.set_size_inches(15,5)
    plt.show()
    # plt.savefig(fn,bbox_inches='tight')

def plot_patient_ticker_time(patient_id, fn, dilate_radius, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):

    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    fig, ax = plt.subplots()
    times = np.arange(-300.0,600,window)*1.0/60

    iter_id = 0
    for clip_id in all_cres[patient_id].keys():
        cres = all_cres[patient_id][clip_id]
        if('abs' in feature and 'no_abs' not in feature):
            cres = np.abs(cres)
        box = ax.boxplot(np.reshape(cres,(cres.shape[0]/window,window)).T+iter_id*0.5);
        for patch in box['whiskers']:
            patch.set_c('black')

        for patch, time in zip(box['boxes'], times):
            if(time < 0):
                patch.set_c('black')
                patch.set_fillstyle('full')
            else:
                patch.set_c('red')
                patch.set_fillstyle('full')
        plt.hold(True)
        iter_id += 1
    plt.xticks(map(lambda x: int(x+1), range(times.shape[0]))[::2],times[::2])
    plt.yticks([])
    plt.xlabel('Time (min)')
    plt.ylabel('c_res(t)')
    plt.title('Patient %s (Poor Outcome) - Cres(t) in broadband cross-correlation connectivity \nover 30 second windows'%(patient_id))
    plt.show()


def plot_group_ticker_time(fn, dilate_radius, skip_chop = True, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):

    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    good_cres = []
    bad_cres = []

    for patient_id in all_cres.keys():
        if('CHOP' in patient_id and skip_chop):
            continue
        # Load patient outcome
        outcome = get_outcome(data['PATIENTS'][patient_id]['Outcome'])
        for clip_id in all_cres[patient_id].keys():
            cres = all_cres[patient_id][clip_id]
            if(cres.shape[0] == 901):
                cres = cres[:900]
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
        good = np.nanmean(np.abs(good_cres),axis=0)
        poor = np.nanmean(np.abs(bad_cres),axis=0)
    elif(feature == 'cres'):
        good = np.nanmean(good_cres,axis=0)
        poor = np.nanmean(bad_cres,axis=0)
    elif(feature == 'cres_only_after_EEC'):
        good = good_cres
        poor = bad_cres
    elif(feature == 'normalized_cres_ai_abs'):
        good = np.nanmean(np.abs(good_cres),axis=0)
        poor = np.nanmean(np.abs(bad_cres),axis=0)
    elif(feature == 'normalized_cres_ai_no_abs'):
        good = good_cres
        poor = bad_cres

    good = good[~np.isnan(good)]
    poor = poor[~np.isnan(poor)]


    fig, ax = plt.subplots()
    times = np.arange(-300.0,600,window)*1.0/60

    cres = poor
    box = ax.boxplot(np.reshape(cres,(cres.shape[0]/window,window)).T);
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')
    plt.hold(True)

    cres = good + 0
    box = ax.boxplot(np.reshape(cres,(cres.shape[0]/window,window)).T);
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')


    plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('c_res(t)')
    plt.title('')
    plt.show()
    # plt.savefig(fn,bbox_inches='tight')


def plot_group_event_based_ticker_time(fn, dilate_radius, skip_chop = True, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):
    '''
    Event-based boxplots
    '''

    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    good_cres = []
    bad_cres = []

    for patient_id in all_cres.keys():
        if('CHOP' in patient_id and skip_chop):
            continue
        # Load patient outcome
        outcome = get_outcome(data['PATIENTS'][patient_id]['Outcome'])
        for clip_id in all_cres[patient_id].keys():
            cres = all_cres[patient_id][clip_id]
            if(cres.shape[0] == 901):
                cres = cres[:900]
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

    fig, ax = plt.subplots()
    times = np.arange(-300.0,600,window)*1.0/60

    cres = []
    for k in poor:
        cres.append(k[~np.isnan(k)])
    box = ax.boxplot(cres)
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')
    plt.hold(True)

    cres = []
    for k in good:
        cres.append(k[~np.isnan(k)])
    box = ax.boxplot(cres)
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')


    plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('c_res(t)')
    plt.title('')
    plt.show()
    # plt.savefig(fn,bbox_inches='tight')

def plot_group_patient_based_ticker_time(fn, dilate_radius, skip_chop = True, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):
    '''
    IntraPatient-averaged boxplots
    '''

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
                avg_data = np.hstack((avg_data,np.reshape(clip_data,(900,1))))
            except Exception:
                avg_data = np.reshape(clip_data,(900,1))
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

    fig, ax = plt.subplots()
    times = np.arange(-300.0,600,window)*1.0/60

    cres = []
    for k in poor:
        cres.append(k[~np.isnan(k)])
    box = ax.boxplot(cres)
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')
    plt.hold(True)

    cres = []
    for k in good:
        cres.append(k[~np.isnan(k)])
    box = ax.boxplot(cres)
    for patch in box['whiskers']:
        patch.set_c('black')

    for patch, time in zip(box['boxes'], times):
        if(time < 0):
            patch.set_c('black')
            patch.set_fillstyle('full')
        else:
            patch.set_c('red')
            patch.set_fillstyle('full')


    plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('c_res(t)')
    plt.title('')
    plt.show()
    # plt.savefig(fn,bbox_inches='tight')

def plot_group_event_based_line(fn, dilate_radius, skip_chop = True, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):
    '''
    Event-based line plot
    '''

    # All cres
    all_cres = gather_results(dilate_radius, fconn)

    good_cres = []
    bad_cres = []

    for patient_id in all_cres.keys():
        if('CHOP' in patient_id and skip_chop):
            continue
        # Load patient outcome
        outcome = get_outcome(data['PATIENTS'][patient_id]['Outcome'])
        for clip_id in all_cres[patient_id].keys():
            cres = all_cres[patient_id][clip_id]
            if(cres.shape[0] == 901):
                cres = cres[:900]
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
    for k in poor:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(cres)
    plt.plot(times,cres,'r-')
    plt.fill_between(times,cres-error,cres+error)
    plt.hold(True)

    cres = []
    for k in good:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(cres)
    plt.plot(times,cres,'b-')
    plt.fill_between(times,cres-error,cres+error)

    # plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('c_res(t)')
    plt.title('')
    plt.show()
    # plt.savefig(fn,bbox_inches='tight')


def plot_group_patient_based_line(fn, dilate_radius, skip_chop = True, zscore = False, feature = 'cres', fconn = 'highgamma', window = 30):
    '''
    IntraPatient-averaged line plots
    '''

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
                avg_data = np.hstack((avg_data,np.reshape(clip_data,(900,1))))
            except Exception:
                avg_data = np.reshape(clip_data,(900,1))
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
    for k in poor:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(cres)
    plt.plot(times,cres,'r-')
    plt.fill_between(times,cres-error,cres+error)
    plt.hold(True)

    cres = []
    for k in good:
        cres.append(np.nanmean(k[~np.isnan(k)]))
    cres = cres[:-1]
    error = scipy.stats.sem(cres)
    plt.plot(times,cres,'b-')
    plt.fill_between(times,cres-error,cres+error)

    # plt.xticks(map(lambda x: x+1, range(times.shape[0]))[::2],times[::2])
    # plt.yticks([])
    plt.xlabel('Time (min.)')
    plt.ylabel('c_res(t)')
    plt.title('')
    plt.show()
    # plt.savefig(fn,bbox_inches='tight')

def run_all_FDA_test(dilate_radius, fconn):
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

    g = all_good_data.shape[1]
    b = all_poor_data.shape[1]
    print curve_test(np.hstack((all_good_data,all_poor_data)),np.arange(0,g),np.arange(g,g+b))



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
    g = all_good_data.shape[1]
    b = all_poor_data.shape[1]
    print curve_test(np.hstack((all_good_data,all_poor_data)),np.arange(0,g),np.arange(g,g+b))


def gather_nodal_results(fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all node virtual resection results.

    Parameters
    ----------
        fconn: str,
            Connectivity metric

    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with nodal cres in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue

        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')

        # Open all adj
        try:
            for fn in os.listdir(comp_dir):
                if('noderes' in fn):
                    if(patient_id not in results.keys()):
                        results[patient_id] = {}
                    print fn
                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    # if('CPS' not in seizure_type):
                    #     continue
                    results[patient_id][fn.split('.')[2]] = np.load('%s/%s'%(comp_dir,fn))['control_centrality_%s'%fconn]
        except:
            continue

    return results


def gather_null_nodal_results(fconn = 'highgamma'):
    '''
    Utility function to output a dictionary of all node virtual resection results.

    Parameters
    ----------
        fconn: str,
            Connectivity metric

    Returns
    -------
        results: dict,
            Results that contain as key patient id, and has as value another dictionary with nodal cres in each clip.
    '''
    # All cres
    results = {}

    for patient_id in os.listdir(os.path.expanduser(data['COMP_DIR'])):
        if(patient_id == 'TEST1'):
            continue

        comp_dir = os.path.join(os.path.expanduser(data['COMP_DIR']),patient_id,'aim3')

        # Open all adj
        try:
            for fn in os.listdir(comp_dir):
                if('nodenull' in fn):
                    if(patient_id not in results.keys()):
                        results[patient_id] = {}
                    print fn
                    clip_id = fn.split('.')[2]
                    seizure_type = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type):
                        continue
                    if clip_id not in results[patient_id].keys():
                        results[patient_id][clip_id] = []
                    results[patient_id][clip_id].append(np.load('%s/%s'%(comp_dir,fn))['control_centrality_%s'%fconn])
        except:
            continue

    for patient_id, events in results.items():
        for clip_id, clip_data in events.items():
            results[patient_id][clip_id] = np.array(clip_data)
    return results

def plot_CC_variability(fig_fn,nodal_control_centrality, null_nodal_control_centrality=None):
    '''
    Utility function to plot control centrality variability.

    Parameters
    ----------
        fig_fn: str,
            Full path filename to save figure

        nodal_control_centrality: nadarray,
            Control centrality of shape N x T where N is the number of nodes, and T are the total time epochs.

        null_nodal_control_centrality: ndarray,
            Control centrality null models of shape P x N x T where P is the number of null permutations, N is the number of nodes, and T are the total time epochs.

    Returns
    -------
        None
    '''

    nodal_sem = scipy.stats.sem(nodal_control_centrality, axis=1)
    nodal_mean = np.mean(nodal_control_centrality, axis=1)

    plt.figure(dpi=1200)
    plt.errorbar(np.arange(len(nodal_mean))+1,nodal_mean[np.argsort(nodal_mean)],yerr=nodal_sem[np.argsort(nodal_mean)],fmt='--')

    if null_nodal_control_centrality is not None:
        null_nodal_low = np.percentile(null_nodal_control_centrality.flatten(),2.5)
        null_nodal_high = np.percentile(null_nodal_control_centrality.flatten(),97.5)
        null_nodal_mean = np.mean(null_nodal_control_centrality.flatten())

        ax = plt.gca()
        ax.fill_between(np.arange(len(nodal_mean))+1,null_nodal_low,null_nodal_high, facecolor='gray',alpha=0.5)

    # plt.show()
    plt.savefig(fig_fn,bbox_inches='tight')


def plot_all_CC_variability():
    '''
    Utility function to plot all control centrality variability.

    Parameters
    ----------
        None

    Returns
    -------
        None
    '''

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    for fconn in ['highgamma','broadband_CC']:
        all_nodal_control_centrality = gather_nodal_results(fconn)
        all_null_nodal_control_centrality = gather_null_nodal_results(fconn
                )
        for patient_id, events in all_nodal_control_centrality.items():
            for clip_id, nodal_control_centrality in all_nodal_control_centrality[patient_id].items():
                try:
                    # Get null
                    null_nodal_control_centrality = all_null_nodal_control_centrality[patient_id][clip_id]
                except KeyError:
                    null_nodal_control_centrality = None

                # Get pre-seizure and seizure epoch
                pre_nodal_control_centrality = nodal_control_centrality[:,:nodal_control_centrality.shape[1]/2]
                post_nodal_control_centrality = nodal_control_centrality[:,nodal_control_centrality.shape[1]/2:]

                if null_nodal_control_centrality is not None:
                    pre_null_nodal_control_centrality = null_nodal_control_centrality[:,:,:nodal_control_centrality.shape[1]/2]
                    post_null_nodal_control_centrality = null_nodal_control_centrality[:,:,nodal_control_centrality.shape[1]/2:]
                else:
                    pre_null_nodal_control_centrality = None
                    post_null_nodal_control_centrality = None

                print pre_nodal_control_centrality.shape, post_nodal_control_centrality.shape
                if null_nodal_control_centrality is not None:
                    print pre_null_nodal_control_centrality.shape, post_null_nodal_control_centrality.shape

                fig_fn = '%s/../fig/%s.Ictal.%s.%s.preseizure.variability.png'%(comp_dir,patient_id,clip_id,fconn)
                plot_CC_variability(fig_fn, pre_nodal_control_centrality, pre_null_nodal_control_centrality)

                fig_fn = '%s/../fig/%s.Ictal.%s.%s.seizure.variability.png'%(comp_dir,patient_id,clip_id,fconn)
                plot_CC_variability(fig_fn, post_nodal_control_centrality, post_null_nodal_control_centrality)


def write_hub_coord_csv(ele_csv_fn, patient_id, clip_idx, clip_idx_label, epoch='post', fconn='highgamma'):
    '''
    Utility function to write csv file given electrode locations all node-based control centrality variability.

    Parameters
    ----------
        ele_csv_fn: str,
            Full path filename to load soz electrode labels and locations

        patient_id: str,
            Patient ID

        clip_idx: nd.array or list,
            List of clip IDs that belong to the same seizure type

    Returns
    -------
        None
    '''

    # Load event
    for clip_id in clip_idx:
        event = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]
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


        # Load channels from
        fn = os.path.join(data_dir, patient_id, 'eeg', event['FILE'])
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
        # Create final list of ordered cartoon electrode labels
        labels = np.array(sorted(labels.keys(), key=lambda x: labels[x][0]))
        labels = np.delete(labels, ignored_node_idx, axis=0)

        # All nodal results
        all_nodal_results = gather_nodal_results(fconn=fconn)
        nodal_control_centrality = all_nodal_results[patient_id][clip_id]
        if(epoch == 'pre'):
            nodal_control_centrality = nodal_control_centrality[:,:nodal_control_centrality.shape[1]/2]
        if(epoch == 'post'):
            nodal_control_centrality = nodal_control_centrality[:,nodal_control_centrality.shape[1]/2:]
        mean_nodal_control_centrality = np.mean(nodal_control_centrality, axis=1)

        # Compute dictionary of electrode label to nodal control centrality mean
        electrode_control_centrality = {}
        for kk, label in enumerate(labels):
            if label not in electrode_control_centrality.keys():
                electrode_control_centrality[label] = []
            electrode_control_centrality[label].append(mean_nodal_control_centrality[kk])

    medianmean_nodal_control_centrality = {}
    for label in electrode_control_centrality.keys():
        medianmean_nodal_control_centrality[label] = np.median(electrode_control_centrality[label])

    # # Fit exponential
    # nodal_control_centrality = medianmean_nodal_control_centrality.values()
    # min_nodal_control_centrality = min(nodal_control_centrality)
    # for label, val in medianmean_nodal_control_centrality.items():
    #     corr_val = np.log(val + 1 - min_nodal_control_centrality)
    #     medianmean_nodal_control_centrality[label] = corr_val

    # Compute mean and normalize nodal_control_centrality for each node
    nodal_control_centrality = medianmean_nodal_control_centrality.values()
    mean_nodal_control_centrality = np.mean(nodal_control_centrality)
    min_nodal_control_centrality = min(nodal_control_centrality)
    max_nodal_control_centrality = max(nodal_control_centrality)

    # Normalize with mean, min, max
    for label, val in medianmean_nodal_control_centrality.items():
        corr_val = (val- min_nodal_control_centrality)/(max_nodal_control_centrality-min_nodal_control_centrality)
        medianmean_nodal_control_centrality[label] = corr_val

    # Load electrode labels csv
    lines = open(ele_csv_fn,'r').readlines()
    out_txt = ''

    # Write normalized nodal_control_centrality
    out_txt = ''
    for line in lines:
        label = line.split(',')[3].replace('\n','').replace('\r','')
        try:
            out_txt += '%s,%s,%0.8f\n'%(','.join(line.split(',')[:3]),label,medianmean_nodal_control_centrality[label])
        except KeyError:
            continue

    open(os.path.join(comp_dir, patient_id, 'aim3', ele_csv_fn.split('/')[-1].replace('_soz_coord.csv','_Ictal_%s_%s_%s_hub_coord.csv'%(epoch,clip_idx_label,fconn))),'w').write(out_txt)


def plot_figure1C(patient_id, clip_idx, clip_idx_label, dilate_radius, fconn='highgamma'):

    comp_dir = os.path.expanduser(data['COMP_DIR'])

    all_cres = gather_results(dilate_radius, fconn)
    all_base_sync = gather_sync_results(dilate_radius, fconn)

    min_seizure_len = 1E100

    for clip_id in clip_idx:
        cres = all_cres[patient_id][clip_id]
        if cres.shape[0] < min_seizure_len:
            min_seizure_len = cres.shape[0]

    for clip_id in clip_idx:
        cres = all_cres[patient_id][clip_id]
        base_sync = all_base_sync[patient_id][clip_id]

        cres = np.interp(np.linspace(-1.0,1.0,min_seizure_len),np.linspace(-1.0,1.0,cres.shape[0]),cres.flatten())
        base_sync = np.interp(np.linspace(-1.0,1.0,min_seizure_len),np.linspace(-1.0,1.0,base_sync.shape[0]),base_sync.flatten())

        try:
            avg_cres_data = np.hstack((avg_cres_data,np.reshape(cres,(cres.shape[0],1))))
        except Exception:
            avg_cres_data = np.reshape(cres,(cres.shape[0],1))
        try:
            avg_base_sync_data = np.hstack((avg_base_sync_data,np.reshape(base_sync,(base_sync.shape[0],1))))
        except Exception:
            avg_base_sync_data = np.reshape(base_sync,(base_sync.shape[0],1))

    avg_cres_error = scipy.stats.sem(avg_cres_data,axis=1,nan_policy='omit')
    avg_base_sync_error = scipy.stats.sem(avg_base_sync_data,axis
        =1,nan_policy='omit')

    avg_cres_data = np.nanmedian(avg_cres_data,axis=1)
    avg_base_sync_data = np.nanmedian(avg_base_sync_data,axis=1)

    plt.figure(dpi=1200)
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_cres_data,'g-',alpha=0.5)
    ax1.plot(np.linspace(-1.0,1.0,avg_base_sync_data.shape[0]),avg_base_sync_data,'b-',alpha=0.5)


    ax1.fill_between(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_cres_data-avg_cres_error,avg_cres_data+avg_cres_error,facecolor='green',alpha=0.25)
    ax1.fill_between(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_base_sync_data-avg_base_sync_error,avg_base_sync_data+avg_base_sync_error,facecolor='blue',alpha=0.25)

    ax1.set_xlabel('Normalized Time')
    ax1.set_ylabel('$cc_{res}(t)$', color='g')
    ax2.set_ylabel('$s(t)$', color='b')
    ax1.set_ylim([-0.6,0.8])
    ax2.set_ylim([0.0,1.0])
    ax1.grid(False)
    ax2.grid(False)
    # plt.show()
    plt.savefig('%s/../fig/Figure1C.%s.%s.%s.png'%(comp_dir,patient_id,clip_idx_label,fconn),bbox_inches='tight')


def plot_figure5(patient_id, clip_idx, clip_idx_label):
    '''
    '''

    comp_dir = os.path.expanduser(data['COMP_DIR'])
    # Create figure
    outcome = get_outcome(patient_id)
    fig,axs = plt.subplots(1,5,sharey=True)
    fig.set_size_inches((12,4))
    fig.suptitle('Patient %s - Subtype %s - %s Outcome'%(patient_id,clip_idx_label,outcome))

    # Figure font options
    font1 = {'family':'raleway',
            'color': 'lightblue',
            'weight':'bold',
            'size':8,
            }
    font2 = {'family':'raleway',
            'color': 'darkred',
            'weight':'bold',
            'size':8,
            }
    max_y = 0.02
    min_y = -0.02

    for ax, fconn in zip(axs,['alphatheta','beta','lowgamma','highgamma','broadband_CC']):
        if(fconn == 'alphatheta'):
            title = 'Alpha/Theta:\n 5-15 Hz'
        if(fconn == 'beta'):
            title = 'Beta:\n 15-25 Hz'
        if(fconn == 'lowgamma'):
            title = 'Low Gamma:\n 30-40 Hz'
        if(fconn == 'highgamma'):
            title = 'High Gamma:\n 95-105 Hz'
        if(fconn == 'broadband_CC'):
            title = 'Broadband\n Cross-Correlation'

        # All nodal results
        all_nodal_results = gather_nodal_results(fconn=fconn)

        # plt.figure(dpi=600)

        # Set up colors
        colors = ['k','r','g','b','c','m','y','w']
        color_iter = 0

        # Load event
        for clip_id in clip_idx:
            event = data['PATIENTS'][patient_id]['Events']['Ictal'][clip_id]
            nodal_control_centrality = all_nodal_results[patient_id][clip_id]

            resected_node_idx, channels = get_resected_node_dx(patient_id)
            non_resected_node_idx = []
            for k in range(nodal_control_centrality.shape[0]):
                if k not in resected_node_idx:
                    non_resected_node_idx.append(k)

            nodal_control_centrality = nodal_control_centrality[:,:nodal_control_centrality.shape[1]/2]
            mean_nodal_control_centrality = np.mean(nodal_control_centrality, axis=1)

            ax.scatter(np.zeros((len(resected_node_idx),)),mean_nodal_control_centrality[resected_node_idx],color=colors[color_iter],alpha=0.1)
            ax.hold(True)
            ax.scatter(np.ones((len(non_resected_node_idx),)),mean_nodal_control_centrality[non_resected_node_idx],color=colors[color_iter],alpha=0.1)

            nodal_control_centrality = nodal_control_centrality[:,nodal_control_centrality.shape[1]/2:]
            mean_nodal_control_centrality = np.mean(nodal_control_centrality, axis=1)

            ax.scatter(2*np.ones((len(resected_node_idx),)),mean_nodal_control_centrality[resected_node_idx],color=colors[color_iter],alpha=0.1)
            ax.hold(True)
            ax.scatter(3*np.ones((len(non_resected_node_idx),)),mean_nodal_control_centrality[non_resected_node_idx],color=colors[color_iter],alpha=0.1)

            ax.set_ylim([min_y,max_y])
            ax.set_xticks([0,1,2,3])
            ax.set_xticklabels(['Resected','Non-Resected','Resected','Non-Resected'],fontdict={'size':8,'weight':'bold'},rotation='vertical')
            color_iter += 1
        ax.text(0.05,max_y-(max_y-min_y)*0.04,'Pre-Seizure',fontdict=font1)
        ax.text(2.05,max_y-(max_y-min_y)*0.04,'Seizure',fontdict=font2)
        ax.text(0.0, min_y+(max_y-min_y)*0.1,title,fontdict={'family':'raleway','size':12,'color':'black'})
    plt.tight_layout()
    fig.savefig('%s/../fig/Figure3.%s.%s.png'%(comp_dir,patient_id,clip_idx_label))


def hack():
    # import sys
    # import glob
    # import json
    # import time

    # from util import *
    # from util_connectivity import *
    # from util_virtual_resection import *
    # from util_plot import *

    patient_id = 'HUP074'
    epoch_length = 1
    dilate_radius = 0
    event_id = '2'

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


    events = data['PATIENTS'][patient_id]['Events']['Ictal']
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

    # Correspond lable names
    labels = correspond_label_names(channels, labels)

    # Load electrodes to ignore
    ignored_node_idx = map(lambda x: labels[x][0],ignored_node_labels)
    for ii,node_id in enumerate(ignored_node_idx):
        print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
    evData = np.delete(evData, ignored_node_idx, axis=1)



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
        print 'Virtually resecting node label: %s because label %s is in the resection zone'%(channels[node_id],resected_node_labels[ii])



    # Generate list of SOZ electrodes and write to CSV file
    soz_node_labels = events[event_id]["SEIZURE_ONSET_ELECTRODES"]

    # Map the resected electrodes to channels
    clean_soz_node_labels = []
    for soz_node_label in soz_node_labels:
        if soz_node_label in ignored_node_labels:
            continue
        else:
            clean_soz_node_labels.append(soz_node_label)
    soz_node_idx = map(lambda x: labels_dict[x][0], clean_soz_node_labels)
    for ii,node_id in enumerate(soz_node_idx):
        print 'Virtually resecting SOZ node label: %s because label %s is in the SOZ zone'%(channels[node_id],soz_node_labels[ii])


    eec = evData.shape[0]/2 # in samples
    ueo = Fs * (events[event_id]['SeizureUEO']-events[event_id]['SeizureEEC']) + eec
    end = Fs * (events[event_id]['SeizureEnd']-events[event_id]['SeizureEEC']) + eec

    # amt = 80.0
    # evData = evData[eec-(end-eec)*amt:eec+(end-eec)*amt:100,:]
    evData = evData[::100,:]

    T = np.float(evData.shape[0])*100
    sep = 0.5
    lw = 2

    # Parameter set
    param = {}
    param['Notch_60Hz'] = {'wpass': [58.0, 62.0],
                           'wstop': [59.0, 61.0],
                           'gpass': 0.1,
                           'gstop': 60.0}
    param['HPF_5Hz'] = {'wpass': [5.0],
                        'wstop': [4.0],
                        'gpass': 0.1,
                        'gstop': 60.0}
    param['LPF_115Hz'] = {'wpass': [115.0],
                          'wstop': [120.0],
                          'gpass': 0.1,
                          'gstop': 60.0}
    param['LPF_50Hz'] = {'wpass': [50.0],
                          'wstop': [55.0],
                          'gpass': 0.1,
                          'gstop': 60.0}
    param['XCorr'] = {'tau': 0.25}

    # Build pipeline
    data_hat = reref.common_avg_ref(evData)
    data_hat = prewhiten.ar_one(data_hat)
    data_hat = filters.elliptic(data_hat, Fs, **param['Notch_60Hz'])
    data_hat = filters.elliptic(data_hat, Fs, **param['HPF_5Hz'])
    if(Fs > 230):
        data_hat = filters.elliptic(data_hat, Fs, **param['LPF_115Hz'])
    else:
        data_hat = filters.elliptic(data_hat, Fs, **param['LPF_50Hz'])

    num_channels = data_hat.shape[1]

    ev_min = np.min(data_hat[:])
    ev_max = np.max(data_hat[:])
    ev_range = ev_max-ev_min

    ev_iter = 0
    count_iter = 0
    plt.figure(1,dpi=1200)
    for channel in range(num_channels):
        if channel not in soz_node_idx and channel not in resected_node_idx:
            count_iter += 1
            if np.mod(count_iter,2) == 0:
                plt.plot(np.linspace(-T/(2*Fs),T/(2*Fs),data_hat.shape[0]),data_hat[:,channel]+sep*ev_iter*ev_range, color='k', linewidth=lw)
                ev_iter += 1
                plt.hold(True)
        else:
            if channel in resected_node_idx:
                if channel in soz_node_idx:
                    plt.plot(np.linspace(-T/(2*Fs),0,data_hat[:data_hat.shape[0]/2,channel].shape[0]),data_hat[:data_hat.shape[0]/2,channel]+sep*ev_iter*ev_range, color='r', linewidth=lw)
                    plt.plot(np.linspace(0,T/(2*Fs),data_hat[data_hat.shape[0]/2:,channel].shape[0]),data_hat[data_hat.shape[0]/2:,channel]+sep*ev_iter*ev_range, color='m', linewidth=lw)
                    ev_iter += 1
                    plt.hold(True)
                else:
                    plt.plot(np.linspace(-T/(2*Fs),T/(2*Fs),data_hat.shape[0]),data_hat[:,channel]+sep*ev_iter*ev_range, color='r', linewidth=lw)
                    ev_iter += 1
                    plt.hold(True)
            else:
                if channel in soz_node_idx:
                    plt.plot(np.linspace(-T/(2*Fs),0,data_hat[:data_hat.shape[0]/2,channel].shape[0]),data_hat[:data_hat.shape[0]/2,channel]+sep*ev_iter*ev_range, color='k', linewidth=lw)
                    plt.plot(np.linspace(0,T/(2*Fs),data_hat[data_hat.shape[0]/2:,channel].shape[0]),data_hat[data_hat.shape[0]/2:,channel]+sep*ev_iter*ev_range, color='b', linewidth=lw)
                    ev_iter += 1
                    plt.hold(True)
                else:
                    blah
    plt.grid(True)
    frame1 = plt.gcf()
    frame1.axes[0].get_yaxis().set_visible(False)
    plt.xlim([-T/(2*Fs),T/(2*Fs)])
    plt.ylim([ev_min,ev_max+sep*ev_iter*ev_range])
    plt.plot([0,0],[ev_min,ev_max+sep*ev_iter*ev_range],color='#444444')
    plt.plot([np.float(ueo-eec)/Fs+1,np.float(ueo-eec)/Fs+1],[ev_min,ev_max+sep*ev_iter*ev_range],color='#666666')
    plt.plot([np.float(end-eec)/Fs+1,np.float(end-eec)/Fs+1],[ev_min,ev_max+sep*ev_iter*ev_range],color='#aaaaaa')
    plt.xlabel('Time (s)')
    plt.ylabel('EEG (uV)')
    plt.savefig('../../fig/Figure4_1.svg',bbox_inches='tight')


def write_all_nodal_mean_csv_individual_seizure():
    '''
    This writes the nodal mean csv for a given patient id and seizure id.
    '''
    for patient_id in ['HUP064','HUP065','HUP068','HUP070','HUP073','HUP074','HUP075','HUP078','HUP080','HUP082','HUP083','HUP086','HUP087','HUP088','HUP094','HUP105','HUP106','HUP107','HUP111A','HUP111B','HUP116','Study012','Study016','Study017','Study019','Study020','Study022','Study028','Study029']:
        # for each band
        for fconn in ['alphatheta','beta','lowgamma','highgamma','broadband_CC']:
            res = gather_nodal_results(fconn)
            # for each nodal clip  - gather_nodal_results
            for event_id in res[patient_id].keys():
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
                fn = os.path.join(data_dir, patient_id, 'eeg', data['PATIENTS'][patient_id]['Events']['Ictal'][event_id]['FILE'])
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


                channel_names = map(lambda x: x[0], sorted(labels_dict.items(), key=lambda x: x[1][0]))


                # Create PREICTAL patient_id.ictal.event_id.nodres.frequency.preictal.mean.csv
                noderes = res[patient_id][event_id]
                df = pd.DataFrame(zip(channel_names,scipy.stats.zscore(np.mean(noderes[:,:noderes.shape[1]/2.0],axis=1))))
                df.to_csv('%s/%s/aim3/%s.Ictal.%s.noderes.%s.preictal.mean.csv'%(comp_dir,patient_id,patient_id,event_id,fconn),header=False,index=False)
                df = pd.DataFrame(zip(channel_names,scipy.stats.zscore(np.mean(noderes[:,noderes.shape[1]/2.0:],axis=1))))
                df.to_csv('%s/%s/aim3/%s.Ictal.%s.noderes.%s.ictal.mean.csv'%(comp_dir,patient_id,patient_id,event_id,fconn),header=False,index=False)

def plot_network_measures_individual_seizure():
    '''
    Duplication of plot_figure1C but modified to run for all patients and individual seizures
    '''

    dilate_radius = 0
    for patient_id in ['HUP064','HUP065','HUP068','HUP070','HUP073','HUP074','HUP075','HUP078','HUP080','HUP082','HUP083','HUP086','HUP087','HUP088','HUP094','HUP105','HUP106','HUP107','HUP111A','HUP111B','HUP116','Study012','Study016','Study017','Study019','Study020','Study022','Study028','Study029']:
        # for each band
        for fconn in ['alphatheta','beta','lowgamma','highgamma','broadband_CC']:
            all_cres = gather_results(dilate_radius, fconn)
            all_base_sync = gather_sync_results(dilate_radius, fconn)

            # for each nodal clip  - gather_nodal_results
            for event_id in all_cres[patient_id].keys():
                clip_idx = [event_id]
                comp_dir = os.path.expanduser(data['COMP_DIR'])
                min_seizure_len = 1E100
                for clip_id in clip_idx:
                    cres = all_cres[patient_id][clip_id]
                    if cres.shape[0] < min_seizure_len:
                        min_seizure_len = cres.shape[0]

                for clip_id in clip_idx:
                    cres = all_cres[patient_id][clip_id]
                    base_sync = all_base_sync[patient_id][clip_id]

                    cres = np.interp(np.linspace(-1.0,1.0,min_seizure_len),np.linspace(-1.0,1.0,cres.shape[0]),cres.flatten())
                    base_sync = np.interp(np.linspace(-1.0,1.0,min_seizure_len),np.linspace(-1.0,1.0,base_sync.shape[0]),base_sync.flatten())

                    try:
                        avg_cres_data = np.hstack((avg_cres_data,np.reshape(cres,(cres.shape[0],1))))
                    except Exception:
                        avg_cres_data = np.reshape(cres,(cres.shape[0],1))
                    try:
                        avg_base_sync_data = np.hstack((avg_base_sync_data,np.reshape(base_sync,(base_sync.shape[0],1))))
                    except Exception:
                        avg_base_sync_data = np.reshape(base_sync,(base_sync.shape[0],1))

                avg_cres_error = scipy.stats.sem(avg_cres_data,axis=1,nan_policy='omit')
                avg_base_sync_error = scipy.stats.sem(avg_base_sync_data,axis
                    =1,nan_policy='omit')

                avg_cres_data = np.nanmedian(avg_cres_data,axis=1)
                avg_base_sync_data = np.nanmedian(avg_base_sync_data,axis=1)

                plt.figure(dpi=1200)
                fig, ax1 = plt.subplots()

                ax2 = ax1.twinx()
                ax1.plot(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_cres_data,'g-',alpha=0.5)
                ax1.plot(np.linspace(-1.0,1.0,avg_base_sync_data.shape[0]),avg_base_sync_data,'b-',alpha=0.5)


                ax1.fill_between(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_cres_data-avg_cres_error,avg_cres_data+avg_cres_error,facecolor='green',alpha=0.25)
                ax1.fill_between(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_base_sync_data-avg_base_sync_error,avg_base_sync_data+avg_base_sync_error,facecolor='blue',alpha=0.25)

                ax1.set_xlabel('Normalized Time')
                ax1.set_ylabel('$cc_{res}(t)$', color='g')
                ax2.set_ylabel('$s(t)$', color='b')
                ax1.set_ylim([-0.6,0.8])
                ax2.set_ylim([0.0,1.0])
                ax1.grid(False)
                ax2.grid(False)
                # plt.show()
                plt.savefig('%s/../fig/pt/%s.Ictal.%s.cres_and_sync.%s.png'%(comp_dir,patient_id,event_id,fconn),bbox_inches='tight')

def plot_average_network_measures_all_seizures():
    '''
    Averages seizures across PATIENTS
    '''
    dilate_radius = 0
    comp_dir = data['COMP_DIR']
    for patient_id in ['HUP064','HUP065','HUP068','HUP070','HUP073','HUP074','HUP075','HUP078','HUP080','HUP082','HUP083','HUP086','HUP087','HUP088','HUP094','HUP105','HUP106','HUP107','HUP111A','HUP111B','HUP116','Study012','Study016','Study017','Study019','Study020','Study022','Study028','Study029']:
        # for each band
        for fconn in ['alphatheta','beta','lowgamma','highgamma','broadband_CC']:
            all_cres = gather_results(dilate_radius, fconn)
            all_base_sync = gather_sync_results(dilate_radius, fconn)

            # Get minimum length
            min_event_len  = 1E100
            for event_id in all_cres[patient_id].keys():
                event_len = all_cres[patient_id][event_id].shape[0]
                if event_len < min_event_len:
                    min_event_len = event_len

            # Norm all seizures to min_event_len
            width = min_event_len
            dictionary_values_norm = {}
            dictionary_values_norm[patient_id] = {}
            for event_id, clip in all_cres[patient_id].items():
                xp = np.linspace(-1.0,1.0,clip.shape[0])
                dictionary_values_norm[patient_id][event_id] = np.interp(np.linspace(-1.0,1,width),xp,clip)
            all_cres = dictionary_values_norm
            dictionary_values_norm = {}
            dictionary_values_norm[patient_id] = {}
            for event_id, clip in all_base_sync[patient_id].items():
                xp = np.linspace(-1.0,1.0,clip.shape[0])
                dictionary_values_norm[patient_id][event_id] = np.interp(np.linspace(-1.0,1,width),xp,clip)
            all_base_sync = dictionary_values_norm

            clip_idx = all_cres[patient_id].keys()

            for clip_id in clip_idx:
                    cres = all_cres[patient_id][clip_id]
                    base_sync = all_base_sync[patient_id][clip_id]

                    cres = np.interp(np.linspace(-1.0,1.0,min_event_len),np.linspace(-1.0,1.0,cres.shape[0]),cres.flatten())
                    base_sync = np.interp(np.linspace(-1.0,1.0,min_event_len),np.linspace(-1.0,1.0,base_sync.shape[0]),base_sync.flatten())

                    try:
                        avg_cres_data = np.hstack((avg_cres_data,np.reshape(cres,(cres.shape[0],1))))
                    except Exception:
                        avg_cres_data = np.reshape(cres,(cres.shape[0],1))
                    try:
                        avg_base_sync_data = np.hstack((avg_base_sync_data,np.reshape(base_sync,(base_sync.shape[0],1))))
                    except Exception:
                        avg_base_sync_data = np.reshape(base_sync,(base_sync.shape[0],1))

            avg_cres_error = scipy.stats.sem(avg_cres_data,axis=1,nan_policy='omit')
            avg_base_sync_error = scipy.stats.sem(avg_base_sync_data,axis
                =1,nan_policy='omit')

            avg_cres_data = np.nanmedian(avg_cres_data,axis=1)
            avg_base_sync_data = np.nanmedian(avg_base_sync_data,axis=1)

            plt.figure(dpi=1200)
            fig, ax1 = plt.subplots()

            ax2 = ax1.twinx()
            ax1.plot(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_cres_data,'g-',alpha=0.5)
            ax1.plot(np.linspace(-1.0,1.0,avg_base_sync_data.shape[0]),avg_base_sync_data,'b-',alpha=0.5)


            ax1.fill_between(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_cres_data-avg_cres_error,avg_cres_data+avg_cres_error,facecolor='green',alpha=0.25)
            ax1.fill_between(np.linspace(-1.0,1.0,avg_cres_data.shape[0]),avg_base_sync_data-avg_base_sync_error,avg_base_sync_data+avg_base_sync_error,facecolor='blue',alpha=0.25)

            ax1.set_xlabel('Normalized Time')
            ax1.set_ylabel('$cc_{res}(t)$', color='g')
            ax2.set_ylabel('$s(t)$', color='b')
            ax1.set_ylim([-0.6,0.8])
            ax2.set_ylim([0.0,1.0])
            ax1.grid(False)
            ax2.grid(False)
            # plt.show()
            plt.savefig('%s/../fig/pt/%s.Ictal.average.cres_and_sync.%s.png'%(comp_dir,patient_id,fconn),bbox_inches='tight')

def write_all_nodal_csv_individual_seizure(width=120):
    '''
    This writes the nodal csv for a given patient id and seizure id, stretched to width after z scoring.
    '''
    for patient_id in ['HUP106','HUP107','HUP111A','HUP111B','Study012','Study016','Study017','Study019','Study020','Study022','Study028','Study029']:
        # for each band
        for fconn in ['alphatheta','beta','lowgamma','highgamma','broadband_CC']:
            res = gather_nodal_results(fconn)
            # for each nodal clip  - gather_nodal_results
            for event_id in res[patient_id].keys():
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
                fn = os.path.join(data_dir, patient_id, 'eeg', data['PATIENTS'][patient_id]['Events']['Ictal'][event_id]['FILE'])
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


                channel_names = map(lambda x: x[0], sorted(labels_dict.items(), key=lambda x: x[1][0]))

                # Create PREICTAL patient_id.ictal.event_id.noderes.frequency.preictal.mean.csv
                noderes = scipy.stats.zscore(res[patient_id][event_id],axis=0)
                xp = np.linspace(-1.0,1.0,noderes.shape[1])
                f = scipy.interpolate.interp1d(xp,noderes)
                noderes_norm = f(np.linspace(-1.0,1.0,width))

                noderes_out = np.concatenate((np.array(channel_names).reshape([noderes_norm.shape[0],1]),noderes_norm),axis=1)

                df = pd.DataFrame(noderes_out)
                df.to_csv('%s/%s/aim3/%s.Ictal.%s.noderes.%s.Movie.csv'%(comp_dir,patient_id,patient_id,event_id,fconn),header=False,index=False)

def make_html(patient_id):
    comp_dir = DATA['COMPD']
    # Compute all event idx
    event_idx = map(str,sorted(map(int,data['PATIENTS'][patient_id]['Events']['Ictal'].keys())))

    # Generate html subcode
    subcode_html = ''
    for event_id in event_idx:
        subcode_html += '''
                <div class="tab">
                    <button class="tablinks" onclick="openType(event, '%s')" id="defaultOpen">%s</button>
                </div>
'''%(event_id,event_id)
    for event_id in event_idx:
        subcode_html += '''
                <div id="%s" class="tabcontent">
                    <h2 style="font-family:Raleway;"> Clip %s </h2>
                    <h2 style="font-family:Raleway;"> Network Measures </h2>
                    <h3 style="font-family:Raleway;"> High Gamma </h3>
                    <img src="pt/%s.Ictal.%s.cres_and_sync.highgamma.png">
                    <h3 style="font-family:Raleway;"> Broadband </p>
                    <img src="pt/%s.Ictal.%s.cres_and_sync.broadband_CC.png">
                    <h3 style="font-family:Raleway;"> Alpha/Theta </p>
                    <img src="pt/%s.Ictal.%s.cres_and_sync.alphatheta.png">
                    <h3 style="font-family:Raleway;"> Beta </p>
                    <img src="pt/%s.Ictal.%s.cres_and_sync.beta.png">
                    <h3 style="font-family:Raleway;"> Low Gamma </p>
                    <img src="pt/%s.Ictal.%s.cres_and_sync.lowgamma.png">
                    <h2 style="font-family:Raleway;"> Animation of node-level </h2>
                    <h3 style="font-family:Raleway;"> Broadband </h3>
                    <img src="pt/%s.Ictal.%s.noderes.broadband.mp4">
                </div>
'''%(event_id,event_id,patient_id,event_id,patient_id,event_id,patient_id,event_id,patient_id,event_id,patient_id,event_id,patient_id,event_id)

    # Print final html code
    final_html = '''
        <!DOCTYPE html>
        <html lang="en">

        <head>

            <meta charset="utf-8">
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <meta name="description" content="">
            <meta name="author" content="">

            <title>Virtual Resection Project </title>

            <!-- Bootstrap Core CSS -->
            <link href="static/css/bootstrap-social.css" rel="stylesheet">
            <link href="static/css/tab.css" rel="stylesheet">
            <link href="static/css/bootstrap.min.css" rel="stylesheet">

            <!-- Custom CSS -->
            <link href="static/css/business-frontpage.css" rel="stylesheet">
            <style>
                img {
                    filter: gray; /* IE6-9 */
                    filter: grayscale(1); /* Microsoft Edge and Firefox 35+ */
                    -webkit-filter: grayscale(1);  Google Chrome, Safari 6+ & Opera 15+
                }

                /* Disable grayscale on hover */
                img:hover {
                    filter: none;
                    -webkit-filter: grayscale(0);
                }
                /*body{font-family:Garamond}*/
            </style>
            <link href='https://fonts.googleapis.com/css?family=Raleway:300,500|Bitter:400,400italic,700,Lato:300' rel='stylesheet' type='text/css'>
            <link href="https://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css" rel="stylesheet">

            <!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
            <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
            <!--[if lt IE 9]>
                <script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
                <script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
            <![endif]-->

        </head>

        <body>
            <!-- jQuery -->
            <script src="static/js/jquery-1.11.2.min.js"></script>

            <!-- Bootstrap Core JavaScript -->
            <script src="static/bootstrap/dist/js/bootstrap.js"></script>

            <!-- Metis Menu Plugin JavaScript -->
            <script src="static/js/metisMenu.js"></script>

            <!-- SB Admin Custom Theme JavaScript -->
            <script src="static/sb-admin/dist/js/sb-admin-2.js"></script>

            <!-- Chosen JavaScript -->
            <script type="text/javascript" src="static/chosen/chosen.jquery.min.js"></script>
            <script type="text/javascript" src="static/chosen/chosen.proto.min.js"></script>

            <!-- Twitter Typeahead -->
            <script src="static/js/typeahead.bundle.js"></script>
            <script src="static/js/handlebars-v3.0.0.js"></script>

            <!-- Custom Scripts -->
            <script src="static/js/dynamic-number-coloring.js"></script>

            <link href="static/css/dc.css" rel="stylesheet" type="text/css">
            <script type="text/javascript" src="static/js/d3.v3.js"></script>
            <script type="text/javascript" src="static/js/crossfilter.js"></script>
            <script type="text/javascript" src="static/js/dc.js"></script>
            <script type="text/javascript" src="static/js/colorbrewer.js"></script>

            <!-- Navigation -->

            <!-- Page Content -->
            <div class="container">

                <hr>
                <h1 style="font-family:Raleway;font-size:300%%;"> Virtual Resection - %s </h1>

                <hr style="border: none;height: 3px;color: #333;background-color: #333;">
                <h2 style="font-family:Raleway;"> Overview </h2>
                <img src="pt/%s_Figure1A.png">
                <h3>Outcome: </h3>
                <h4>Lesion Status: </h4>
                <h4>Gender: </h4>
                <h4>Age at Onset: </h4>
                <h4>Age at Surgery: </h4>
                <h4>Seizure Onset: </h4>
                <img src="pt/%s_Figure1A.png">

%s


                <hr style="border: none;height: 3px;color: #333;background-color: #333;">
                <h2 style="font-family:Raleway;"> All Seizures </h2>
                <img src="pt/%s.average.cres_and_sync.broadband_CC.png">
                <img src="pt/%s.average.cres_and_sync.highgamma.png">
                <img src="pt/%s.average.cres_and_sync.lowgamma.png">
                <img src="pt/%s.average.cres_and_sync.beta.png">
                <img src="pt/%s.average.cres_and_sync.alphatheta.png">

            </div>

            <!-- /.container -->
            <script type="text/javascript">
            function openType(evt, evtName) {
                // Declare all variables
                var i, tabcontent, tablinks;

                // Get all elements with class="tabcontent" and hide them
                tabcontent = document.getElementsByClassName("tabcontent");
                for (i = 0; i < tabcontent.length; i++) {
                    tabcontent[i].style.display = "none";
                }

                // Get all elements with class="tablinks" and remove the class "active"
                tablinks = document.getElementsByClassName("tablinks");
                for (i = 0; i < tablinks.length; i++) {
                    tablinks[i].className = tablinks[i].className.replace(" active", "");
                }

                // Show the current tab, and add an "active" class to the button that opened the tab
                document.getElementById(evtName).style.display = "block";
                evt.currentTarget.className += " active";
            }
            </script>
            <script type="text/javascript">
            // Get the element with id="defaultOpen" and click on it
            document.getElementById("defaultOpen").click();
            </script>

            <!-- jQuery -->
            <script src="static/js/jquery.js"></script>

            <!-- Bootstrap Core JavaScript -->
            <script src="static/js/bootstrap.min.js"></script>

        </body>

        </html>

        '''%(patient_id,patient_id,patient_id,subcode_html,patient_id,patient_id,patient_id,patient_id,patient_id)
        # pass
    open('%s/../fig/demo/%s.html'%(comp_dir,patient_id),'w').write(final_html)
    return




