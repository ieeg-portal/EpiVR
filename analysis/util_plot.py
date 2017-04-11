#!/usr/bin/python
'''
Utility module to plot results from the Virtual Resection project.
'''

from util import *
from Echobase.Statistics.FDA.fda import *

np.random.seed(sum(map(ord, "aesthetics")))

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

warnings.filterwarnings('ignore')

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
                    seizure_type = data['PATIENTS'][patient_id]['EVENTS']['Ictal'][clip_id]['SeizureType']
                    if('CPS' not in seizure_type):
                        continue
                    results[patient_id][clip_id] = np.load('%s/%s'%(comp_dir,fn))['control_centrality_%s'%fconn]
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
