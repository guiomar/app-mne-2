# Copyright (c) 2020 brainlife.io
#
# This file is a MNE python-based brainlife.io App
#
# Author: Guiomar Niso
# Indiana University

# Required libraries
# pip install mne-bids coloredlogs tqdm pandas scikit-learn json_tricks fire

# set up environment
#import mne-study-template
import os
import json
from shutil import copyfile

# Current path
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

# Path to mne-study-template 
#mnest_path = '/Users/guiomar/Documents/GitHub/mne-bids-pipeline'
mnest_path = '/mne-bids-pipeline'

# Populate mne_config.py file with brainlife config.json
with open(__location__+'/config.json') as config_json:
    config = json.load(config_json)


bids_root = str(config['output']) 
deriv_root = 'out_dir'

'''
#study_name = 'ds000246'
bids_root = str(config['output']) # '/Users/guiomar/Projects/ds000246'
deriv_root = 'out_dir'
subjects = ['0001']

l_freq = .3
h_freq = 100.
decim = 10 #4
reject = dict(mag=4e-12, eog=250e-6)
conditions = ['standard', 'deviant', 'button']


contrasts = [('deviant', 'standard')]
decode = True
daysback = -365 * 110
on_error = 'debug'

'''


'''
#study_name = 'ds000248'
subjects = ['01']
rename_events = {'Smiley': 'Emoji','Button': 'Switch'}
conditions = ['Auditory', 'Visual', 'Auditory/Left', 'Auditory/Right']
ch_types = ['meg']
#mf_reference_run = '1' ##bl2bids just 1 not 01
find_flat_channels_meg = True
find_noisy_channels_meg = True
use_maxwell_filter = True
process_er = True
reject = dict(mag=4e-12, eog=250e-6)

contrasts = [('Visual', 'Auditory'),('Auditory/Right', 'Auditory/Left')]
bem_mri_images = 'FLASH'
recreate_bem = True
noise_cov = 'emptyroom'


'''

# Create new MNE config .py file

fname = 'mne_config.py'

with open(fname, 'w') as f: 

    # == GENERAL SETTINGS ==

    f.write("bids_root = '{}'".format(bids_root)+'\n')
    f.write("deriv_root = '{}'".format(deriv_root)+'\n')
    #For freesurfer
    '''
    if config['subjects_dir']:      f.write('subjects_dir = {}'.format(config['subjects_dir'])+'\n')

    if config['study_name']:        f.write('study_name = {}'.format(config['study_name'])+'\n')
    if config['interactive']:       f.write('interactive = {}'.format(config['interactive'])+'\n')
    if config['crop']:              f.write('crop = {}'.format(config['crop'])+'\n')
 
    if config['sessions']:          f.write('sessions = {}'.format(config['sessions'])+'\n')
    if config['task']:              f.write('task = {}'.format(config['task'])+'\n')
    if config['runs']:              f.write('runs = {}'.format(config['runs'])+'\n')
    if config['acq']:               f.write('acq = {}'.format(config['acq'])+'\n')
    if config['proc']:              f.write('proc = {}'.format(config['proc'])+'\n')
    if config['rec']:               f.write('rec = {}'.format(config['rec'])+'\n')
    if config['space']:             f.write('space = {}'.format(config['space'])+'\n')
    if config['subjects']:          f.write('subjects = {}'.format(config['subjects'])+'\n')
    if config['exclude_subjects']:  f.write('exclude_subjects = {}'.format(config['exclude_subjects'])+'\n')
    '''
    if config['process_er']:        f.write('process_er = {}'.format(config['process_er'])+'\n')
    if config['ch_types']:          f.write("ch_types = {}".format(config['ch_types'])+'\n')
    if config['data_type']:         f.write("data_type = {}".format(config['data_type'])+'\n')
    if config['eog_channels']:      f.write('eog_channels = {}'.format(config['eog_channels'])+'\n')
    if config['eeg_bipolar_channels']:  f.write('eeg_bipolar_channels = {}'.format(config['eeg_bipolar_channels'])+'\n')
    if config['eeg_reference']:     f.write("eeg_reference = '{}'".format(config['eeg_reference'])+'\n')
    if config['eeg_template_montage']:  f.write('eeg_template_montage = {}'.format(config['eeg_template_montage'])+'\n')
    if config['drop_channels']:     f.write('drop_channels = {}'.format(config['drop_channels'])+'\n')
    if config['analyze_channels']:  f.write('analyze_channels = {}'.format(config['analyze_channels'])+'\n')
 
       
    # == MAXFLTER (for fif) ==
    # Bad channels
    f.write('find_flat_channels_meg = {}'.format(config['find_flat_channels_meg'])+'\n')
    f.write('find_noisy_channels_meg = {}'.format(config['find_noisy_channels_meg'])+'\n')
  
    f.write('use_maxwell_filter = {}'.format(config['use_maxwell_filter'])+'\n')
    if config['mf_st_duration']:    f.write('mf_st_duration = {}'.format(config['mf_st_duration'])+'\n')
    if config['mf_head_origin']:    f.write('mf_head_origin = {}'.format(config['mf_st_dmf_head_originuration'])+'\n')
    if config['mf_reference_run']:  f.write('mf_reference_run = {}'.format(config['mf_reference_run'])+'\n')
    if config['mf_cal_fname']:      f.write('mf_cal_fname = {}'.format(config['mf_cal_fname'])+'\n')
    if config['mf_ctc_fname']:      f.write('mf_ctc_fname = {}'.format(config['mf_ctc_fname'])+'\n')
       
    # == FILTER & RESAMPLING ==
    # Filter
    if config['l_freq']:    f.write('l_freq = {}'.format(config['l_freq'])+'\n')
    if config['h_freq']:    f.write('h_freq = {}'.format(config['h_freq'])+'\n')

    # Resampling
    if config['resample_sfreq']:    f.write('resample_sfreq = {}'.format(config['resample_sfreq'])+'\n')
    if config['decim']:             f.write('decim = {}'.format(config['decim'])+'\n')   

    # == EPOCHING ==
    if config['rename_events']:     f.write("rename_events = {}".format(config['rename_events'])+'\n')
    #?? if config['on_rename_missing_events']:  f.write("on_rename_missing_events = '{}'".format(config['on_rename_missing_events'])+'\n')
    if config['event_repeated']:    f.write("event_repeated = '{}'".format(config['event_repeated'])+'\n')
    if config['conditions']:        f.write("conditions = {}".format(config['conditions'])+'\n')
    if config['epochs_tmin']:       f.write("epochs_tmin = {}".format(config['epochs_tmin'])+'\n')
    if config['epochs_tmax']:       f.write("epochs_tmax = {}".format(config['epochs_tmax'])+'\n')
    if config['baseline']:          f.write("baseline = {}".format(config['baseline'])+'\n')
  
    if config['epochs_metadata_tmin']:       f.write("epochs_metadata_tmin = {}".format(config['epochs_metadata_tmin'])+'\n')
    if config['epochs_metadata_tmax']:       f.write("epochs_metadata_tmax = {}".format(config['epochs_metadata_tmax'])+'\n')
    if config['epochs_metadata_keep_first']: f.write("epochs_metadata_keep_first = {}".format(config['epochs_metadata_keep_first'])+'\n')
    if config['epochs_metadata_keep_last']:  f.write("epochs_metadata_keep_last = {}".format(config['epochs_metadata_keep_last'])+'\n')


    # == ARTIFACT REMOVAL ==

    # Stimulation Artifact    
    f.write('fix_stim_artifact = {}'.format(config['fix_stim_artifact'])+'\n')
    if config['stim_artifact_tmin']:  f.write('stim_artifact_tmin = {}'.format(config['stim_artifact_tmin'])+'\n')
    if config['stim_artifact_tmax']:  f.write('stim_artifact_tmax = {}'.format(config['stim_artifact_tmax'])+'\n')
 
    # SSP & ICA
    if config['spatial_filter']:    f.write("spatial_filter = '{}'".format(config['spatial_filter'])+'\n')
    # add: ica_reject
    if config['ica_algorithm']:     f.write("ica_algorithm = '{}'".format(config['ica_algorithm'])+'\n')
    if config['ica_l_freq']:        f.write("ica_l_freq = {}".format(config['ica_l_freq'])+'\n')
    if config['ica_max_iterations']: f.write("ica_max_iterations = {}".format(config['ica_max_iterations'])+'\n')
    if config['ica_n_components']:  f.write("ica_n_components = {}".format(config['ica_n_components'])+'\n')
    if config['ica_decim']:         f.write("ica_decim = {}".format(config['ica_decim'])+'\n')
    if config['ica_ctps_ecg_threshold']: f.write("ica_ctps_ecg_threshold = {}".format(config['ica_ctps_ecg_threshold'])+'\n')
    if config['ica_eog_threshold']: f.write("ica_eog_threshold = {}".format(config['ica_eog_threshold'])+'\n')
 
    # Amplitude-Based
    if config['reject']:            f.write("reject = {}".format(config['reject'])+'\n') 
    if config['reject_tmin']:       f.write("reject_tmin = '{}'".format(config['reject_tmin'])+'\n')
    if config['reject_tmax']:       f.write("reject_tmax = '{}'".format(config['reject_tmax'])+'\n')

    f.close() 

'''
    # SENSOR

    # == STATS ==
    if config['contrasts']:         f.write("contrasts = {}".format(config['contrasts'])+'\n')
    if config['decode']:            f.write("decode = {}".format(config['decode'])+'\n')
    if config['decoding_metric']:   f.write("decoding_metric = {}".format(config['decoding_metric'])+'\n')
    if config['decoding_n_splits']: f.write("decoding_n_splits = {}".format(config['decoding_n_splits'])+'\n')
    if config['n_boot']:            f.write("n_boot = {}".format(config['n_boot'])+'\n')

    # == TF ==
    if config['time_frequency_conditions']:  f.write("time_frequency_conditions = {}".format(config['time_frequency_conditions'])+'\n')
 

    # SOURCES

    # General Settings
    if config['run_source_estimation']: f.write("run_source_estimation = {}".format(config['run_source_estimation'])+'\n')

    # BEM surface
    if config['bem_mri_images']:        f.write("bem_mri_images = {}".format(config['bem_mri_images'])+'\n')
    if config['recreate_bem']:          f.write("recreate_bem = {}".format(config['recreate_bem'])+'\n')

    # Source space & forward solution
    if config['mri_t1_path_generator']: f.write("mri_t1_path_generator = {}".format(config['mri_t1_path_generator'])+'\n')
    if config['spacing']:               f.write("spacing = {}".format(config['spacing'])+'\n')
    if config['mindist']:               f.write("mindist = {}".format(config['mindist'])+'\n')

    # Inverse solution
    if config['inverse_method']:    f.write("inverse_method = {}".format(config['inverse_method'])+'\n')
    if config['noise_cov']:         f.write("noise_cov = {}".format(config['noise_cov'])+'\n')
'''

# Run mne-study-template python script
os.system( mnest_path + '/run.py --config=' + __location__+'/mne_config.py \
    --steps=preprocessing,report/make_reports.py')


# Find the reports and make a copy in out_html folder
for dirpath, dirnames, filenames in os.walk(__location__+"/out_dir"):
    for filename in [f for f in filenames if f.endswith(".html")]:
        if not "sub-average" in filename:
            print(filename)
            copyfile(os.path.join(__location__,"out_dir", dirpath,filename), os.path.join(__location__,"html_report",filename))