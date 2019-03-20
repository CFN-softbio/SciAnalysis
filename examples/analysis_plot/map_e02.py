#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from fun_map import *

mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.titlesize'] = 12
mpl.rcParams['lines.linewidth'] = 2 
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12


# =============================================================================
# Input files
# 1) Do runXS.py to generate necessary files (eg qr_image, circular_average), make sure mask and beam center are correct
# 2) Exclude bad scans
# =============================================================================
dir_path = '/home/etsai/BNL/Research/KY_platelets/saxs/analysis/'
dir_path = '/home/etsai/BNL/Users/SMI/CMurray/2018C3_CMurray_data/saxs/'
feature_args = {
                'source_dir': dir_path,
                'filename'  : 'medium_as-synth_highC_f*10.00s', 
                'filename'  : 'medium_G1_13mgml_*5.00s', 
                #'filename'  : 'medium_G2-2G1_highC_m*10.00s',  
                #'filename'  : 'medium_G2-3G1_20mgml_x*_y*5.00s', 
                #'filename'  : '14_As-synthesized_DEG_Grid',  #x-0.350_y0.20 #14_As-synthesized_DEG_Grid',
                'filename': '*076288*', 
                'exclude': ['072641', '072729', '079729', '081511',  
                            '0698', '06996', '06997',
                            '074831', '074833'], 
                'direct': 1, ###
                'feature_id': 1, # ignore
                'map_type': 'xy',
                'log10'  : 0, 
                'verbose': 1,
                'plot_interp':  [None, 0.001], #None, 'linear'(recommended), 'cubic', 'nearest', pixel in mm
                'subplot': 1
               } 

feature_1_args = {'source_dir' : dir_path+'analysis/qr_image/', 'ext' : '.npz',
             'protocols': [Protocols.qr_image()], 
             'targets' : [],  # [[190, 43], [123, 123]]
             'roi': [1, 'mean'],    # [*] Choose +/- n pixels to include
             }

feature_2_args = {'source_dir' : dir_path+'analysis/circular_average/', 'ext' : '.dat', 'data_col' : [0, 2],
             'protocols': [Protocols.circular_average()], 
             'targets' : [0.036], #0.053  # [*] Choose q0 or q0,q1
             'roi': [3, 'mean'],    # [*] Choose the half-width (data points) of the peak q
             'N_peaks_find': 5
             }
                   
feature_3_args = {'source_dir' : dir_path+'analysis/linecut_angle092/', 'ext' : '.dat', 'data_col' : [0, 1],
             'protocols': [Protocols.linecut_angle(q0=feature_2_args['targets'], dq=0.002)], 
              'angle_roi': [6, 'mean'], #[-61,  1], # range [0, 60] or N_fold [6, 'mean']
             'targets': ['argmax','var'], #['argmax','var'], #[6.6, 9.5, 12.6,32, 35, 49], #, 'var', 10, 26, 36, 42 , 57, 59, 69], #'max', #[21] # 'max', 'var', or specify angle 
             'normalize': True, # normalize by sum(I)
             'N_peaks_find': 15,
             }

feature_4_args = {'source_dir' : dir_path+'analysis/circular_average/', 'ext' : '.dat', 'data_col' : [0, 2],
             'protocols': [Protocols.circular_average()], 
              'fit_range': [0.085, 0.097],  
             'chi2_thr': 0.001, # only consider pixels above this threshold
             'targets': ['b', 'prefactor1', 'x_center1', 'd_spacing_nm', 'grain_size_nm', 'chi2'] #b, prefactor1, x_center1, sigma1, chi2
             }

feature_args.update(feature_1_args=feature_1_args, feature_2_args=feature_2_args, feature_3_args=feature_3_args, feature_4_args=feature_4_args)


# =============================================================================
# The usual: calibration, mask, and process args
# 
# The typical process.run([infile], protocols, output_dir='./', force=True)
# is in function get_map -> get_feature
# =============================================================================
calibration = Calibration(wavelength_A=0.770088) # 16.1 keV
calibration.set_image_size(981, height=1043) # Pilatus1M
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_distance(5.300)

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus1M_main_gaps-mask.png')
mask.load(dir_path+'analysis/Pilatus1M_bad_pixel.png') 
filename = feature_args['filename']
if ('medium_G2-2G1' in filename) or ('medium_as-synth' in filename):
    calibration.set_distance(5.300)
    calibration.set_beam_position(452.0, 566.0) # medium_G2; medium_as-synth (round2 samples)
    mask.load(dir_path+'analysis/Pilatus1M_current-mask.png') 
else:
    calibration.set_distance(2.300)
    calibration.set_beam_position(493.0, 560.0) # medium_G1_13mgml_; medium_G2-3G1_
    mask.load(dir_path+'analysis/Pilatus1M_current-mask2.png')
    
load_args = { 'calibration' : calibration,
             'mask' : mask,
             }
run_args = { 'verbosity' : 0,
            'plot_save': False,
            'threshold': 65000
            }
process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)
if feature_args['direct']:
    feature_args['process'] = process
        
    
# =============================================================================
# Feature maps
# Get maps, plot, apply math, overlay
# =============================================================================
features_map_list = []; t0 = time.time()

feature_ids = [2]
for idx in feature_ids:
    feature_args['feature_id'] = idx; 
    
    ## Find matching files   
    infiles, match_re = get_filematch(feature_args)  
    
    ## Get map
    features_map = get_map(infiles, match_re, feature_args)
    features_map_list.append(features_map)
    
    
    ## Plot one data 
    if True:
        fig = plt.figure(200+feature_args['feature_id'], figsize=[8,8]); plt.clf()
        cmap = plt.get_cmap('jet');  feature_args.update(cmap=cmap)    
        #feature_args.update(filename='*074852');   infiles, match_re = get_filematch(feature_args)
        feature_args.update(log10=1)
        #feature_args.update(val_stat = [0, 3])
        plot_data(infiles[0], **feature_args)
        
    
    t1 = time.time()-t0
    print('----------------------')
    print('Total of {} maps'.format(count_maps(features_map_list)))
    print('Time = {:.1f} s = {:.1f} min'.format(t1, t1/60))
    print('----------------------')
   

    
## Plot one measurement
if 1:
    feature_args['feature_id'] = 1;
    fig = plt.figure(1, figsize=[8,8]); plt.clf()
    cmap = plt.get_cmap('jet');  feature_args.update(cmap=cmap) 
    #feature_args.update(filename='*082070');   
    infiles, match_re = get_filematch(feature_args)
    feature_args.update(log10=1)
    plot_data(infiles[0], **feature_args)    
    
    

## Plot all maps
fig = plt.figure(300, figsize=[16,8]); plt.clf()  
features_map_all = extract_maps(features_map_list) 
feature_args.update(log10=0)
plot_map(features_map_all, **feature_args)


gc.collect()

