
import time, os, sys, re, glob, random, copy
import numpy as np
import matplotlib.pyplot as plt
import PIL.Image as Image
from skimage import color
from skimage import io

from fun_map import *
import matplotlib as mpl
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
dir_path = '/home/etsai/BNL/Users/SMI/CMurray/2018C3_CMurray_data/saxs/analysis/'
feature_args = {#'filename'  : 'large_G1_15mgml_finegrid2*5.00s', # [*] Specify
                'filename'  : 'medium_as-synth_highC_f*10.00s', #Round 2 Sample1
                'filename'  : 'medium_G1_13mgml_*5.00s*72990', # m*y-7*5
                #'filename'  : 'medium_G2-3G1_20mgml_x*_y*5.00s', 
                'filename'  : 'medium_G2-2G1_highC_m*10.00s*82818',  #81484, 082969
                #'filename'  : dir_path+'large_G2-2G1_2_med*10.00s',
                #'filename'  : '14_As-synthesized_DEG_Grid',  #x-0.350_y0.20 #14_As-synthesized_DEG_Grid',
                'exclude': ['072641', '081511',  
                            '0698', '06996', '06997',
                            '074831', '074833'], 
                'feature_id': 1, # ignore
                'map_type': 'xy',
                'log10'  : 0, 
                'verbose': 1,
                'plot_interp':  [ 'linear', 0.001], #None, 'linear'(recommended), 'cubic', 'nearest', pixel in mm
                'subplot': 1
               } 

# =============================================================================
# Specify features
# For each feature_id, one can specify what specific feature maps to extract.
# Example: 
#   If feature_id == 1, for each specified pixel, we get a feature map
#   If feature_id == 2, feature maps for I(q=q0), I(q=q1), ...
#   If feature_id == 3, feature maps for I(chi=chi0), argmax(I(chi|q0)), ...
#   If feature_id == 4, feature maps based on fitting results, d_spacing, grain_size, ...
# =============================================================================
feature_1_args = {'source_dir' : dir_path+'qr_image/', #thumbnails2/
             'ext' : '.npz',
             'targets' : [ ],  # [190, 43], [*] Choose pixels
             'roi': [1, 'mean'],    # [*] Choose +/- n pixels to include
             }

feature_2_args = {'source_dir' : dir_path+'circular_average/', #'../circular_average/',
             'ext' : '.dat',
             'data_col' : [0, 2],
             'targets' : [0.074], #0.053  # [*] Choose q0 or q0,q1
             'roi': [3, 'mean'],    # [*] Choose the half-width (data points) of the peak q
             'N_peaks_find': 7
             }
                   
feature_3_args = {'source_dir' : dir_path+'linecut_angle092/',
             'ext' : '.dat',
             'data_col' : [0, 1],
             'angle_roi': [6, 'mean'], #[-61,  1], # range [0, 60] or N_fold [6, 'mean']
             'targets': [], #['argmax','var'], #[6.6, 9.5, 12.6,32, 35, 49], #, 'var', 10, 26, 36, 42 , 57, 59, 69], #'max', #[21] # 'max', 'var', or specify angle 
             'normalize': True, # normalize by sum(I)
             'N_peaks_find': 15,
             }

feature_4_args = {'source_dir' : dir_path+'circular_average/',
             'ext' : '.dat',
             'data_col' : [0, 2],
             'fit_range': [0.085, 0.097],  
             'chi2_thr': 0.001, # only consider pixels above this threshold
             'targets': ['b', 'prefactor1', 'x_center1', 'd_spacing_nm', 'grain_size_nm', 'chi2'] #b, prefactor1, x_center1, sigma1, chi2
             }

feature_args.update(feature_1_args=feature_1_args, feature_2_args=feature_2_args, feature_3_args=feature_3_args, feature_4_args=feature_4_args)

    
## Plot one measurement
if 1:
    feature_args['feature_id'] = 1;
    fig = plt.figure(1, figsize=[8,8]); plt.clf()
    cmap = plt.get_cmap('jet');  feature_args.update(cmap=cmap) 
    #feature_args.update(filename='*082070');   
    infiles, match_re = get_filematch(feature_args)
    feature_args.update(log10=1)
    plot_data(infiles[0], **feature_args)

    
# =============================================================================
# Feature maps
# Get maps, plot, apply math, overlay
# =============================================================================
features_map_list = []; 
t0 = time.time()

## Get maps for each feature_ids
feature_ids = [2]
for idx in feature_ids:
    feature_args['feature_id'] = idx; 
    
    ## Find matching files   
    infiles, match_re = get_filematch(feature_args)  
    
    ## Get map   
    features_map = get_map(infiles, match_re, feature_args)
    features_map_list.append(features_map)
    
    ## Plot map
    if False:
        fig = plt.figure(100+feature_args['feature_id'], figsize=[16,8]); plt.clf()  
        cmap = plt.get_cmap('jet');  feature_args.update(cmap=cmap)    
        feature_args.update(log10=0)
        plot_map(features_map, **feature_args)
    
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


## Plot sum
if feature_args['feature_id'] ==2:
    log10 = 1
    line_x = np.asarray(features_map['info_map'][-1][0].x)
    line_y = np.zeros(line_x.shape)
    for idx, line in enumerate(features_map['info_map']):        
        line_y += np.asarray(line[0].y)
    if log10: line_y = np.log10(line_y)
    
    line = DataLine(x=line_x, y=line_y)
    plt.figure(98); plt.clf()
    plot_peaks(line, 0.02)
    plt.title('feature_id {}'.format(feature_args['feature_id'])+'\n'+'{}'.format(feature_args['filename']))
    plt.ylabel('sum over all {} measurements'.format((idx+1)))
        
## Plot all maps
fig = plt.figure(300, figsize=[16,8]); plt.clf()  
features_map_all = extract_maps(features_map_list) 
feature_args.update(log10=0)
plot_map(features_map_all, **feature_args)


## Apply math to selected maps
if False:  
    # Manually:
    # feature_a = features_map_all['features'][0]
    # feature_b = features_map_all['features'][1]
    # feature_c = feature_a + feature_b
    # Or:
    feature_args['math_ab'] = [1, 0, 'multiply']
    feature_c = math_features(features_map_list, **feature_args)
    print('Total of {} maps'.format(count_maps(features_map_list)))
    
    feature_args['math_ab'] = [4, 0, 'divide']
    feature_c = math_features(features_map_list, **feature_args)
    print('Total of {} maps'.format(count_maps(features_map_list)))
    
    features_map_all = extract_maps(features_map_list)
    
    fig = plt.figure(200, figsize=[16,5]); plt.clf()  
    plot_map(features_map_all, **feature_args)

## Plot overlay of three maps (RGB)  
if True:
    feature_args['overlay_rgb'] = [0,1] # starts from 0
    feature_args['normalize_each'] = 1
    feature_args.update(log10=1)
    overlay = plot_overlay(features_map_list, **feature_args)    

            
## Plot histogram
if False:
    x = features_map_all['features'][1]
    x = np.asarray(x)
    plt.figure(99); plt.clf(); 
    plt.hist(x, bins=50)
    plt.grid()
    
## Show FOV
if 0:
    roi_x = [np.min(features_map_list[0]['x_pos']), np.max(features_map_list[0]['x_pos'])]
    roi_y = [np.min(features_map_list[0]['y_pos']), np.max(features_map_list[0]['y_pos'])]
    print(roi_x)
    print(roi_y)






