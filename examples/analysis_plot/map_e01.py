
import time
import os, sys
import re
import glob
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
import random
import PIL.Image as Image
from skimage import color
from skimage import io

from fun_map import *

########## Input 
dir_path = '/home/etsai/BNL/Research/KY_platelets/saxs/analysis/'
dir_path = '/home/etsai/BNL/Users/SMI/CMurray/2018C3_CMurray_data/saxs/analysis/'
feature_args = {#'filename'  : 'large_G1_15mgml_finegrid2*5.00s', # [*] Specify
                'filename'  : 'medium_G1_13mgml_f*5.00s', # m*y-7*5
                #'filename'  : 'medium_G2-3G1_20mgml_*x-2*5.00s', 
                #'filename'  : 'medium_as-synth_highC_fine*10.00s', Round 2 Sample1
                #'filename'  : 'medium_G2-2G1_highC_med*10.00s', 
                #'filename'  : dir_path+'large_G2-2G1_2_med*10.00s',
                #'filename'  : '14_As-synthesized_DEG_Grid',  #x-0.350_y0.20 #14_As-synthesized_DEG_Grid',
                'feature_id': 1,
                'map_type': 'xy',
                'log10'  : 0,
                'verbose': 1,
                'plot_interp':  ['linear', 0.001], #'none', 'linear'(recommended), 'cubic', 'nearest', pixel in mm
               } 

feature_1_args = {'source_dir' : dir_path+'thumbnails2/', #thumbnails2/
             'ext' : '.jpg',
             #'pixels' : [[89, 122], [102, 201], [152, 73]],  # [*] Choose pixels
             'targets' : [[228, 61], [181, 39]],  # [190, 43], [*] Choose pixels
             'roi': [1, 'mean'],    # [*] Choose +/- n pixels to include
             }

feature_2_args = {'source_dir' : dir_path+'circular_average/', #'../circular_average/',
             'ext' : '.dat',
             'data_col' : [0, 2],
             'targets' : [0.038], #, 0.059], #0.053  # [*] Choose q0 or q0,q1
             'roi': [1, 'mean'],    # [*] Choose the half-width (data points) of the peak q
             }
                   
feature_3_args = {'source_dir' : dir_path+'linecut_angle059/',
             'ext' : '.dat',
             'data_col' : [0, 1],
             'targets': ['max'], #'max', #[21] # 'max', 'var', or specify angle 
             'angle_roi': [5,  65], # range to consider for max or var 
             }

feature_4_args = {'source_dir' : dir_path+'results/',
             'ext' : '.xml',
             'targets': 'fit_peaks_grain_size'
             }

feature_args.update(feature_1_args=feature_1_args, feature_2_args=feature_2_args, feature_3_args=feature_3_args, feature_4_args=feature_4_args)

########## Feature map
features_map_list = []
for idx in [4]:
    feature_args['feature_id'] = idx; 
    
    ## Find matching files   
    infiles, match_re = get_filematch(feature_args)  
    
    ## Get map
    #scans, x_pos, y_pos, feature = get_map(infiles, match_re, feature_args)
    features_map = get_map(infiles, match_re, feature_args)
    features_map_list.append(features_map)
    N_maps = len(features_map['info'][1])
    
    ## Plot map
    fig = plt.figure(100+feature_args['feature_id'], figsize=[16,4]); plt.clf()  
    plot_map(features_map, **feature_args)
    
    ## Plot one data 
    ax2 = plt.subplot2grid((1, N_maps+1), (0, 0), colspan=1); ax2.cla()
    cmap = plt.get_cmap('magma');  feature_args.update(cmap=cmap)    
    #feature_args.update(filename='*74852') # Sample 1 70526
    #infiles, match_re = get_filematch(feature_args)
    img = plot_data(infiles[0], **feature_args)


overlay = plot_overlay(features_map_list, **feature_args) 











