
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
feature_args = {#'filename'  : 'large_G1_15mgml_finegrid2*5.00s', # [*] Specify
                #'filename'  : 'medium_G1_13mgml_*5.00s', # [*] Specify
                #'filename'  : 'medium_G2-3G1_20mgml_*x-2*5.00s', 
                #'filename'  : 'medium_as-synth_highC_fine*10.00s', Round 2 Sample1
                #'filename'  : 'medium_G2-2G1_highC_med*10.00s', 
                #'filename'  : dir_path+'large_G2-2G1_2_med*10.00s',
                'filename'  : '1_G1_Grid',  #14_As-synthesized_DEG_Grid',
                'feature_id': 2,  # [*] Specify
                'map_type': 'xy',
                'log10'  : 0,
                'verbose': 1,
               } 

feature_1_args = {'source_dir' : dir_path+'thumbnails/',
             'ext' : '.jpg',
             #'pixels' : [[89, 122], [102, 201], [152, 73]],  # [*] Choose pixels
             'pixels' : [[317, 191]],  # [190, 43], [*] Choose pixels
             'pixels_stat' : 0,     # [*] 'mean', 'max', 'var', or idx
             }

feature_2_args = {'source_dir' : dir_path+'circular_average/', #'../circular_average/',
             'ext' : '.dat',
             'data_col' : [0, 2],
             'q_targets' : [0.027], #0.053  # [*] Choose q0 or q0,q1
             'n' : 5     # [*] Choose the half-width (data points) of the peak q
             }

feature_3_args = {'source_dir' : dir_path+'linecut_angle/',
             'ext' : '.dat',
             'data_col' : [0, 1],
             'angle_targets': [30] #[21] # 'max' or specify angle 
             }
    
feature_args.update(feature_1_args=feature_1_args, feature_2_args=feature_2_args, feature_3_args=feature_3_args)

########## Feature map
for idx in [2]:
    feature_args.update(feature_id=idx); 
    
    ## Find matching files
    infiles, match_re = get_filematch(feature_args)  
    
    ## Get map
    scans, x_pos, y_pos, feature = get_map(infiles, match_re, feature_args) 
    
    ## Plot map
    fig = plt.figure(100+feature_args['feature_id'], figsize=[10,4]); plt.clf()
    ax1 = plt.subplot2grid((1, 5), (0, 3), colspan=2); ax1.cla()
    cmap = plt.get_cmap('viridis');    feature_args.update(cmap=cmap)
    #feature_args.update(val_stat = [0, 0.1])
    plot_map(x_pos, y_pos, feature, feature_args)
    
    #ax2 = plt.subplot2grid((1, 5), (0, 0), colspan=2); ax2.cla()
    #feature_args.update(val_stat = [0, 3])
    #plot_map(x_pos, y_pos, feature, feature_args)
    
    ## Plot one data 
    ax2 = plt.subplot2grid((1, 5), (0, 0), colspan=2); ax2.cla()
    cmap = plt.get_cmap('magma');  feature_args.update(cmap=cmap)
    
    #feature_args.update(filename='medium_G1_13mgml_*5.00s*74852') # sample 4
    #feature_args.update(filename='*83261') # Sample 1 70526
    #infiles, match_re = get_filematch(feature_args)
    plot_data(infiles[0], feature_args)


