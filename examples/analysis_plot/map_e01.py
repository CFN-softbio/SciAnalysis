
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
feature_args = {'filename'  : 'large_G1_15mgml_x*5.00s', # [*] Specify
                'feature_id': 1,  # [*] Specify
                'log10'  : 1,
                'verbose': 0,
               } 

feature_1_args = {'source_dir' : '../thumbnails2/',
             'ext' : '.jpg',
             'pixels' : [[89, 122], [102, 201], [152, 73]],  # [*] Choose pixels
             #'pixels' : [[186, 111], [191, 115]],  # [*] Choose pixels
             'pixels_stat' : 2,     # [*] 'mean', 'max', 'var', or idx
             }

feature_2_args = {'source_dir' : '../circular_average/',
             'ext' : '.dat',
             'data_col' : [0, 2],
             'q_targets' : [0.0275], #0.053  # [*] Choose q0 or q0,q1
             'n' : 5     # [*] Choose the half-width (data points) of the peak q
             }

feature_3_args = {'source_dir' : '../linecut_angle/',
             'ext' : '.dat',
             'data_col' : [0, 1],
             'angle_targets': 'max' # 'max' or specify angle 
             }
    
feature_args.update(feature_1_args=feature_1_args, feature_2_args=feature_2_args, feature_3_args=feature_3_args)

########## Feature map
for idx in [1,2]:
    feature_args.update(feature_id=idx); 
    
    ## Find matching files
    infiles, match_re = get_filematch(feature_args)  
    
    ## Get map
    x_pos, y_pos, feature = get_map(infiles, match_re, feature_args) 
    feature_args.update(val_stat=[np.min(feature), np.max(feature)])
    
    ## Plot map
    fig = plt.figure(100+feature_args['feature_id'], figsize=[15,5]); plt.clf()
    ax1 = plt.subplot2grid((1, 5), (0, 2), colspan=3); ax1.cla()
    cmap = plt.get_cmap('viridis');    feature_args.update(cmap=cmap)
    #feature_args.update(val_stat = [0, 20])
    plot_map(x_pos, y_pos, feature, feature_args)
    
    ## Plot one data 
    ax2 = plt.subplot2grid((1, 5), (0, 0), colspan=2); ax2.cla()
    cmap = plt.get_cmap('magma');  feature_args.update(cmap=cmap)
    plot_data(infiles[0], feature_args)



