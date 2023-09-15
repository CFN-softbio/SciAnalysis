import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
import glob, os, time
#import databroker
import imageio, random, math


#from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols
#from SciAnalysis.Result import *

### Temporary place for putting together various plotting functions

from scipy import interpolate
from scipy.interpolate import griddata
from scipy.signal import find_peaks

## Calculate interpolated map
def interp_map(x_pos, y_pos, feature, plot_interp): 
    x_ax_fine = np.arange(np.min(x_pos), np.max(x_pos), plot_interp[1]) 
    y_ax_fine = np.arange(np.min(y_pos), np.max(y_pos), plot_interp[1])
    x_pos_fine, y_pos_fine = np.meshgrid(x_ax_fine, y_ax_fine)
    feature_fine = griddata((x_pos, y_pos), feature, (x_pos_fine, y_pos_fine), method=plot_interp[0])
    feature_fine = np.asarray(feature_fine)
    return x_pos_fine, y_pos_fine, feature_fine


## Plot interpolated map (2D color image)
# Specify 3 lists: x, y, intensity
# Usage: 
#	plot_interp = ['linear', 0.01] 
#	x_pos_fine, y_pos_fine, feature = interp_map(x_pos, y_pos, intensity, plot_interp) 
# 	extent = (np.nanmin(x_pos_fine), np.nanmax(x_pos_fine), np.nanmin(y_pos_fine), np.nanmax(y_pos_fine))
# 	plt.imshow(feature, extent=extent, origin='lower', cmap='viridis') 


# =============================================================================
# Generate a random color, each channel with range (a,b), 0 dark
# =============================================================================
def rand_color(a, b):
    r = b-a
    color = (random.random()*r+a, random.random()*r+a, random.random()*r+a)
    return color

# =============================================================================
# Load a Line and find/return peaks
# =============================================================================
def plot_peaks(line, N_peaks_find=2, fit_param = [0,1,0.001], verbose=0, line_color = 'k', label_color='r', flag_log=[1, 0], roundup=3, gridline=True):
    if verbose>0:
        print('fit_param = [height, width, prominence]')

    if flag_log[0]==1:
        line_y = np.log(line.y)
    else:
        line_y = line.y

    if flag_log[1] == 1:
        line_x = np.log(line.x)
    else:
        line_x = line.x

    plt.plot(line_x, line_y, color = line_color); 
    peaks, _ = find_peaks(line_y, height=fit_param[0], width=fit_param[1], prominence=(fit_param[2], None))
    while len(peaks)>N_peaks_find:
        #print('  N_peaks = {}, increase fit_prom to reduce N_peaks'.format(len(peaks)))
        fit_param[2] = fit_param[2]*1.01
        peaks, _ = find_peaks(line_y, height=fit_param[0], width=fit_param[1], prominence=(fit_param[2], None)) 
    print('{} peaks found: {}'.format(len(peaks), np.round(line_x[peaks],4)) +' for fit_prom {:.5f}'.format(fit_param[2]))
    
    ylim = [np.nanmin(line_y[line_y != -np.inf]), np.nanmax(line_y)]
    yrange = ylim[1]-ylim[0]
    for idx, peak in enumerate(peaks):
        plt.plot([line_x[peak], line_x[peak]], ylim, '--', color=label_color)
        if verbose > 0:
            plt.text(line_x[peak], ylim[1]-np.mod(idx,4)*yrange*0.04, str(np.round(line_x[peak],roundup)),fontweight='bold')
    if gridline:
        plt.grid(True, which='major', color='k', linestyle='-', alpha=0.3) 
    
    return line_x[peaks]


# =============================================================================
# Label peaks
# =============================================================================
def label_peaks(q_labels, flag_log=[1, 0], y_range=[0, 1], color='g', roundup=3, verbose=1):
    for idx, q in enumerate(q_labels):
        if flag_log[0]==1:
            plot_range = [np.log(np.max([1, y_range[0]])), np.log(y_range[1])]
        else:
            plot_range = [y_range[0], y_range[1]]
            
        if math.isnan(q)==False:
            plt.plot([q, q], plot_range, color=color,  alpha=0.5)
            if verbose>0:
                plt.text(q, plot_range[1]/2+np.mod(idx,4)*(plot_range[1]-plot_range[0])*0.04, str(np.round(q, roundup)), color=color)
            









