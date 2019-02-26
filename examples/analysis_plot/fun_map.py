#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 13:51:10 2019

@author: etsai
"""

import time
import os, sys
import re
import glob
from scipy import ndimage
from scipy import interpolate
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import random
import PIL.Image as Image
from skimage import color
from skimage import io

# =============================================================================
# Load data from .dat 
# - Extract columns col[0] and col[1]
# =============================================================================   
def extract_data(filename, col):
    infile = open(filename, 'r')
    infile.readline() # skip the first line
    q = []
    I = []
    for line in infile:
        data = line.split()
        q.append(float(data[col[0]]))
        I.append(float(data[col[1]]))
    infile.close()
    return q, I

# =============================================================================
# Get files with matching dir, filename, ext
# - Depending on feature_id, it loads the corresponding args and files
# - sort: scan number (better implementation?)
# ============================================================================= 
def get_filematch(feature_args):
    filename = feature_args['filename']    
    feature_id = feature_args['feature_id']
    verbose = feature_args['verbose']
    if feature_id == 1:
        kwargs = feature_args['feature_1_args']
    elif feature_id == 2:
        kwargs = feature_args['feature_2_args']
    elif feature_id == 3:
        kwargs = feature_args['feature_3_args']
         
    source_dir = kwargs['source_dir']
    ext = kwargs['ext']

    pattern = filename+'*'+ext
    print(pattern)
    infiles = glob.glob(os.path.join(source_dir, pattern))
    infiles.sort()
    #infiles.sort(key=lambda name: int(name[-15:-9]))  #key=lambda x:float(re.findall("(\d+)",x)[0])
    
    #parse_re = '^.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_SAXS{}$'.format(ext)
    if feature_args['map_type']=='xy':
        parse_re = '^.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)_\w+{}$'.format(ext)
    elif feature_args['map_type']=='xT':
        parse_re = '^.+_x(-?\d+\.\d+)_T(-?\d+\.\d+)_.+_(\d+)_\w+{}$'.format(ext)
    else:
        print('Specify map type (eg. xy, T)!')
        match_re = [];
        return infilles, match_re
        
    match_re = re.compile(parse_re)    
    if verbose>0:
        print(pattern)
        print('Considering {} files...'.format(len(infiles)))    
    return infiles, match_re

# =============================================================================
# Get files with matching pattern
# ============================================================================= 
def get_filematch_s(pattern):
    infiles = glob.glob(pattern)
    infiles.sort()  
    parse_re = '^.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)$'
    match_re = re.compile(parse_re)      
    return infiles

# =============================================================================
# Given x y position, find the file
# =============================================================================
def find_file(xf, yf, feature_args):
    filename = feature_args['filename']
    feature_id = feature_args['feature_id']
    if feature_id == 1:
        kwargs = feature_args['feature_1_args']
    elif feature_id == 2:
        kwargs = feature_args['feature_2_args']
    elif feature_id == 3:
        kwargs = feature_args['feature_3_args']
    source_dir = kwargs['source_dir']
    
    n = filename.find('*') # assume before this is the sample name
    
    temp = '*x{:.3f}*_y{:.3f}*'.format(xf, yf) 
    temp = filename[0:n-1]+temp # ignore char filename[n]
    pattern = os.path.join(source_dir, temp) 
    infiles = get_filematch_s(pattern)
    return infiles

# =============================================================================
# Given x,y and a list of data positions, find the closest point with data
# =============================================================================
def get_closest(pos, post_list):# pos_list is 2 by N
    r_min = 1e10;
    for idx, item in enumerate(post_list[0]):
        x = post_list[0][idx]
        y = post_list[1][idx]
        r = calc_distance(pos, [x, y])
        if r<r_min:
            r_min = r
            xf = x; yf = y
            #idxf = int(idx)
    return xf, yf

def calc_distance(p0, p1):
    r =  math.hypot(p0[0]-p1[0], p0[1]-p1[1])
    return r
# =============================================================================
#
# Define features! 
#
# =============================================================================
def get_feature(infile, feature_args):
    log10 = feature_args['log10']
    feature_id = feature_args['feature_id']
    if feature_id == 1:
        kwargs = feature_args['feature_1_args']
    elif feature_id == 2:
        kwargs = feature_args['feature_2_args']
    elif feature_id == 3:
        kwargs = feature_args['feature_3_args']
    
    val = []    
    if feature_id == 1:
        pixels = kwargs['pixels']
        n = roi[0]
        im = color.rgb2gray(io.imread(infile))
        imarray = np.array(im)
        for pixel in pixels:
            temp_roi = imarray[pixel[1]-n:pixel[1]+n+1,pixel[0]-n:pixel[0]+n+1] #TEMP
            if roi[1]=='mean':
                temp = np.mean(temp_roi)
            elif roi[1]=='max':
                temp = np.max(temp_roi) 
            else:
                temp = imarray[pixel[1], pixel[0]]
            if log10: temp = np.log10(temp)
            val.append(temp) 
        
    elif feature_id == 2:  
        data_col = kwargs['data_col']
        q_targets = kwargs['q_targets']
        roi = kwargs['roi']
        n = roi[0]
        q, I = extract_data(infile, data_col)
        for q_target in q_targets:
            cen = get_target_idx(q, q_target)
            if roi[1]=='mean':
                temp = np.mean(I[cen-n:cen+n+1]) 
            elif roi[1]=='max':
                temp = np.max(I[cen-n:cen+n+1]) 
            else:
                temp = I[cen]                
            if log10: temp = np.log10(temp)

    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['angle_targets']
        angle_roi = kwargs['angle_roi']
        angle, I = extract_data(infile, data_col)
        temp = np.nan
        i0 = get_target_idx(angle, angle_roi[0])
        i1 = get_target_idx(angle, angle_roi[1])
        I_crop = I[i0:i1+1]
        for angle_target in angle_targets:            
            if angle_target =='max':
                if np.var(I_crop) > 0:
                    val.append(angle[i0+np.nanargmax(I_crop)])
            elif angle_target =='var':
                val.append(np.nanvar(I_crop))
            else: 
                try:
                    temp = I[get_target_idx(angle, angle_target)]
                except:
                    temp = np.nan
                val.append(temp)      # I(chi0)

    return val

# =============================================================================
# Get the index (for array q) closest to q_target
# =============================================================================  
def get_target_idx(q, target):
    q = np.array(q)
    idx = np.argmin(np.abs(q-target))
    return idx

# =============================================================================
# Fill the map: coordinates x, y, and the feature
# =============================================================================
def get_map(infiles, match_re, feature_args):
    scans = []
    x_pos = []
    y_pos = []
    feature = []
    for idx, infile in enumerate(infiles):
           
        filebase, filename = os.path.split(infile)
        m = match_re.match(filename)
        
        if m!=None:
            x = float(m.groups()[0]) 
            y = float(m.groups()[1]) # note: y is sometimes off by 0.5um because filename has only 3 decimal
            scan = int(m.groups()[2]) # scan number
            x_pos.append(x)
            y_pos.append(y)
            scans.append(scan)
    
            val = get_feature(infile, feature_args)
            feature.append(val)

    feature_args.update(val_stat=[np.min(feature), np.max(feature)])
    feature = np.asarray(feature)
    
    return scans, x_pos, y_pos, feature
    
    
# =============================================================================
# Plot one data
# =============================================================================        
def plot_data(infile, feature_args):
    log10 = feature_args['log10']
    feature_id = feature_args['feature_id']
    if feature_id == 1:
        kwargs = feature_args['feature_1_args']
    elif feature_id == 2:
        kwargs = feature_args['feature_2_args']
    elif feature_id == 3:
        kwargs = feature_args['feature_3_args']
    if 'cmap' in feature_args and feature_args['cmap']:
        cmap = feature_args['cmap']
    else:
        cmap = 'viridis'
    
    if feature_id == 1:
        pixels = kwargs['pixels']
        im = color.rgb2gray(io.imread(infile))
        if log10:
            im = np.log10(im)
        plt.imshow(im, cmap=cmap)
        plt.colorbar(shrink=0.8, aspect=24)
        for pixel in pixels:
            plt.plot(pixel[0],pixel[1], 'o', markersize=8, markeredgewidth=1, markeredgecolor='w', markerfacecolor='None')
        plt.title(infile)
        return im
    elif feature_id == 2: 
        q_targets = kwargs['q_targets']
        data_col = kwargs['data_col']
        n = kwargs['n']
        q, I = extract_data(infile, data_col)        
        I = np.log10(I)
        plt.plot(q, I)     
        for idx, q_target in enumerate(q_targets):
            # plot q_target 
            plt.plot([q_target, q_target], [-1, 4])
            plt.text(q_target, -0.9+idx*0.5, '('+str(q_target)+')')
            # plot integration region
            cen = get_target_idx(q, q_target)
            plt.plot([q[cen-n], q[cen+n]], [-1, -1]) 
        plt.ylabel('log10(I)')
        plt.xlabel('q ($\AA$^-1)')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)      
    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['angle_targets']
        angle_roi = kwargs['angle_roi']
        angle, I0 = extract_data(infile, data_col)
        I = np.log10(I0)
        plt.plot(angle, I)     
        if angle_targets =='max':
            i0 = get_target_idx(angle, angle_roi[0])
            i1 = get_target_idx(angle, angle_roi[1])
            plt.plot([angle[i0], angle[i1]], [0, 0])
            I_crop = I[i0:i1+1]
            val = angle[i0+np.argmax(I_crop)]
            plt.plot([val, val], [0, 3])
            plt.text(val, 0.1, str(np.round(val,3)))
        elif angle_targets =='var':
            val = np.var(I0)
        else: 
            for idx, angle_target in enumerate(angle_targets):
                plt.plot([angle_target, angle_target], [0, 0])
                plt.plot([angle_target, angle_target], [0, 3])
                plt.text(angle_target, 0.1+idx*0.1, '('+str(angle_target)+')')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
        plt.xlabel('$\chi$ (degree)')
        
    plt.title(infile)
 
# =============================================================================
# Plot map based on feature
# =============================================================================       
def plot_map(x_pos, y_pos, feature, feature_args):
    filename = feature_args['filename']
    val_stat = feature_args['val_stat']
    feature_id = feature_args['feature_id']
    if 'plot_interp' in feature_args:
        plot_interp = feature_args['plot_interp']
    else:
        plot_interp = ['none', 1]
    if 'cmap' in feature_args and feature_args['cmap']:
        cmap = feature_args['cmap'];
    else:
        cmap = plt.get_cmap('viridis')
    if feature_id == 1:
        kwargs = feature_args['feature_1_args']
    elif feature_id == 2:
        kwargs = feature_args['feature_2_args']
    elif feature_id == 3:
        kwargs = feature_args['feature_3_args']
    source_dir = kwargs['source_dir']
    
    if plot_interp[0]!='none':
        x_pos_fine, y_pos_fine, feature_fine = interp_map(x_pos, y_pos, feature, plot_interp)      
        plt.pcolormesh(x_pos_fine, y_pos_fine, feature_fine, vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 
        #plt.pcolormesh(x_pos, y_pos, feature) 
    else:
        plt.scatter(x_pos, y_pos, c=feature, marker="s", vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 
        
    plt.colorbar(shrink=1, pad=0.02, aspect=24);
    plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
    #plt.title(source_dir+filename)
    plt.axis('equal')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    
# =============================================================================
# Give interpolated map with finer discretization
# note - griddata works better than interpolate.interp2d
# =============================================================================          
def interp_map(x_pos, y_pos, feature, plot_interp): 
    x_ax_fine = np.arange(np.min(x_pos), np.max(x_pos), plot_interp[1]) 
    y_ax_fine = np.arange(np.min(y_pos), np.max(y_pos), plot_interp[1])
    x_pos_fine, y_pos_fine = np.meshgrid(x_ax_fine, y_ax_fine)
    feature_fine = griddata((x_pos, y_pos), feature, (x_pos_fine, y_pos_fine), method=plot_interp[0])
    return x_pos_fine, y_pos_fine, feature_fine

# =============================================================================
# Overlay three features
# =============================================================================       
def plot_overlay(x_pos, y_pos, feature_array, feature_args):
    fig = plt.figure(200, figsize=[8,8]); plt.clf()
    ax = fig.add_subplot(1, 1, 1)
    #ax.set_facecolor((0, 0, 0))
    ax.set_facecolor((1, 1, 1))
    
    if 'plot_interp' in feature_args:
        plot_interp = feature_args['plot_interp']
    else:
        plot_interp = ['none', 1]    
    overlay = []
    for idx, feature in enumerate(feature_array):
        x_pos_fine, y_pos_fine, feature_fine = interp_map(x_pos, y_pos, feature[0], plot_interp) 
        feature_fine = (feature_fine-np.nanmin(feature_fine)) / (np.nanmax(feature_fine)-np.nanmin(feature_fine))
        feature_fine[np.isnan(feature_fine)] = 1 #np.nanmean(feature_fine)
        if idx<=2:
            overlay.append(feature_fine)
        else:
            print('More then 3 features, only use the first three for RGB')
 
    while idx<2:
        overlay.append(feature_fine*0.0)
        idx = idx+1
    overlay = np.asarray(overlay)
    overlay = np.transpose(overlay, (1,2,0))
    extent = (np.nanmin(x_pos_fine), np.nanmax(x_pos_fine), np.nanmin(y_pos_fine), np.nanmax(y_pos_fine))
    plt.imshow(overlay, extent=extent,origin='lower')
    
    #plt.colorbar(shrink=1, pad=0.02, aspect=24);
    plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
    plt.axis('equal')
    #plt.axis('tight')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')

    
    