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
import numpy as np
import matplotlib.pyplot as plt
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
    infiles = glob.glob(os.path.join(source_dir, pattern))
    infiles.sort(key=lambda name: int(name[-15:-9]))  #key=lambda x:float(re.findall("(\d+)",x)[0])
    
    #parse_re = '^.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_SAXS{}$'.format(ext)
    parse_re = '^.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)_SAXS{}$'.format(ext)
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
    temp = '*x{:.2f}*_y{:.2f}*'.format(xf, yf) # ignoring decimal position 3
    temp = filename[0:n-1]+temp # ignore char filename[n]
    pattern = os.path.join(source_dir, temp) 
    infiles = get_filematch_s(pattern)
    return infiles

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
        
    if feature_id == 1:
        pixels = kwargs['pixels']
        pixels_stat = kwargs['pixels_stat']
        im = color.rgb2gray(io.imread(infile))
        imarray = np.array(im)
        val_list = []
        for idx, pixel in enumerate(pixels):
            temp = imarray[pixel[0],pixel[1]]
            if log10: temp = np.log10(temp)
            val_list.extend([temp]) 
        if pixels_stat=='mean':
            val = np.mean(val_list)
        elif pixels_stat=='max':
            val = np.max(val_list)
        elif pixels_stat=='var':
            val = np.var(val_list)
        else:
            val = val_list[pixels_stat]
        
    elif feature_id == 2:  
        data_col = kwargs['data_col']
        q_targets = kwargs['q_targets']
        n = kwargs['n']
        q, I = extract_data(infile, data_col)
        val_list = []
        for idx, q_target in enumerate(q_targets):
            cen = get_idx_q(q, q_target)
            temp = np.mean(I[cen-n:cen+n+1])  
            if log10: temp = np.log10(temp)
            if idx==0:
                val = temp      # ring intensity I(q0)
            else:
                val = val/(temp+1e-5)  # I(q0)/I(q1)

    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['angle_targets']
        angle, I = extract_data(infile, data_col)
        val = 0
        if angle_targets =='max':
            i0 = get_idx_q(angle, 5)
            i1 = get_idx_q(angle, 65)
            I_crop = I[i0:i1+1]
            if np.var(I_crop) > 0:
                val = angle[i0+np.nanargmax(I_crop)]
        else: 
            val_list = []
            for idx, angle_target in enumerate(angle_targets):
                temp = I[get_idx_q(angle, angle_target)]
                if idx==0:
                    val = temp      # I(chi0)
                else:
                    val = val/(temp+1e-5)  # I(chi0)/I(chi1)

    return val

# =============================================================================
# Get the index (for array q) closest to q_target
# =============================================================================  
def get_idx_q(q, target):
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
    #print('Done mapping')
    feature_args.update(val_stat=[np.min(feature), np.max(feature)])
    
    return scans, x_pos, y_pos, feature
    
    
# =============================================================================
# Plot one data
# =============================================================================        
def plot_data(infile, feature_args):
    feature_id = feature_args['feature_id']
    if feature_id == 1:
        kwargs = feature_args['feature_1_args']
    elif feature_id == 2:
        kwargs = feature_args['feature_2_args']
    elif feature_id == 3:
        kwargs = feature_args['feature_3_args']
    if feature_id == 1:
        pixels = kwargs['pixels']
        if 'cmap' in kwargs and kwargs['cmap']:
            cmap = kwargs['cmap']
        else:
            cmap = plt.get_cmap('viridis')
        im = color.rgb2gray(io.imread(infile))
        plt.imshow(im, cmap=cmap)
        plt.colorbar(shrink=0.8, aspect=24)
        for pixel in pixels:
            plt.plot(pixel[0],pixel[1], 'o', markersize=8, markeredgewidth=1, markeredgecolor='w', markerfacecolor='None')
        plt.title(infile)
        return im
    elif feature_id == 2: 
        q_targets = kwargs['q_targets']
        data_col = kwargs['data_col']
        q, I = extract_data(infile, data_col)
        I_log10 = np.log10(I)
        plt.plot(q, I_log10)     
        for idx, q_target in enumerate(q_targets):
            plt.plot([q_target, q_target], [-1, 4])
            plt.text(q_target,-1+idx*0.2, '('+str(q_target)+')')
        plt.ylabel('log10(I)')
        plt.xlabel('q ($\AA$^-1)')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)      
    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['angle_targets']
        angle, I = extract_data(infile, data_col)
        I_log10 = np.log10(I)
        plt.plot(angle, I_log10)     
        if angle_targets =='max':
            i0 = get_idx_q(angle, 5)
            i1 = get_idx_q(angle, 65)
            plt.plot([angle[i0], angle[i1]], [0, 0])
            I_crop = I[i0:i1+1]
            val = angle[i0+np.argmax(I_crop)]
            plt.plot([val, val], [0, 3])
        else: 
            for idx, angle_target in enumerate(angle_targets):
                plt.plot([angle_target, angle_target], [0, 0])
                plt.plot([angle_target, angle_target], [0, 3])
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
        
    plt.title(infile)
 
# =============================================================================
# Plot map based on feature
# =============================================================================       
def plot_map(x_pos, y_pos, feature, feature_args):
    filename = feature_args['filename']
    val_stat = feature_args['val_stat']
    feature_id = feature_args['feature_id']
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
    plt.scatter(x_pos, y_pos, c=feature, marker="s", vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 
    plt.colorbar(shrink=1, pad=0.02, aspect=24);
    plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
    plt.title(source_dir+filename)
    plt.axis('equal')
    plt.xlabel('x (mm)')
    #plt.ylabel('y (mm)')
    
    
    
    