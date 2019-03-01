#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 13:51:10 2019

@author: etsai
"""

import time, os, sys, re, glob, random, copy
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap
import PIL.Image as Image
from skimage import color
from skimage import io

from fun_ky import *
import lmfit
        
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
    exclude = feature_args['exclude']
    feature_id = feature_args['feature_id']
    verbose = feature_args['verbose']
    kwargs = feature_args['feature_{}_args'.format(feature_id)]        
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
    
    # Exclude some files
    for idx, infile in enumerate(infiles):
        for file in exclude:
            if infile.find(file)>-1:
                infiles.pop(idx)
    if verbose>0:
        print('  - Now considering {} files...'.format(len(infiles)))   
            
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
    kwargs = feature_args['feature_{}_args'.format(feature_id)] 
    source_dir = kwargs['source_dir']
    ext = kwargs['ext']
    
    n = filename.find('*') # assume before this is the sample name
    
    temp = '*x{:.3f}*_y{:.3f}*'.format(xf, yf) 
    temp = filename[0:n-3]+temp+ext  # ignore some char
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

# =============================================================================
# Calculate distance between 2 Cartersian points
# =============================================================================
def calc_distance(p0, p1):
    r =  math.hypot(p0[0]-p1[0], p0[1]-p1[1])
    return r

# =============================================================================
#
# Define features! 
#
# =============================================================================
def get_feature(infile, feature_args):
    log10 = feature_args['log10'][0]
    feature_id = feature_args['feature_id']
    kwargs = feature_args['feature_{}_args'.format(feature_id)] 
    
    val = []; info = [] # additional infor to store 
    if feature_id == 1:
        pixels = kwargs['targets']
        roi = kwargs['roi']
        n = roi[0]
        im = color.rgb2gray(io.imread(infile))
        imarray = np.array(im)
        if log10: imarray = np.log10(imarray)
        for pixel in pixels:
            temp_roi = imarray[pixel[1]-n:pixel[1]+n+1,pixel[0]-n:pixel[0]+n+1] #TEMP
            if roi[1]=='mean':
                temp = np.nanmean(temp_roi)
            elif roi[1]=='max':
                temp = np.nanmax(temp_roi) 
            else:
                temp = imarray[pixel[1], pixel[0]]
            val.append(temp)
        
    elif feature_id == 2:  
        data_col = kwargs['data_col']
        q_targets = kwargs['targets']
        roi = kwargs['roi']
        n = roi[0]
        #t1 = time.time()
        q, I = extract_data(infile, data_col)
        if log10: I = np.log10(I)
        #print('time.f2 = {}'.format(time.time()-t1))
        for q_target in q_targets:
            cen = get_target_idx(q, q_target)
            if roi[1]=='mean':
                temp = np.nanmean(I[cen-n:cen+n+1]) 
            elif roi[1]=='max':
                temp = np.nanmax(I[cen-n:cen+n+1]) 
            else:
                temp = I[cen] 
            val.append(temp)
        
    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['targets']
        angle_roi = kwargs['angle_roi']
        angle, I = extract_data(infile, data_col)
        if log10: I = np.log10(I)
        i0 = get_target_idx(angle, angle_roi[0])
        i1 = get_target_idx(angle, angle_roi[1])
        I_crop = I[i0:i1+1]
        for angle_target in angle_targets:  
            temp = np.nan
            if angle_target =='max':
                if np.var(I_crop) > 0: # TEMP 
                    temp = angle[i0+np.nanargmax(I_crop)]
            elif angle_target =='var':
                temp = np.nanvar(I_crop)
            else: 
                try:
                    temp = I[get_target_idx(angle, angle_target)]
                except:
                    print('Cannot find I[get_target_idx(angle, angle_target)] \n')
            val.append(temp)

    elif feature_id == 4:  
        data_col = kwargs['data_col']
        feats = kwargs['targets']
        q, I = extract_data(infile, data_col) 
        if log10: I = np.log10(I)
        line = DataLine(x=q, y=I)
        run_args = {'fit_range': [0.02, 0.06], 'sigma': 0.001, 'verbosity': 0}
        lm_result, fit_line, fit_line_extended = Protocols.circular_average_q2I_fit()._fit_peaks(line=line, q0=None, vary=True, **run_args)
        for feat in feats:
            if feat == 'd_spacing_nm':
                temp = 0.1*2.*np.pi/lm_result.params['x_center1'] #d in nm
            elif feat == 'grain_size_nm':
                temp = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/lm_result.params['sigma1'] #nm
            elif feat == 'chi2':
                temp = lm_result.chisqr/lm_result.nfree
            else:
                temp = lm_result.params[feat]  
            val.append(temp)
        info.append(line)
        info.append(fit_line)
            
    return val, info


# =============================================================================
# Get the index (for array q) closest to q_target
# =============================================================================  
def get_target_idx(q, target):
    q = np.array(q)
    idx = np.argmin(np.abs(q-target))
    return idx


# =============================================================================
# get_map(infiles, match_re, feature_args)
#
# Fill the map on coordinates x, y, with the feature
# Calls "get_feature" for each file
#
# Inputs: 
#   infiles: file list
#   match_re: extract positions x, y, and scan ID
#   feature_args: which feature, what specifics
# Outputs:
#   features_map:
#       scans: list of scans (not sure how useful yet)
#       x_pos
#       y_pos
#       tag: [id, feature_names]
#       features: for each feature_id, there can be several features
# =============================================================================
def get_map(infiles, match_re, feature_args):
    filename = feature_args['filename']
    feature_id = feature_args['feature_id']
    kwargs = feature_args['feature_{}_args'.format(feature_id)] 
    targets = kwargs['targets']
    ids = []
    tags = []
    scans = []
    x_pos = []
    y_pos = []
    features = [] 
    info_map = []
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
            ids.append(feature_id)
            tags.append(targets)
    
            val, info = get_feature(infile, feature_args) # val can be an array
            features.append(val)
            if info: info_map.append(info)

    features = np.asarray(features)
    features = (features.T).tolist()
    features_map = {'ids': ids, 'tags': tags, 'scans': scans, 'x_pos':x_pos, 'y_pos':y_pos, 'features':features, 
                   'info_map':info_map, 'filename': filename}
    
    return features_map
    
    
# =============================================================================
# Plot one data
# =============================================================================        
def plot_data(infile, **feature_args):
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
    
    if 'log10' in feature_args:
        log10 = feature_args['log10'][1]
    else:
        log10 = 0
    if 'feature_id' in feature_args:
        feature_id = feature_args['feature_id']
    else:
        feature_id = 1
    verbose = feature_args['verbose']        
    kwargs = feature_args['feature_{}_args'.format(feature_id)] 

    if 'cmap' in feature_args and feature_args['cmap']:
        cmap = feature_args['cmap']
    else:
        cmap = 'viridis'
       
    result = []
    if feature_id == 1:
        pixels = kwargs['targets']
        im = color.rgb2gray(io.imread(infile))
        if log10:
            im = np.log10(im)
        plt.imshow(im, cmap=cmap)
        plt.colorbar(shrink=0.8, aspect=24)
        for pixel in pixels:
            plt.plot(pixel[0],pixel[1], 'o', markersize=8, markeredgewidth=3, markeredgecolor='w', markerfacecolor='None')
        result = im
        
    elif feature_id == 2: 
        q_targets = kwargs['targets']
        data_col = kwargs['data_col']
        if 'roi' in kwargs:
            n = kwargs['roi'][0]
        else:
            n = 0
        q, I = extract_data(infile, data_col)        
        if log10: I = np.log10(I)
        plt.plot(q, I)     
        for idx, q_target in enumerate(q_targets):
            if type(q_target) is not str:
                plt.plot([q_target, q_target], [-1, 4])
                plt.text(q_target, -0.9+idx*0.5, '('+str(q_target)+')')
                # plot integration region
                cen = get_target_idx(q, q_target)
            plt.plot([q[cen-n], q[cen+n]], [-1, -1]) 
        plt.xlabel('q ($\AA$^-1)')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)  
        if log10: 
            plt.ylabel('log10(I)')
        else:
            plt.ylabel('Intensity (a.u.)')
   
    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['targets']
        angle_roi = kwargs['angle_roi']
        angle, I = extract_data(infile, data_col)
        if log10: I = np.log10(I)
        plt.plot(angle, I)  
        y_lim = [np.nanmin(I), np.nanmax(I)]
        for idx, angle_target in enumerate(angle_targets):
            if angle_target =='max':
                i0 = get_target_idx(angle, angle_roi[0])
                i1 = get_target_idx(angle, angle_roi[1])
                plt.plot([angle[i0], angle[i1]], [y_lim[0], y_lim[0]])
                I_crop = I[i0:i1+1]
                val = angle[i0+np.argmax(I_crop)]
                plt.plot([val, val], y_lim)
                plt.text(val, np.max(I_crop)*0.95, 'argmax='+str(np.round(val,2)))
            elif type(angle_target) is not str:
                plt.plot([angle_target, angle_target], [y_lim[0], y_lim[0]])
                plt.plot([angle_target, angle_target], y_lim)
                plt.text(angle_target, y_lim[0]+idx*0.1, '('+str(angle_target)+')')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
        plt.xlabel('$\chi$ (degree)')
        if log10: 
            plt.ylabel('log10(I)')
        else:
            plt.ylabel('Intensity (a.u.)')
            
    elif feature_id == 4:
        feats = kwargs['targets']
        val, info = get_feature(infile, feature_args)
        for line in info:
            I = np.log10(line.y)
            plt.plot(line.x, I) 
        
        ys = np.nanmin(I)
        yf = np.nanmax(I)
        space = yf*0.11
        xs = np.min(line.x)
        xf = np.max(line.x)*0.9
        for idx, feat in enumerate(feats):
            temp = np.asarray(val[idx]) #Note - type(val[0])=lmfit.parameter.Parameter; val[idx].value
            plt.text(xf, yf-space*idx, feat+'={:.3f}'.format(temp), **font)
        plt.ylim([ys*0.3, yf*1.2])
        plt.xlim([xs*0.5, xf*1.5])
        plt.xlabel('q ($\AA$^-1)')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
        
    if verbose: plt.title(infile)
    return result
 
# =============================================================================
# Plot map based on feature
# =============================================================================       
def plot_map(features_map, **kwargs):
    if 'filename' in features_map:
        filename = features_map['filename']
    else:
        filename = ''
    if 'ids' in features_map:
        ids = features_map['ids']
    else:
        print('WHY no ids?')
        return
    if 'tags' in features_map:
        tags = features_map['tags']
    else:
        print('WHY no tags?')
        return
    x_pos = features_map['x_pos']
    y_pos = features_map['y_pos']
    features = features_map['features']
    if 'log10' in kwargs:
        log10 = kwargs['log10'][1]
    else:
        log10 = 0
    if 'val_stat' in kwargs:
        val_stat = kwargs['val_stat']
    if 'cmap' in kwargs and kwargs['cmap']:
        cmap = kwargs['cmap'];
    else:
        cmap = plt.get_cmap('viridis')
    if 'plot_interp' in kwargs:
        plot_interp = kwargs['plot_interp']
    else:
        plot_interp = [None, 1]
    
    N_maps = len(features)
    for idx, feature in enumerate(features):
        ax = plt.subplot2grid((1, N_maps+1), (0, idx+1), colspan=1); 
        feature = np.asarray(feature)
        if log10:
            feature = np.log10(feature)
        if 'val_stat' not in kwargs:
            #val_stat = [np.nanmin(feature), np.mean([np.nanmedian(feature), np.nanmax(feature)]) ]
            val_stat = [np.nanmin(feature), np.nanmax(feature)]
        if plot_interp[0] is not None:
            #print('Plotting map using imshow')
            x_pos_fine, y_pos_fine, feature_fine = interp_map(x_pos, y_pos, feature, plot_interp) 
            #plt.pcolormesh(x_pos_fine, y_pos_fine, feature_fine, vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 
            extent = (np.nanmin(x_pos_fine), np.nanmax(x_pos_fine), np.nanmin(y_pos_fine), np.nanmax(y_pos_fine))
            plt.imshow(feature_fine, vmin=val_stat[0], vmax=val_stat[1], extent=extent, origin='lower')
        else:
            #print('Plotting map using scatter')
            plt.scatter(x_pos, y_pos, c=feature, marker="s", vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 
            
        plt.colorbar(shrink=1, pad=0.02, aspect=24);
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
        #plt.title(source_dir+filename)
        plt.title('Map {}'.format(idx))
        plt.axis('equal')
        plt.xlabel('x (mm)  [feature_id '+ str(ids[idx]) + ',  ' + str(tags[idx])+']')
        plt.ylabel('y (mm)')
    
# =============================================================================
# Give interpolated map with finer discretization
# Note - griddata works better than interpolate.interp2d
# =============================================================================          
def interp_map(x_pos, y_pos, feature, plot_interp): 
    x_ax_fine = np.arange(np.min(x_pos), np.max(x_pos), plot_interp[1]) 
    y_ax_fine = np.arange(np.min(y_pos), np.max(y_pos), plot_interp[1])
    x_pos_fine, y_pos_fine = np.meshgrid(x_ax_fine, y_ax_fine)
    feature_fine = griddata((x_pos, y_pos), feature, (x_pos_fine, y_pos_fine), method=plot_interp[0])
    return x_pos_fine, y_pos_fine, feature_fine

# =============================================================================
# Overlay three features (RGB)
#   features_map_list: list of feautres_maps, with len = # of feature ids
#   features = feature_map['features']: 1D or 2D array, axes are [postision, feature]
#   feature in features: 1D array, axis is [position]
#   feature_array: list of 1D or 2D arrays, from all the feature_ids, [postision, feature]
#
#   Example of features_map_list[0]['tag']:  
#       [4, ['b', 'prefactor1', 'd_spacing_nm', 'grain_size_nm', 'chi2']]
# =============================================================================       
def plot_overlay(features_map_list, **kwargs):
    if 'overlay_rgb' in kwargs:
        overlay_rgb = kwargs['overlay_rgb']
    else:
        overlay_rgb = [0, 1, 2]
    if 'log10' in kwargs:
        log10 = kwargs['log10'][1]    
    else:
        log10 = 0
    if 'plot_interp' in kwargs:
        plot_interp = kwargs['plot_interp']
        if plot_interp[0] is None:
            plot_interp[0] = 'linear' 
    else:
        plot_interp = ['linear', 1] 
         
    ## Get all the maps into one 2D array, feature_array
    features_map, legends = extract_maps(features_map_list)
    x_pos = features_map['x_pos']
    y_pos = features_map['y_pos']
    feature_array = features_map['features']
    
    ## Take three channels for plotting
    overlay = []; overlay_legend = []    
    if feature_array!=[]:
        fig = plt.figure(500, figsize=[10, 8]); plt.clf()
        ax = plt.subplot2grid((1, 5), (0, 0), colspan=4); ax.cla()
        ax.set_facecolor('k')
          
        ## Take three channels, interpolate to fine grid
        if len(feature_array)>3: 
            print('More then 3 features, using only {} for RGB'.format(overlay_rgb))
        for ii, feature in enumerate(feature_array):
            feature = np.asarray(feature)
            if log10: feature = np.log10(feature)
            x_pos_fine, y_pos_fine, feature_fine = interp_map(x_pos, y_pos, feature, plot_interp) 
            feature_fine = (feature_fine-np.nanmin(feature_fine)) / (np.nanmax(feature_fine)-np.nanmin(feature_fine)) # Normalize each channel
            feature_fine[np.isnan(feature_fine)] = 0  # Replace nan 
            if ii in overlay_rgb:
                overlay.append(feature_fine)
                overlay_legend.append(legends[ii])
     
        ## Fill empty channels
        nc = len(overlay) # number of channels (RGB)
        while nc<3:
            print('Less than 3 features, filling channel with 0')
            overlay.append(feature_fine*0.0)
            overlay_legend.append('empty')
            nc = nc+1
        
        ## Plot with imshow
        overlay = np.asarray(overlay)
        overlay = np.transpose(overlay, (1,2,0))
        extent = (np.nanmin(x_pos_fine), np.nanmax(x_pos_fine), np.nanmin(y_pos_fine), np.nanmax(y_pos_fine))
        plt.imshow(overlay, extent=extent,origin='lower')        
        plt.title('(R){}, (G){}, (B){}'.format(overlay_legend[0],overlay_legend[1],overlay_legend[2]))
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
        plt.axis('tight')
        plt.axis('equal')
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        
        ## Plot the colorcone
        ax2 = plt.subplot2grid((3, 5), (0, 4), colspan=1); ax2.cla()
        colorbar = Image.open('hsl_cone_graphic.jpg')
        plt.imshow(colorbar)
        plt.axis('off')
    else:
        print('feature_array is empty!\n')
    return overlay


# =============================================================================
# Extract maps from all feature_ids
# Example:
#    features_map, legends = extract_maps(features_map_list)
# Input:
#   features_map_list: list of feautres_maps, with len = # of feature ids
#       features = feature_map['features']: 1D or 2D array, axes are [postision, feature]
#       feature in features: 1D array, axis is [position]
# Output: 
#   features_map (see output of get_map)
#       x_pos, x_pos, tag
#       feature_array: list of 1D or 2D arrays, from all the feature_ids, [postision, feature]
#   legends: id, tag (e.g. grain_size_nm)
# =============================================================================
def extract_maps(features_map_list):
    feature_array = []; 
    ids = []       # feature_id
    tags = [] # feature name
    legends = []
    for ii, feature_map in enumerate(features_map_list): # ii the index for feature_ids
        if ii==0:
            x_pos = feature_map['x_pos']
            y_pos = feature_map['y_pos']
        features = feature_map['features']  # 2D map
        for jj, feature in enumerate(features):  # jj the index for each features within each feature_id
            feature_array.append(feature)
            id_here = features_map_list[ii]['ids'][jj]
            ids.append(id_here)
            for tag in features_map_list[ii]['tags'][jj]:
                tags.append(tag)
                legends.append('id={}, {}'.format(id_here, tag))
    
   
    # Repack into features_map (good/bad?)
    features_map = {}
    features_map.update(x_pos=x_pos, y_pos=y_pos, features=feature_array)
    features_map.update(ids=ids, tags=tags)
    
    return features_map, legends


# =============================================================================
# Do math on two selected features
# Example (feature_c = feature_a/feature_b):
#   feature_args['math_ab'] = [0, 1, 'divide'] 
#   feature_c = math_features(features_map_list, **feature_args)
# Output:
#   feature_c: new feature map
#   features_map_list: updated (appended) with a new feature_id (100+ii) contianing the new feature_c map
# =============================================================================
def math_features(features_map_list, **kwargs):
    print('Current features_map_list len = {}'.format(len(features_map_list)))
    feature_array = []; legends = []
    if 'math_ab' in kwargs:
        math_ab = kwargs['math_ab']
    else:
        math_ab = [1, 2, 'divide']
    if 'log10' in kwargs:
        log10 = kwargs['log10'][1]    
    else:
        log10 = 0
    if 'plot_interp' in kwargs:
        plot_interp = kwargs['plot_interp']
        if plot_interp[0] is None:
            plot_interp[0] = 'linear' 
    else:
        plot_interp = ['linear', 1] 
         
    ## Get all the maps into one 2D array, feature_array
    features_map, legends = extract_maps(features_map_list)
    feature_array = features_map['features']

    feature_a = np.asarray(feature_array[math_ab[0]])
    feature_b = np.asarray(feature_array[math_ab[1]])
    if math_ab[2] == 'divide':
        feature_c = feature_a / feature_b
    elif math_ab[2] == 'substract':
        feature_c = feature_a - feature_b   
    elif math_ab[2] == 'multiply':
        feature_c = feature_a * feature_b
    elif math_ab[2] == 'correlation':
        feature_c = np.corrcoef(feature_a, feature_b)  

    idx = len(features_map_list)-1
    current_id = int(np.asarray(features_map_list[idx]['ids'][0]))
    if current_id<100:
        math_id = 100
    else:
        math_id = current_id+1
    original_list = copy.deepcopy(features_map_list[idx])
    features_map_list.append(original_list)
    features_map_list[idx+1]['features'] = [feature_c] # see def get_map for features_map structure
    features_map_list[idx+1]['ids'].append(math_id)
    features_map_list[idx+1]['tags'].append(math_ab[2])
    
    print('Current features_map_list len = {}'.format(len(features_map_list)))
    
    return feature_c


# =============================================================================
# Count # of feature maps
# =============================================================================
def count_maps(features_map_list):
    N_maps = 0
    for ii, feature_map in enumerate(features_map_list): 
        N_maps += len(feature_map['features'])
        
    return N_maps









