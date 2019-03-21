#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from fun_ky import * # SciAnalysis_PATH defined here, TEMP solution

import time, io, os, sys, re, glob, random, copy, gc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy import ndimage
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap
import PIL.Image as PILImage
from skimage import color
import skimage.io as skiio

import lmfit
from scipy.signal import find_peaks

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
        
from bqplot import *
import bqplot.pyplot as bplt
#import ipywidgets as widgets

mpl.rcParams['font.size'] = 15
#mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['figure.titlesize'] = 12
mpl.rcParams['lines.linewidth'] = 2 
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12

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
    q = np.asarray(q)
    I = np.asarray(I)
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
    verbose = feature_args['verbose']  if 'verbose' in feature_args else 0
    
    direct = feature_args['direct']
    if direct:
        source_dir = feature_args['source_dir']
        ext = '.tif'
    else:
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
        print('Considering {} files...'.format(len(infiles)))   
    
    # Exclude some files
    infiles_new = []
    for idx, infile in enumerate(infiles):
        include = True
        for file in exclude:
            if infile.find(file)>-1:
                include = False
        if include: 
            infiles_new.append(infiles[idx])
                
    if verbose>0:
        print('  - Now considering {} files...'.format(len(infiles_new)))   
            
    return infiles_new, match_re

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

    direct = feature_args['direct']
    if direct:
        source_dir = feature_args['source_dir']
        ext = '.tif'
    else:
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
# For file 'infile', get features and store as lists in val
# Input: 
#   infile, feature_id (in feature_args)
# Output: 
#   lists: val and info
# =============================================================================
def get_feature(infile, feature_args):
    direct = feature_args['direct']  # direct processing from .tif
    verbose = feature_args['verbose']
    feature_id = feature_args['feature_id']
    verbose = feature_args['verbose']
    kwargs = feature_args['feature_{}_args'.format(feature_id)] 
 
    if direct==1:
        if verbose:
            print('Processing {}'.format(infile))
        process = feature_args['process']
        protocols = kwargs['protocols']
        if verbose<=1: # Suppress output, temporary solution
            text_trap = io.StringIO()
            sys.stdout = text_trap
        result = process.run([infile], protocols, output_dir='./', force=True, store=False)
        if verbose: 
            sys.stdout = sys.__stdout__       
    
    val = []; info = [] # additional info to store 
    if feature_id == 1:
        pixels = kwargs['targets']
        roi = kwargs['roi']
        n = roi[0]
        
        if direct:
            imarray = np.array(result[0].data)  
            x_axis = result[0].x_axis
            y_axis = result[0].y_axis
            x_scale = result[0].x_scale
            y_scale = result[0].y_scale
        else:
            #im = color.rgb2gray(skiio.imread(infile)); imarray = np.array(im)
            im = np.load(infile).items()
            imarray = np.array(im[0][1])   
            x_axis = im[1][1] 
            y_axis = im[2][1] 
            x_scale = im[3][1] 
            y_scale = im[4][1]            
        image_dict = {'image': imarray, 'x_axis': x_axis, 'y_axis': y_axis, 'x_scale': x_scale, 'y_scale':y_scale}
        
        for pixel in pixels:
            temp_roi = imarray[pixel[1]-n:pixel[1]+n+1,pixel[0]-n:pixel[0]+n+1] #TEMP
            if roi[1]=='mean':
                temp = np.nanmean(temp_roi)
            elif roi[1]=='max':
                temp = np.nanmax(temp_roi) 
            else:
                temp = imarray[pixel[1], pixel[0]]
            val.append(temp)
        info.append(image_dict)
        
    elif feature_id == 2:  
        data_col = kwargs['data_col']
        q_targets = kwargs['targets']
        roi = kwargs['roi']
        n = roi[0]      
        
        if direct:
            line = result[0]
            q = line.x; I = line.y
        else:
            q, I = extract_data(infile, data_col)
            line = DataLine(x=q, y=I)
       
        for q_target in q_targets:
            cen = get_target_idx(q, q_target)
            if roi[1]=='mean':
                temp = np.nanmean(I[cen-n:cen+n+1]) 
            elif roi[1]=='max':
                temp = np.nanmax(I[cen-n:cen+n+1]) 
            else:
                temp = I[cen] 
            val.append(temp)
        info.append(line)
        
    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['targets']
        angle_roi = kwargs['angle_roi']
        normalize = kwargs['normalize']
        if type(angle_roi[1]) is str:
            N_fold = int(np.asarray(angle_roi[0]))
            stat = angle_roi[1]
        else:
            N_fold = 0           
        
        if direct:
            line_orig = result[0]
            angle = line_orig.x; I = line_orig.y
        else:
            angle, I = extract_data(infile, data_col)        
            line_orig = DataLine(x=angle, y=I)
        info.append(line_orig)
        
        if N_fold:
            angle_fold, I_fold_stat = fold_line(angle, I, N_fold, verbose)
            if stat=='mean':
                I_fold = I_fold_stat[:,1] 
            elif stat=='max':
                I_fold = I_fold_stat[:,2]   
            line_fold = DataLine(x=angle_fold, y=I_fold)
            info.append(line_fold)            
        else:
            i0 = get_target_idx(angle, angle_roi[0])
            i1 = get_target_idx(angle, angle_roi[1])
            I_crop = I[i0:i1+1]
        
        for angle_target in angle_targets:  
            temp = np.nan
            if angle_target =='argmax':
                if N_fold: 
                    temp = angle_fold[np.nanargmax(I_fold)]
                else:
                    temp = angle[i0+np.nanargmax(I_crop)]
            elif angle_target =='var':
                temp = np.nanvar(I)
            else: 
                if N_fold:
                    temp = I_fold[get_target_idx(angle_fold, angle_target)] 
                else:
                    temp = I[get_target_idx(angle, angle_target)] 
                if normalize:
                    temp = temp / np.sum(I)
            val.append(temp)

    elif feature_id == 4:  
        data_col = kwargs['data_col']
        feats = kwargs['targets']
        fit_range = kwargs['fit_range']
        chi2_thr = kwargs['chi2_thr']
        
        if direct:
            line = result[0]
            q = line.x; I = line.y
        else:
            q, I = extract_data(infile, data_col)
            line = DataLine(x=q, y=I)
        
        run_args = {'fit_range': fit_range, 'sigma': 0.001, 'verbosity': 0}
        lm_result, fit_line, fit_line_extended = Protocols.circular_average_q2I_fit()._fit_peaks(line=line, q0=None, vary=True, **run_args)
        chi2 = lm_result.chisqr/lm_result.nfree
        for feat in feats:
            if feat == 'd_spacing_nm':
                temp = 0.1*2.*np.pi/lm_result.params['x_center1'] #d in nm
            elif feat == 'grain_size_nm':
                temp = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/lm_result.params['sigma1'] #nm
            elif feat == 'chi2':
                temp = chi2
            else:
                temp = lm_result.params[feat]  
            ## Threshold 
            if chi2> chi2_thr or (feat is 'b') or (feat is 'prefactor1'): 
                temp = np.asarray(temp)
            else:
                temp = np.nan
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
    ids = []
    tags = []
    scans = []
    x_pos = []
    y_pos = []
    features = [] 
    info_map = []
    
    for target in kwargs['targets']:
        ids.append(feature_id) 
        tags.append(target) 
 
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
    
            val, info = get_feature(infile, feature_args) # val can be an array
            features.append(val)
            if info: 
                info_map.append(info)

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
    direct = feature_args['direct']  
    val_stat = feature_args['val_stat'] if 'val_stat' in feature_args else 'auto'
    log10 = feature_args['log10'] if 'log10' in feature_args else 0
    feature_id = feature_args['feature_id'] if 'feature_id' in feature_args else 1
    subplot = feature_args['subplot'] if 'subplot' in feature_args else 0
    verbose = feature_args['verbose']   if 'verbose' in feature_args else 0
        
    kwargs = feature_args['feature_{}_args'.format(feature_id)] 
    if 'cmap' in feature_args and feature_args['cmap']:
        cmap = feature_args['cmap']
    else:
        cmap = 'jet'
       
    val, info = get_feature(infile, feature_args)  
    if feature_id == 1:
        pixels = kwargs['targets']
        imarray = info[0]['image']
        if log10: 
            imarray[imarray<=0] = 1e-1
            imarray = np.log10(imarray)

        if val_stat=='auto':
            val_stat = [np.nanmin(imarray), np.nanmax(imarray)]   

        x_axis = info[0]['x_axis']
        y_axis = info[0]['y_axis']
        extent = (np.nanmin(x_axis), np.nanmax(x_axis), np.nanmin(y_axis), np.nanmax(y_axis)) 
        if subplot==1: # Plot also the q-axis
            host = host_subplot(111,axes_class=AA.Axes)
            plt.subplots_adjust(right=0.8)       
            host.imshow(imarray, cmap=cmap, origin='lower', vmin=val_stat[0], vmax=val_stat[1])  
            host.set_facecolor('k')
            host.colorbar()
            par2 = host.twinx()
            new_fixed_axis = par2.get_grid_helper().new_fixed_axis
            par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(10, 0))
            par2.axis["right"].toggle(all=True)
            par2.set_ylabel("q ($\AA^{-1}$)", color='k')
            par2.set_ylim(extent[2],extent[3]) 
        else:
            plt.imshow(imarray, cmap=cmap, extent=extent, origin='lower', vmin=val_stat[0], vmax=val_stat[1])   
            ax = plt.gca()
            ax.set_facecolor('k')
            plt.colorbar()
        for pixel in pixels: # Label pixels of interest
            plt.plot(pixel[0],pixel[1], 'o', markersize=8, markeredgewidth=1, markeredgecolor='w', markerfacecolor='None')
        if verbose>1: # Extra plot with q-axes
            plt.figure(161)
            plt.imshow(imarray, cmap=cmap, origin='lower') 
    
    elif feature_id == 2: 
        q_targets = kwargs['targets']
        data_col = kwargs['data_col']
        N_peaks_find = kwargs['N_peaks_find']
        if 'roi' in kwargs:
            n = kwargs['roi'][0]
        else:
            n = 0
        q = info[0].x
        I = info[0].y
        if log10: I = np.log10(I)
        y_lim = [np.nanmin(I[I!=-np.inf]), np.nanmax(I)]
        y_range = y_lim[1] - y_lim[0]
        #plt.plot(q, I)     
        plot_peaks(DataLine(x=q,y=I), N_peaks_find, 1)
        if subplot==1:
            for idx, q_target in enumerate(q_targets):
                if type(q_target) is not str:
                    plt.plot([q_target, q_target], y_lim)
                    plt.text(q_target, y_lim[0]-y_range*0.05, '('+str(q_target)+')')
                    # plot integration region
                    cen = get_target_idx(q, q_target)
                plt.plot([q[cen-n], q[cen+n]], [y_lim[0], y_lim[0]])         
        plt.xlabel('q ($\AA$^-1)')
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)  
        if subplot==1: 
            if log10: 
                plt.ylabel('log10(I)')
            else:
                plt.ylabel('Intensity (a.u.)')
   
    elif feature_id == 3:  
        data_col = kwargs['data_col']
        angle_targets = kwargs['targets']
        angle_roi = kwargs['angle_roi']
        N_peaks_find = kwargs['N_peaks_find']
        protocols = kwargs['protocols']
        if type(angle_roi[1]) is str:
            N_fold = int(np.asarray(angle_roi[0]))
            stat = angle_roi[1]
        else:
            N_fold = 0
        angle = info[0].x   # the original curve
        I = info[0].y       # the original curve
        if log10: I = np.log10(I)
        
        if N_fold:
            angle_fold = info[1].x  # the folded curve
            I_fold = info[1].y      # the folded curve
            if log10: 
                I_fold = np.log10(I_fold)
            if subplot==1: ax2 = plt.subplot2grid((2, 1), (1, 0), colspan=1) 
            if subplot>=0:
                for nn in np.arange(-N_fold/2,N_fold/2):
                    plt.plot(angle_fold+nn*360/N_fold, I_fold)  
                    #if nn==0:
                        #plot_peaks(DataLine(x=angle_fold, y=I_fold), N_peaks_find, 0)
                plt.xlim([np.min(angle), np.max(angle)])
            elif subplot==-1: # TEMP SOL
                #plt.plot(angle_fold, I_fold) 
                plot_peaks(DataLine(x=angle_fold, y=I_fold), N_peaks_find, 1)
                plt.xlim([np.min(angle_fold), np.max(angle_fold)])
            plt.xlabel('$\chi$ (degree), q={:.3f}, dq={:.3g}'.format(protocols[0].run_args['q0'], protocols[0].run_args['dq']))
            plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)

        y_lim = [np.nanmin(I), np.nanmax(I)]
        if subplot==1: 
            ax1 = plt.subplot2grid((2, 1), (0, 0), colspan=1) 
        
        if subplot==1 or N_fold==0:
            plt.plot(angle, I) 
            #plot_peaks(DataLine(x=angle,y=I), N_peaks_find, 0)
            plt.xlim([np.min(angle), np.max(angle)])
            for idx, angle_target in enumerate(angle_targets):
                if angle_target =='argmax':
                    if N_fold==0:                    
                        ax1.plot([val[idx], val[idx]], y_lim,'--')
                        ax1.plot(angle_roi, [y_lim[0], y_lim[0]])
                        ax1.text(val[idx], y_lim[1], 'argmax=('+str(np.round(val[idx],2))+')')
                    else:
                        ax2.plot([val[idx], val[idx]], y_lim,'--')
                        ax2.text(val[idx], y_lim[1], 'argmax='+str(np.round(val,2)))
                elif type(angle_target) is not str:
                    plt.plot([angle_target, angle_target], [y_lim[0], y_lim[0]])
                    plt.plot([angle_target, angle_target], y_lim,'--')
                    plt.text(angle_target, y_lim[0]+idx*0.1, '('+str(angle_target)+')')
            plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
            plt.xlabel('$\chi$ (degree) at q={:.3f}, dq={:.3g}'.format(protocols[0].run_args['q0'], protocols[0].run_args['dq']))
            if log10: 
                plt.ylabel('log10(I)')
            else:
                plt.ylabel('Intensity (a.u.)')
            
    elif feature_id == 4:
        feats = kwargs['targets']
        val, info = get_feature(infile, feature_args)
        for ii, line in enumerate(info):
            I = line.y
            if log10: I = np.log10(line.y)
            if ii==0: 
                plt.plot(line.x, I) 
            else:            
                plt.plot(line.x, I, '--') 
        
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
        
    if verbose: 
        ll = len(infile)
        l0 = int(np.round(ll/3))
        plt.title(infile[0:l0]+'\n'+infile[l0:l0*2]+'\n'+infile[l0*2:])
                
# =============================================================================
# Plot map based on feature
# =============================================================================       
def plot_map(features_map, **feature_args):
    val_stat = feature_args['val_stat'] if 'val_stat' in feature_args else 'auto'
    log10 = feature_args['log10'] if 'log10' in feature_args else 0
    subplot = feature_args['subplot'] if 'subplot' in feature_args else 0
    verbose = feature_args['verbose']   if 'verbose' in feature_args else 0
    if 'cmap' in feature_args and feature_args['cmap']:
        cmap = feature_args['cmap']
    else:
        cmap = plt.get_cmap('viridis')
    
    ids = features_map['ids'] 
    tags = features_map['tags']
    x_pos = features_map['x_pos']
    y_pos = features_map['y_pos']
    features = features_map['features']
    if 'plot_interp' in feature_args and len(x_pos)>1:
        plot_interp = feature_args['plot_interp'] 
    else:
        plot_interp = [None, 1]
        
    N_maps = len(features)
    for idx, feature in enumerate(features):
        if subplot==1: ax = plt.subplot2grid((2, N_maps), (0, idx), colspan=1); 
        
        feature = np.asarray(feature)
        if log10:
            feature = np.log10(feature)
        #if val_stat == 'auto':
        val_stat = [0.1*np.nanmedian(feature)+0.9*np.nanmin(feature), 0.1*np.nanmedian(feature)+0.9*np.nanmax(feature) ]
        #val_stat = [np.nanmin(feature), np.nanmax(feature)]
        if plot_interp[0] is not None:
            #print('Plotting map using imshow')
            x_pos_fine, y_pos_fine, feature = interp_map(x_pos, y_pos, feature, plot_interp) 
            #plt.pcolormesh(x_pos_fine, y_pos_fine, feature, vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 
            extent = (np.nanmin(x_pos_fine), np.nanmax(x_pos_fine), np.nanmin(y_pos_fine), np.nanmax(y_pos_fine))
            plt.imshow(feature, vmin=val_stat[0], vmax=val_stat[1], extent=extent, origin='lower', cmap=cmap) 
        else:
            #print('Plotting map using scatter')
            plt.scatter(x_pos, y_pos, c=feature, marker="s", vmin=val_stat[0], vmax=val_stat[1], cmap=cmap) 

        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
        plt.axis('tight')
        plt.axis('equal')
        plt.xlabel('[feature_id '+ str(ids[idx]) + ',  ' + str(tags[idx])+']')            
        if idx==0: plt.ylabel('y (mm)')
        if subplot==1: 
            plt.colorbar(shrink=1, pad=0.02, aspect=24);
            plt.title('Map {}'.format(idx))
        
        if subplot==1: # Plot histogram
            ax = plt.subplot2grid((2, N_maps), (1, idx), colspan=1); 
            feature = feature[~np.isnan(feature)]
            N_bins = 50;
            plt.hist(feature, bins=N_bins, color=rand_color(0.3,0.9))
            #plt.xlabel('val (bins={})'.format(N_bins))
            if idx==0: plt.ylabel('count')
            plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.25)
            

    
# =============================================================================
# Give interpolated map with finer discretization
# Note - griddata works better than interpolate.interp2d
# =============================================================================          
def interp_map(x_pos, y_pos, feature, plot_interp): 
    x_ax_fine = np.arange(np.min(x_pos), np.max(x_pos), plot_interp[1]) 
    y_ax_fine = np.arange(np.min(y_pos), np.max(y_pos), plot_interp[1])
    x_pos_fine, y_pos_fine = np.meshgrid(x_ax_fine, y_ax_fine)
    feature_fine = griddata((x_pos, y_pos), feature, (x_pos_fine, y_pos_fine), method=plot_interp[0])
    feature_fine = np.asarray(feature_fine)
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
    if 'rgb' in kwargs:
        rgb = kwargs['rgb'] 
    else:
        rgb = 'RGB'       
    if 'normalize_each' in kwargs:
        normalize_each = kwargs['normalize_each'] 
    else:
        normalize_each = True
    log10 = kwargs['log10'] 
    if 'plot_interp' in kwargs:
        plot_interp = kwargs['plot_interp']
        if plot_interp[0] is None:
            plot_interp[0] = 'linear' 
    else:
        plot_interp = ['linear', 1] 
         
    ## Get all the maps into one 2D array, feature_array
    features_map = extract_maps(features_map_list)
    x_pos = features_map['x_pos']
    y_pos = features_map['y_pos']
    feature_array = features_map['features']
    ids = features_map['ids']
    tags = features_map['tags']
    
    ## Take three channels for plotting
    overlay = []; overlay_legend = [] ; channel = 0; #rgb = 'RGB'
    if feature_array!=[] and len(feature_array)>=len(rgb):
        fig = plt.figure(500, figsize=[10, 8]); plt.clf()
        
        ## Get max and min for normalization 
        max_val = -1; min_val = 1e15;
        for ii in overlay_rgb:
            if np.nanmax(feature_array[ii]) > max_val:
                max_val = np.nanmax(feature_array[ii])
            if np.nanmin(feature_array[ii]) < min_val:
                min_val = np.nanmin(feature_array[ii])
        if log10: 
            min_val = np.log10(min_val)         
            max_val = np.log10(max_val) 

        ## Take three channels, interpolate to fine grid
        if len(feature_array)>3: 
            print('More then 3 features available, using only {} for RGB'.format(overlay_rgb))
        for ii in overlay_rgb:
            feature = np.asarray(feature_array[ii])
            if log10: feature = np.log10(feature)
            x_pos_fine, y_pos_fine, feature_fine = interp_map(x_pos, y_pos, feature, plot_interp) 
            if normalize_each:
                feature_fine = (feature_fine-np.nanmin(feature_fine)) / (np.nanmax(feature_fine)-np.nanmin(feature_fine)) # Normalize each channel
            else:
                feature_fine = (feature_fine-min_val) / (max_val-min_val) # Normalize wrt max_val 
            feature_fine[np.isnan(feature_fine)] = 0  # Replace nan 

            ## Plot each channel                
            ax = plt.subplot2grid((3, 7), (channel, 0), colspan=2); 
            image_channel = np.asarray(image_RGB(feature_fine, rgb[channel]))
            if overlay==[]:
                overlay = image_channel
                extent = (np.nanmin(x_pos_fine), np.nanmax(x_pos_fine), np.nanmin(y_pos_fine), np.nanmax(y_pos_fine))
            else: 
                overlay += image_channel
            plt.imshow(image_channel, extent=extent, origin='lower') 
            plt.title('({}) id={}, {}, [min, max]=[{:.3f}, {:.3f}]'.format(rgb[channel], ids[ii], tags[ii], np.nanmin(feature), np.nanmax(feature)))
            channel += 1
        
        ## Plot with imshow
        ax = plt.subplot2grid((3, 7), (0, 2), rowspan=3, colspan=4); ax.cla()
        ax.set_facecolor('k')
        plt.imshow(overlay, extent=extent,origin='lower')        
        plt.title('normalize_each={}, log10={}'.format(normalize_each,log10))
        #plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
        plt.axis('tight')
        plt.axis('equal')
        plt.xlabel('x (mm)')
        #plt.ylabel('y (mm)')
        if 0:
            plt.plot([-3.26, -3.26], [-7.46, -7.07],color='w')
            plt.plot([-3.26, -2.87], [-7.46, -7.46],color='w')
            
            plt.plot([-3.16, -3.16], [-6.96, -6.57],color='w')
            plt.plot([-3.16, -2.77], [-6.96, -6.96],color='w')
            # medium_G1_13mgml finegrid
            plt.plot([-3.116, -3.116], [-6.715, -6.568],color='w')
            plt.plot([-3.116, -2.98], [-6.715, -6.715],color='w')
            # [-1.375, -0.655] [-6.375, -5.655]
        if 0:
            plt.plot([-1.375, -1.375], [-6.375, -5.655],color='w')
            plt.plot([-1.375, -0.655], [-6.375, -6.375],color='w')
            plt.plot([-1.375, -0.655], [-5.655, -5.655],color='w')
            plt.plot([-0.655, -0.655], [-6.375, -5.655],color='w')
            
        ## Plot the colorcone
        ax2 = plt.subplot2grid((3, 7), (0, 6), colspan=1); ax2.cla()
        colorcone = PILImage.open(SciAnalysis_PATH+'examples/analysis_plot/hsl_cone_graphic.jpg')
        plt.imshow(colorcone)
        plt.axis('off')
    else:
        print('feature_array is empty or too short! No overlay plotted\n')
    return overlay

# =============================================================================
# 
# =============================================================================
def image_RGB(image, rgb):
    dim = image.shape
    image_stack = np.zeros([dim[0], dim[1], 3])    
    if 'R' in rgb:
        image_stack[:,:,0] = image
    if 'G' in rgb:
        image_stack[:,:,1] = image     
    if 'B' in rgb:
        image_stack[:,:,2] = image
    
    return image_stack


# =============================================================================
# Extract maps from all feature_ids
# Example:
#    features_map = extract_maps(features_map_list)
# Input:
#   features_map_list: list of feautres_maps, with len = # of feature ids
#       features = feature_map['features']: 1D or 2D array, axes are [postision, feature]
#       feature in features: 1D array, axis is [position]
# Output: 
#   features_map (see output of get_map)
#       x_pos, x_pos, tag
#       feature_array: list of 1D or 2D arrays, from all the feature_ids, [postision, feature]
# =============================================================================
def extract_maps(features_map_list):
    feature_array = []; 
    ids = []  # feature_id
    tags = [] # feature name
    for ii, feature_map in enumerate(features_map_list): # ii the index for feature_ids
        if ii==0:
            x_pos = feature_map['x_pos']
            y_pos = feature_map['y_pos']
        features = feature_map['features']  # 2D map
        for jj, feature in enumerate(features):  # jj the index for each features within each feature_id
            feature_array.append(feature)
            ids.append(features_map_list[ii]['ids'][jj])
            tags.append(features_map_list[ii]['tags'][jj])
   
    # Repack into features_map (good/bad?)
    features_map = {}
    features_map.update(x_pos=x_pos, y_pos=y_pos, features=feature_array)
    features_map.update(ids=ids, tags=tags)
    
    return features_map


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
    print('Apply math to features...')
    print('  - Current features_map_list len = {}'.format(len(features_map_list)))
    print('  - Current N_maps = {}'.format(count_maps(features_map_list)))
    feature_array = []; 
    if 'math_ab' in kwargs:
        math_ab = kwargs['math_ab'] 
    else: 
        math_ab = [1, 2, 'divide']
    if 'plot_interp' in kwargs:
        plot_interp = kwargs['plot_interp']
        if plot_interp[0] is None:
            plot_interp[0] = 'linear' 
    else:
        plot_interp = ['linear', 1] 
         
    ## Get all the maps into one 2D array, feature_array
    features_map = extract_maps(features_map_list)
    feature_array = features_map['features']
    tags = features_map['tags']

    feature_a = np.asarray(feature_array[math_ab[0]])
    feature_b = np.asarray(feature_array[math_ab[1]])
    if math_ab[2] == 'divide':
        feature_c = feature_a / feature_b
    elif math_ab[2] == 'substract':
        feature_c = feature_a - feature_b   
    elif math_ab[2] == 'multiply':
        feature_c = feature_a * feature_b
    elif math_ab[2] == 'correlation': ## change
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
    features_map_list[idx+1].update(ids=[math_id])
    temp = '({})({}){}, {}'.format(math_ab[0],math_ab[1],math_ab[2], tags[math_ab[0]])
    features_map_list[idx+1].update(tags=[temp])
    
    print('  - Current features_map_list len = {}'.format(len(features_map_list)))
    print('  - Current N_maps = {}'.format(count_maps(features_map_list)))
    
    return feature_c

# =============================================================================
# Plot peaks
# =============================================================================
def plot_peaks(line, N_peaks_find, verbose):
    plt.plot(line.x, line.y); 
    fit_prom = 0.01
    peaks, _ = find_peaks(line.y, height=0, width=2, prominence=(fit_prom, None))
    while len(peaks)>N_peaks_find:
        #print('  N_peaks = {}, increase fit_prom to reduce N_peaks'.format(len(peaks)))
        fit_prom = fit_prom*1.1
        peaks, _ = find_peaks(line.y, height=0, width=2, prominence=(fit_prom, None)) 
    print('Peaks found at {}'.format(np.round(line.x[peaks],3)) +' for fit_prom {:.2f}'.format(fit_prom))
    
    ylim = [np.nanmin(line.y[line.y != -np.inf]), np.nanmax(line.y)]
    yrange = ylim[1]-ylim[0]
    for idx, peak in enumerate(peaks):
        plt.plot([line.x[peak], line.x[peak]], ylim, '--', color=rand_color(0.3, 0.9))
        if verbose and idx<15:
            plt.text(line.x[peak], ylim[0]+idx*yrange*0.08, str(np.round(line.x[peak],3)),fontweight='bold')
    plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3) 
    
    return line.x[peak]
        
# =============================================================================
# Count # of feature maps
# =============================================================================
def count_maps(features_map_list):
    N_maps = 0
    for ii, feature_map in enumerate(features_map_list): 
        N_maps += len(feature_map['features'])
        
    return N_maps


# =============================================================================
# Assume N symmetry
# =============================================================================
def fold_line(x, y, N, verbose):     
    step = 0.05
    period = 360/N/step
    x_fine = np.arange(np.min(x), np.max(x)+step, step)
    y_fine = np.interp(x_fine, x, y)  
    length = len(x_fine)
   
    idx0 = get_target_idx(x_fine, 0)
    idx1 = get_target_idx(x_fine, 360/N)
    
    #print(idx1-idx0+1)
    x_fold = x_fine[idx0:idx1+1]
    y_fold = np.zeros([len(x_fold),4])
    
    for ii in np.arange(idx0,idx1+1):
        val = []
        cen = get_target_idx(x_fine, x_fine[ii])
        for nn in np.arange(-N/2,N/2):
            #idx = get_target_idx(x_fine, x_fine[ii]+nn*360/N)
            idx = int(cen + nn*period)
            if idx<length: val.append(y_fine[idx])
        jj = ii-idx0
        y_fold[jj,0] = np.min(val)
        y_fold[jj,1] = np.mean(val)
        y_fold[jj,2] = np.max(val)
        y_fold[jj,3] = np.var(val)
    
    if verbose>1:
        plt.figure(24); plt.clf()
        for jj in [0,1,2,3]:
            ax = plt.subplot2grid((4,1), (jj,0))
            plt.plot(x_fold, y_fold[:,jj])
            plt.grid()
    
    return x_fold, y_fold



# =============================================================================
# Generate a random color, each channel with range (a,b), 0 dark
# =============================================================================
def rand_color(a, b):
    r = b-a
    color = (random.random()*r+a, random.random()*r+a, random.random()*r+a)
    return color


