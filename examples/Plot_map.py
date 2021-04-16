#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, io, random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import find_peaks

############
dn = '/home/etsai/BNL/Users/CMS/ETsai_shell/2021C1/waxs/analysis/maps/'
sample = 'FN_B2_run1bdry_offline'

protocals = ['circular_average_002', 'circular_average_4-30']
namelist = ['prefactor1', 'd01', 'grain_size1', 'strain1' ]
N = len(namelist)

FS = 12
plt.rcParams.update({'font.size': FS})
plt.rc('xtick',labelsize=FS)
plt.rc('ytick',labelsize=FS)

for p, protocal in enumerate(protocals):
    plt.figure(10+p); plt.clf()
    for ii, name in enumerate(namelist):
        fn = '/{}__fit_peaks_{}/FN_B2_run1bdry-{}__fit_peaks_{}-final.npz'.format(protocal, name, protocal, name)
        
        data = np.load(dn+sample+fn)
        image = data['image'].astype(float)
        x_axis = data['x_axis'].astype(float)
        y_axis = data['y_axis'].astype(float)
        X, Y = np.meshgrid(x_axis, y_axis)
        corner = [np.nanmin(x_axis)+0.01, np.nanmax(y_axis)+0.01]
        
        plt.subplot(N,1,ii+1)
        vmin = np.nanquantile(image, 0.1)
        vmax = np.nanquantile(image, 0.9)
        if ii==0 or ii==1:
            cmap = 'plasma'
        elif ii==2 or ii==3:
            cmap = 'plasma'
        else:
            cmap = 'magma'        
            
        plt.pcolormesh(X,Y, (image), vmin=vmin, vmax=vmax, cmap=cmap, alpha = 1); 
        plt.colorbar(); plt.axis('off')
        plt.text(corner[0], corner[1], '{}'.format(name), fontsize=10, fontweight='bold', color='k')
        if ii==0:
            plt.title("{}\n{}".format(sample, protocal))
        elif ii==N-1:
            plt.axis('on')
            plt.xlabel('x (mm)',fontsize=FS)


'''
fn = 'FN_B2_run1bdry-extracted.txt'

datalist = np.loadtxt(dn+fn, skiprows=1, dtype='str') 

strain = []
for data in datalist:
   name = data[0]
   x = float(data[1])
   y = float(data[2])
   seq = int(float(data[3]))
   
   p = float(data[4])
   d = float(data[5])
   grain = float(data[6])
   strain = float(data[7])
'''   





