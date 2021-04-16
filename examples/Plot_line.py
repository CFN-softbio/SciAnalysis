#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, io, random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import find_peaks

fnlist = [] 
#dn = '/home/etsai/BNL/Users/CMS/RVerduzco_2020_3/waxs/analysis/circular_average/'
#fn = 'DM_dcm-1-134-3_pos1_th0.100_x0.000_10.00s_130151_waxs_stitched.dat'
dn = '/home/etsai/BNL/Users/CMS/ETsai_shell/2021C1/waxs/analysis/circular_average_2/'
fnlist.append('FN_B2_run1bdry_x-0.164_y0.600_RH3.522_2.00s_172255_waxs_q2I.dat')

'''
dn = '/home/etsai/BNL/Users/CMS/2019C3/ETsai2/waxs/analysis/circular_average/' 
fnlist.append('Eli_sample_Co_Si_13.5kev_th0.100_x3.500_15.00s_2666489_waxs.dat')
fnlist.append('Eli_sample_Co_Si_13.5kev_th0.100_x7.500_15.00s_2666497_waxs.dat')
'''
  
flag_log10 = 1
if 0:
    wavelength = 0.9184
else:
    wavelength = 1.5418 

plt.figure(1); plt.clf()
for idx, fn in enumerate(fnlist):
    x  = np.loadtxt(dn+fn)    
    q = x[:,0]
    intensity = x[:,2]
    if flag_log10 ==1:
        intensity = np.log10(intensity)
    twotheta = np.arcsin(q*wavelength/4/np.pi)/np.pi*180*2
    
    color = rand_color(0.3, 0.6)
    plt.subplot(2,1,1); 
    plt.plot(q, intensity, label=fn); plt.grid('on')
    plt.xlabel('q (A^-1)',fontweight='bold',fontsize=10)

    plt.legend()    
    #plt.ylabel('log10',fontweight='bold')
    plt.title(dn,fontweight='bold')
    
    fit_prom = 0.0001
    peaks, _ = find_peaks(intensity, height=0.6, width=4.5, prominence=(fit_prom, None))
    ylim = [np.nanmin(intensity[intensity != -np.inf]), np.nanmax(intensity)]
    yrange = ylim[1]-ylim[0]
    for idx, peak in enumerate(peaks):
        plt.plot([q[peak], q[peak]], ylim, '--', color= color) #rand_color(0.3, 0.9))
        plt.text(q[peak], ylim[0]+idx*yrange*0.06, str(np.round(q[peak],2)),fontweight='bold', color=color)
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3) 
        
    
    plt.subplot(3,1,3); 
    plt.plot(twotheta, intensity); plt.grid('on')
    plt.xlabel('twotheta, lambda={} A'.format(wavelength),fontweight='bold')
    
    for idx, peak in enumerate(peaks):
        plt.plot([twotheta[peak], twotheta[peak]], ylim, '--', color=color)
        plt.text(twotheta[peak], ylim[0]+idx*yrange*0.06, str(np.round(twotheta[peak],2)),fontweight='bold', fontsize=8, color=color)
        plt.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3) 

# =============================================================================
# Generate a random color, each channel with range (a,b), 0 dark
# =============================================================================
def rand_color(a, b):
    r = b-a
    color = (random.random()*r+a, random.random()*r+a, random.random()*r+a)
    return color


    
