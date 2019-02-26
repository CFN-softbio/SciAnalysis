#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 17:03:24 2019

@author: etsai
"""

import numpy as np
import matplotlib.pyplot as plt


if False:
    #make some sample data
    r, g = np.meshgrid(np.linspace(0,255,100),np.linspace(0,255,100))
    b=255-r
    
    #this is now an RGB array, 100x100x3 that I want to display
    rgb = np.array([r,g,b]).T
    
    color_tuple = rgb.transpose((1,0,2)).reshape((rgb.shape[0]*rgb.shape[1],rgb.shape[2]))/255.0
    
    fig = plt.figure(222); fig.clf()
    plt.pcolormesh(r, g, g*0, color=color_tuple, linewidth=4)
    
    plt.colorbar()


# https://stackoverflow.com/questions/29232439/plotting-an-irregularly-spaced-rgb-image-in-python
# https://stackoverflow.com/questions/41389335/how-to-plot-geolocated-rgb-data-faster-using-python-basemap
# https://github.com/matplotlib/matplotlib/issues/4277
# NOTE: color changes the grid color, the dimension is N-1

x, y = np.meshgrid([0, 1, 2, 3],[0, 1, 2, 3]) # 4 by 4

color_tuple = ([0, 0, 1],[0, 0, 1],[0, 0, 1],[0, 0, 1],[0, 0, 1],[0, 0, 1],[0, 0, 1],[0, 0, 1],[0, 1, 1]) # (4-1) by (4-1)

plt.figure(200); 
#m = plt.pcolormesh(x, y, (x+y)*0, facecolor=[0,0,0], edgecolors=[0,0,0], color=color_tuple, linewidth=100)







extent = (np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(x))
zz = np.asarray([(x)/x, x/3, y/3]).T
plt.imshow(zz, extent=extent)
plt.colorbar


