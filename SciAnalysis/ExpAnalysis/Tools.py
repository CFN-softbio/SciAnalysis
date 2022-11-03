import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
import glob, os, time
import databroker
import imageio

#from SciAnalysis import tools
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols
#from SciAnalysis.Result import *

### Temporary place for putting together various plotting functions

from scipy import interpolate
from scipy.interpolate import griddata

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

