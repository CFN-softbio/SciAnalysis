#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols


# Files to analyze
########################################
source_dir = './merge_images/'
output_dir = './'

pattern = '*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.npy'))
infiles.sort()


# Experimental parameters
########################################
from SciAnalysis.XSAnalysis.DataQ import *

axes_file = glob.glob(os.path.join(source_dir, pattern+'axes.npz'))
calibration = CalibrationQ(axes_file[0])

mask_file = glob.glob(os.path.join(source_dir, pattern+'mask.png'))
mask = Mask(mask_file[0])


            

# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.005], zmin=0) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 5.5, None, None]) ,
    Protocols.q_image(name='q_images_merged', cmap=cmap_vge, plot_range=[-0.1, 3, 0, 3], xticks=[0, 2, 4, 6], yticks=[0, 2, 4], plot_buffers=[0.2, 0.05, 0.15, 0.05], ztrim=[0.3, 0.001])
    Protocols.thumbnails(name='thumbnails_merged', crop=None, resize=1, blur=None, cmap=cmap_vge_hdr, ztrim=[0.1, 0.002])
    ]

# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=1)


# Loop
########################################
#process.monitor_loop(source_dir=source_dir, pattern='*.tif', protocols=protocols, output_dir=output_dir, force=False)
 
