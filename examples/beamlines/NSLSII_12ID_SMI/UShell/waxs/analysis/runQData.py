#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
#SciAnalysis_PATH='/GPFS/xf12id1/analysis/CFN/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols




# Experimental parameters
########################################
from SciAnalysis.XSAnalysis.DataQ import *
calibration = CalibrationQ(infile='q_axes.npz')
mask = Mask('q_mask.png')


# Files to analyze
########################################
source_dir = '../'
output_dir = './'


pattern = '*'
#pattern = 'ET_bar1_C1a_50nm_sam1_16100.0eV_0.08deg*'


#source_dir = './sum_images/'
source_dir = './merge_images/'
infiles = glob.glob(os.path.join(source_dir, pattern+'.npy'))
infiles.sort()




            

# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


thumb_merged = Protocols.thumbnails(crop=None, resize=1, blur=None, cmap=cmap_vge_hdr, ztrim=[0.1, 0.002])
thumb_merged.name = 'thumbnails_merged'
q_images_merged = Protocols.q_image(cmap=cmap_vge, plot_range=[-0.1, 7.0, 0, 5.5], xticks=[0, 2, 4, 6], yticks=[0, 2, 4], plot_buffers=[0.2, 0.05, 0.15, 0.05], ztrim=[0.3, 0.001])
q_images_merged.name = 'q_images_merged'

protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.005], zmin=0) ,
    Protocols.circular_average(ylog=True, plot_range=[0, 5.5, None, None]) ,
    thumb_merged,
    q_images_merged,
    ]

# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=False)


# Loop
########################################
#process.monitor_loop(source_dir=source_dir, pattern='*.tif', protocols=protocols, output_dir=output_dir, force=False)
