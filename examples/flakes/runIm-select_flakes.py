#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys # To get commandline arguments
import glob

SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)
from SciAnalysis import tools
from SciAnalysis.ImAnalysis.Data import *
#from SciAnalysis.ImAnalysis import Protocols
from SciAnalysis.ImAnalysis.Flakes import Protocols
from SciAnalysis.ImAnalysis.Flakes import Multiple





# Files to analyze
########################################
source_dir = '../'
output_dir = './'

pattern = 'tile*'
infiles = glob.glob(os.path.join(source_dir, pattern+'.tif'))
infiles.sort()



# Bright
image_contrast1 = (0, 1)
image_contrast2 = (0.2, 0.8)
image_contrast3 = (0.3, 0.7)
image_contrast4 = (0.32, 0.6)
image_contrast5 = (0.33, 0.55)
image_contrast6 = (0.35, 0.48)

# Dark
image_contrastA = (0.45, 0.7)
image_contrastB = (0.5, 0.6)

# Using
image_contrast = image_contrastA

scale_bar_pix, scale_bar_um = 284, 50.0
pixel_size_um = scale_bar_um/scale_bar_pix # 0.176 um/pixel



process = Protocols.ProcessorImRGB()

run_args = { 
    'verbosity' : 5,
    'num_jobs' : 10, # Parallel processing
    }
load_args = {
    'defer_load' : False ,
    'scale' : pixel_size_um # um/pixel
    }


if True:
    # Classification
    ########################################
    infiles = glob.glob(output_dir+'find_flakes/tile_x???_y???.pkl')
    load_args['defer_load'] = True
    
    selection = {
        'radius_um' : [2, 45],
        'flake_contrast' : [-0.02, -0.01],
        #'gray std __relative' : [0, 0.1], # Relative standard deviation of gray channel
        #'R' : [0, 40], # Red channel (in raw units; i.e. 0-255 for this channel)
        #'B __rescaled' : [0.2, 1.5], # Blue channel (in rescaled unit where avg=0 and std=1)
        #'entropy_inner __rescaled' : [-1, 3],
        }
    
    from SciAnalysis.ImAnalysis.Flakes import cluster
    protocols = [ 
        cluster.cluster(image_contrast=image_contrast, overlays=3) ,
        cluster.select_flakes(image_contrast=image_contrast, selection=selection, overlays=3, output_all=True) ,
        ]
    process.run_multiple_all(basename='cluster', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)    
    
    
    



