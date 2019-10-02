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

# Background
infiles = [ 'x003_y010', 'x004_y010' ]
infiles = [ '{}tile_{}.tif'.format(source_dir, infile) for infile in infiles ]


pattern = 'tile*'
#pattern = 'tile_x001_y005'
#pattern = 'tile_x006_y010'
#pattern = 'tile_x001_y00*'
#pattern = 'tile_x00*_y00*'
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


if False:
    # Combine images into maps
    ########################################
    load_args['defer_load'] = True
    protocols = [ 
        Multiple.tile_img(image_contrast=image_contrast, overlap=0.0) ,
        Multiple.tile_svg(_subdir='../thumbnails/', overlap=0.0) ,
        Multiple.average_image(basename='tile', file_extension='.tiff') ,
        ]
    process.run_multiple_all(basename='tile', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)            


if False:
    # Analyze flakes
    ########################################
    load_args['defer_load'] = False
    protocols = [ 
                    #Protocols.thumbnails_contrast(resize=0.5, image_contrast=image_contrast, ),
                    #Protocols.find_flakes(image_contrast=image_contrast, background='../background.tif', size_threshold=50, overlays=4),
                    Protocols.flake_images(image_contrast=image_contrast3, image_contrast2=image_contrastB),
                    ]

    print('Processing {} infiles...'.format(len(infiles)))
    process.run(infiles, protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=False)
    #process.run_parallel(infiles, protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=False)



if False:
    # Maps of tagged images
    ########################################
    infiles = glob.glob(output_dir+'find_flakes/tile_x???_y???.png')
    load_args['defer_load'] = True
    protocols = [ 
        Multiple.tile_img(image_contrast=image_contrast1, overlap=0.0) ,
        Multiple.tile_svg(subdir='../find_flakes/', subdir_ext='.png', overlap=0.0) ,
        ]
    process.run_multiple_all(basename='tile-flakes', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)    


if False:
    # Combine results
    ########################################
    infiles = glob.glob(output_dir+'find_flakes/tile_x???_y???.pkl')
    load_args['defer_load'] = True
    protocols = [ 
        Multiple.histogram(min_area_pixels=80, interact=True) ,
        ]
    process.run_multiple_all(basename='aggregate', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)    





if False:
    # Machine-learning classification
    ########################################
    infiles = glob.glob(output_dir+'find_flakes/tile_x???_y???.pkl')
    load_args['defer_load'] = True
    
    from SciAnalysis.ImAnalysis.Flakes import cluster
    protocols = [ 
        #cluster.cluster(image_contrast=image_contrast, overlays=3) ,
        cluster.select_flakes(image_contrast=image_contrast, overlays=3) ,
        #cluster.classify() , # TODO
        ]
    process.run_multiple_all(basename='cluster', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)    
    
    
    



