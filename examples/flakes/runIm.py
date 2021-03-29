#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Import libraries
########################################

import os
import sys # To get commandline arguments
import glob

#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/home/qpress/current/code/SciAnalysis/main/'
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

# Background image
background_img = '../background.tiff' # An image we collected from a blank sample area
background_img = './average_image/average_median.tiff' # A computed image obtained by averaging the tile scan

# Tile-scan images
name_convention = 'ixiy'
#pattern = '*_tile__ix*_iy*.tiff' # All the images
pattern = '*_tile__ix000_iy023.tiff' # A specific image
infiles = glob.glob(os.path.join(source_dir, pattern))
infiles.sort()



# Settings
########################################

# image_contrast is used to rescale the images for better display

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
image_contrast = image_contrast2


# pixel_size_um provides the calibration for the microscope
microscopes = {
    'Nikon' : {
        '50X' : { 'scale_bar_pix': 284, 'scale_bar_um': 50.0 }, # 0.176 µm/pixel
        },
    'Witec' : {
        '5X' : { 'scale_bar_pix': 1920, 'scale_bar_um': 2150.2 }, # 1.120 µm/pixel
        '10X' : { 'scale_bar_pix': 1920, 'scale_bar_um': 1087.0 }, # 0.566 µm/pixel
        '20X' : { 'scale_bar_pix': 1920, 'scale_bar_um': 534.7 }, # 0.278 µm/pixel
        '50X' : { 'scale_bar_pix': 1920, 'scale_bar_um': 213.2 }, # 0.111 µm/pixel
        '100X' : { 'scale_bar_pix': 1920, 'scale_bar_um': 107.6 }, # 0.056 µm/pixel
        },
    }
for microscope, md in microscopes.items():
    for magnification, calibration in md.items():
        if 'pixel_size_um' not in calibration:
            calibration['pixel_size_um'] = calibration['scale_bar_um']/calibration['scale_bar_pix']
        #print('{} {}: {:.3f} µm/pixel'.format(microscope, magnification, calibration['pixel_size_um']))



# Select the appropriate microscope and magnification
pixel_size_um = microscopes['Witec']['20X']['pixel_size_um']



process = Protocols.ProcessorImRGB()

run_args = { 
    'verbosity' : 10,
    'num_jobs' : 10, # Parallel processing
    }
load_args = {
    'defer_load' : False ,
    'scale' : pixel_size_um # µm/pixel
    }





# Run analysis
########################################



if False:
    # Generate maps, estimate background
    ########################################
    load_args['defer_load'] = True
    protocols = [ 
        Multiple.tile_img(image_contrast=image_contrast, overlap=0.0, name_convention=name_convention) ,
        #Multiple.tile_svg(_subdir='../thumbnails/', overlap=0.0) ,
        #Multiple.mean_image(outname='average', file_extension='.tiff') ,
        Multiple.average_image(outname='average', file_extension='.tiff', average_type='median') ,
        ]
    process.run_multiple_all(basename='tile', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)
    




if False:
    # Analyze flakes
    ########################################
    load_args['defer_load'] = False
    protocols = [ 
                    #Protocols.thumbnails_contrast(resize=0.5, image_contrast=image_contrast, ),
                    Protocols.find_flakes(image_contrast=image_contrast, background=background_img, size_threshold=50, overlays=4, image_threshold=10),
                    Protocols.flake_images(image_contrast=image_contrast2, image_contrast2=image_contrast3),
                    ]

    print('Processing {} infiles...'.format(len(infiles)))
    process.run(infiles, protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)
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
        Multiple.histogram(min_area_pixels=80, interact=False, ztrim=[0, 0.001], name_convention=name_convention) ,
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
    
    
    



