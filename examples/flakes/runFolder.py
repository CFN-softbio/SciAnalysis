#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys # To get commandline arguments
import glob



PARALLELLIZE = True
if PARALLELLIZE:
    # Set the backend to 'Agg' to avoid bugs related to parralelization (joblib)
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams['mathtext.fontset'] = 'cm'
    import pylab as plt


SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)
from SciAnalysis import tools
from SciAnalysis.ImAnalysis.Data import *
#from SciAnalysis.ImAnalysis import Protocols
from SciAnalysis.ImAnalysis.Flakes import Protocols
from SciAnalysis.ImAnalysis.Flakes import Multiple



import re
def multi_replace(string, replacements, ignore_case=False):
    """
    Given a string and a dict, replaces occurrences of the dict keys found in the 
    string, with their corresponding values. The replacements will occur in "one pass", 
    i.e. there should be no clashes.
    :param str string: string to perform replacements on
    :param dict replacements: replacement dictionary {str_to_find: str_to_replace_with}
    :param bool ignore_case: whether to ignore case when looking for matches
    :rtype: str the replaced string
    """
    # Based on:
    #  https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
    rep_sorted = sorted(replacements, key=lambda s: len(s[0]), reverse=True)
    if ignore_case:
        replacements = dict((pair[0].lower(), pair[1]) for pair in sorted(replacements.items()))
    rep_sorted = sorted(replacements, key=lambda s: (len(s), s), reverse=True)
    rep_escaped = [re.escape(replacement) for replacement in rep_sorted]
    pattern = re.compile("|".join(rep_escaped), re.I if ignore_case else 0)
    return pattern.sub(lambda match: replacements[match.group(0)], string)

def analysis(folder, step=0):

    # Files to analyze
    ########################################
    source_dir = './'+folder+'/'
    output_dir = source_dir+'analysis/'

    tools.make_dir(output_dir)


    # Directory names with '[' or ']' will confused "glob"...
    source_dir_e = multi_replace(source_dir, {'[':'[[]', ']':'[]]'})
    output_dir_e = multi_replace(output_dir, {'[':'[[]', ']':'[]]'})

    # Background
    #infiles = [ 'x003_y010', 'x004_y010' ]
    #infiles = [ '{}tile_{}.tif'.format(source_dir, infile) for infile in infiles ]

    pattern = 'tile*'
    #pattern = 'tile_x001_y005'
    #pattern = 'tile_x006_y010'
    #pattern = 'tile_x001_y00*'
    #pattern = 'tile_x00*_y00*'
    infiles = glob.glob(os.path.join(source_dir_e, pattern+'.tif'))
    #infiles = glob.glob(os.path.join(source_dir_e, pattern+'.jpg'))
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
        'verbosity' : 3,
        'num_jobs' : 5, # Parallel processing
        }
    load_args = {
        'defer_load' : False ,
        'scale' : pixel_size_um # um/pixel
        }


    if step<=1:
        # Combine images into maps
        ########################################
        load_args['defer_load'] = True
        protocols = [ 
            Multiple.tile_img(image_contrast=image_contrast, overlap=0.0) ,
            Multiple.tile_svg(subdir='../../', subdir_ext='.tif', overlap=0.0) ,
            Multiple.average_image(basename='tile', file_extension='.tiff') ,
            ]
        process.run_multiple_all(basename='tile', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=False)            


    if step<=5:
        # Analyze flakes
        ########################################
        load_args['defer_load'] = False
        protocols = [ 
                        Protocols.thumbnails_contrast(resize=0.5, image_contrast=image_contrast, ),
                        #Protocols.find_flakes(image_contrast=image_contrast, background='./background.tif', size_threshold=50, overlays=4),
                        #Protocols.flake_images(image_contrast=image_contrast3, image_contrast2=image_contrastB),
                        ]

        print('Processing {} infiles...'.format(len(infiles)))
        if PARALLELLIZE:
            process.run_parallel(infiles, protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=False)
        else:
            process.run(infiles, protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=True)



    if step<=7:
        # Maps of tagged images
        ########################################
        infiles = glob.glob(output_dir_e+'find_flakes/tile_x???_y???.png')
        load_args['defer_load'] = True
        protocols = [ 
            Multiple.tile_img(image_contrast=image_contrast1, overlap=0.0) ,
            Multiple.tile_svg(subdir='../find_flakes/', subdir_ext='.png', overlap=0.0) ,
            ]
        process.run_multiple_all(basename='tile-flakes', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=False)    


    if step<=10:
        # Combine results
        ########################################
        infiles = glob.glob(output_dir_e+'find_flakes/tile_x???_y???.pkl')
        load_args['defer_load'] = True
        protocols = [ 
            Multiple.histogram(min_area_pixels=80, interact=False) ,
            ]
        process.run_multiple_all(basename='aggregate', infiles=infiles, protocols=protocols, output_dir=output_dir, load_args=load_args, run_args=run_args, force=False)    




ignore = [ 'collate',
        'hand',
        'depth2time1N1shear[0.4,10]velocity5',
        'depth2time1N5shear0velocity5',
        'depth2time1N10shear0velocity5',
        'depth1time1N1shear[0.4,40]velocity5',
        'depth0.3time1N1shear0velocity5',
        'depth2time1N1shear0veloocity5',
        'depth2time1N1shear[0.2,10]velocity5',
        'depth1time1N1shear[0.4,10]velocity5',
        'depth1time1N1shear[0.1,10]velocity5',
        'depth0.1time1N1shear0velocity5',
        ]

for folder in glob.glob('depth*'):
#for folder in ['./depth2time1N1shear[0.4,10]velocity5']:

    
    if folder not in ignore:
        
        print('########################################')
        print('Analysis for folder:')
        print('    {}'.format(folder))
        print('########################################')
        analysis(folder, step=0)
    
