#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/etsai/BNL/Users/software/SciAnalysis/'
SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols


# Experimental parameters
########################################

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'

if False:
    # SAXS detector on CMS
    calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
    calibration.set_image_size(1475, height=1679) # Pilatus2M
    calibration.set_pixel_size(pixel_size_um=172.0)    
    calibration.set_beam_position(754, 1679-604) # SAXSx = -65, SAXSy = -73
    calibration.set_distance(5.05)


    
    mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
    mask.load('./Pilatus2M_current-mask.png')
    
    
    #mask = Mask(mask_dir+'Pilatus2M_generic-mask.png')
    #mask.load('./Pilatus2M_current-mask_post-July5-1930.png')
    #mask.load('./Pilatus2M_current-mask_SAXSy_5.png')
    

else:
    # Pilatus800 WAXS detector on CMS
    calibration = Calibration(wavelength_A=0.9184) # 17 keV
    calibration.set_image_size(981, height=1043) # Pilatus2M
    calibration.set_pixel_size(pixel_size_um=172.0)    
    
    
    calibration.set_beam_position(456, 1043-392) # default
    calibration.set_beam_position(456, 1043-392, named='pos1') # pos1
    calibration.set_beam_position(456, 1043-362, named='pos2') # pos2
    calibration.set_distance(0.257)

    
    mask = Mask(mask_dir+'Dectris/Pilatus800k_gaps-mask.png')
    #mask.load('./Pilatus1M_current-mask.png')
    





# Files to analyze
########################################
source_dir = '../raw/'
output_dir = './'


pattern = '*_pos2_*' # We search for pos2 images since they are created last (corresponding pos1 should exist)

infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 5,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)




from SciAnalysis.XSAnalysis import Multiple
#pattern_re = '^(.+_th\d\.\d\d\d_)_waxs.+\.tif$'
#protocols = [ Multiple.sum_images() ]
protocols = [ Multiple.stitch_images_position(positions=['pos1', 'pos2']) ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
#process.run(infiles, protocols, output_dir=output_dir, force=1)
#process.run_multiple(pattern_re=pattern_re, infiles=infiles, protocols=protocols, output_dir=output_dir, force=True)


# Log un-paired files
log_filename = output_dir+'stitch_skipped.txt'
def log_file(filename, textstring):
    with open(filename, 'a') as fout:
        fout.write('{}\n'.format(textstring))


import re
rename_re = re.compile('^(.+)(_pos2_)(.+_)(\d+_waxs\.tiff)$')
allfiles = glob.glob(os.path.join(source_dir, '*.tiff'))
for infile in infiles:

    print(infile)
    
    # Find corresponding "non-pos2" file
    pos1 = None
    m = rename_re.match(infile)
    if m:
        els = m.groups()
        search_for = '{}{}{}'.format(els[0], '_pos1_', els[2])
        try:
            pos1 = [s for s in allfiles if search_for in s][0]
            #outname = pos1
            outname = pos1[1:-6] # Exclude .tiff
        except:
            log_file(log_filename, infile)
            continue
        
    else:
        print('WARNING: No re match for {}'.format(infile))
        log_file(log_filename, infile)
    
            
    if pos1 is None:
        print('No pos1 file found for:')
        print('    {}'.format(infile))
        log_file(log_filename, infile)

    else:
        print('Will combine:')
        print('    {}'.format(pos1))
        print('    {}'.format(infile))

        process.run_multiple_all(basename=outname, infiles=[pos1,infile], protocols=protocols, output_dir=output_dir, force=False)
    
    print('\n')





# Loop
########################################
# This code is typically only used at the beamline (it loops forever, watching for new files).
import time
donefiles = []
while False:

    infiles = glob.glob(os.path.join(source_dir, '*.tiff'))

    for infile in infiles:
        if infile in donefiles:
            pass

        else:
            process.run_multiple(pattern_re=pattern_re, infiles=[infile], protocols=protocols, output_dir=output_dir, force=False)

            donefiles.append(infile)

    time.sleep(4)


