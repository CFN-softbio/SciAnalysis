#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/etsai/BNL/Users/software/SciAnalysis/'
SciAnalysis_PATH='/nsls2/data/cms/legacy/xf11bm/software/SciAnalysis/'
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
    calibration.set_distance(3.0)


    
    mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
    mask.load('./Pilatus2M_current-mask.png')
    
    
    #mask = Mask(mask_dir+'Pilatus2M_generic-mask.png')
    #mask.load('./Pilatus2M_current-mask_post-July5-1930.png')
    #mask.load('./Pilatus2M_current-mask_SAXSy_5.png')
    

else:
    # Pilatus800 WAXS detector on CMS
    calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
    calibration.set_image_size(981, height=1043) # Pilatus800k
    calibration.set_pixel_size(pixel_size_um=172.0)    
    
    
    # calibration.set_beam_position(460, 1043-404) # default
    # calibration.set_beam_position(430, 1043-404, named='inner_normal') # pos1
    # calibration.set_beam_position(254, 1043-288, named='outer_normal') # pos2

    # calibration.set_beam_position(445, 1043-404) # default
    # calibration.set_beam_position(445, 1043-404, named='inner-normal') # pos1
    # calibration.set_beam_position(269, 1043-228, named='outer_normal') # pos2


    calibration.set_beam_position(475, 1043-402) # default
    calibration.set_beam_position(475, 1043-402, named='inner_normal') # pos1
    calibration.set_beam_position(201, 1043-267, named='outer_normal') # pos2   # this is for leftover samples
    # calibration.set_beam_position(258, 1043-273, named='outer_normal') # pos2 # this is for Jom and Zhenxing's samples
    # calibration.set_beam_position(201, 1043-274, named='outer_normal') # pos2   # this is for Qin

    calibration.set_distance(0.260)

    
   
    # mask = Mask(mask_dir+'Dectris/Pilatus800k_gaps-mask.png')
    mask = Mask(mask_dir+'Dectris/Pilatus800k_vertical_gaps-mask.png')
    # mask.load('./Pilatus800k_current-mask.png') 
    





# Files to analyze
########################################
# source_dir = '../raw/'
source_dir = '../stitched/DKL_DKL/th0.05/Temperature/'
output_dir = './'

# pattern = 'AgBH_*outer*'
# pattern = 'Qin_*_outer_normal*' # We search for pos2 images since they are created last (corresponding pos1 should exist)
# pattern = '*leftover*_outer_normal*' # We search for pos2 images since they are created last (corresponding pos1 should exist)
pattern = '*'

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
protocols = [ Multiple.stitch_images_position(positions=['inner_normal', 'outer_normal']) ]
    
    
    
# Log un-paired files
def log_file(filename, textstring):
    file1 = open(filename,'a')
    file1.write('{}\n'.format(textstring))
    file1.close() 
    
log_filename = output_dir+'stitch_skipped.txt'



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
# rename_re = re.compile('^(.+)(_outer_normal_)(.+_)(\d+_waxs_stitched\.tiff)$')
# rename_re = re.compile('^(.+)(_outer_normal_)(.+_)(\d+_waxs_stitched\.tiff)$')
# rename_re = re.compile('^(.+)(_outer_normal_pos1_T)(\d+.+\d)(C.+_)(\d+_waxs_stitched\.tiff)$')
## allfiles = glob.glob(os.path.join(source_dir, '*T26.549*th0.5*.tiff'))
# infiles = glob.glob(os.path.join(source_dir, 'DKL_DKL-T39*1683476*.tiff'))
# allsearchfiles = glob.glob(os.path.join(source_dir, 'DKL_DKL-T39*_inner_normal*.tiff'))


### Input outer files
rename_re = re.compile('^(.+)(_outer_normal_pos1_T)(\d+.+\d)(C.+_th)(\d+.+\d)(_.+_)(\d+_waxs_stitched\.tiff)$')
infiles = glob.glob(os.path.join(source_dir, '*outer_normal*C_*.tiff'))
# ../stitched/DKL_DKL/th0.05/DKL_DKL-T39_outer_normal_pos1_T150.006C_x0.000_th0.050_10.00s_1683476_waxs_stitched.tiff
print(len(infiles))
### For each outer, search the following inner files to matchs
allsearchfiles = glob.glob(os.path.join(source_dir, '*_inner_normal_pos1_T*.tiff'))
rename_re_inner = re.compile('^(.+)(_inner_normal_pos1_T)(\d+.+\d)(C.+_th)(\d+.+\d)(_.+_)(\d+_waxs_stitched\.tiff)$')


for infile in infiles:

    print(infile)
    
    # Find corresponding "non-pos2" file
    pos1 = None
    m = rename_re.match(infile)
    if m:
        els = m.groups()
        sample_name = els[0]
        T = float(els[2])
        th = float(els[4])
        search_for = '{}{}'.format(sample_name, '_inner_normal_')
        searchfiles = [s for s in allsearchfiles if search_for in s]
        #search_for = '{}{}{}'.format(els[0], '_inner_normal_pos1_T', els[2][0:3])
        try: 
            #print(search_for)
            for ii, searchfile in enumerate(searchfiles):
                #print('[{}] {}'.format(ii, searchfile))
                searchm = rename_re_inner.match(searchfile)
                searchels = searchm.groups()
                searchT = float(searchels[2]) ##Change for T or theta etc
                #print(searchT)
                if abs(T - searchT) < 1.0:
                    #print('##{} found'.format(T))
                    # print(searchfile)
                    pos1 = searchfile
                    outname = searchfile[1:-5] # exclude .tiff 

            print('outname = {}'.format(outname))
        # try:           
        #     pos1 = [s for s in allfiles if search_for in s][0]
            # outname = pos1[1:-6] # exclude .tiff 
        except:
            print('FAILED {}'.format(infile))
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

        process.run_multiple_all(basename=outname, infiles=[pos1,infile], protocols=protocols, output_dir=output_dir, force=1)

        
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


