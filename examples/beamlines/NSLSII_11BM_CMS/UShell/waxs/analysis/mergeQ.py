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
# For pilatus2M
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
    #calibration.set_beam_position(380, 1043-151) # SAXSx = -65, SAXSy = -73
    

    calibration.set_beam_position(456, 1043-392) # SAXSx = -65, SAXSy = -73
    #calibration.set_beam_position(456, 1043-362) # pos2 $379
    calibration.set_distance(0.257)

    
    mask = Mask(mask_dir+'Dectris/Pilatus1M_main_gaps-mask.png')
    #mask.load('./Pilatus1M_current-mask.png')
    





# Files to analyze
########################################

#root_dir = '/GPFS/xf11bm/Pilatus300/'
#root_dir = '/GPFS/xf11bm/Pilatus300/2016-3/CFN_aligned-BCP/'
#source_dir = os.path.join(root_dir, '')

source_dir = '../raw/'



#output_dir = os.path.join(source_dir, 'analysis/')
output_dir = './'


pattern = '*pos2*'
#pattern = 'AgBH_5m_13.5kev_*'
#pattern = 'HEA*x-32.0*10.00s*'
#pattern = 'JKS_AuPS_53k_mono_LL_UVO_th0.150_x0.000_10.00s_850784_saxs'
#pattern = 'SC_bulk4_x0.000_y0.000_1.00s_750568_saxs.tiff'
#pattern = 'AgBH_*'
#pattern = 'SImount_A1-1_cube30_2_x-0.390_y-0.400_30.00s_438588_saxs'
#pattern = 'SImount_AA4_x0.400_y-1.000_1.00s_438601_saxs'
#pattern = 'JKS_AuPS_53k*'
#pattern = 'TA*854826*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
#infiles = glob.glob(os.path.join(source_dir, 'AgBH*.tiff'))
#infiles = glob.glob(os.path.join(source_dir, 'AgBH_2m_cali_053018_5.00s_335521_saxs.tiff'))
#infiles = glob.glob(os.path.join(source_dir, 'AgBH_cali_5m_20180529_3879.1s_T-273.150C_1.00s_334237_saxs.tiff'))
#infiles = glob.glob(os.path.join(source_dir, 'AgBH_5m_th0.000_10.00s_20323_saxs.tiff'))

infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)




from SciAnalysis.XSAnalysis import Multiple
pattern_re = '^(.+_th\d\.\d\d\d_)_waxs.+\.tif$'
#protocols = [ Multiple.sum_images() ]
protocols = [ Multiple.merge_images_position(q_range=[-2.0, 2.5, -1.0, 2.5], dq=0.008, save_axes=True, save_mask=True, save_maps=True) ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
#process.run(infiles, protocols, output_dir=output_dir, force=1)
#process.run_multiple(pattern_re=pattern_re, infiles=infiles, protocols=protocols, output_dir=output_dir, force=True)

import re
rename_re = re.compile('^(.+)(_pos2)(_th\d\.\d\d\d_\d+\.\d+s)(_.+)')
allfiles = glob.glob(os.path.join(source_dir, '*.tiff'))
for infile in infiles:

    print(infile)
    
    # Find corresponding "non-pos2" file
    pos1 = None
    m = rename_re.match(infile)
    if m:
        els = m.groups()
        for curfile in allfiles:
            if els[0] in curfile and els[2] in curfile and 'pos1' in curfile:
                pos1 = curfile
        
        outname = els[0]+els[2]
        
    else:
        print('WARNING: No re match for {}'.format(infile))
        if infile=='../AgBH_5m_13.5kev_pos2_x0.200_y0.000_5.00s_2281002_waxs.tiff':
            pos1 = '../AgBH_5m_13.5kev_x0.000_y0.000_5.00s_2280695_waxs.tiff'
            outname = '../AgBH_5m_13.5kev'
            
            
    
            
    if pos1 is None:
        print('No pos1 file found for:')
        print('    {}'.format(infile))
        
    else:
            
        print('Will combine:')
        print('    {}'.format(pos1))
        print('    {}'.format(infile))
        
        process.run_multiple_all(basename=outname, infiles=[pos1,infile], protocols=protocols, output_dir=output_dir, force=True)
    
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


 
