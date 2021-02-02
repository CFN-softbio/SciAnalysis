#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH='/GPFS/xf12id1/analysis/CFN/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols




# Experimental parameters
########################################

# WAXS detector on SMI
from SciAnalysis.XSAnalysis.DataGonio import *
calibration = CalibrationGonio(wavelength_A=0.770088) # 16.1 keV
#calibration = CalibrationGonio(wavelength_A=0.619920987) # 20.0 keV
calibration.set_image_size(195, height=1475) # Pilatus300kW vertical
#calibration.set_image_size(1475, height=195) # Pilatus300kW horizontal
calibration.set_pixel_size(pixel_size_um=172.0)


#calibration.set_beam_position(97.0, 1414.0)
#calibration.set_distance(0.275)
calibration.set_beam_position(96.0, 1388.0)
calibration.set_distance(0.273900)

calibration.set_angles(det_phi_g=0, det_theta_g=0)
print('ratio Dw = {:.3f}'.format(calibration.get_ratioDw()))

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus300kWv_main_gaps-mask.png')
mask.load('./mask.png')




# Files to analyze
########################################
source_dir = '../'
output_dir = './'


pattern = '*'
#pattern = 'ET_bar1_C1a_50nm_sam1_16100.0eV_0.08deg_*'
infiles = glob.glob(os.path.join(source_dir, pattern+'.tif'))
    
infiles.sort()




# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             'rotCCW' : True,
             #'flip' : True,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


from SciAnalysis.XSAnalysis import Multiple
pattern_re = '^(.+deg)_waxs.+\.tif$'
#protocols = [ Multiple.sum_images() ]
protocols = [ Multiple.merge_images_gonio_phi(q_range=[-.10, 7.5, -1.0, 5.5], dq=0.015, phi_re='.+_waxs(-?\d+\.\d+)_.+', save_axes=False, save_mask=False, save_maps=False) ]





# Run
########################################
process.run_multiple(pattern_re=pattern_re, infiles=infiles, protocols=protocols, output_dir=output_dir, force=True)


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






