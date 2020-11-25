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
#calibration = CalibrationGonio(wavelength_A=0.770088) # 16.1 keV
calibration = CalibrationGonio(wavelength_A=0.619920987) # 20.0 keV
calibration.set_image_size(195, height=1475) # Pilatus300kW vertical
#calibration.set_image_size(1475, height=195) # Pilatus300kW horizontal
calibration.set_pixel_size(pixel_size_um=172.0)

#calibration.set_beam_position(97.0, 1414.0)
#calibration.set_distance(0.275)

calibration.set_beam_position(96.0, 1388.0) # det_phi_g=0.0
#calibration.set_beam_position(106-7.0, 1475.0-88.0) # det_phi_g=-3.5
calibration.set_distance(0.273900)


calibration.set_angles(det_phi_g=0., det_theta_g=0.)
print('ratio Dw = {:.3f}'.format(calibration.get_ratioDw()))

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus300kWv_main_gaps-mask.png')
#mask.load(mask_dir+'NSLSII_12ID/SMI/Pilatus300kWv_SMI_badpixels-mask.png')
#mask.load('./Pilatus300kWv_current-mask.png')




# Files to analyze
########################################
source_dir = '../'
output_dir = './'


pattern = '*'
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

protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.005], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge_hdr, ztrim=[0.01, 0.001]) ,
    ]
    

# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=False)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
process.monitor_loop(source_dir=source_dir, pattern='*.tif', protocols=protocols, output_dir=output_dir, force=False)
