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


calibration = Calibration(wavelength_A=0.770088) # 16.1 keV
#calibration = Calibration(wavelength_A=0.619920987) # 20.0 keV
binning = 2
calibration.set_image_size(3840/binning) # RayonixMAXS3840
calibration.set_pixel_size(width_mm=210.0)


calibration.set_beam_position( (1920-23)/binning, (1920+46)/binning)
calibration.set_distance(0.795)



mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
#mask = Mask(mask_dir+'/NSLSII_12ID_SMI/RayonixMAXS3840ccw.png')
mask = Mask(mask_dir+'/NSLSII_12ID_SMI/RayonixMAXS1920ccw.png')
#mask.load('./current-mask.png')







# Files to analyze
########################################

source_dir = '../'
output_dir = './'

infiles = glob.glob(os.path.join(source_dir, '*.tif'))
#infiles = glob.glob(os.path.join(source_dir, 'test_sample_01_th0.521_2.10s_068624_SAXS.tif'))

infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration,
             'mask' : mask,
             'rotCCW' : True,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)

# Examples:


protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.005], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    Protocols.thumbnails(crop=None, resize=0.5, blur=None, cmap=cmap_vge_hdr, ztrim=[0.1, 0.002]) ,
    ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=False)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
process.monitor_loop(source_dir=source_dir, pattern='*.tif', protocols=protocols, output_dir=output_dir, force=False)
