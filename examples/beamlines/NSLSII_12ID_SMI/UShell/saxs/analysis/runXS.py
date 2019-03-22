#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH='/GPFS/xf12id1/analysis/CFN/SciAnalysis/SciAnalysis2018C3/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols





# Experimental parameters
########################################



calibration = Calibration(wavelength_A=0.770088) # 16.1 keV
calibration.set_image_size(981, height=1043) # Pilatus1M
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_distance(5.300)
calibration.set_beam_position(457.0, 569.0)


mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus1M_main_gaps-mask.png')
mask.load('./Pilatus1M_current-mask.png')




# Files to analyze
########################################

source_dir = '../'
output_dir = './'

infiles = glob.glob(os.path.join(source_dir, '*.tif'))
infiles = glob.glob(os.path.join(source_dir, 'test_sample_01_th0.521_2.10s_068624_SAXS.tif'))

infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration,
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)

# Examples:


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
# This code is typically only used at the beamline (it loops forever, watching for new files).
import time
donefiles = []
while True:

    infiles = glob.glob(os.path.join(source_dir, '*.tif'))

    for infile in infiles:
        if infile in donefiles:
            pass

        else:
            process.run([infile], protocols, output_dir=output_dir, force=False)

        donefiles.append(infile)

    time.sleep(4)






