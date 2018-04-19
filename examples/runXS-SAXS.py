#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols





# Experimental parameters
########################################



calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
calibration.set_image_size(487, height=619) # Pilatus300k
calibration.set_pixel_size(pixel_size_um=172.0)


#cms.SAXS.setCalibration([402, 443], 5.0, [25.00, 16.00]) # 2017-04-01, 13.5 keV, 5m, GISAXS
calibration.set_beam_position(402.0, 443.0)
calibration.set_distance(5.038)

#cms.SAXS.setCalibration([286, 448], 5.038, [5.00, 17.00]) # 2017-04-01, 13.5 keV, 5m, GISAXS
calibration.set_beam_position(286.0, 448.0)


mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus300k_main_gaps-mask.png')
#mask.load('./Pilatus300k_current-mask.png')




# Files to analyze
########################################

#root_dir = '/GPFS/xf11bm/Pilatus300/'
#root_dir = '/GPFS/xf11bm/Pilatus300/2016-3/CFN_aligned-BCP/'
#source_dir = os.path.join(root_dir, '')

source_dir = '../'


#output_dir = os.path.join(source_dir, 'analysis/')
output_dir = './'

infiles = glob.glob(os.path.join(source_dir, '*.tiff'))
#infiles = glob.glob(os.path.join(source_dir, 'Ag*.tiff'))
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

# Examples:


protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.05], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge, ztrim=[0.0, 0.01]) ,
    ]
    



# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=False)


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
            process.run([infile], protocols, output_dir=output_dir, force=False)

        donefiles.append(infile)

    time.sleep(4)






