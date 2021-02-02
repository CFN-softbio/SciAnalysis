#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
SciAnalysis_PATH='/home/yager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols



# Define some custom analysis routines
########################################
# TBD



# Experimental parameters
########################################

calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
#calibration = Calibration(wavelength_A=0.8856) # 14.0 keV
calibration.set_image_size(1024, height=1024) # MarCCD
calibration.set_pixel_size(pixel_size_um=161.0)
calibration.set_beam_position(507.0, 653.0)

detector_width = 165 # mm = 1024*0.161
Dw = 32.40 # Dw = Distance/detector_width
distance = (Dw*detector_width)/1000
print('Detector distance = {:.3f} m'.format(distance))
calibration.set_distance(distance) # ~5m

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'CCD/MarCCD1024-mask.png')
mask.load('./MarCCD1024-mask.png')


if False:
    # Pilatus 1M
    detector_width = 168.732 # mm = 981*0.172
    calibration.set_distance(Dw*detector_width/1000)
    calibration.set_image_size(981, height=1043) # Pilatus1M
    calibration.set_pixel_size(pixel_size_um=172.0) # Pilatus1M
    mask = Mask(mask_dir+'Dectris/Pilatus1M_main_gaps-mask.png')
    


# Files to analyze
########################################
source_dir = '../'
output_dir = './'



pattern = 'box1_5_saumil_th015_phi000.00_spot1_10sec_SAXS'
pattern = '*th015*_SAXS'
#pattern = '*_SAXS'

infiles = glob.glob(os.path.join(source_dir, pattern+''))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             'format' : 'tiff',
             'full_name' : True, # Don't truncate at a decimal point
             }
run_args = { 'verbosity' : 3,
            #'save_results' : ['xml', 'plots', 'txt', 'hdf5'],
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


# Examples:
#protocols = [ Protocols.q_image(q_max=0.14, blur=2.0, bins_relative=0.25, xticks=[-.1, 0, .1], ztrim=[0.01, 0.001])]
#protocols = [ Protocols.linecut_angle(q0=0.01687, dq=0.00455*1.5, show_region=False) ]
q0, dq0, qp = 0.0136, 0.005, 0.05
protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.005], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    Protocols.q_image(plot_range=[-qp,+qp,0,2*qp], plot_buffers=[0.28,0.07,0.25,0.05]),
    Protocols.linecut_qr_fit(qz=0.033, dq=0.005, show_region=False, q0=q0, plot_range=[0, qp/2, 0, None], fit_range=[q0-dq0, q0+dq0]),
    #Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge_hdr, ztrim=[0.01, 0.001]) ,
    ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=True)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
#process.monitor_loop(source_dir=source_dir, pattern='*.tiff', protocols=protocols, output_dir=output_dir, force=False)
