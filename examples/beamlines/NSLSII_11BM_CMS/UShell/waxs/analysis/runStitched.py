#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
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


# Expanded detector defined by a stitched image
offset = [0, 30] # offset between pos1 and pos2
calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
calibration.set_image_size(981+offset[0], height=1043+offset[1]) # Pilatus1M expanded to accomodate offset
calibration.set_pixel_size(pixel_size_um=172.0)

#calibration.set_beam_position(456, 1043-392, named='pos1') # pos1
#calibration.set_beam_position(456-offset[0], 1043-392-offset[1], named='pos2') # pos2
calibration.set_beam_position( 456+max(0,offset[0]) , 1043-392+max(0,offset[1]) )

calibration.set_distance(0.257)

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'NSLSII_11BM_CMS/Pilatus800k_stitch_x0_y30.png')
mask = Mask('./mask-stitched.png')
    
    




# Files to analyze
########################################
source_dir = './stitched/'
output_dir = './'

pattern = '*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
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
#protocols = [ Protocols.circular_average_q2I(plot_range=[0, 0.2, 0, None]) ]
#protocols = [ Protocols.linecut_angle(q0=0.01687, dq=0.00455*1.5, show_region=False) ]
#protocols = [ Protocols.q_image(blur=1.0, bins_relative=0.5, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.2, 0.01]) ]
#protocols = [ Protocols.qr_image(blur=1.0, bins_relative=0.5, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], zmin=1010., ztrim=[None, 0.01]) ]
#protocols = [ Protocols.qr_image(blur=None, bins_relative=0.8, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.38, 0.002], dezing_fill=True) ]
#protocols = [ Protocols.q_phi_image(bins_relative=0.25, plot_range=[0, 3.0, 0, +90]) ]

protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=1.369*0.25, dq=0.002, num_rings=10, ztrim=[0.2, 0.01], dpi=300) ,
    #Protocols.circular_average(ylog=False, plot_range=[0, 3.5, None, None]) ,
    Protocols.thumbnails(crop=None, blur=None, resize=1.0, cmap=cmap_vge_hdr, vmin=0, ztrim=[0.02, 0.001]) ,
    ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=False)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
#process.monitor_loop(source_dir=source_dir, pattern='*.tiff', protocols=protocols, output_dir=output_dir, force=False)
