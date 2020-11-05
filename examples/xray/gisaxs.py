#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
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
calibration.set_image_size(1475, height=1679) # Pilatus2M
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_beam_position(780.0, 1680-604) # SAXSx -60, SAXSy -71

calibration.set_distance(5.038) # 5m

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
#mask.load('./Pilatus2M_current-mask.png')



# Files to analyze
########################################
source_dir = './'
output_dir = './'

pattern = 'BCP*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 4,
            #'save_results' : ['xml', 'plots', 'txt', 'npz', 'hdf5'],
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


# Examples:
#protocols = [ Protocols.q_image(q_max=0.14, blur=2.0, bins_relative=0.25, xticks=[-.1, 0, .1], ztrim=[0.01, 0.001])]
#protocols = [ Protocols.linecut_angle(q0=0.01687, dq=0.00455*1.5, show_region=False) ]

q0, qw = 0.0259, 0.0015
#q0, qw = 0.0313, 0.0008
protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.05], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    #Protocols.thumbnails(name='thumbnails_hdr_goldish', cmap=cmap_hdr_goldish,  blur=None, crop=None, resize=1),
    #Protocols.thumbnails(name='thumbnails_jet', crop=None,  blur=None, resize=1) , 
    #Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge, ztrim=[0.01, 0.001]) ,
    #Protocols.q_image(blur=None, cmap='jet', plot_range=[-0.05, 0.05, 0, 0.1], ztrim=[0.05, 0.002], colorbar=False, plot_buffers=[0.25, 0.06, 0.2, 0.05], xticks=[-0.04, 0, 0.04]),
    Protocols.linecut_qz_fit(qr=0.0216, dq=0.0027, show_region=False, label_filename=False, trim_range=[0, 0.1], fit_range=[q0-qw, q0+qw], plot_range=[0, 0.06, 0, None], critical_angle_substrate=0.132, critical_angle_film=0.094, ),
    ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=True)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
#process.monitor_loop(source_dir=source_dir, pattern='*.tiff', protocols=protocols, output_dir=output_dir, force=False)
