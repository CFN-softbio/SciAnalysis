#!/usr/bin/python3
# -*- coding: utf-8 -*-


# conda activate analysis-2019-3.0.1-smi 
# python runXS.py


# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
#SciAnalysis_PATH='/nsls2/xf12id2/analysis/CFN/SciAnalysis/'
# 2023-May
SciAnalysis_PATH = '/nsls2/data/smi/legacy/gpfs-smi/analysis/CFN/SciAnalysis/'
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
if 0:
    calibration = Calibration(wavelength_A=0.770088) # 16.1 keV
    # #calibration = Calibration(wavelength_A=0.619920987) # 20.0 keV
    calibration.set_image_size(1475, height=619) # Pilatus1M
    calibration.set_pixel_size(pixel_size_um=172.0)

    calibration.set_beam_position(220.5-1.5, 310)  
    calibration.set_distance(0.281-0.012) #sm7200

    mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
    mask = Mask('./mask_pil900kw.png')
    # mask.load('./mask_waxs_badpixel2.png')
    # mask.load('./mask_waxs_bs.png') #if waxs00

# WAXS detector on SMI
if 1:
    from SciAnalysis.XSAnalysis.DataGonioMulti import *
    calibration = CalibrationGonioMulti(wavelength_A=0.770088) # 16.1 keV
    calibration.set_image_size(1475, height=619) # Pilatus1M
    # calibration.set_image_size(619, height=1475) # Pilatus1M
    # calibration.set_image_size(195, height=1475) 
    calibration.set_pixel_size(pixel_size_um=172.0)

    # calibration.set_beam_position(96, 1388) # det_phi_g=0.0
    # calibration.set_beam_position(106-7, 1475.0-88.0) # det_phi_g=-3.5
    # calibration.set_beam_position(310-195, 1475.0-219.0) # det_phi_g=-3.5
    # calibration.set_beam_position(219, 310-195-17)  
    calibration.set_beam_position(219, 310-195-17)  
    # calibration.set_beam_position(310-195, 1475.0-219.0) # det_phi_g=-3.5
    # calibration.set_distance(0.273900)
    calibration.set_distance(0.268)

 
    calibration.set_angles(det_phi_g=0, det_theta_g=0.)

    print('ratio Dw = {:.3f}'.format(calibration.get_ratioDw()))

    mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
    #mask = Mask(mask_dir+'Dectris/Pilatus300kWv_main_gaps-mask.png')
    #mask.load(mask_dir+'NSLSII_12ID/SMI/Pilatus300kWv_SMI_badpixels-mask.png')
    #mask.load('./Pilatus300kWv_current-mask.png')
    mask = Mask('./mask_waxs_auto.png')   
    mask.load('./mask_waxs_bs.png') 




# Files to analyze
########################################
#source_dir = '/nsls2/xf12id2/data/images/users/2022_2/309389_Subramanian/900KW/'  ###CHANGE THIS
source_dir = '/nsls2/data/smi/legacy/results/data/2024_1/312283_Subramanian/900KW/'  ###CHANGE THIS
# source_dir = '/nsls2/data/smi/legacy/results/data/2024_1/310000_ETsai/900KW/'  ###CHANGE THIS
output_dir = './WAXS00/'  ###CHANGE THIS

#pattern = 'StaticT_*Ag*waxs00*084*'
# pattern = '*KaptonTape*smz7200*waxs00*'
# pattern = '*Ag*waxs00*630*'
# pattern = '*test*waxs0*'
#pattern = '*T_KL13*waxs0*9640*'
pattern = '*waxs0*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tif'))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
              'mask' : mask,
            #   'rotCCW' : True,
            # 'flip' : True,
             }
run_args = { 'verbosity' : 3,        
        'rcParams': {'axes.labelsize': 15,
                        'xtick.labelsize': 15,
                        'ytick.labelsize': 15,
                        'xtick.major.pad': 10,
                        'ytick.major.pad': 10,
                        },
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


# Examples:
#protocols = [ Protocols.q_image(q_max=0.14, blur=2.0, bins_relative=0.25, xticks=[-.1, 0, .1], ztrim=[0.01, 0.001])]
#protocols = [ Protocols.linecut_angle(q0=0.01687, dq=0.00455*1.5, show_region=False) ]

protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=1.076, dq=0.008, num_rings=10, ztrim=[0.001, 0.001],label_filename=1, transparent=0) ,
    #Protocols.calibration_check(show=0, AgBH=True, q0=0.01076, num_rings=4, ztrim=[0.005, 0.005]) ,
    Protocols.circular_average(ylog=True, plot_range=[0, 5.3, None, None], gridlines=True, label_filename=1, transparent=0, colorbar=True) ,
    Protocols.q_image(blur=None, cmap=cmap_vge, ztrim=[0.0005, 0.0005], label_filename=1, transparent=0, colorbar=1, plot_buffers = [0.12, 0.12, 0.12, 0.12], dpi=200) ,
    
    #Protocols.thumbnails(crop=0.5, cmap=cmap_vge, make_square=True, resize=0.5, blur=0.8, ztrim=[0.05, 0.004]) ,
    #Protocols.thumbnails(show=False, cmap=cmap_vge, resize=0.5, ztrim=[0.05, 0.004]) ,
    #Protocols.qr_image(show=False, blur=None, ztrim=[0.05, 0.004], label_filename=1, transparent=0, save_data=True, colorbar=True) ,
    ]
    


# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force = 1)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
#process.monitor_loop(source_dir=source_dir, pattern=pattern, protocols=protocols, output_dir=output_dir, force=0)

# /GPFS/xf12id1/analysis/2019_3/305435_Murray/





