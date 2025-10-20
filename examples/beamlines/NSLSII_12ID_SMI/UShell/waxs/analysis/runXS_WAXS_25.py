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
from SciAnalysis.XSAnalysis.DataGonioMulti import *


# Define some custom analysis routines
########################################
# TBD



# Experimental parameters
########################################
# #calibration = Calibration(wavelength_A=0.770088) # 16.1 keV
# calibration = Calibration(wavelength_A=1.102044) # 11.25
# #calibration = Calibration(wavelength_A=0.619920987) # 20.0 keV
# calibration.set_image_size(1475, height=619) # Pilatus1M
# calibration.set_pixel_size(pixel_size_um=172.0)

# WAXS detector on SMI
#beam_center = [97, 1254]  for waxs=0,
#then set waxs = 25, i get  [ -655, 1254 ] 

if 1:
    #calibration = CalibrationGonioMulti(wavelength_A=0.770088) # 16.1 keV
    calibration = CalibrationGonioMulti(wavelength_A=1.102044) # 11.25
    # calibration.set_image_size(1475, height=619) # Pilatus1M
    calibration.set_image_size(619, height=1475) # Pilatus1M
    calibration.set_pixel_size(pixel_size_um=172.0)

    # calibration.set_beam_position(96, 1388) # det_phi_g=0.0
    # calibration.set_beam_position(106-7, 1475.0-88.0) # det_phi_g=-3.5
    # calibration.set_beam_position(310-195, 1475.0-219.0) # det_phi_g=-3.5
    # calibration.set_beam_position(219, 310-195-17)  
    # calibration.set_beam_position(219, 310-195-17)  

    # calibration.set_beam_position(221, 97) #-655+23, 97
    # calibration.set_beam_position(97-8, 1475-221-25)  #n45*04361*
    calibration.set_beam_position(97+5, 1475-221)
    # calibration.set_beam_position(315, 1475-221)
 
    #calibration.set_beam_position(310-195, 1475.0-219.0) # det_phi_g=-3.5
    # calibration.set_distance(0.273900)
    # calibration.set_distance(0.268)
    calibration.set_distance(0.276)


    calibration.set_angles(det_phi_g=-25, det_theta_g=0)

    print('ratio Dw = {:.3f}'.format(calibration.get_ratioDw()))

    mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
    #mask = Mask(mask_dir+'Dectris/Pilatus300kWv_main_gaps-mask.png')
    #mask.load(mask_dir+'NSLSII_12ID/SMI/Pilatus300kWv_SMI_badpixels-mask.png')
    #mask.load('./Pilatus300kWv_current-mask.png')
    mask = Mask('./mask/mask_pil900kw_v.png')   
    #mask.load('./mask_waxs_bs.png') 


# Files to analyze
########################################
#source_dir = '/nsls2/xf12id2/data/images/users/2022_2/309389_Subramanian/900KW/'  ###CHANGE THIS
source_dir = '/nsls2/data/smi/legacy/results/data/2024_1/312437_Jones/900KW/'  ###CHANGE THIS
source_dir = '/nsls2/data/smi/proposals/2025-2/pass-318110/projects/static/user_data/900KW/'

output_dir = '/nsls2/data/smi/proposals/2025-2/pass-318110/projects/insitu/analysis/WAXS25/'  

# pattern = '*n45*04361*'
pattern = '*AgBH*3597*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tif'))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
              'mask' : mask,
              'rotCCW' : True,
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
    Protocols.calibration_check(show=False, temp=True, q0=1.076, dq=0.005, num_rings=50, ztrim=[0.0001, 0.005], transparent=0) ,
    #Protocols.calibration_check(show=0, AgBH=True, q0=0.01076, num_rings=4, ztrim=[0.005, 0.005]) ,
    Protocols.circular_average(name='circular_average_', ylog=True, plot_range=[0, None, None, None], gridlines=True, label_filename=1) ,
    # Protocols.linecut_angle(name='linecut_angle_2.7', q0=2.70, dq=0.1, show_region='save'), 
    # Protocols.linecut_angle(name='linecut_angle_3.1',q0=3.10, dq=0.1, show_region='save'), 
    Protocols.linecut_angle(name='linecut_angle_wide',q0=3.0, dq=0.4, show_region='save', plot_range=[50, 180, 0, None], gridlines=True), 
    Protocols.circular_average_q2I_fit(plot_range=[2, 4, 0, None], qn_power=0.0, trim_range=[2, 4], fit_range=[2.6, 3.5], num_curves=2, q0=[2.67, 3.08], sigma=0.1, show_curves=1, label_filename=True),
    #Protocols.qr_image(blur=None, cmap=cmap_vge, ztrim=[0.0005, 0.0005], label_filename=1, transparent=0, colorbar=1, plot_buffers = [0.12, 0.12, 0.12, 0.12], dpi=200) ,
    
    #Protocols.thumbnails(crop=0.5, cmap=cmap_vge, make_square=True, resize=0.5, blur=0.8, ztrim=[0.05, 0.004]) ,
    #Protocols.thumbnails(show=False, cmap=cmap_vge, resize=0.5, ztrim=[0.05, 0.004]) ,
    # Protocols.linecut_qr(qz=-2.5, dq=0.1, ylog=True, show_region='save', gridlines=True), 
    Protocols.qr_image(show=False, blur=None, ztrim=[0.05, 0.004], save_data=True, colorbar=True, transparent=0) ,
    # Protocols.q_image(show=False, blur=None, ztrim=[0.05, 0.004], save_data=True, colorbar=True, transparent=0) ,
    ]
    


# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force = 1)


# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
# process.monitor_loop(source_dir=source_dir, pattern=pattern, protocols=protocols, output_dir=output_dir, force=0)

# /GPFS/xf12id1/analysis/2019_3/305435_Murray/





