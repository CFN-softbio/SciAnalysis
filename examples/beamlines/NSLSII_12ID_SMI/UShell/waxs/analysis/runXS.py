#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols



# Experimental parameters
########################################


mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'

# WAXS detector on SMI
from SciAnalysis.XSAnalysis.DataGonio import *
calibration = CalibrationGonio(wavelength_A=0.619920987) # 20.0 keV
calibration.set_image_size(195, height=1475) # Pilatus300kW vertical
#calibration.set_image_size(1475, height=195) # Pilatus300kW horizontal
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_beam_position(97.0, 1314.0)
calibration.set_distance(0.275)

calibration.set_angles(det_phi_g=0., det_theta_g=0.)
print('ratio Dw = {:.3f}'.format(calibration.get_ratioDw()))

mask = Mask(mask_dir+'Pilatus300kWh_main_gaps-mask.png')
#mask.load('./Pilatus300kWh_current-mask.png')




# Files to analyze
########################################
source_dir = '../'
output_dir = './'


if True:
    # Working with raw data
    pattern = 'NY_sample5_*'
    #pattern = 'NY_sample5_0*'
    #pattern = 'NY_zero_angle_-4.649_001'
    #pattern = '*'
    infiles = glob.glob(os.path.join(source_dir, pattern+'.tif'))
    
else:
    # Working with merged data
    #source_dir = './sum_images/'
    source_dir = './merge_images/'
    pattern = '*'
    infiles = glob.glob(os.path.join(source_dir, pattern+'.npy'))
    mask.data = np.rot90(mask.data) # rotate CCW

infiles.sort()



def get_phi(filename, phi_offset=4.649, phi_start=1.0, phi_spacing=4.0, polarity=-1):
    
    pattern_re='^.+\/?([a-zA-Z0-9_]+_)(\d\d\d)(\.tif)$'
    phi_re = re.compile(pattern_re)
    phi_offset = 4.649
    m = phi_re.match(filename)
    
    if m:
        idx = float(m.groups()[1])
        phi_c = polarity*( phi_offset + phi_start + (idx-1)*phi_spacing )
        
    else:
        print("ERROR: File {} doesn't match phi_re".format(filename))
        phi_c = 0.0
        
    return phi_c
            

# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             'rot' : True,
             #'flip' : True,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)


    

if False:
    protocols = [ q_image_phi(calibration=calibration, blur=None, bins_relative=0.5, plot_range=[0.0, 10.0, -1.0, 6.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.01, 0.001], cmap=cmap_vge_hdr) ]


if True:
    phis = [get_phi(infile) for infile in infiles]

    #protocols = [ Protocols.sum_images(pattern_re='^.+\/([a-zA-Z0-9_]+_)(\d\d\d)(\.tif)$', infiles=infiles, processor=process, load_args=load_args, force=True) ]
    protocols = [ Protocols.merge_images_gonio_phi(q_range=[0, 10.0, -1.0, 6.0], dq=0.015, pattern_re='^.+\/([a-zA-Z0-9_]+_)(\d\d\d)(\.tif)$', infiles=infiles, phis=phis, processor=process, load_args=load_args, force=False) ]



# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=False)


# Loop
########################################
import time
donefiles = []
while False:

    infiles = glob.glob(os.path.join(source_dir, '*.tif'))

    for infile in infiles:
        if infile in donefiles:
            pass

        else:
            process.run([infile], protocols, output_dir=output_dir, force=False)

        donefiles.append(infile)

    time.sleep(4)






