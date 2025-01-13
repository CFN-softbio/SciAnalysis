#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/nsls2/data/cms/legacy/xf11bm/software/SciAnalysis/'
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

stitched = 0

calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
calibration.set_image_size(1475, height=1679) # Pilatus2M
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_beam_position(765.0, 1680-579) # SAXSx -60, SAXSy -71

calibration.set_distance(5.038) # 5m
#calibration.set_distance(2.001) # 2m

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
#mask.load('./Pilatus2M_current-mask.png')




# Files to analyze
########################################
# conda activate nsls2-analysis-2021-1.2

if stitched==0:
    source_dir = '../raw/'
    output_dir = './'
    pattern = 'Ag*10.0*'
    pattern = '*pos1*'
    pattern= '*'
else:
    source_dir = '../stitched/'
    output_dir = '../stitched_analysis/'
    pattern = '*'


infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))

infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             #'rot180' : False,
             #'flip' : True, # PSCCD
             }
#run_args = { 'verbosity' : 3,
            #'save_results' : ['xml', 'plots', 'txt', 'hdf5'],
#            }
run_args = { 'verbosity' : 3,
            'rcParams': {'axes.labelsize': 25,
                            'xtick.labelsize': 18,
                            'ytick.labelsize': 18,
                            'xtick.major.pad': 5,
                            'ytick.major.pad': 5,
                            },
            }


process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)
#process.connect_databroker('cms') # Access databroker metadata

# Examples:
# Protocols.circular_average_q2I(plot_range=[0, 0.2, 0, None])
# Protocols.sector_average(angle=-70, dangle=25, show_region=False) 
# Protocols.linecut_q(chi0= 90+70, dq= .5, gridlines=True, label_filename=True, save_results = [ 'hdf5' ] )
# Protocols.linecut_angle(q0=0.01687, dq=0.00455*1.5, show_region=False)
# Protocols.q_image(blur=1.0, bins_relative=0.5, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.2, 0.01])
# Protocols.qr_image(blur=1.0, bins_relative=0.5, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], zmin=1010., ztrim=[None, 0.01])
# Protocols.qr_image(blur=None, bins_relative=0.8, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.38, 0.002], dezing_fill=True)
# Protocols.qr_image(blur=None, colorbar=True, save_data=False, transparent=False, label_filename=True)
# Protocols.q_phi_image(bins_relative=0.25, plot_range=[0, 3.0, 0, +90])


# Run
########################################
def main ():
    def get_protocols(protocol):
        if protocol == "qr_image":
            protocols = [
                Protocols.qr_image(blur=None, colorbar=True, save_data=False, transparent=False, label_filename=True, plot_buffers = [0.2, 0.1, 0.2, 0.1]),
            ]
        elif protocol == "q_image":
            protocols = [
                Protocols.q_image(blur=None, colorbar=True, save_data=False, transparent=False, label_filename=True, plot_buffers = [0.2, 0.1, 0.2, 0.1]),
            ]
        elif protocol == "circular_average":
            protocols = [
                Protocols.circular_average(ylog=False, plot_range=[0, None, 0, None], gridlines=True, dezing=True, label_filename=True, transparent=False) ,
            ]
        else:
            protocols = [
                Protocols.thumbnails(crop=None, cmap=cmap_vge, ztrim=[0.02, 0.001]) , # Pilatus800k
            ]
                        
        return protocols
        
        
    # if len(sys.argv) < 3:
    #     print("Usage: python runXS.py <protocal> ")
        #print("nfiles: 0 for all, -1 for latest, -10 for latest 10")
    
    protocol = sys.argv[1]  
    #Nfile = int(sys.argv[2]) 
    protocols = get_protocols(protocol)

    print('Processing latest infile...')
    process.run(infiles[-1:], protocols, output_dir=output_dir, force = 0)



if __name__ == "__main__":
    main()



# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
# process.monitor_loop(source_dir=source_dir, pattern=pattern+'*.tiff', protocols=protocols, output_dir=output_dir, force=False)



