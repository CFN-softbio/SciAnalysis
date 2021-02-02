#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys # To get commandline arguments
import glob

SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)
from SciAnalysis import tools
#from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.ImAnalysis.Data import *
from SciAnalysis.ImAnalysis import Protocols




#L0 = 27 # nm 'C48'
#L0 = 32 # nm 'C48'
#L0 = 44 # nm 'C67'
#L0 = 43.6 # nm 'C99'
#L0 = 29 # nm 'L36'
#L0 = 51 # nm 'L104'
#L0 = 38.5 # nm 'L75'
#L0 = 76 # nm 'C177'
#L0 = 79 # nm 'O184'
#L0 = 65 # nm 'L176'
#L0 = 128 # nm 'L570'
#L0 = 32 # nm 'SEO30'
#L0 = 30 # nm 'S2VP45'


# layering distance
L0 = 44 # nm 'C67'
q0 = 2*np.pi/(L0)

# cyl-cyl distance
d_cc = L0/(np.sqrt(3.0)/2.0)




process = Protocols.ProcessorIm()

run_args = { 'verbosity' : 3,
                'local_partition_image_size' : 30, # pixels
                'local_partition_step' : 5.0, # relative to image_size
                #'blur' : 1.0, # pixels
                'q0' : q0, # nm^-1
                'dq' : q0*0.6, # nm^-1
                'symmetry' : 2,
                #'rcParams': {'axes.labelsize': 45,},
                'NN_cutoff_distance_nm' : L0*1.36,
                'area_min' : L0*0.35,
                #'area_max' : L0*0.35*20, # Not yet implemented
                }

load_args = { 'format' : 'custom',
                #'scale' : 1000.0/1008, # nm/pixel
                'crop_edges' : [0, 0, 248, 0],
                'load_param_file' : True,
                }
            
            
protocols = [ 
                Protocols.fft(blur=0.6, **run_args),
                #Protocols.grain_size_hex(**run_args),
                #Protocols.grain_size(**run_args),
                #Protocols.thumbnails(resize=0.5, crop=0.5),
                ]




root_dir = './'


if False:
    # Single directory
    pattern = '*'
    data_dir = './Dow 30/'
    #load_args['scale'] = 1000.0/1008 # nm/pixel

    source_dir = os.path.join( root_dir, '', data_dir )
    output_dir = os.path.join( root_dir, '', data_dir, 'analysis' )
    tools.make_dir(output_dir)

    infiles = glob.glob(source_dir + '/'+pattern+'.jpg')
    infiles.sort()
    
    print('{} infiles'.format(len(infiles)))
    process.run(infiles, protocols, output_dir=output_dir, force=True, load_args=load_args, run_args=run_args)
    
    
    


if True:
    # Walk directories
    for cur_root, subdirs, files in os.walk(root_dir):
        
        if '/analysis/' in cur_root:
            # Don't analyze the analysis sub-folders!
            pass
        else:
        
            print('Considering: {} ({:d} subdirs, {:d} files)'.format(cur_root, len(subdirs), len(files)))
            infiles = glob.glob(cur_root + '/*.jpg')
            
            if len(infiles)>0:

                source_dir = cur_root
                output_dir = os.path.join( cur_root, 'analysis' )

                tools.make_dir(output_dir)

                infiles.sort()

                process.run(infiles, protocols, output_dir=output_dir, force=True, load_args=load_args, run_args=run_args)

