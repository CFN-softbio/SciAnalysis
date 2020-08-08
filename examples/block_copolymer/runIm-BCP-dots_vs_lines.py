#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys # To get commandline arguments
import glob

#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/home/yager/current/code/SciAnalysis/main/'
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
                'q0' : q0, # nm^-1
                'dq' : q0*0.6, # nm^-1
                'NN_cutoff_distance_nm' : L0*1.36,
                'correlation_step_size_points' : 40, # Ignore image edges for grain size analysis
                'correlation_edge_exclusion' : 50, # Step for grain size analysis (speeds up code)
                'radius_min' : L0*0.08, # nm
                'dot_size_cutoff_nm' : L0*0.3, # size cutoff for distinguishing a dot vs. line
                }

load_args = { 'format' : 'custom',
                'scale' : 500.0/403, # nm/pixel
                'crop_edges' : [0, 0, 124, 0],
                'load_param_file' : False,
                }
            
            
protocols = [ 
                #Protocols.fft(blur=0.6, **run_args),
                #Protocols.thumbnails(resize=0.5, crop=0.5),
                #Protocols.particles(threshold=190, invert=False, preprocess='highloweq', **run_args),
                Protocols.dots_vs_lines(threshold=190, invert=False, **run_args),
                Protocols.grain_size_hex(name='grain_size_hex_dots', threshold=190, invert=False, symmetry=6, mask='dots', **run_args),
                Protocols.grain_size(name='grain_size_lines', symmetry=2, mask='lines', **run_args),
                Protocols.fft(name='fft_dots', mask='dots', blur=0.6, **run_args),
                Protocols.fft(name='fft_lines', mask='lines', blur=0.6, **run_args),
                ]




source_dir = '../'
output_dir = './'
pattern = '*'

infiles = glob.glob(source_dir + '/'+pattern+'.tif')
infiles.sort()

print('{} infiles'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=False, load_args=load_args, run_args=run_args)




