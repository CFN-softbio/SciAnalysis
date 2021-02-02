#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################


import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
#from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.ImAnalysis.Data import *
from SciAnalysis.ImAnalysis import Protocols


# Define some custom analysis routines
########################################
class particles(Protocols.particles):
    
    def preprocess_custom(self, data, **run_args):
    
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(20):
            data.blur(2)
                            
        data.maximize_intensity_spread()
        
        return data
    


        


# Files to analyze
########################################

#filename, scale = '1min_proc', 50.0/182
#filename, scale = '45min_proc', 100.0/168
filename, scale = '322min_proc', 100.0/152

source_dir = '../'
output_dir = './'
pattern = '*'
pattern = filename

infiles = glob.glob(os.path.join(source_dir, pattern+'.jpg'))
infiles.sort()



# Analysis to perform
########################################


load_args = { 'format' : 'custom',
                'scale' : scale, # nm/pixel
                'load_param_file' : False,
                }
            



run_args = { 
                'verbosity' : 4,
                'preprocess' : 'custom',
                'radius_min' : 5, # nm
                'radius_max' : 100.0, # nm
                
                }



#process = Protocols.ProcessorIm(load_args=load_args, run_args=run_args)
process = Protocols.ProcessorImRGB(load_args=load_args, run_args=run_args)
            
protocols = [ 
                #particles(q0=0.02, threshold=40),
                Protocols.particles_annotated(),
                ]




# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=True)
