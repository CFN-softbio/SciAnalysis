#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
import numpy as np
import re

from SciAnalysis import tools
from SciAnalysis.Data import *
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols


# Settings
########################################
verbosity = 3
pattern = 'AE_sample' # Files to consider
source_dir = '../' # The location of the SciAnalysis outputs
output_dir = './{}/'.format(pattern)
tools.make_dir(output_dir)

 


# Extract results from xml files
########################################
from SciAnalysis.Result import * # Results() object
def extract_results(infiles, extractions, outfile, verbosity=3):
    if verbosity>=3:
        print("Extracting results for {} infiles...".format(len(infiles)))
    
    results = Results().extract_multi_save_txt(outfile, infiles, extractions, verbosity=verbosity)
    
    return results

def load_file(infile, verbosity=3):
    if verbosity>=3:
        print(" Loading data from file: {}".format(infile))
    
    with open(infile, 'r') as fin:
        names = fin.readline().split()
        lines = fin.readlines()

    if verbosity>=4:
        print('  Saved data has {} columns:'.format(len(names)))
        print(names)
        
    return names, lines


# Extract results from XML files
results_dir = source_dir + '/results/' # Location of xml files
infiles = glob.glob(os.path.join(results_dir, '{}*.xml'.format(pattern)))
outfile = os.path.join(output_dir, '{}-extracted.txt'.format(pattern))



extractions = [ [ 'metadata_extract', ['x_position', 'y_position', 'sequence_ID', 'anneal_time'] ] ,
            ['circular_average_q2I_fit', ['fit_peaks_prefactor1', 'fit_peaks_x_center1', 'fit_peaks_sigma1', 'fit_peaks_chi_squared', 'fit_peaks_d0', 'fit_peaks_grain_size' ] ],
            ]    

results = extract_results(infiles, extractions, outfile=outfile, verbosity=verbosity)

names, lines = load_file(outfile, verbosity=verbosity)

columns = [
    'metadata_extract__x_position', 
    'metadata_extract__y_position', 
    'metadata_extract__anneal_time',
    'circular_average_q2I_fit__fit_peaks_chi_squared'
    ]

indices = [names.index(col) for col in columns]

data = []
for i, line in enumerate(lines):
    line = line.split()
    if len(line)>=len(indices):
        row = [ float(line[i]) for i in indices ]
        data.append(row)

if verbosity>=3:
    print('Extracted {} rows'.format(len(data)))    
if verbosity>=5:
    print(np.asarray(data))

outfile = os.path.join(output_dir, '{}-{:d}.npy'.format(pattern, len(data)))
np.save(outfile, data)

