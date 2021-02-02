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

import time


# Define some custom analysis routines
########################################
from SciAnalysis.Result import * # The Results() class allows one to extract data from XML files.


# DEPRECATED: Use Results().extract_dict instead
def extractions_keys(extractions):
    keys = []
    for protocol, signals in extractions:
        for signal in signals:
            keys.append('{}__{}'.format(protocol, signal))
    return keys
def extractions_dictionary(extractions, results):
    keys = extractions_keys(extractions)
    
    results_dict = {'infile' : results[0]}
    for key, result in zip(keys, results[1:]):
        results_dict[key] = result
    
    return results_dict

def autonomous_result(xml_file, method='dict', clean=True, verbosity=3):
    
    time.sleep(1) # Kludge to avoid corrupting XML file?
    
    extractions = [ 
        [ 'metadata_extract', 
            [
            'x_position', 
            'y_position',
            #'sequence_ID',
            ] 
        ],
        ['circular_average_q2I_fit', 
            [
            'fit_peaks_prefactor1', 
            'fit_peaks_x_center1', 
            'fit_peaks_sigma1', 
            'fit_peaks_chi_squared', 
            'fit_peaks_prefactor1_error', 
            'fit_peaks_x_center1_error', 
            'fit_peaks_sigma1_error', 
            ] 
        ],
        #['linecut_angle_fit', 
            #[
            #'fit_eta_eta', 
            #'orientation_factor', 
            #'orientation_angle', 
            #'fit_eta_span_prefactor',
            #] 
        #],
        ]


    infiles = [xml_file]
    if method=='dict':
        results_dict = Results().extract_dict(infiles, extractions, verbosity=verbosity)[0]
        
    else:
        extracted_results = Results().extract_multi(infiles, extractions, verbosity=verbosity)
        
        results = extracted_results[0]
        for i, result in enumerate(results):
            if i>0: # Skip results[0], since this is just infile (string)
                results[i] = np.nan_to_num(float(result))
        
        if method=='map_to_dict':
            # DEPRECATED: Use Results().extract_dict instead
            results_dict = extractions_dictionary(extractions, results)
            
        else:
            # DEPRECATED: This code relies on manually assigning keys and results
            results_dict = {
                    #'x_position' : results[1],
                    #'y_position' : results[2],
                    'prefactor' : results[3],
                    'prefactor variance' : results[7],
                    'q' : results[4],
                    'q variance' : results[8],
                    'sigma' : results[5],
                    'sigma variance' : results[9],
                    'chi_squared' : results[6],
                    'chi_squared variance' : 0,
                    #'eta' : results[10],
                    #'orientation_factor' : results[11],
                    #'orientation_angle' : results[12],
                    #'eta_prefactor' : results[13],
                }


    if clean:
        for key, result in results_dict.items():
            if key!='filename':
                results_dict[key] = np.nan_to_num(float(result))
        

    return results_dict     


# Experimental parameters
########################################

calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
calibration.set_image_size(1475, height=1679) # Pilatus2M
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_beam_position(754, 1075)

calibration.set_distance(5.03) # 5m


mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
mask.load('./mask.png')



# Files to analyze
########################################
source_dir = '../raw/'
output_dir = './'

pattern = '*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             #'background' : source_dir+'empty*saxs.tiff',
             #'transmission_int': '../../data/Transmission_output.csv', # Can also specify an float value.
             }
run_args = { 'verbosity' : 3,
            #'save_results' : ['xml', 'plots', 'txt', 'hdf5'],
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)



patterns = [
            ['theta', '.+_th(\d+\.\d+)_.+'] ,
            ['x_position', '.+_x(-?\d+\.\d+)_.+'] ,
            ['y_position', '.+_yy(-?\d+\.\d+)_.+'] ,
            ['anneal_time', '.+_anneal(\d+)_.+'] ,
            #['cost', '.+_Cost(\d+\.\d+)_.+'] ,
            #['annealing_temperature', '.+_T(\d+\.\d\d\d)C_.+'] ,
            #['annealing_time', '.+_(\d+\.\d)s_T.+'] ,
            #['annealing_temperature', '.+_localT(\d+\.\d)_.+'] ,
            #['annealing_time', '.+_clock(\d+\.\d\d)_.+'] ,
            #['o_position', '.+_opos(\d+\.\d+)_.+'] ,
            #['l_position', '.+_lpos(\d+\.\d+)_.+'] ,
            ['exposure_time', '.+_(\d+\.\d+)s_\d+_saxs.+'] ,
            ['sequence_ID', '.+_(\d+).+'] ,
            ]

protocols = [
    #Protocols.HDF5(save_results=['hdf5'])
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, num_rings=4, ztrim=[0.05, 0.05], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None], label_filename=True) ,
    #Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge, ztrim=[0.01, 0.001]) ,
    
    #Protocols.circular_average_q2I_fit(show=False, q0=0.0140, qn_power=2.5, sigma=0.0008, plot_range=[0, 0.06, 0, None], fit_range=[0.008, 0.022]) ,
    #Protocols.circular_average_q2I_fit(qn_power=3.5, trim_range=[0.005, 0.03], fit_range=[0.007, 0.019], q0=0.0120, sigma=0.0008) ,
    Protocols.circular_average_q2I_fit(qn_power=3.0, trim_range=[0.005, 0.035], fit_range=[0.008, 0.03], q0=0.0180, sigma=0.001) ,
    
    Protocols.metadata_extract(patterns=patterns) ,
    ]
    



# Helpers
########################################
def print_d(d, i=4):
    '''Simple helper to print a dictionary.'''
    for k, v in d.items():
        if isinstance(v,dict):
            print('{}{} : <dict>'.format(' '*i,k))
            print_d(v, i=i+4)
        elif isinstance(v,(np.ndarray)):
            print('{}{} : Ar{}: {}'.format(' '*i,k,v.shape,v))
        elif isinstance(v,(list,tuple)):
            print('{}{} : L{}: {}'.format(' '*i,k,len(v),v))
        else:
            print('{}{} : {}'.format(' '*i,k,v))

def print_results(results):
    '''Simple helper to print out a list of dictionaries.'''
    for i, result in enumerate(results):
        print(i)
        print_d(result)

def print_n(d):
    '''Simple helper to print nested arrays/dicts'''
    if isinstance(d, (list,tuple,np.ndarray)):
        print_results(d)
    elif isinstance(d, dict):
        print_d(d)
    else:
        print(d)

def val_stats(values, name='z'):
    span = np.max(values)-np.min(values)
    print("  {} = {:.2g} ± {:.2g} (span {:.2g}, from {:.3g} to {:.3g})".format(name, np.average(values), np.std(values), span, np.min(values), np.max(values)))


# Inspect .npy files
########################################
# This code can be pasted into a new ipython shell
# to help inspect the .npy files that are passed
# in the AE loop.
'''

import numpy as np

def print_d(d, i=4):
    for k, v in d.items():
        if isinstance(v,dict):
            print('{}{} : <dict>'.format(' '*i,k))
            print_d(v, i=i+4)
        else:
            print('{}{} : {}'.format(' '*i,k,v))

def print_results(results):
    for i, result in enumerate(results):
        print(i)
        print_d(result)

        
results = np.load('analyze-received.npy', allow_pickle=True); print_results(results)
results = np.load('analyze-sent.npy', allow_pickle=True); print_results(results)

results = np.load('../../measure-received.npy', allow_pickle=True); print_results(results)
results = np.load('../../measure-sent.npy', allow_pickle=True); print_results(results)

results = np.load('../../gpcamv4and5/scripts/decision-received.npy', allow_pickle=True); print_results(results)
results = np.load('../../gpcamv4and5/scripts/decision-sent.npy', allow_pickle=True); print_results(results)

print_results(results)

'''


# Run autonomous loop
########################################
def run_autonomous_loop(protocols, clear=False, verbosity=3, simulate=False):
    
    # IMPORTANT NOTE: Search for "# TOCHANGE" in the code below for
    # beamline-specific and experiment-specific assumptions that need
    # to be adjusted.


    # Connect to queue to receive the next analysis command
    #queue_PATH='/nsls2/xf12id2/data/CFN/2020_3/MFukuto/'
    queue_PATH='../../../'
    queue_PATH in sys.path or sys.path.append(queue_PATH)

    from CustomQueue import Queue_analyze as queue
    q = queue()

    if clear:
        q.clear()
    
    
    if verbosity>=3:
        print('\n\n\n')
        print('=============================')
        print('==  Autonomous Experiment  ==')
        print('=============================')
        print('\n')
        
        
    filename_re = re.compile('.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)_saxs.+') # TOCHANGE
    
    # Loop forever
    while True:
        
        results = q.get() # Get analysis command
        
        num_to_analyze = int( sum( 1.0 for result in results if result['analyzed'] is False ) )
        
        if verbosity>=3:
            print('Analysis requested for {} results (total {} results)'.format(num_to_analyze, len(results)))
        if verbosity>=10:
            print_results(results)
        
        ianalyze = 0
        for i, result in enumerate(results):
            
            if 'analyzed' in result and result['analyzed'] is False and 'filename' in result:
                ianalyze += 1
                infile = source_dir + result['filename']
                infile = infile + '_saxs.tiff' # TOCHANGE


                # Code to handle bug where saved filename doesn't exactly match what is specified in metadata:
                if not os.path.exists(infile):
                    if verbosity>=1:
                        print('Specified infile is missing. We will attempt to locate the right file based on sequence_ID.')
                    if verbosity>=5:
                        print('  infile: {}'.format(infile))
                    m = filename_re.match(infile)
                    if m:
                        sID = int(m.groups()[2])
                        if verbosity>=2:
                            print('    sequence ID: {:d}'.format(sID))
                        mfiles = glob.glob('{}*_{:d}_saxs.tiff'.format(source_dir, sID)) # TOCHANGE
                        if len(mfiles)<1:
                            if verbosity>=1:
                                print('    No file matches sequence ID {}.'.format(sID))
                        elif len(mfiles)==1:
                            infile = mfiles[0]
                            if verbosity>=1:
                                print('    Using filename: {}'.format(infile))
                        else:
                            if verbosity>=1:
                                print('    {} files match sequence ID {}'.format(len(mfiles), sID))
                                print('    Aborting.')
                            return
                    else:
                        if verbosity>=1:
                            print("    RE did not match. Aborting.")
                        return
                    

                if verbosity>=3:
                    print('        Analysis for result {}/{}, file: {}'.format(ianalyze, num_to_analyze, infile))


                if simulate:
                    val, var = np.random.random()*10, 1.0
    
                else:
                    process.run([infile], protocols, output_dir=output_dir, force=True)
                
                    # Get the result of this analysis
                    xml_file = '{}{}{}'.format( './results/', infile[len(source_dir):-5], '.xml' ) # TOCHANGE
                    new_result = autonomous_result(xml_file)
                    #result.update(new_result)
                    if 'metadata' not in result:
                        result['metadata'] = None
                    if result['metadata'] is None:
                        result['metadata'] = { 'SciAnalysis': new_result }
                    else:
                        result['metadata'].update({ 'SciAnalysis': new_result })
                

                    # TOCHANGE                
                    #val, var = new_result['prefactor'], np.square(new_result['prefactor variance'])
                    val, var = new_result['chi_squared']*1e9, 0

                
                result['measurement values'] = {
                    'values' : np.asarray([val]) ,
                    'variances' : np.asarray([var]) ,
                    'value positions' : np.asarray([[0.]]) # Positions/indices for multi-task GP
                    }
                result['analyzed'] = True


        if verbosity>=3:
            print('Analyzed {} results'.format(ianalyze))
        if verbosity>=1 and ianalyze<1:
            print('WARNING: No results were analyzed.')
        
        q.publish(results)
    

#process.run(infiles, protocols, output_dir=output_dir, force=True)
run_autonomous_loop(protocols, clear=False, verbosity=3)

