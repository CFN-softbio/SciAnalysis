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



def autonomous_result(xml_file):
    
    time.sleep(1) # Kludge to avoid corrupting XML file?
    
    extractions = [ [ 'metadata_extract', ['x_position', 'y_position'] ] , # 1, 2
                ['circular_average_q2I_fit', ['fit_peaks_prefactor1', 'fit_peaks_x_center1', 'fit_peaks_sigma1', 'fit_peaks_chi_squared', 'fit_peaks_prefactor1_error', 'fit_peaks_x_center1_error', 'fit_peaks_sigma1_error' ] ], # 3, ... 9
                #['linecut_angle_fit', ['fit_eta_eta', 'orientation_factor', 'orientation_angle', 'fit_eta_span_prefactor'] ], # 10, ... 13
                ]

    infiles = [xml_file]
    extracted_results = Results().extract_multi(infiles, extractions)
    
    result = extracted_results[0]
    #print(result)
    data_vector = {
            #'x_position' : float(result[1]),
            #'y_position' : float(result[2]),
            'prefactor' : float(result[3]),
            'prefactor variance' : np.nan_to_num(float(result[7])),
            'q' : float(result[4]),
            'q variance' : np.nan_to_num(float(result[8])),
            'sigma' : float(result[5]),
            'sigma variance' : np.nan_to_num(float(result[9])),
            'chi_squared' : float(result[6]),
            'chi_squared variance' : 0,
            #'eta' : float(result[10]),
            #'orientation_factor' : float(result[11]),
            #'orientation_angle' : float(result[12]),
            #'eta_prefactor' : float(result[13]),
        }
    

    return data_vector     


# Experimental parameters
########################################

calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
calibration.set_image_size(1475, height=1679) # Pilatus2M
calibration.set_pixel_size(pixel_size_um=172.0)
calibration.set_beam_position(753.5, 1074)

calibration.set_distance(5.025) # 5m


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
    


def load_retry(infile, max_retries=10, wait_time=0.1):
    attempts = 0
    while attempts < max_retries:
        try:
            #results = np.load(infile, allow_pickle=True).item() # SMART style
            results = np.load(infile, allow_pickle=True) # gpcam style
            break
        except Exception as ex:
            attempts += 1
            print("    Attempt {}: Exception {} occurred while loading {}".format(attempts+1, ex.__class__.__name__, infile))
            print(ex)
            
            time.sleep(wait_time)
            wait_time *= 2
            
            if attempts>=max_retries:
                raise ex
            
    return results    

def run_autonomous_loop(protocols, communication_directory='../../gpcam/data/', verbosity=3):

    if verbosity>=3:
        print('\n\n\n')
        print('=============================')
        print('==  Autonomous Experiment  ==')
        print('=============================')
        print('\n')
        
    input_file = communication_directory+'command/analysis_command.npy'
    output_file = communication_directory+'result/result.npy'

    #filename_re = re.compile('.+_(\d+)_saxs.+') # In case we need the sequence_ID
    filename_re = re.compile('.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)_saxs.+')
    
    # Loop forever
    while True:
    
        # Wait for analysis command to appear
        start_time = time.time()
        while not os.path.isfile(input_file):
            if verbosity>=3:
                print(' Waiting for command file (analysis_command.npy) [{:.1f} minutes]'.format( (time.time()-start_time)/60 ))
            time.sleep(10)
            
        if verbosity>=2:
            print('New command file found (waited {:.1f} minutes)'.format((time.time()-start_time)/60) )
            
        # Open command file
        time.sleep(0.2)
        results = load_retry(input_file)
        
        
        if verbosity>=3:
            num_items = sum(1 for v in results if 'analyzed' in v and v['analyzed'] is False)
            print('    SciAnalysis processing {} files...'.format(num_items))

        if verbosity>=10:
            print_dictionary(results)
        
        for i, result in enumerate(results): # gpcam style

            infile = source_dir + result['filename']
            #infile = infile[:-4] + '_SAXS.tif'
            if infile=='AE_sample1_anneal5_run2_x12.745_yy45.943_5.00s_30851_saxs.tiff':
                skip = False
            else:
                skip = False
            
            #if 'filename' in result: # Update all analysis
            if 'analyzed' in result and result['analyzed'] is False and 'filename' in result and not skip:


                if verbosity>=3:
                    #print('        Analysis for index {}, file: {}'.format(index, infile))
                    print('        Analysis for result {}/{}, file: {}'.format(i, num_items, infile))
                    #print('\n')
                
                process.run([infile], protocols, output_dir=output_dir, force=True)
                
                #if verbosity>=3:
                    #print('\n')
                
                xml_file = '{}{}{}'.format( './results/', infile[len(source_dir):-5], '.xml' )
                new_result = autonomous_result(xml_file)

                val, var = new_result['q'], np.square(new_result['q variance']) # TOCHANGE

                result.update(new_result)
                
                if verbosity>=5:
                    print(result)
                
                result['measurement values'] = {}
                result['measurement values']['values'] = np.asarray([val])
                result['measurement values']['variances'] = np.asarray([var])
                #result['measurement values']['variances'] = np.asarray([0.0]) # Force variances to zero
                result['measurement values']['value positions'] = np.array([[0.0]])

                result['analyzed'] = True

            else:
                # Modify old records
                #result['measurement values']['variances'] = np.asarray([0.0]) # Set all variances to zero
                pass

            
            
        
        # Update the file for gpCAM to analyze
        if verbosity>=3:
            print('    Removing command file (analysis_command.npy)')
        
        #os.remove(input_file)
        #os.replace(input_file, communication_directory+'/new_experiment_command/analysis_command_bak.npy')
        os.replace(input_file, communication_directory+'/command/analysis_command_bak.npy')
        
        if os.path.isfile(output_file):
            time.sleep(10)
            
        if os.path.isfile(output_file):
            if verbosity>=3:
                print('    Saving experiment_result.npy (file exists).')
        else:
            if verbosity>=3:
                print('    Saving experiment_result.npy (file missing).')
        np.save(output_file, results)
            
        time.sleep(1)





# Run
########################################
#print('Processing {} infiles...'.format(len(infiles)))
#infiles = [ 'AE_sample1_anneal5_run2_x12.745_yy45.943_5.00s_30851_saxs.tiff' ]
#process.run(infiles, protocols, output_dir=output_dir, force=True)

# Loop
########################################
# This is typically only used at the beamline (it loops forever, watching for new files).
#process.monitor_loop(source_dir=source_dir, pattern='*.tiff', protocols=protocols, output_dir=output_dir, force=False)

# Autonomous Experiment
########################################
run_autonomous_loop(protocols)

