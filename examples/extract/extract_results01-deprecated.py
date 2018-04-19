#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################################
# analysis.py
#    version 0.1.0
################################################################################

import numpy as np
import glob
import re
 
 
class BatchProcessor(object):
    
    def __init__(self, source_dir='./', outfile='results.dat', output_dir='./', results_file='results.dat', delimiter='\t' ):
         
        self.source_dir = source_dir
        self.output_dir = output_dir
        self.outfile = outfile
         
        self.results_file = results_file
        
        self.delimiter = delimiter
        
        self.define_res()

        
    def define_res(self):
        
        self.results_file_protocol_re = re.compile( '^\[\[.+\]\]$' )
        
        self.bracketed_num_re = re.compile( '^\[ ?(\d+\.?\d*) ?]$' )
     
     
    def extract(self, parameters, subdir_match=None, show_header=True):
         
        fout = open( self.output_dir + self.outfile , 'w' )

        sample_dirs = glob.glob( self.source_dir + '/*/' )
        sample_dirs.sort()
         
        for main_sample_dir in sample_dirs:
            
            
            # Identify subdirectories (if any)
            if subdir_match==None:
                subdirs = [ main_sample_dir ]
            else:
                subdirs = glob.glob( main_sample_dir + subdir_match + '*/' )
                subdirs.sort()
            
            
            for sample_dir in subdirs:
                
                sample_name = sample_dir[len(self.source_dir):-1]
                
                fout.write( sample_name + '\n' )
                fout.write( sample_dir + '\n' )
                fout.write( 'Comments:'+self.delimiter+' \n' )
                fout.write( 'Morphology:'+self.delimiter+' \n' )
                
                image_dirs = glob.glob( sample_dir + '/*/' )
                image_dirs.sort()
                
                if show_header:
                    # Output a header row
                    fout.write( 'image_name' )
                    for parameter in parameters:
                        protocol, par_name = parameter
                        fout.write( self.delimiter+par_name+' ('+protocol+')' )
                    fout.write( '\n' )
                    
                
                running_vals = []
                for parameter in parameters:
                    running_vals.append( [] )
                    
                for image_dir in image_dirs:
                    image_name = image_dir[len(sample_dir):-1]
                    
                    fout.write( image_name )
                    
                    results_file = image_dir + self.results_file
                    
                    for ip, parameter in enumerate(parameters):
                        protocol, par_name = parameter
                        
                        results_dict = self.retrieve_protocol_output( results_file, protocol )
                        
                        try:
                            result = results_dict[par_name]
                        except KeyError:
                            result = '?'
                        
                        if protocol=='current_custom' and par_name=='f_NN6':
                            r = results_dict = self.retrieve_protocol_output( results_file, 'hex_grain_size' )
                            result = (r['NN_equal6']*1.0)/(r['particle_count']*1.0)
                            
                        
                        fout.write( self.delimiter+str(result) )
                        
                        running_vals[ip].append( result )
                        
                    fout.write('\n')

                    
                fout.write('avg:')
                for ip, parameter in enumerate(parameters):
                    try:
                        avg = np.average( running_vals[ip] )
                    except TypeError:
                        avg = '?'
                    fout.write(self.delimiter+str(avg))
                fout.write('\n')
                fout.write('std:')
                for ip, parameter in enumerate(parameters):
                    try:
                        std = np.std( running_vals[ip] )
                    except TypeError:
                        std = '?'
                    fout.write(self.delimiter+str(std))
                fout.write('\n')
                        
                        
                fout.write( '\n\n' )
         
         
        fout.close()
         
        # END OF: def extract(self, parameters)
         
         
    def retrieve_protocol_output(self, results_file, protocol_name ):

        fin = open( results_file, 'r')
        results_lines = fin.readlines()
        fin.close()
        ipos = -1
        self.results_file_protocol_re = re.compile( '^\[\[(.+) \(.+\.pyc?\) \d\d\d\d-\d\d?-\d\d? \d\d?:\d\d:\d\d\]\]$' )
        for i, line in enumerate( results_lines ):
            m = self.results_file_protocol_re.match(line)
            if m and m.groups()[0]==protocol_name:
                ipos = i
        # ipos now points to last (most recent) instance of 'protocol_name'
                
        if ipos<0:
            print( '  Warning: protocol %s not found in %s' % (protocol_name, results_file) )
            return {}
            
        
        
        data_dictionary = {}
        # Go line by line (starting after the identified protocol header) until we reach another protocol header
        ipos += 1
        while ipos<len(results_lines):
            
            line = results_lines[ipos]
            
            m = self.results_file_protocol_re.match(line)
            if m:
                # We're done (we found the next protocol header)
                ipos = len(results_lines)
            else:
                els = line.split()
                key = els[0]
                # Remove trailing ":" (if any)
                if key[-1]==':':
                    key = key[:-1]
                
                try:
                    val = float( els[1] )
                except ValueError:
                    # Treat this as a string rather than a number
                    val = line[ len(els[0]) : ].strip()
                    
                    # Check if it is of the form [ ##.##]
                    m = self.bracketed_num_re.match(val)
                    if m:
                        val = float(m.groups()[0])

                        
                data_dictionary[key] = val
                
            ipos += 1
        
        return data_dictionary         
        
        # END OF: def retrieve_protocol_output(...)
        
        
         
    # END OF: class BatchProcessor(object)
    ########################################

     
     
     
     
#main_directory = '/home/kyager/BNL/Laser_ZA/research/2013_08-success_air_vac/SEM_raw/analysis/results/20131016_OVEN_SERIES__Glass_1mm_Ge_100nm_C48_170nm/'

main_directory = '/home/kyager/BNL/Laser_ZA/research/2013_08-success_air_vac/SEM_raw/analysis_redo2/results/MISSING_OVEN_SERIES__Glass_1mm_Ge_100nm_C48_170nm/'


subdir_match = None

#main_directory = './results/20131016_LZAVAC__6TURN_OD_VELOCITY_SERIES_Glass_Ge100_C48_170nm_05x/'
#subdir_match = 'smx'

parameters_to_extract = [
                            [ 'dots_vs_lines' , 'dot_fractional_area' ] ,
                            [ 'dots_vs_lines' , 'lines_fractional_area' ] ,
                        ]

process = BatchProcessor(main_directory, 'results_extracted.dat')
process.extract(parameters_to_extract, subdir_match=subdir_match)