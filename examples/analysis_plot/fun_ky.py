#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os, time
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
#SciAnalysis_PATH='/GPFS/xf12id1/analysis/CFN/SciAnalysis/SciAnalysis2018C3/'
#SciAnalysis_PATH='/home/etsai/BNL/Users/software/SciAnalysis2018C3/'
SciAnalysis_PATH='/home/etsai/BNL/Users/software/SciAnalysis_pr'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols


# =============================================================================
# fit_result
# =============================================================================
def fit_result(xml_file):
    
    time.sleep(1) # Kludge to avoid corrupting XML file?
    
    extractions = [ ['circular_average_q2I', ['fit_peaks_prefactor1', 'fit_peaks_x_center1', 'fit_peaks_sigma1', 'fit_peaks_chi_squared' ] ],
                ]

    infiles = [xml_file]
    extracted_results = Results().extract_multi(infiles, extractions)
    
    result = extracted_results[0]
    data_vector = {
            'fit_peaks_b' : float(result[1]),
            'fit_peaks_prefactor1' : float(result[1]),
            'fit_peaks_x_center1' : float(result[2]),
            'fit_peaks_sigma1' : float(result[3]),
            'fit_peaks_grain_size' : float(result[4]),
            'fit_peaks_chi_squared' : float(result[4]),            
        }
    
    return data_vector


# =============================================================================
# Results
# =============================================================================
################################################################################
class Results(object):
    '''Simple object to help extract result values from a bunch of xml files.'''
    
    def __init__(self):
        
        #import xml.etree.ElementTree as etree
        from lxml import etree
        #import xml.dom.minidom as minidom
        
        self.etree = etree
        
        
    def extract_save_txt(self, outfile, infiles, protocol, result_names):
        
        results = self.extract(infiles, protocol, result_names)
        
        with open(outfile, 'w') as fout:
            
            header = '#filename\t{}\n'.format('\t'.join(result_names))
            fout.write(header)
            
            for result in results:
                line = '{}\n'.format('\t'.join(str(r) for r in result))
                fout.write(line)
                
        return results
                
    
    def extract(self, infiles, protocol, result_names):
        '''Extract the specified results-values (for the given protocol), from
        the specified files. The most recent run of the protocol is used.'''
        
        results = []
        
        for i, infile in enumerate(infiles):
            
            if i%100==0:
                print( 'Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            line = [infile]
            
            result_names_e, results_e = self.extract_results_from_xml(infile, protocol)
            
            for result_name in result_names:
                idx = result_names_e.index(result_name)
                line.append(results_e[idx])
                
            results.append(line)
                
                
        return results
            
            
    def extract_multi_save_txt(self, outfile, infiles, extractions, delimeter='__'):
        
        results = self.extract_multi(infiles, extractions)
        print('Generated {} results.'.format(len(results)))
        
        result_names_all = []
        for protocol, result_names in extractions:
            for result_name in result_names:
                result_names_all.append( '{}{}{}'.format(protocol, delimeter, result_name) )
        
        with open(outfile, 'w') as fout:
            
            header = '#filename\t{}\n'.format('\t'.join(result_names_all))
            fout.write(header)
            
            for result in results:
                line = '{}\n'.format('\t'.join(str(r) for r in result))
                fout.write(line)
                
                
        return results
    
    
            
    def extract_multi(self, infiles, extractions):
        
        results = []
        for i, infile in enumerate(infiles):
            results.append( [infile] )
        
        for i, infile in enumerate(infiles):
            
            if i%100==0:
                print( 'Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            
            for protocol, result_names in extractions:
                
                result_names_e, results_e = self.extract_results_from_xml(infile, protocol)
                
                for result_name in result_names:
                    if result_name in result_names_e:
                        idx = result_names_e.index(result_name)
                        results[i].append(results_e[idx])
                    else:
                        results[i].append('-')
                    
                
        return results                  
            
        
    def extract_results_from_xml(self, infile, protocol):
        
        parser = self.etree.XMLParser(remove_blank_text=True)
        root = self.etree.parse(infile, parser).getroot()
        
        # Get the latest protocol
        element = root
        children = [child for child in element if child.tag=='protocol' and child.get('name')==protocol]
        children_v = [float(child.get('end_timestamp')) for child in element if child.tag=='protocol' and child.get('name')==protocol]
        
        idx = np.argmax(children_v)
        protocol = children[idx]
        
        # In this protocol, get all the results (in order)
        element = protocol
        children = [child for child in element if child.tag=='result']
        children_v = [child.get('name') for child in element if child.tag=='result']
        
        idx = np.argsort(children_v)
        #result_elements = np.asarray(children)[idx]
        result_elements = [children[i] for i in idx]
        
        result_names = []
        results = []
        for element in result_elements:
            
            #print( element.get('name') )
            
            if element.get('value') is not None:
                result_names.append(element.get('name'))
                try:
                    results.append(float(element.get('value')))
                except ValueError:
                    results.append( element.get('value') )
                
                if element.get('error') is not None:
                    result_names.append(element.get('name')+'_error')
                    results.append(float(element.get('error')))
                
            elif element.get('type') is not None and element.get('type')=='list':
                
                # Elements of the list
                children = [child for child in element if child.tag=='element']
                children_v = [int(child.get('index')) for child in element if child.tag=='element']
                #print(children_v)
                
                # Sorted
                idx = np.argsort(children_v)
                children = [children[i] for i in idx]
                
                # Append values
                for child in children:
                    #print( child.get('index') )
                    result_names.append('{}_{}'.format(element.get('name'), child.get('index')))
                    results.append(float(child.get('value')))
                    
            
            else:
                print('    Errror: result has no usable data ({})'.format(element))
            
        
        return result_names, results
                
    
    # End class Results(object)
    ########################################