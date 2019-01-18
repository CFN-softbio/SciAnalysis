#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
import numpy as np
from SciAnalysis import tools
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols



# Results
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
            
            #fout.write('
            
            header = '#filename\t{}\n'.format('\t'.join(result_names))
            fout.write(header)
            
            for result in results:
                line = '{}\n'.format('\t'.join(str(r) for r in result))
                fout.write(line)
                
    
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
                results.append(float(element.get('value')))
                
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
    




# Files to analyze
########################################

root_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/2015_10Oct_18-SAXS_all/'
source_dir = os.path.join(root_dir, 'CNC_825_1p2T_kin/')
source_dir = os.path.join(root_dir, 'CNC_825_0p56T_kin/')
source_dir = os.path.join(root_dir, 'CNC_825_0T_kin/')
source_dir = os.path.join(root_dir, 'CNC_413_1p2T_kin/')
source_dir = os.path.join(root_dir, 'CNC_165_1p2T_kin/')
source_dir = os.path.join(root_dir, 'CNC_165_0T_kin/')
source_dir = os.path.join(root_dir, 'CNC_1p65_0p56T_1_kin/')


results_dir = os.path.join(source_dir, 'SciAnalysis/', 'results/')

infiles = glob.glob(os.path.join(results_dir, '*.xml'))
infiles.sort()
#infiles = [infiles[0]]


# Extract
########################################

output_dir = os.path.join(source_dir, 'SciAnalysis', 'trend')
tools.make_dir(output_dir)

if False:
    outfile = os.path.join(source_dir, 'SciAnalysis', 'trend', 'q2I.txt')
    results = Results().extract_save_txt(outfile, infiles, protocol='circular_average_q2I', result_names=['fit_peaks_prefactor1', 'fit_peaks_prefactor1_error', 'fit_peaks_x_center1', 'fit_peaks_x_center1_error', 'fit_peaks_sigma1', 'fit_peaks_sigma1_error'])

if False:
    outfile = os.path.join(source_dir, 'SciAnalysis', 'trend', 'angle_fit.txt')
    results = Results().extract_save_txt(outfile, infiles, protocol='linecut_angle', result_names=['fit_eta_eta', 'fit_eta_eta_error', 'fit_library_eta_S', 'fit_MaierSaupe_m', 'fit_MaierSaupe_m_error', 'fit_library_MaierSaupe_S'])

if True:
    outfile = os.path.join(source_dir, 'SciAnalysis', 'trend', 'angle_fit_qm.txt')
    results = Results().extract_save_txt(outfile, infiles, protocol='linecut_angle_qm', result_names=['fit_eta_eta', 'fit_eta_eta_error', 'fit_library_eta_S', 'fit_MaierSaupe_m', 'fit_MaierSaupe_m_error', 'fit_library_MaierSaupe_S'])
