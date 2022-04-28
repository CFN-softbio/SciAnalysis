#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.Result - Class for extracting results from xml files
================================================
.. module:: SciAnalysis.XSAnalysis.Result
   :synopsis: Extracting results from xml files
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  Class for extracting results previously saved by SciAnalysis into xml
# files (by default into the directory "./analysis/results/").
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################
 
import numpy as np




# ResultsXML
################################################################################
try:
    # 'Fancy' xml library
    from lxml import etree
    USE_LXML = True
except:
    # 'Regular' xml library
    import xml.etree.ElementTree as etree # XML read/write
    USE_LXML = False
#import xml.dom.minidom as minidom
class ResultsXML(object):
    '''Simple object to help extract result values from a bunch of xml files.'''
    
    def __init__(self):
        
        #import xml.etree.ElementTree as etree
        #from lxml import etree
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
            
            if len(infiles)>250 and i%100==0:
                print( '    Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            line = [infile]
            
            result_names_e, results_e = self.extract_results_from_xml(infile, protocol)
            
            for result_name in result_names:
                idx = result_names_e.index(result_name)
                line.append(results_e[idx])
                
            results.append(line)
                
                
        return results
            
            
    def extract_multi_save_txt(self, outfile, infiles, extractions, delimeter='__', verbosity=3):
        
        results = self.extract_multi(infiles, extractions, verbosity=verbosity)
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
    
    
    def extract_dict(self, infiles, extractions, verbosity=3):
        '''Extracts results from the specified infiles, organizing them into
        a list of dictionaries.'''
        
        # TODO: kwarg to extract all possible results?
        
        results = [ {'filename': infile} for infile in infiles ]
        for i, infile in enumerate(infiles):
            
            if verbosity>=5 or (verbosity>=3 and len(infiles)>250 and i%100==0):
                print( '    Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            if verbosity>=5:
                print( '     filename: {}'.format(infile))

            try:
                for protocol, result_names in extractions:
                    if verbosity>=5:
                        print('      Procotol {} (looking for {} results)'.format(protocol, len(result_names)))
                    
                    result_names_e, results_e = self.extract_results_from_xml(infile, protocol, verbosity=verbosity)

                    if verbosity>=5:
                        print('        (found {} results)'.format(len(result_names_e)))
                    
                    for result_name in result_names:
                        if result_name in result_names_e:
                            idx = result_names_e.index(result_name)
                            
                            key = '{}__{}'.format(protocol, result_name)
                            results[i][key] = results_e[idx]
                            
            except Exception as e:
                if verbosity>=1:
                    print( '    ERROR: Extraction failed for {}'.format(infile))
                if verbosity>=5:
                    print(e)
        
        return results
    
            
    def extract_multi(self, infiles, extractions, verbosity=3):
        
        results = []
        for i, infile in enumerate(infiles):
            results.append( [infile] )
        
        ifailed = 0
        for i, infile in enumerate(infiles):
            
            if verbosity>=5 or (verbosity>=3 and len(infiles)>250 and i%100==0):
                print( '    Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            if verbosity>=5:
                print( '     filename: {}'.format(infile))

            try:
                for protocol, result_names in extractions:
                    if verbosity>=5:
                        print('      Procotol {} (looking for {} results)'.format(protocol, len(result_names)))
                    
                    result_names_e, results_e = self.extract_results_from_xml(infile, protocol, verbosity=verbosity)

                    if verbosity>=5:
                        print('        (found {} results)'.format(len(result_names_e)))
                    
                    for result_name in result_names:
                        if result_name in result_names_e:
                            idx = result_names_e.index(result_name)
                            results[i].append(results_e[idx])
                        else:
                            results[i].append('-')
                            
            except Exception as e:
                ifailed += 1
                if verbosity>=1:
                    print( '    ERROR: Extraction failed for {}'.format(infile))
                if verbosity>=5:
                    print(e)

        if verbosity>=2 and len(infiles)>0:
            print( '  Extracted {} results (failed on {}/{} = {:.1f}%)'.format(len(results)-ifailed, ifailed, len(infiles), 100.0*ifailed/len(infiles)) )

                
        return results                  
            
        
    def extract_results_from_xml(self, infile, protocol, verbosity=3):
        
        if USE_LXML:
            parser = self.etree.XMLParser(remove_blank_text=True)
        else:
            parser = self.etree.XMLParser()
        root = self.etree.parse(str(infile), parser).getroot()
        
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
                
                if element.get('error') is not None and element.get('error')!='None':
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
                if verbosity>=1:
                    print('    Errror: result has no usable data ({})'.format(element))
            
        
        return result_names, results
                
    
    # End class ResultsXML(object)
    ########################################
    
    
    
# ResultsDB
################################################################################
from pathlib import Path
import pickle
class ResultsDB():
    
    def __init__(self, source_dir='./', results_dir='results', db_file='results.db'):
        import sqlite3
        
        infile = Path(source_dir, results_dir, db_file)
        
        self.db_connection = sqlite3.connect(str(infile))
        self.db_connection.row_factory = sqlite3.Row
        self.db_cursor = self.db_connection.cursor()


    def __del__(self):
        '''Destructor for the class, called when this object
        is no longer needed.'''
        
        # Close our connection to the results database
        if self.db_connection is not None:
            if self.db_connection:
                self.db_connection.close()


    def extract_single(self, infile, remove_ext=True, verbosity=3):
        
        if remove_ext:
            infile = Path(infile).stem # Just the important filename part
        else:
            infile = Path(infile).name
        
        sql = '''-- Retrieve the most recent analyses for a given filename
        SELECT analysis_id, protocol, MAX(save_timestamp) 
        FROM analyses 
        WHERE filename=?
        GROUP BY protocol
        '''
        self.db_cursor.execute(sql, ( str(infile), ))
        analysis_rows = self.db_cursor.fetchall()
        
        results = {}
        for analysis_row in analysis_rows:
            analysis_row = dict(analysis_row)
            
            sql = '''-- Retrieve the results associated with a protocol
            SELECT *
            FROM results 
            WHERE analysis_id=?
            '''
            self.db_cursor.execute(sql, ( analysis_row['analysis_id'], ))
            result_rows = self.db_cursor.fetchall()
            
            entries = {}
            for result_row in result_rows:
                result_row = dict(result_row)
                #if result_row['value_blob'] is not None: result_row['value_blob'] = pickle.loads(result_row['value_blob'])
                
                if result_row['value'] is not None:
                    value = result_row['value']
                    if result_row['units'] is not None:
                        key = '{}_units'.format(result_row['result_name'])
                        entries[key] = result_row['units']
                    if result_row['error'] is not None:
                        key = '{}_error'.format(result_row['result_name'])
                        entries[key] = result_row['error']
                elif result_row['value_text'] is not None:
                    value = result_row['value_text']
                elif result_row['value_blob'] is not None:
                    value = pickle.loads(result_row['value_blob'])
                else:
                    value = None

                entries[ result_row['result_name'] ] = value
            
            results[ analysis_row['protocol'] ] = entries
            
        
        return results


    def extract(self, infiles, remove_ext=True, print_every=20, verbosity=3):
        
        results = {}
        
        import time
        start_time = time.time()
        for i, infile in enumerate(infiles):
            
            if verbosity>=3 and i%print_every==0:
                took = time.time()-start_time
                if i>0 and i<len(infiles):
                    estimate = (len(infiles)-i)*took/i
                    estimate = '; done in ~{:.1f}s'.format(estimate)
                else:
                    estimate = ''
                print('        extracting file {}/{} = {:.1f}% ({:.1f} s){}'.format(i, len(infiles), 100.*i/len(infiles), took, estimate))
                
            
            result = self.extract_single(infile, remove_ext=remove_ext, verbosity=verbosity)
            results[infile] = result
            
        return results
        

    def extract_pattern(self, pattern='', any_before=True, any_after=True, verbosity=3):
        
        pattern = pattern.replace('_', r'\_') # SQL treats underscore as "any single character" wildcard
        
        if any_before:
            pattern = '%{}'.format(pattern)
        if any_after:
            pattern = '{}%'.format(pattern)
        
        sql = '''-- Retrieve list of files
        SELECT DISTINCT filename 
        FROM analyses
        WHERE filename like ? ESCAPE ? 
        ORDER BY filename
        '''
        self.db_cursor.execute(sql, ( pattern, "\\", ))
        result_rows = self.db_cursor.fetchall()
        
        #for result_row in result_rows:
            #result_row = dict(result_row)
            #print(result_row)
        
        infiles = [result_row['filename'] for result_row in result_rows]
        
        return self.extract(infiles, remove_ext=False, verbosity=verbosity)
        

    # End class ResultsDB()
    ########################################


# TODO: Switch default to DB/SQL
Results = ResultsXML
#Results = ResultsDB
