#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.tools` - Helpful tools
================================================
.. module:: SciAnalysis.tools
   :synopsis: Provides small tools helpful in different contexts
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  Small tools.
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
#  Things like "load_args={}" in function definition can lead to side-effects.
################################################################################




import os
import time
import numpy as np


from pathlib import Path
import sqlite3, pickle # Used for storing/retrieving from database


from SciAnalysis.settings import * #from .settings import *

try:
    # 'Fancy' xml library
    from lxml import etree
    USE_LXML = True
except:
    # 'Regular' xml library
    import xml.etree.ElementTree as etree # XML read/write
    USE_LXML = False
import xml.dom.minidom as minidom



# Helpers
########################################

def make_dir(directory):
    if not os.path.isdir(directory):
        #os.mkdir( directory )
        os.makedirs( directory )
 
def timestamp(filepath):
    statinfo = os.stat(filepath)
    filetimestamp = statinfo.st_mtime
    return filetimestamp


# Printout helpers
########################################

def print_array(data, name='array', verbosity=3):
    '''Helper code for inspecting arrays (e.g. for debugging).'''
    span = np.max(data)-np.min(data)
    if verbosity>=3:
        print('print_array for: {} (shape: {})'.format(name, data.shape))
    if verbosity>=1:
        print('    values: {:.4g} ± {:.4g} (span {:.3g}, from {:.3g} to {:.3g})'.format(np.average(data), np.std(data), span, np.min(data), np.max(data)))
    if verbosity>=4:
        print(data)

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

def val_stats(values, name='z', sizing=False):
    span = np.max(values)-np.min(values)
    if sizing:
        sizing = " [{} = {} elements]".format(values.shape, values.size)
    else:
        sizing = ""
    print("  {} = {:.2g} ± {:.2g} (span {:.2g}, from {:.3g} to {:.3g}){}".format(name, np.average(values), np.std(values), span, np.min(values), np.max(values), sizing))



# Filename
################################################################################
# TODO: Modernize by using Python3's new capabilities:
#    from pathlib import Path
class Filename(object):
    '''Parses a filename into pieces following the desired pattern.'''
    
    def __init__(self, filepath):
        '''Creates a new Filename object, to make it easy to separate the filename
        into its pieces (path, file, extension).'''
        
        self.full_filepath = filepath
        self._update()
        
        
    def _update(self):
    
        path, filename, filebase, ext = self.file_split(self.full_filepath)
        
        self.path = path            # filesystem path to file
        self.filename = filename    # filename (including extension)
        self.filebase = filebase    # filename (without extension)
        self.ext = ext              # extension
        

    def file_split(self, filepath):
        
        
        filepath, filename = os.path.split(filepath)
        filebase, ext = os.path.splitext(filename)
        
        return filepath, filename, filebase, ext


    def split(self):
        # filepath, filename, filebase, ext = f.split()
        return self.path, self.filename, self.filebase, self.ext

    def get_filepath(self):
        return self.full_filepath

    def get_path(self):
        return self.path+'/'
        
    def get_filename(self):
        return self.filename
        
    def get_filebase(self):
        return self.filebase
        
    def get_ext(self):
        return self.ext

    def timestamp(self):
        statinfo = os.stat(self.full_filepath)
        filetimestamp = statinfo.st_mtime
        return filetimestamp

    def matches_basename(self, filepath):
        path, filename, filebase, ext = self.file_split(filepath)
        
        return self.filebase==filebase
    
    def append(self, text):
        self.full_filepath = os.path.join(self.path, self.filebase + text + self.ext)
        self._update()
        return self.get_filepath()
        
    def path_append(self, path):
        self.full_filepath = os.path.join(self.path, path, self.filename)
        self._update()
        return self.get_filepath()

    def get_best_match(self, df): 
        import re, difflib

        found = 0
		# Find the scan
        m = re.findall('[0-9]*_[a-z]axs', self.filename)
        if len(m)>0:
            scan_id = int(m[0].split('_')[0])
            x = df['b_scanID'].isin([scan_id])
            if x.any():
                samplename = scan_id
                df_matched = df[x]
                found = 1
                        
        if found==0:
            #('Did not find matching scan id.')
            # Find a string from the filelist that matches the filebase the most
            length = -1
            while length > (-len(self.filebase)):
                samplename = difflib.get_close_matches(self.filebase[0:length], df['a_filename'], cutoff=0.01)[0] 
                # get the best match
                length = length-5 # TODO: Explain why a -5 step is used
                if samplename in self.filebase: break # break if the best mathch is actually in the filename
            
            # TODO: Replace above code with for loop, along the lines of:
            #samplename = None
            # for length in range(-1, -len(self.filebase), -5):
                #if (samplename is not None) and samplename in self.filebase: break
                #samplename = difflib.get_close_matches(self.filebase[0:length], filelist, cutoff=0.01)[0] 

            df_matched = df[df['a_filename'].str.contains(samplename)]

        return samplename, df_matched

    def get_best_match_old_version(self, filelist): 
        import difflib
        # Find a string from the filelist that matches the filebase the most
        length = -1
        while length > (-len(self.filebase)):
            samplename = difflib.get_close_matches(self.filebase[0:length], filelist, cutoff=0.01)[0] 
            # get the best match
            length = length-5 # TODO: Explain why a -5 step is used
            if samplename in self.filebase: break # break if the best mathch is actually in the filename
        
        # TODO: Replace above code with for loop, along the lines of:
        #samplename = None
        # for length in range(-1, -len(self.filebase), -5):
            #if (samplename is not None) and samplename in self.filebase: break
            #samplename = difflib.get_close_matches(self.filebase[0:length], filelist, cutoff=0.01)[0] 
            
            
        return samplename	

    # End class Filename(object)
    ########################################





# Processor
################################################################################
class Processor(object):
    '''Base class for processing a bunch of data files.'''
    
    def __init__(self, load_args={}, run_args={}, **kwargs):
        
        self.load_args = load_args
        self.run_args = run_args
        
        self.db_connection = None
        self.db_cursor = None


    def __del__(self):
        '''Destructor for the Processor class, called when this object
        is no longer needed.'''
        
        # Close our connection to the results database
        if self.db_connection is not None:
            if self.db_connection:
                self.db_connection.close()
    
    
    def set_files(self, infiles):
        
        self.infiles = infiles
        
        
    def set_protocols(self, protocols):
        
        self.protocols = protocols
        
        
    def set_output_dir(self, output_dir):
        
        self.output_dir = output_dir


    def access_dir(self, base, extra=''):
        '''Returns a string which is the desired output directory.
        Creates the directory if it doesn't exist.'''

        output_dir = os.path.join(base, extra)
        make_dir(output_dir)

        return output_dir

        
    def run(self, infiles=None, protocols=None, output_dir=None, force=False, ignore_errors=False, sort=False, load_args={}, run_args={}, verbosity=3, **kwargs):
        '''Process the specified files using the specified protocols.'''
        
        l_args = self.load_args.copy()
        l_args.update(load_args)
        r_args = self.run_args.copy()
        r_args.update(run_args)
        
        if infiles is None:
            infiles = self.infiles
        if sort:
            infiles.sort()
                
        if protocols is None:
            protocols = self.protocols
        for protocol in protocols:
            protocol._processor = self # Allow a protocol to access global connections
            
        if output_dir is None:
            output_dir = self.output_dir
            
            
        for infile in infiles:
            
            try:
                data = self.load(infile, **l_args)
            
                for protocol in protocols:
                    
                    output_dir_current = self.access_dir(output_dir, protocol.name)
                    
                    if not force and protocol.output_exists(data.name, output_dir_current):
                        # Data already exists
                        if verbosity>=2:
                            print(' Skipping {} for {}'.format(protocol.name, data.name))
                        
                    else:
                        if verbosity>=2:
                            print('Running {} for {}'.format(protocol.name, data.name))
                        
                        results = protocol.run(data, output_dir_current, **r_args)
                        
                        md = {}
                        md['infile'] = data.infile
                        if 'full_name' in l_args:
                            md['full_name'] = l_args['full_name']
                        if 'save_results' in r_args:
                            md['save_results'] = r_args['save_results']
                            
                        self.store_results(results, output_dir, infile, protocol, **md)
                        

            except Exception as exception:
                if SUPPRESS_EXCEPTIONS or ignore_errors:
                    # Ignore errors, so that execution doesn't get stuck on a single bad file
                    if verbosity>=1:
                        print('  ERROR ({}) with file {}.'.format(exception.__class__.__name__, infile))
                else:
                    raise


    def run_parallel(self, infiles=None, protocols=None, output_dir=None, force=False, ignore_errors=False, sort=False, load_args={}, run_args={}, verbosity=3, **kwargs):
        '''Process the specified files using the specified protocols.'''
        
        #from multiprocessing import Pool
        from joblib import Parallel, delayed
        
        l_args = self.load_args.copy()
        l_args.update(load_args)
        r_args = self.run_args.copy()
        r_args.update(run_args)
        
        if infiles is None:
            infiles = self.infiles
        if sort:
            infiles.sort()
                
        if protocols is None:
            protocols = self.protocols
        for protocol in protocols:
            protocol._processor = self # Allow a protocol to access global connections
            
        if output_dir is None:
            output_dir = self.output_dir
            
        n_jobs = r_args['num_jobs'] if 'num_jobs' in r_args else 5
        with Parallel(n_jobs=n_jobs) as parallel:
            ret = parallel( delayed(self.run_parallel_file)(infile, protocols, output_dir, force, ignore_errors, l_args, r_args, verbosity) for infile in infiles )

            
    def run_parallel_file(self, infile, protocols, output_dir, force, ignore_errors, l_args, r_args, verbosity):
            
        try:
            data = self.load(infile, **l_args)
        
            for protocol in protocols:
                
                output_dir_current = self.access_dir(output_dir, protocol.name)
                
                if not force and protocol.output_exists(data.name, output_dir_current):
                    # Data already exists
                    if verbosity>=2:
                        print(' Skipping {} for {}'.format(protocol.name, data.name))
                    
                else:
                    if verbosity>=2:
                        print('Queueing {} for {}'.format(protocol.name, data.name))
                    
                    results = protocol.run(data, output_dir_current, **r_args)
                    
                    md = {}
                    md['infile'] = data.infile
                    if 'full_name' in l_args:
                        md['full_name'] = l_args['full_name']
                    if 'save_results' in r_args:
                        md['save_results'] = r_args['save_results']
                    self.store_results(results, output_dir, infile, protocol, **md)
                    

        except Exception as exception:
            if SUPPRESS_EXCEPTIONS or ignore_errors:
                # Ignore errors, so that execution doesn't get stuck on a single bad file
                if verbosity>=1:
                    print('  ERROR ({}) with file {}.'.format(exception.__class__.__name__, infile))
            else:
                raise
            
        return 'done'



    def load(self, infile, **kwargs):
        
        data = Data2D(infile, **kwargs)
        data.infile = infile
        
        return data
        
    
    def store_results(self, results, output_dir, name, protocol, **md):
        
        if 'save_results' not in md:
            md['save_results'] = ['xml', 'sql']
        
        output_dir = self.access_dir(output_dir, 'results')
        
        if 'xml' in md['save_results']:
            # Save the results to an XML file with the name of the input file
            
            if 'full_name' in md and md['full_name']:
                outfile = os.path.join( output_dir, Filename(name).get_filename()+'.xml' )
            else:
                outfile = os.path.join( output_dir, Filename(name).get_filebase()+'.xml' )

            if os.path.isfile(outfile):
                # Result XML file already exists
                if USE_LXML:
                    parser = etree.XMLParser(remove_blank_text=True)
                else:
                    parser = etree.XMLParser()
                root = etree.parse(outfile, parser).getroot()

            else:
                # Create new XML file        
                # TODO: Add characteristics of outfile
                root = etree.Element('DataFile', name=name)


            attributes = {}
            attributes['name'] = protocol.name
            attributes['start_timestamp'] = protocol.start_timestamp
            attributes['end_timestamp'] = protocol.end_timestamp
            attributes['runtime'] = protocol.end_timestamp - protocol.start_timestamp
            attributes['save_timestamp'] = time.time()
            attributes['output_dir'] = output_dir
            attributes['outfile'] = outfile
            
            attributes.update(md)
            
            attributes = dict([k, str(v)] for k, v in attributes.items())
            prot = etree.SubElement(root, 'protocol', **attributes)

            for result_name, content in results.items():
                import numpy as np

                if isinstance(content, dict):
                    content = dict([k, str(v)] for k, v in content.items())
                    etree.SubElement(prot, 'result', name=result_name, **content)
                    
                elif isinstance(content, list) or isinstance(content, np.ndarray):
                    
                    res = etree.SubElement(prot, 'result', name=result_name, type='list')
                    for i, element in enumerate(content):
                        etree.SubElement(res, 'element', index=str(i), value=str(element))
                        
                else:
                    etree.SubElement(prot, 'result', name=result_name, value=str(content))

            tree = etree.ElementTree(root)
            if USE_LXML:
                tree.write(outfile, pretty_print=True)
            else:
                tree.write(outfile)
                
                
        if 'sql' in md['save_results']:
            # Save the results to an SQLite database
            
            #outfile = os.path.join( output_dir, 'results.db' )
            #if not os.path.isfile(outfile):
            outfile = Path(output_dir, 'results.db')
            if not outfile.is_file():
                # Create database
                self.db_connection = sqlite3.connect(str(outfile))
                self.db_cursor = self.db_connection.cursor()
                
                # Create db structure
                sql = '''-- analyses table (for each run of a Protocol)
                CREATE TABLE IF NOT EXISTS analyses (
                        analysis_id integer PRIMARY KEY,
                        protocol text NOT NULL,
                        infile text NOT NULL,
                        filename text NOT NULL,
                        infile_resolved text NOT NULL,
                        start_timestamp timestamp NOT NULL,
                        end_timestamp timestamp NOT NULL,
                        runtime real NOT NULL,
                        save_timestamp timestamp NOT NULL
                );
                '''
                self.db_cursor.execute(sql)
                sql = '''-- results table (for each result/value)
                CREATE TABLE IF NOT EXISTS results (
                        result_id integer PRIMARY KEY,
                        analysis_id integer NOT NULL,
                        protocol text NOT NULL,
                        result_name text NOT NULL,
                        value real,
                        units text,
                        error real,
                        value_text text,
                        value_blob blob
                );
                '''
                self.db_cursor.execute(sql)
                
                
            if self.db_connection is None:
                # Connect to database
                self.db_connection = sqlite3.connect(str(outfile), timeout=30)
                self.db_cursor = self.db_connection.cursor()
            
            # Store data in SQLite database
            
            # First store the analysis run
            sql = '''-- Add a new analysis run
            INSERT INTO analyses (protocol, infile, filename, infile_resolved, start_timestamp, end_timestamp, runtime, save_timestamp)
            VALUES(?,?,?,?,?,?,?,CURRENT_TIMESTAMP);
            '''
            
            insert_tuple = (protocol.name, name, str(Path(name).stem), str(Path(name).resolve()), protocol.start_timestamp, protocol.end_timestamp, protocol.end_timestamp-protocol.start_timestamp )
            
            itry = 0
            while True:
                # If multiple SciAnalysis codes are running,
                # you may get "database is locked" errors.
                # So, we just keep trying until db is free.
                try:
                    self.db_cursor.execute(sql, insert_tuple)
                except sqlite3.OperationalError as exc:
                    print('  DB locked (attempt {:d}); retrying...'.format(itry+1))
                    itry += 1
                    time.sleep(2)
                finally:
                    break
            
            
            analysis_id = self.db_cursor.lastrowid
            
            # Then store each result value
            for result_name, content in results.items():
                insert_tuple = [analysis_id, protocol.name, result_name, None, None, None, None, None]
                if isinstance(content, (int,float)):
                    insert_tuple[3] = content
                elif isinstance(content, str):
                    insert_tuple[6] = content
                elif isinstance(content, dict):
                    if 'value' in content:
                        insert_tuple[3] = content['value']
                    if 'units' in content:
                        insert_tuple[4] = content['units']
                    if 'error' in content:
                        insert_tuple[5] = content['error']
                    insert_tuple[7] = pickle.dumps(content)
                else:
                    insert_tuple[7] = pickle.dumps(content)
                
                sql = '''-- Add a result
                INSERT INTO results (analysis_id, protocol, result_name, value, units, error, value_text, value_blob) 
                VALUES(?,?,?,?,?,?,?,?);
                '''
                self.db_cursor.execute(sql, insert_tuple)
                
                
                

    def rundirs(self, indir, pattern='*', protocols=None, output_dir=None, force=False, check_timestamp=False, ignore_errors=False, sort=True, load_args={}, run_args={}, verbosity=3, **kwargs):
        
        import glob
        
        dirs = [name for name in os.listdir(indir) if os.path.isdir(os.path.join(indir, name))]
        
        if sort:
            dirs.sort()
        
        
        for directory in dirs:
            if verbosity>=2:
                print('Running directory {}'.format(directory))
            
            #infiles = glob.glob('{}/{}'.format(os.path.join(indir, directory), pattern))
            infiles = glob.glob(os.path.join(indir, directory, pattern))
        
            output_dir_current = os.path.join(output_dir, directory)
        
            self.run(infiles=infiles, protocols=protocols, output_dir=output_dir_current, force=force, ignore_errors=ignore_errors, sort=sort, load_args=load_args, run_args=run_args, **kwargs)



    def run_alternate_inner(self, infiles=None, protocols=None, output_dir=None, force=False, ignore_errors=False, sort=False, load_args={}, run_args={}, verbosity=3, **kwargs):
        '''Process the specified files using the specified protocols.
        This version defers loading data until necessary. If running multiple
        protocols, the data is reloaded many times (inefficient), but if
        running on a directory with most data already processed, this
        avoids useless loads.'''
        
        l_args = self.load_args.copy()
        l_args.update(load_args)
        r_args = self.run_args.copy()
        r_args.update(run_args)
        
        if infiles is None:
            infiles = self.infiles
        if sort:
            infiles.sort()
                
        if protocols is None:
            protocols = self.protocols
        for protocol in protocols:
            protocol._processor = self # Allow a protocol to access global connections
            
        if output_dir is None:
            output_dir = self.output_dir
            
            
        for infile in infiles:
            
            try:
                
                data_name = Filename(infile).get_filebase()
            
                for protocol in protocols:
                    
                    output_dir_current = self.access_dir(output_dir, protocol.name)
                    
                    if not force and protocol.output_exists(data_name, output_dir_current):
                        # Data already exists
                        if verbosity>=2:
                            print(' Skipping {} for {}'.format(protocol.name, data_name))
                        
                    else:
                        data = self.load(infile, **l_args)
                        
                        if verbosity>=2:
                            print('Running {} for {}'.format(protocol.name, data.name))
                        
                        results = protocol.run(data, output_dir_current, **r_args)
                        
                        md = {}
                        md['infile'] = data.infile
                        if 'full_name' in l_args:
                            md['full_name'] = l_args['full_name']
                        self.store_results(results, output_dir, infile, protocol, **md)


            except (OSError, ValueError):
                if SUPPRESS_EXCEPTIONS:
                    if verbosity>=1:
                        print('  ERROR with file {}.'.format(infile))
                else:
                    raise
                
                
    def monitor_loop(self, source_dir, pattern, protocols, output_dir=None, force=False, sleep_time=4, load_args={}, run_args={}, **kwargs):
        
        import glob

        if protocols is None:
            protocols = self.protocols
        for protocol in protocols:
            protocol._processor = self # Allow a protocol to access global connections

        if output_dir is None:
            output_dir = self.output_dir
        
        donefiles = []
        while True:

            infiles = glob.glob(os.path.join(source_dir, pattern))

            for infile in infiles:
                if infile in donefiles:
                    pass

                else:
                    self.run([infile], protocols, output_dir=output_dir, force=force, load_args=load_args, run_args=run_args, **kwargs)

                    donefiles.append(infile)

            time.sleep(sleep_time)
            
            
            
    def run_multiple_all(self, basename, infiles=None, protocols=None, output_dir=None, minimum_number=None, force=False, ignore_errors=False, sort=False, load_args={}, run_args={}, verbosity=3, **kwargs):
        '''Process the specified file sets using the specified protocols. The protocols must be able to operate on sets of datas.'''
        # This version runs on all the supplied infiles (i.e. they are all assumed to be part of the group/set).
        
        l_args = self.load_args.copy()
        l_args.update(load_args)
        r_args = self.run_args.copy()
        r_args.update(run_args)
        
        if infiles is None:
            infiles = self.infiles
        if sort:
            infiles.sort()
                
        if protocols is None:
            protocols = self.protocols
        for protocol in protocols:
            protocol._processor = self # Allow a protocol to access global connections
            
        if output_dir is None:
            output_dir = self.output_dir
            
            
        basename = Filename(basename).get_filename()
        setfiles = infiles
        try:
        
            # Load all the files into data-objects
            datas = []
            for setfile in setfiles:
                datas.append( self.load(setfile, **l_args) )
                
                
            for protocol in protocols:
                
                #outfile = protocol.get_outfile(basename, output_dir)
                output_dir_current = self.access_dir(output_dir, protocol.name)

                if not force and protocol.output_exists(basename, output_dir_current):
                    # Data already exists
                    if verbosity>=2:
                        print(' Skipping {} for {}'.format(protocol.name, basename))
                else:
                    if verbosity>=2:
                        print('Running {} for {}'.format(protocol.name, basename))
                        
                    results = protocol.run(datas, output_dir_current, basename=basename, **r_args)

                    md = {}
                    self.store_results(results, output_dir, infiles[0], protocol, **md)
                    
                
        except Exception as exception:
            if SUPPRESS_EXCEPTIONS or ignore_errors:
                # Ignore errors, so that execution doesn't get stuck on a single bad file
                if verbosity>=1:
                    print('  ERROR ({}) with file {}.'.format(exception.__class__.__name__, infile))
            else:
                raise                        
                            
                    

            
            
    def run_multiple(self, pattern_re, infiles=None, protocols=None, output_dir=None, minimum_number=None, force=False, ignore_errors=False, sort=False, load_args={}, run_args={}, verbosity=3, **kwargs):
        '''Process the specified file sets using the specified protocols. The protocols must be able to operate on sets of datas.'''
        
        l_args = self.load_args.copy()
        l_args.update(load_args)
        r_args = self.run_args.copy()
        r_args.update(run_args)
        
        if infiles is None:
            infiles = self.infiles
        if sort:
            infiles.sort()
                
        if protocols is None:
            protocols = self.protocols
        for protocol in protocols:
            protocol._processor = self # Allow a protocol to access global connections
            
        if output_dir is None:
            output_dir = self.output_dir
            
            
        import re
        pattern_re = re.compile(pattern_re)
            
        donefiles = []
            
        for infile in infiles:

            # Determine basename, from which we identify other files in the group/set
            m = pattern_re.match(infile)
            if m:
                if verbosity>=5:
                    print('    RE match for: {}'.format(infile))
                
                basename = m.groups()[0]
                basename = Filename(basename).get_filename()
                
                if basename in donefiles:
                    # Processing for this set of files already completed
                    pass
                
                else:
                    
                    setfiles = [s for s in infiles if basename in s]
                    if minimum_number is None or len(setfiles)>=minimum_number:
                        
                        try:
                        
                            # Load all the files into data-objects
                            datas = []
                            for setfile in setfiles:
                                datas.append( self.load(setfile, **l_args) )
                                
                                
                            for protocol in protocols:
                                
                                #outfile = protocol.get_outfile(basename, output_dir)
                                output_dir_current = self.access_dir(output_dir, protocol.name)

                                if not force and protocol.output_exists(basename, output_dir_current):
                                    # Data already exists
                                    if verbosity>=2:
                                        print(' Skipping {} for {}'.format(protocol.name, basename))
                                else:
                                    if verbosity>=2:
                                        print('Running {} for {}'.format(protocol.name, basename))
                                        
                                    results = protocol.run(datas, output_dir_current, basename=basename, **r_args)

                                    md = {}
                                    md['basename'] = basename
                                    self.store_results(results, output_dir, infile, protocol, **md)
                                    
                            donefiles.append(basename)
                                    
                                    
                        except Exception as exception:
                            if SUPPRESS_EXCEPTIONS or ignore_errors:
                                # Ignore errors, so that execution doesn't get stuck on a single bad file
                                if verbosity>=1:
                                    print('  ERROR ({}) with file {}.'.format(exception.__class__.__name__, infile))
                            else:
                                raise                        
                            
                    
            else:
                if verbosity>=3:
                    print('    RE did not match for: {}'.format(infile))



                

    # class Processor(object)
    ########################################



# Protocol
################################################################################
def run_default(inner_function):
    '''Standard book-keeping required for the 'run' method of any protocol.'''
    def _run_default(self, data, output_dir, **kwargs):
        
        run_args = self.run_args.copy()
        run_args.update(kwargs)
        
        if 'save_results' not in run_args: run_args['save_results'] = DEFAULT_SAVE_RESULTS
        self.extra = run_args['extra'] if 'extra' in run_args else ''
        self.ir = 1
        self.start_timestamp = time.time()

        results = inner_function(self, data, output_dir, **run_args)

        self.end_timestamp = time.time()

        return results

    return _run_default

class Protocol(object):
    '''Base class for defining an analysis protocol, which can be applied to data.'''

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.out'
        self.run_args = {}
        self.run_args.update(kwargs)

    
    def get_outfile(self, name, output_dir, extra=None, ext=None, ir=False):
        
        if ext is None:
            ext = self.default_ext
            
        if extra is None:
            if hasattr(self, 'extra'):
                extra = self.extra
            else:
                extra = ''
            
        if ir:
            name = '{:02d}_{}{}{}'.format(self.ir, name, extra, ext)
            self.ir += 1
        else:
            name = '{}{}{}'.format(name, extra, ext)
            
        return os.path.join(output_dir, name)

        
    def output_exists(self, name, output_dir):
        
        if 'file_extension' in self.run_args:
            ext = self.run_args['file_extension']
        else:
            ext = None
            
        outfile = self.get_outfile(name, output_dir, ext=ext)
        return os.path.isfile(outfile)


    def prepend_keys(self, dictionary, prepend):
        
        new_dictionary = {}
        for key, value in dictionary.items():
            new_dictionary['{}{}'.format(prepend,key)] = value
            
        return new_dictionary


    @run_default
    def run(self, data, output_dir, **run_args):
        '''This class should be over-ridden in the derived class
        used to define a given Protocol. This function will be called
        when the Protocol is run, so the main logic of the analysis
        should be within this function.'''
        
        outfile = self.get_outfile(data.name, output_dir)
        
        results = {}
        
        return results


    def get_outfile_HDF5(self, name, output_dir, extra=None, ext='.h5', save_dir='HDF5'):
        ''' Determine path and filename for the HDF5 file.
        We implicitly assume that we are given an 'output_dir' that is a protocol
        path (such as './circular_average') and so we replace the trailing directory
        with the standard path for HDF5 file saving.'''
        
        head, tail = os.path.split(output_dir)
        output_dir = '/'.join([head, save_dir])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)        
        outfile = self.get_outfile(name, output_dir, extra=extra, ext=ext, ir=False)
        
        return outfile    
    
    def save_Data2D_HDF5(self, data, label, output_dir, results=None, extra=None):
        '''Save analysis results to the HDF5 file 
        Input:
            data: Data2D object, from which data.data will be saved
            label:  a string describing the data
            results: optional  dictionary containing additional results
            output_dir: output directory for protocol (will be used to guess HDF5 path)
            extra: optional addition to filename
        Output:
            None    
        '''             
        outfile = self.get_outfile_HDF5(data.name, output_dir, extra=extra)
        
        from SciAnalysis.IO_HDF import dicttoh5
        to_save = {
            'data': data.data, 
            'label': label, 
            'x_axis': data.x_axis,
            'y_axis': data.y_axis,
            'x_scale': data.x_scale,
            'y_scale': data.y_scale,
            'results':results 
            }
        dicttoh5(to_save, outfile, overwrite_data=True, h5path='/{}'.format(self.name), mode='a')

    def save_DataLine_HDF5(self, line, name, output_dir, results=None, extra=None):
        '''Save analysis results to the HDF5 file 
        Input:
            data: DataLine object, from which data (x and y values) will be saved
            label:  a string describing the data
            results: optional  dictionary containing additional results
            output_dir: output directory for protocol (will be used to guess HDF5 path)
            extra: optional addition to filename
        Output:
            None    
        '''             
        outfile = self.get_outfile_HDF5(name, output_dir, extra=extra)
        
        from SciAnalysis.IO_HDF import dicttoh5
        # TODO: Handle case where there is only an x_err or y_err (but not both)
        if line.x_err is not None and line.y_err is not None:
            label =  [ line.x_label, line.y_label,  line.x_label+ '_err',  line.y_label+ '_err'    ]
            data = np.vstack( [line.x, line.y,  line.x_err, line.y_err  ] ).T    
        else:
            label =  [ line.x_label, line.y_label ]
            data = np.vstack( [line.x, line.y] ).T    
            
        to_save = { 'data': data, 'label': label, 'results':results }
        # TODO: Handle case where results contain arrays.
        dicttoh5(to_save, outfile, overwrite_data=True, h5path='/{}'.format(self.name), mode='a')
        
        
    def label_filename(self, data, datap=None, **run_args):
        '''Handle the common case where we want to apply the data's
        filename as the title for a graph.
        The name of data.name is applied to datap'''
        if 'label_filename' in run_args and run_args['label_filename']:
            if datap is None:
                datap = data
            datap.plot_args.update( {'title': data.name} )
        
    
    # End class Protocol(object)
    ########################################


class ProtocolMultiple(Protocol):
    
    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        outfile = self.get_outfile(basename, output_dir)
        
        results = {}
        
        return results
    
    
    # End class ProtocolMultiple(Protocol)
    ########################################
    


# get_result
################################################################################
# TODO: Rationalize when to use get_result and when to use Results()
def get_result_xml(infile, protocol):
    '''Extracts a list of results for the given protocol, from the specified
    xml file. The most recent run of the protocol is used.'''

    import numpy as np
    
    if USE_LXML:
        from lxml import etree
    else:
        import xml.etree.ElementTree as etree
    #import xml.dom.minidom as minidom

    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(infile, parser).getroot()

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
    
    results = {}
    for element in result_elements:
        
        #print( element.get('name') )
        
        if element.get('value') is not None:
            results[element.get('name')] = float(element.get('value'))
            
            if element.get('error') is not None and element.get('error')!='None':
                results[element.get('name')+'_error'] = float(element.get('error'))
            
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
                name = '{}_{}'.format(element.get('name'), child.get('index'))
                results[name] = float(child.get('value'))
                
        
        else:
            print('    Errror: result has no usable data ({})'.format(element))
        
    
    return results

    # End def get_result_xml()
    ########################################

# TODO: Add get_result_db() for extracting previous results from SQLite results.db





# Notes
################################################################################
# verbosity=0 : Output nothing
# verbosity=1 : Output only final (minimal) result
# verbosity=2 : Output 'regular' amounts of information/data
# verbosity=3 : Output all useful results
# verbosity=4 : Output marginally useful things (e.g. essentially redundant/obvious things)
# verbosity=5 : Output everything (e.g. for testing)



