#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.Multiple` - Data analysis protocols
================================================
.. module:: SciAnalysis.XSAnalysis.Multiple
   :synopsis: Convenient protocols for data analysis of multiple files.
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  Data analysis protocols that operate on file sequences (multiple files).
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################

from .Data import *
from ..tools import *
from .Protocols import *

class ProtocolMultiple(Protocol):
 
    def preliminary(self, infiles, output_dir, **run_args):

        if run_args['processor'] is None:
            print('ERROR: {} needs a reference to the Processor (so that it can load files).'.format(self.name))
            return None
        
        if len(infiles)<1:
            print('ERROR: {} needs 1 or more files ({} files supplied)).'.format(self.name, len(infiles)))
            return None
        
        
        # Identify name_base (for saving output)
        if 'name_base' in run_args:
            name_base = run_args['name_base']
        elif 'name_re' in run_args:
            f = Filename(infiles[0])
            filepath, filename, filebase, ext = f.split()
            
            import re
            match_re = re.compile(run_args['name_re'])
            m = match_re.match(filename)
            if m:
                name_base = m.groups()[0]
                print(name_base)
        else:
            import os
            prefix = os.path.commonprefix(infiles)
            filepath, name_base = os.path.split(prefix)

        if 'append_protocol_name' in run_args:
            output_dir = output_dir + self.name
            make_dir(output_dir)
        outfile = self.get_outfile(name_base, output_dir, ext=run_args['file_extension'])
        
        if (not run_args['force']) and os.path.isfile(outfile):
            print(' Skipping (internal check) {} for {}'.format(self.name, data.name))
            return results

        if run_args['verbosity']>=3:
            print(' Running for {} files; saving to: {}'.format(len(infiles), outfile))
            
        return outfile
        
 
 
class sum_images(ProtocolMultiple):
    
    def __init__(self, name='sum_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : None,
                        'blur' : None,
                        'resize' : None,
                        'ztrim' : [0.05, 0.005],
                        'pattern_re' : '^.+\/([a-zA-Z0-9_]+_)(\d+)(\.+)$',
                        'file_extension' : '-sum.npy',
                        'processor' : None,
                        'append_protocol_name' : True,
                        'force' : False,
                        'verbosity' : 3,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, infiles, output_dir, **run_args):
        
        results = {}
        
        outfile = self.preliminary(infiles, output_dir, **run_args)
        if outfile is None:
            return {}
        
        processor = run_args['processor']
        load_args = processor.load_args
        
        
        # Load first image
        data = processor.load(infiles[0], **load_args)
        data = self.transform(data, **run_args)
        
        # Iterate through remaining images
        for infile in infiles[1:]:
            # Add this new image to the data
            newdata = run_args['processor'].load(infile, **load_args)
            newdata = self.transform(newdata, **run_args)
            data.data += newdata.data
                
            
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'sum of multiple images (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
            
        np.save(outfile, data.data)
        
        
        return results
        
        
    def transform(self, data, **run_args):
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'])
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize'])        
            
        return data
    
    