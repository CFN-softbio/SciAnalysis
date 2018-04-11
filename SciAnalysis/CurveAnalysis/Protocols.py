#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.CurveAnalysis.Protocols` - 1D curve analysis protocols
================================================
.. module:: SciAnalysis.CurveAnalysis.Protocols
   :synopsis: Convenient protocols for data analysis.
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  Data analysis protocols.
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################

import os
from .Data import *
from ..tools import *


class ProcessorCurve(Processor):

    
    def load(self, infile, **kwargs):
        
        load_args = {
            'skiplines' : 1,
            'comment_char' : '#',
            'xindex' : 0,
            'yindex' : -1,
            }
        load_args.update(kwargs)

        data = DataLine(infile, **load_args)
        
        return data
        
        
class ProcessorCurveStructure(ProcessorCurve):

    
    def load(self, infile, **kwargs):
        
        load_args = {
            'skiplines' : 1,
            'comment_char' : '#',
            'xindex' : 0,
            'yindex' : -1,
            }
        load_args.update(kwargs)

        data = DataLineStructured(infile, **load_args)
        
        return data        
        

class plot(Protocol):
    
    def __init__(self, name='plot', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'blur' : None,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['blur'] is not None:
            data.smooth(run_args['blur'])
        
        outfile = self.get_outfile(data.name, output_dir)
        data.plot(outfile)
        
        return results
        
        


class structure(Protocol):
    
    def __init__(self, name='structure', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'blur' : None,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        #output_dir = os.path.join(output_dir, data.name)
        #make_dir(output_dir)
        
        results = {}
        
        if run_args['blur'] is not None:
            data.smooth(run_args['blur'])
        
        
        if run_args['verbosity']>=2:
            plot = True
        else:
            plot = False
        
        
        if True:
            # Analyze the data by 'sorting' the curve's y-values
            outfile = self.get_outfile('sort_{}'.format(data.name), output_dir)
            
            new_results = data.stats(prepend='stats_')
            results.update(new_results)
            
            data_n = DataLineStructuredSort(x=data.x, y=data.y)
            new_results = data_n.analyze(outfile, plot=plot, **run_args)
            new_results = self.prepend_keys(new_results, 'sort_')
            results.update(new_results)
            
            new_results = data_n.stats(prepend='stats_normed_')
            results.update(new_results)
        
        if True:
            # Analyze the variation in variance
            outfile = self.get_outfile('std_{}'.format(data.name), output_dir)

            data_n = DataLineStructuredStd(x=data.x, y=data.y)
            new_results = data_n.analyze(outfile, plot=plot, **run_args)
            new_results = self.prepend_keys(new_results, 'std_')
            results.update(new_results)
            
        if True:
            # Analyze curve spectrally
            outfile = self.get_outfile('fft_{}'.format(data.name), output_dir)

            data_n = DataLineStructuredFFT(x=data.x, y=data.y)
            new_results = data_n.analyze(outfile, plot=plot, **run_args)
            new_results = self.prepend_keys(new_results, 'fft_')
            results.update(new_results)

            new_results = data_n.stats(prepend='fft_stats_')
            results.update(new_results)
        
        
        return results
    
    
    