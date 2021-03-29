#!/usr/bin/python3
# -*- coding: utf-8 -*-

# This generates a map (2D false-color plot or 3D height plot) for a set of
# experiments (that are presumptively defined in some (x,y) space). The code
# assumes you've already used SciAnalysis to process your data, such that you
# have XML files in your "results" sub-folder with the analysis results of
# interest. This code then compiles that data and generates the plot.

# The data can be interpolated in a variety of ways to yield a "surrogate
# model". For instance, a naive linear interpolation can be used, or 
# gpCAM (v6) can be invoked with a physics-aware kernel.

# The code can also be used to generate an animation of the sequence of
# measurements during the experiment.




# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/GPFS/xf12id1/analysis/CFN/SciAnalysis/'
#SciAnalysis_PATH='/nsls2/xf12id2/analysis/CFN/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import time
import numpy as np
#import glob
from pathlib import Path
#import re

from SciAnalysis import tools
from SciAnalysis.Data import *
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols 



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




# Data2DSM(Data2D)
########################################
class Data2DSM(Data2D):
    
    def _plot_extra(self, **plot_args):
        
        xi, xf, yi, yf = self.ax.axis()
        
        # Faded overlay
        rect = mpl.patches.Rectangle((xi,yi), xf-xi, yf-yi, linewidth=1, edgecolor='none', facecolor='white', alpha=plot_args['faded'], zorder=10)
        self.ax.add_patch(rect)
        
        
        # Scatterplot
        cmap = plot_args['cmap'] if 'cmap' in plot_args else 'viridis'
        zmin = plot_args['zmin']
        zmax = plot_args['zmax']
        #self.ax.scatter(self.x_vals, self.y_vals, s=40, c=self.z_vals, cmap=cmap, vmin=zmin, vmax=zmax, edgecolor='k', zorder=100)
        
        
        # Colorbar
        n = 5
        colorbar_labels = [ zmin + i*(zmax-zmin)/(n-1) for i in range(n) ]
        
        tick_positions = self._plot_z_transform(data=colorbar_labels, set_Z=False)
        cbar = self.fig.colorbar(self.im, ax=self.ax, ticks=tick_positions, fraction=0.056, pad=0.02)
        colorbar_labels = ["{:.3g}".format(c) for c in colorbar_labels]
        cbar.ax.set_yticklabels(colorbar_labels, size=18)
        
        # Titles
        size = plot_args['rcParams']['axes.labelsize']
        #size = plot_args['rcParams']['xtick.labelsize']
        plt.figtext(0, 1, '$N = {:,d}$'.format(len(self.z_vals)), size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
        z_label = self.z_rlabel if self.z_rlabel is not None else self.z_label
        if z_label is not None:
            plt.figtext(1, 1, z_label, size=size, weight='bold', verticalalignment='top', horizontalalignment='right')


        #self.ax.axhline(0, linewidth=1, color='g', alpha=0.75)
        #self._plot_guides(**plot_args)
            
        #self.ax.set_aspect('equal', 'datalim')
        # How to set? 'auto', 'equal', num (force ratio)
        # What to adjust? None, 'box', 'datalim'
        
        
    def _plot_guides(self, **plot_args):
        '''Example of overlaying some guide-lines, obtained using the 
        SurrogateModel (self.model).'''
        
        #v_print, y_pos = self.xy_axes() # Index values
        v_print, y_pos = self.x_axis, self.y_axis

        period = self.model.transform_period(v_print)
        phase = self.model.transform_phase(v_print)
        
        miniguide = 0.05
        for i in range(-3, +3+1, 1):
            line = phase + period*i
            self.ax.plot(v_print, line, '-', color='g', dashes=[5,5], linewidth=2, alpha=0.75, zorder=100)

            self.ax.plot(v_print, line+miniguide, '-', color='g', dashes=[5,5], linewidth=1, alpha=0.5, zorder=100)
            self.ax.plot(v_print, line-miniguide, '-', color='g', dashes=[5,5], linewidth=1, alpha=0.5, zorder=100)
            self.ax.plot(v_print, line+2*miniguide, '-', color='g', dashes=[5,5], linewidth=0.5, alpha=0.25, zorder=100)
            self.ax.plot(v_print, line-2*miniguide, '-', color='g', dashes=[5,5], linewidth=0.5, alpha=0.25, zorder=100)
            
        self.ax.plot(self.model.guide_points[:,0], self.model.guide_points[:,1], 'o', color='g', markersize=10, alpha=0.5)
        
        
    def _plot_extra3D(self, **plot_args):
        # Colorbar
        cbar = self.fig.colorbar(self.surf, ax=self.ax, aspect=40, fraction=0.02, pad=0.0)
        cbar.ax.yaxis.set_tick_params(labelsize=15)

        # Titles
        size = plot_args['rcParams']['axes.labelsize']
        #size = plot_args['rcParams']['xtick.labelsize']
        plt.figtext(0, 1, '$N = {:,d}$'.format(len(self.z_vals)), size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
        z_label = self.z_rlabel if self.z_rlabel is not None else self.z_label
        if z_label is not None:
            plt.figtext(1, 1, z_label, size=size, weight='bold', verticalalignment='top', horizontalalignment='right')

    ########################################
    # End: class Data2DSM(Data2D)



class Base():
    # Helpers
    ########################################
    def msg(self, txt, threshold=3, indent=0, indent_txt='  ', verbosity=None, empty_lines=0, **kwargs):
        '''Outputs a status line indicating the current state of execution.'''
        if verbosity is None:
            verbosity = self.verbosity
        if verbosity>=threshold:
            indent = np.clip(indent, 0, 10)
            indent = indent_txt*indent
            for i in range(empty_lines):
                print('')
            print('{}> {}{}'.format(self.name, indent, txt))

    def msgm(self, txt=None, threshold=0, indent=0, indent_txt='  ', verbosity=None, mark='=', nmark=80, empty_lines=1, **kwargs):
        '''Outputs a very noticeable message, demarcated by lines.'''
        if verbosity is None:
            verbosity = self.verbosity
        if verbosity>=threshold:
            for i in range(empty_lines):
                print('')
            
            self.msg(txt=mark*nmark, threshold=threshold, indent=indent, indent_txt=indent_txt, verbosity=verbosity, **kwargs)
            if txt is not None:
                self.msg(txt=txt, threshold=threshold, indent=indent, indent_txt=indent_txt, verbosity=verbosity, **kwargs)
                self.msg(txt=mark*nmark, threshold=threshold, indent=indent, indent_txt=indent_txt, verbosity=verbosity, **kwargs)

    def print_data(self):
        tools.print_d(self.data)
        
    def timing_start(self):
        self.start_time = time.time()
        
    def timing_end(self):
        return time.time() - self.start_time
    
    def timing_end_msg(self, txt='', iterations=None, threshold=3, indent=0, indent_txt='  ', verbosity=None, empty_lines=0, **kwargs):
        took = self.timing_end()
        if iterations is None:
            txt = '{} took {:.1f}s'.format(txt, took)
        else:
            txt = '{} took {:.1f}s for {} iterations ({:.3f} s/iteration)'.format(txt, took, iterations, took/iterations)
        self.msg(txt, threshold=threshold, indent=indent, indent_txt=indent_txt, verbosity=verbosity, empty_lines=empty_lines, **kwargs)

    def timing_progress_msg(self, icurrent, itotal, threshold=4, indent=4, indent_txt='  ', every=50, verbosity=None):
        if verbosity is None:
            verbosity = self.verbosity
        if verbosity>=threshold:
            if icurrent%every==0:
                amt = icurrent/itotal
                took = self.timing_end()
                if icurrent>0 and icurrent<itotal:
                    estimate = (itotal-icurrent)*took/icurrent
                    estimate = '; done in ~{:.1f}s'.format(estimate)
                else:
                    estimate = ''
                print("{}{}/{} = {:.1f}% ({:.1f}s{})".format(indent_txt*indent, icurrent, itotal, 100.*icurrent/itotal, took, estimate))
                                                                                    


# SurrogateModel(Base)
########################################
class SurrogateModel(Base):
    
    def __init__(self, sample='sample', sample_pattern=None, name='SM', source_dir='../', output_dir=None, verbosity=3, **kwargs):
        
        self.kwargs = kwargs
        self.sample = sample
        self.sample_pattern = sample if sample_pattern is None else sample_pattern
        self.name = name
        self.source_dir = source_dir
        self.output_dir = Path('./', self.sample) if output_dir is None else output_dir
        
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        self.verbosity = verbosity
        self.msg("init SurrogateModel named '{}' for sample '{}'".format(self.name, self.sample), 3, 0)
        
        
        # Default associations (assumptions about what keys mean)
        self.associations = {
            'parameters': [
                'metadata_extract__x_position',
                'metadata_extract__y_position',
                ],
            'metadata': [
                '#filename', 
                'filename',
                'metadata_extract__sequence_ID',
                ],
            #'signals': # Everything else will be considered a signal
        }


        # Internal store of all the possible signals to plot
        self.data = { 'parameters': {}, 'signals': {}, 'metadata': {} }
        self.N_limit = None
        
        # Selected values just before interpolating
        self.x_vals, self.y_vals, self.z_vals = None, None, None
        self.x_name, self.y_name, self.z_name = None, None, None
        
        
        # The associations are used to convert a signal codename into a
        # pretty name (suitable for plotting). The list is traversed in 
        # reverse order, so that you can append high-priority matches.
        # Each line should be 5 elements:
        #   [ match, label, zmin, zmax, cmap ]
        self.label_associations = [
            [ 'x', '$x$', None, None, None ] ,
            [ 'y', '$y$', None, None, None ] ,
            [ 'z', '$z$', None, None, None ] ,
            [ 'x_pos', '$x \, (\mathrm{mm})$', None, None, None ] ,
            [ 'y_pos', '$y \, (\mathrm{mm})$', None, None, None ] ,
            [ 'sequence_ID', '$\mathrm{sID}$', 0, None, None ] ,
            [ 'prefactor', r'$p \, (\mathrm{a.u.})$', 0, None, cmap_vge ],
            [ 'chi_squared', '$\chi^2\, (\mathrm{a.u.})$', 0, None, 'plasma' ],            
            [ 'fit_peaks_x_center1', '$q_0 \, (\mathrm{\AA^{-1}})$', None, None, 'viridis' ] ,
            [ 'fit_peaks_d0', '$d_0 \, (\mathrm{nm})$', None, None, 'viridis' ] ,
            [ 'fit_peaks_sigma1', r'$\sigma_0 \, (\mathrm{\AA^{-1}})$', None, None, 'inferno' ] ,
            [ 'fit_peaks_grain_size', r'$\xi \, (\mathrm{nm})$', 0, None, 'inferno' ] ,
            [ 'fit_peaks_prefactor1', '$p \, (\mathrm{a.u.})$', 0, None, cmap_vge ] ,
            [ 'fit_eta_eta', r'$\eta$', 0, 1, 'inferno' ],
            [ 'orientation_factor', r'$f_{\mathrm{ori}}$', -1, 1, 'gray' ],
            [ 'orientation_angle', r'$\chi \, (\mathrm{^{\circ}})$', -90, 90, cmap_cyclic_rb ],
            ]
        
        # User can append overrides for the zscale (zmin, zmax) to be
        # used in plotting.
        self.zscales = [
            #[ 'prefactor1', 0, 1 ],
            ]
        

        
    # Resolvers
    ########################################
    
    def label(self, name, retall=False):
        
        for match, label, zmin, zmax, cmap in reversed(self.label_associations):
            if match in name:
                if retall:
                    return label, zmin, zmax, cmap
                else:
                    return label
            
        if retall:
            return '$x_i$', None, None, None
        else:
            return '$x_i$'
        
        
    def zscale(self, name):
        
        label, zmin, zmax, cmap = self.label(name, retall=True)
        
        if cmap is None: cmap = 'jet'
        
        for match, zmin_c, zmax_c in reversed(self.zscales):
            if match in name:
                if zmin_c is not None: zmin = zmin_c
                if zmax_c is not None: zmax = zmax_c
        
        if zmin is None: zmin = np.min(self.z_vals)
        if zmax is None: zmax = np.max(self.z_vals)
        
        return zmin, zmax, cmap

    
    def add_zscales(self, signals):
        for signal in signals:
            self.zscales.append(signal)
        

    def get_assoc_type(self, name):
        '''Determine the association type for a given signal name.
        Signal names are of the form: protocol__signal'''
        found = False
        for assoc_type in self.associations.keys(): # Consider 'parameters', 'metadata', ...
            if name in self.associations[assoc_type]:
                # The column 'name' is of type 'assoc_type'; e.g. 'x_position' is 'parameters'
                found = True
                break
    
        if not found:
            # Default assumption is that it is a signal
            assoc_type = 'signals'
            
        return assoc_type    

    
    
    # Load data
    ########################################

    def save_data(self, outfile=None):
        '''Save the model data list to a npy file.'''
        if outfile is None:
            outfile = Path(self.output_dir, '{}-{}.npy'.format(self.name, self.sample))
        np.save(outfile, self.data, allow_pickle=True)
        
        
    def load_data(self, infile=None):
        '''Load data from npy file.
        Data should be formatted exactly as SurrogateModel expects it
        (e.g. from model.save_data()).'''
        if infile is None:
            infile = Path(self.output_dir, '{}-{}.npy'.format(self.name, self.sample))
        self.data = np.load(infile, allow_pickle=True).item()


    def load_sql(self, extractions):
        '''Extract results from the SQLite database saved by SciAnalysis.'''
        
        self.msg("Extracting data from SQLite file", 3, 0)
        
        from SciAnalysis.Result import ResultsDB
        results = ResultsDB(source_dir=self.source_dir).extract_pattern(self.sample_pattern)
        
        n = 0
        for protocol, signals in extractions:
            for signal in signals:
                n += 1
                
                name = '{}__{}'.format(protocol, signal)
                assoc_type = self.get_assoc_type(name)
                if name not in self.data[assoc_type]:
                    # Create an empty array for this signal
                    self.data[assoc_type][name] = np.asarray([])


        # Reorganize the results
        self.msg("Reorganizing {} results, assuming {} columns".format(len(results), n), 5, 1)
        
        skips = 0
        for filename, result in results.items():
            # Check that each target key exists
            nc = 0
            for protocol, signals in extractions:
                if protocol in result:
                    for signal in signals:
                        if signal in result[protocol]:
                            nc += 1
                        
            self.msg("n = {}; nc = {}; filename: {}".format(n, nc, filename), 7, 2)
            if nc!=n:
                self.msg("WARNING: nc = {} (expected n = {}) for filename: {}".format(nc, n, filename), 2, 3)
                skips += 1
                
            else:
                # We have all the requested signals; add this result to the data
                for protocol, signals in extractions:
                    for signal in signals:
                        value = result[protocol][signal]
                        name = '{}__{}'.format(protocol, signal)
                        assoc_type = self.get_assoc_type(name)
                        self.data[assoc_type][name] = np.append(self.data[assoc_type][name], value)
                    
            if self.verbosity>=11:
                tools.print_n(result)

        self.msg("Reorganized {} results, skippings {}/{} = {:.1f}%".format(len(results), skips, len(results), 100.*skips/len(results)), 5, 2)
        
    
    def load_xml(self, extractions, save_cache=True, use_cached=False, results_dir='results'):
        '''Extract results from XML files, as saved by SciAnalysis.'''
        
        outfile = Path(self.output_dir, '{}-extracted.npy'.format(self.sample_pattern))
        if use_cached and outfile.is_file():
            self.msg("Using results cached in: {}".format(outfile), 4, 1)
            results = np.load(outfile, allow_pickle=True)
            
        else:
            results_dir = Path(self.source_dir, results_dir)
            self.msg("Extracting data from XML files in {}".format(results_dir), 3, 0)
            
            infiles = list( results_dir.glob('{}*.xml'.format(self.sample_pattern)) )
            self.msg("Found {} XML files".format(len(infiles)), 4, 1)
                
            from SciAnalysis.Result import Results # Results() object
            self.msg("Extracting results for {} infiles".format(len(infiles)), 3, 1)
        
            results = Results().extract_dict(infiles, extractions, verbosity=self.verbosity)
            
            if self.verbosity>=6:
                tools.print_results(results)
            
            if save_cache:
                self.msg("Extracted results saved to: {}".format(outfile), 4, 1)
                np.save(outfile, results, allow_pickle=True)


        # Reorganize the results
        n = len(results[0])
        self.msg("Reorganizing {} results, assuming {} columns".format(len(results), n), 5, 1)
        skips = 0
        for i, result in enumerate(results):
            if len(result.keys())!=n:
                skips += 1
                self.msg("ERROR: Skipping result {}, which has {} keys (expecting {})".format(i, len(result.keys()), n), 2, 1)
            else:
                for name, value in result.items():
                    assoc_type = self.get_assoc_type(name)
                
                    # Put this value where it belongs
                    self.msg('Putting {} into {}'.format(name, assoc_type), 10, 2)
                    if name not in self.data[assoc_type]:
                        self.data[assoc_type][name] = np.asarray([value])
                    else:
                        self.data[assoc_type][name] = np.append(self.data[assoc_type][name], value)
                        
        self.msg("Reorganized {} results, skippings {}/{} = {:.1f}%".format(len(results), skips, len(results), 100.*skips/len(results)), 5, 2)
            
            
    def load_xml_extraction(self, extractions, use_cached=False, results_dir='results'):
        '''Extract results from XML files, as saved by SciAnalysis.'''
        # DEPRECATED: Use model.load_xml() instead (which internally uses Results().extract_dict)
        
        outfile = Path(self.output_dir, '{}-extracted.txt'.format(self.sample_pattern))
        if use_cached and outfile.is_file():
            self.msg("Using results cached in: {}".format(outfile), 4, 1)
        
        else:
            results_dir = Path(self.source_dir, results_dir)
            self.msg("Extracting data from XML files in {}".format(results_dir), 3, 0)
            
            infiles = list( results_dir.glob('{}*.xml'.format(self.sample_pattern)) )
            self.msg("Found {} XML files".format(len(infiles)), 4, 1)
                
            results = self.extract_results_from_xml(infiles, extractions, outfile=outfile)
            
        names, lines = self.load_extraction_file(outfile)
        self.handle_extraction_lines(names, lines)

        
    def extract_results_from_xml(self, infiles, extractions, outfile):
        '''Extract results from a sequence of xml files that
        were saved using SciAnalysis.'''
        from SciAnalysis.Result import Results # Results() object
        self.msg("Extracting results for {} infiles".format(len(infiles)), 3, 1)
    
        results = Results().extract_multi_save_txt(outfile, infiles, extractions, verbosity=self.verbosity)
        self.msg("Extracted results saved to: {}".format(outfile), 4, 1)
        
        return results                                 


    def load_extraction_file(self, infile):
        '''Load results from a simple text file, such as that
        saved from SciAnalysis extraction routines.'''
        
        self.msg("Loading data from file: {}".format(infile), 3, 1)
        
        with open(infile, 'r') as fin:
            names = fin.readline().split()
            lines = fin.readlines()

        self.msg("Saved data has {} columns and {} lines".format(len(names), len(lines)), 4, 2)
        if self.verbosity>=5:
            print(names)
        
        return names, lines
    
    
    def handle_extraction_lines(self, names, lines):
        self.msg("Handling extracted data ({} lines)".format(len(lines)), 3, 1)

        
        # Prepare to place each column into the appropriate part of self.data
        assoc_types = []
        for name in names: # Iterate through the columns
            
            found = False
            for assoc_type in self.associations.keys(): # Consider 'parameters', 'metadata', ...
                if name in self.associations[assoc_type]:
                    # The column 'name' is of type 'assoc_type'; e.g. 'x_position' is 'parameters'
                    assoc_types.append(assoc_type)
                    if name not in self.data[assoc_type]:
                        # Create the array
                        self.data[assoc_type][name] = np.asarray([])
                    found = True
                    break
                
            if not found:
                # Default assumption is that it is a signal
                assoc_types.append('signals')
                if name not in self.data['signals']:
                    self.data['signals'][name] = np.asarray([])
                
            
        skips = 0
        nans = 0
        for i, line in enumerate(lines):
            if len(lines)>250 and i%100==0:
                self.msg("Line {}/{} = {:.1f}% done".format(i, len(lines), 100.*i/len(lines)), 3, 3)

            els = line.split()                
            if len(els)==len(names) and els[0][0]!='#':
                
                for idx, el in enumerate(els):
                    if el=='-':
                        value = np.nan
                        nans += 1
                    else:
                        try:
                            value = float(el)
                        except ValueError:
                            value = el
                    
                    self.data[assoc_types[idx]][names[idx]] = np.append( self.data[assoc_types[idx]][names[idx]], value )
                        
            else:
                skips += 1
                self.msg("Skipping line {}: {}".format(i, line.strip()), 3, 3)
                
        num_signals = len(self.data['signals'].keys())
        self.msg("Handled {}/{} = {:.1f}% lines ({} nans over {} signals ~ {:.1f} nans/signal)".format(len(lines)-skips, len(lines), 100.*(len(lines)-skips)/len(lines), nans, num_signals, nans/num_signals ), 3, 1)



    # Manage/modify data
    ########################################
    
    def sort_data(self, assoc_type='metadata', key='metadata_extract__sequence_ID'):
        indices = np.argsort(self.data[assoc_type][key])
        for data_type, data_list in self.data.items():
            for name, data in data_list.items():
                self.data[data_type][name] = data[indices]
    
    
    def get_N(self):
        if self.N_limit is None:
            return len(self.z_vals)
        else:
            return self.N_limit

    def get_N_max(self):
        signal = list(self.data['signals'].values())[0]
        return len(signal)

    def set_N_limit(self, N_limit=None):
        self.N_limit = N_limit

    def change_N_limit(self, N_limit):
        self.select(x=self.x_name, y=self.y_name, signal=self.z_name, N_limit=N_limit)
        

    def select(self, x='x_pos', y='y_pos', signal=None, N_limit=None):
        
        self.N_limit = N_limit
        
        if signal is None:
            signal = list(self.data['signals'].keys())[0]

        self.msg("Selecting parameters/signals, based on x='{}', y='{}', signal='{}'".format(x, y, signal), 4, 0)
            
        for name, data in self.data['parameters'].items():
            if x in name:
                self.msg("Using '{}' ({} values) as x-coordinate".format(name, len(data)), 4, 1)
                self.x_name = name
                self.x_vals = self.data['parameters'][name]
            if y in name:
                self.msg("Using '{}' ({} values) as y-coordinate".format(name, len(data)), 4, 1)
                self.y_name = name
                self.y_vals = self.data['parameters'][name]
                
        for name, data in self.data['signals'].items():
            if signal in name:
                self.msg("Using '{}' ({} values) as z-values".format(name, len(data)), 4, 1)
                self.z_name = name
                self.z_vals = self.data['signals'][name]

        N = self.get_N()
        if self.N_limit is not None:
            self.x_vals, self.y_vals, self.z_vals = self.x_vals[:N], self.y_vals[:N], self.z_vals[:N]
                
        nans = np.count_nonzero(np.isnan(self.z_vals))
        if nans>0:
            self.msg("NOTE: selected signal includes nan for {}/{} = {:.1f}% of the entries".format(nans, N, 100.*nans/N), 2, 2)


    def select_single(self, find, order=['parameters', 'signals', 'metadata']):
        
        for assoc_type in order:
            for name, data in self.data[assoc_type].items():
                if find in name:
                    return data
                
        return None
        

    # Coordinate transformations
    ########################################
    # These are examples of coordinate transformations that one could apply to
    # data. The intended use is that the coord_transform function adds new
    # entries to self.data['parameters'], which allows the user to then plot
    # using the original/raw coordinates, or the newly-defined coordinates.

    def coord_transform_Th(self):
        
        # Convert x-coordinate to thickness (nm)
        x_vals = self.select_single('x_position')
        h = (x_vals/50)*(200-140) + 140
        self.data['parameters']['thickness'] = h
        self.label_associations.append(['thickness', '$h \, (\mathrm{nm})$', 0, None, None])
        
        # Convert y-coordinate to Temperature (C)
        y_vals = self.select_single('y_position')
        T = ((y_vals+50)/50)*(200-30) + 30 # Temperature
        self.data['parameters']['temperature'] = T
        self.label_associations.append(['temperature', '$T \, (\mathrm{^{\circ}C})$', 25, 500, None])
        
        
    
    def coord_transform_vy(self, acceleration=25.0):
        '''Transform (x,y) mapping data into (v_print, y_corrected) data
        for a 3D-printed material viewed in cross-section.'''
        
        # Convert from x_position into the print velocity at that position
        x_vals = self.select_single('x_position')
        v_print = np.sqrt(2*acceleration*np.abs(x_vals))
        self.data['parameters']['v_print'] = v_print
        self.label_associations.append(['v_print', '$v_{\mathrm{print}} \, (\mathrm{mm/s})$', 0, None, None])
        
        # Undo the undulations in the rows of printed material
        self.guide_points = np.asarray([
            [0, 3.05],
            [10, 3.12],
            [20, 3.194],
            [30, 3.09],
            [37, 2.99],
            [45, 2.95],
            [52, 3.00],
            [55, 3.10],
            [60, 3.294],
            [64, 3.575],
            [67, 4.09],
            ])
        
        y_vals = self.select_single('y_position')
        y_corrected = y_vals - self.transform_phase(v_print)
        self.data['parameters']['y_corrected'] = y_corrected
        self.label_associations.append(['y_corrected', '$y_{\mathrm{corrected}} \, (\mathrm{mm})$', None, None, None])

    def transform_phase(self, v_print):
        # Fit the guide_points to a spline
        import scipy.interpolate
        f = scipy.interpolate.interp1d(self.guide_points[:,0], self.guide_points[:,1], kind='cubic')
        return f(v_print)
    def transform_period(self, v_print, slope=1.59483070e-05, intercept=0.295558503):
        return slope*v_print + intercept
    
    def enforce_periodic(self, axis='y_corrected', period=0.295558503, extend=1):
        '''Assume the system is periodic in a given direction, such that
        we can collapse/repeat the data.'''
        
        coord_vals = self.select_single(axis)
        coord_vals = np.mod(coord_vals, period)
        self.data['parameters'][axis] = coord_vals
        
        if extend>0:
            n_copies = 2*extend+1
            # Repeat the data along the periodic direction
            for data_type, data_list in self.data.items():
                for name, data in data_list.items():
                    if axis in name:
                        # Repeat the data, extending along axis
                        new_data = np.asarray([])
                        for i in range(-extend, extend+1):
                            new_data = np.concatenate( (new_data,data+i*period) )
                        self.data[data_type][name] = new_data
                    else:
                        # Repeat the data n times
                        self.data[data_type][name] = np.tile(data, n_copies)



    # Interpolate
    ########################################
    
    def interpolate(self, interp_mode='griddata', **kwargs):
        self.msg('Interpolating using method: {}'.format(interp_mode), 3, 0)
        if self.verbosity>=4:
            tools.val_stats(self.z_vals, name='z_vals')
        
        getattr(self, 'interpolate_{}'.format(interp_mode))(**kwargs)
    
    
    def interpolate_griddata(self, method='linear', rescale=True, convex_clip=False, **kwargs):
        grid, xi, yi, XI, YI = self.make_grid(**kwargs)
        self.msg("Interpolating {:,} points to {:,}×{:,} = {:,} points (densification = {:.0f})".format(len(self.z_vals), len(xi), len(yi), len(xi)*len(yi), len(xi)*len(yi)/len(self.z_vals)), 3, 1)
        
        import scipy.interpolate
        POINTS = np.column_stack((self.x_vals, self.y_vals))
        VALUES = self.z_vals
    
        if 'interpolate_cyclic' in kwargs and kwargs['interpolate_cyclic'] is not None:
            cycle = kwargs['interpolate_cyclic']
            xlike = np.cos(VALUES*2*np.pi/cycle)
            ylike = np.sin(VALUES*2*np.pi/cycle)
            XLIKE = scipy.interpolate.griddata(POINTS, xlike, (XI, YI), method='linear')
            YLIKE = scipy.interpolate.griddata(POINTS, ylike, (XI, YI), method='linear')
            
            ZI = ( np.arctan2(YLIKE, XLIKE)/(2*np.pi) )*cycle
            
        else:        
            ZI = scipy.interpolate.griddata(POINTS, VALUES, (XI, YI), rescale=rescale, method=method) # method='nearest' 'linear' 'cubic'
        
        if convex_clip:
            # Force the grid to only exist within the convex hull defined by the points
            # This is useful for method='nearest', which otherwise will extrapolate everywhere.
            clip = np.isnan(scipy.interpolate.griddata(POINTS, VALUES, (XI, YI), rescale=rescale, method='linear'))
            ZI = np.ma.masked_where( clip, ZI)
        
        ZI_mask = np.ma.masked_where( np.isnan(ZI), ZI)

        if self.verbosity>=4:
            #tools.val_stats(ZI, name='ZI')
            tools.val_stats(ZI_mask, name='ZI_mask')
            
        self.grid = grid
        self.xi, self.yi = xi, yi
        self.XI, self.YI = XI, YI
        self.ZI = ZI_mask
    
    
    def make_grid(self, grid=None, d_grid=None, n_grid=200, **kwargs):
        
        # Define grid for interpolation
        if grid is None:
            grid = [None, None, None, None]
        if grid[0] is None: grid[0] = np.min(self.x_vals)
        if grid[1] is None: grid[1] = np.max(self.x_vals)
        if grid[2] is None: grid[2] = np.min(self.y_vals)
        if grid[3] is None: grid[3] = np.max(self.y_vals)
                
        if d_grid is None:
            if not isinstance(n_grid, (list, tuple, np.ndarray)):
                n_grid = [n_grid, n_grid]
            d_grid = [ (grid[1]-grid[0])/n_grid[0], (grid[3]-grid[2])/n_grid[1] ]
        elif isinstance(d_grid, float):
            d_grid = [d_grid, d_grid]
        
        xi = np.arange(grid[0], grid[1]+d_grid[0], d_grid[0])
        yi = np.arange(grid[2], grid[3]+d_grid[1], d_grid[1])
        XI, YI = np.meshgrid(xi, yi)
        
        return grid, xi, yi, XI, YI
        
        

    def interpolate_gpcam(self, hps_guess=None, gp_method='global', pre_optimize=False, fill='pixelwise', gpcam_PATH=None, **kwargs):
        
        # TOCHANGE: Make sure gpCAM is available
        #gpcam_PATH='/home/kyager/current/code/gpcam/main/'
        #gpcam_PATH='../../../gpcamv4and5/gpcam/'
        if gpcam_PATH is not None:
            gpcam_PATH in sys.path or sys.path.append(gpcam_PATH)
        from gpcam.gp_optimizer import GPOptimizer
        
        params = np.stack( (self.x_vals, self.y_vals), axis=1 )
        signals = np.stack( (self.z_vals, ), axis=1 )
        
        # index_set_bounds is the valid range for each parameter
        # We set this to the actual min/max range of the parameters.
        index_set_bounds = np.asarray([ [np.min(col), np.max(col)] for col in params.T ])

        # Hyperparameters depend on the kernel definition
        gp_kernel, hps_bounds, hps_guess = self.hps_initialize_default(params, signals, hps_guess=hps_guess)
        #gp_kernel, hps_bounds, hps_guess = self.hps_initialize_periodic2Daniso(params, signals, hps_guess=hps_guess)
        
        gp = GPOptimizer(
                input_space_dimension=params.shape[1],
                output_space_dimension=1,
                output_number=signals.shape[1],
                index_set_bounds=index_set_bounds,
                hyperparameter_bounds=hps_bounds,
                gp_kernel_function=gp_kernel,
                #objective_function=obj_func,
                )
        
        if pre_optimize:
            gp.tell(params, signals, likelihood_optimization_method=None, init_hyperparameters=hps_guess,)
            self.msg('gpCAM will use initial guess hyperparameters: {}'.format(repr(hps_guess)), 4, 1)
            self.msg('with log-likelihood: {:,.1f}'.format(gp.gp.log_likelihood(hps_guess)), 4, 2)
            hps_guess = self.optimize_hps(hps_bounds, hps_guess, gp)

        
        self.timing_start()
        gp.tell(params, signals, 
                likelihood_optimization_method=gp_method, # 'global', 'local', None
                init_hyperparameters=hps_guess,
                likelihood_optimization_max_iter=100,
                likelihood_optimization_pop_size=20,
                likelihood_optimization_tolerance=0.0001,
                dask_client=False)
        self.timing_end_msg('gp.tell {} optimization'.format(gp_method), threshold=4, indent=1)
        
        hps = gp.hyperparameters
        self.msg('gpCAM using hyperparameters: {}'.format(repr(hps)), 4, 1)
        self.msg('with log-likelihood: {:,.1f}'.format(gp.gp.log_likelihood(hps)), 4, 2)
        #res = gp.gp.posterior_mean( np.array([[x[i],y[j]]]) )
        #res = gp.gp.posterior_mean( np.array([[x1,y1], [x2,y2], [x3,y3]]) )
        
        
        grid, xi, yi, XI, YI = self.make_grid(**kwargs)
        self.msg("Interpolating {:,} points to {:,}×{:,} = {:,} points (densification = {:.0f})".format(len(self.z_vals), len(xi), len(yi), len(xi)*len(yi), len(xi)*len(yi)/len(self.z_vals)), 3, 1)

        
        ZI = np.zeros_like(XI)
        
        self.msg('Filling ZI grid ({})'.format(fill), 4, 1)
        self.timing_start()
        if fill=='pixelwise':
            for ix, x in enumerate(xi):
                self.timing_progress_msg(ix, len(xi), 4)
                for iy, y in enumerate(yi):
                    res = gp.gp.posterior_mean( np.array([[x,y]]) )
                    ZI[iy,ix] = res['f(x)'][0]
            # 0.6s/2,601pts = 0.2ms/pt
            # 2.2s/10,201pts = 0.2ms/pt
            # 47.6s/252,252pts = 0.2ms/pt
        
        elif fill=='whole':
            points = np.column_stack((np.ravel(XI), np.ravel(YI)))
            res = gp.gp.posterior_mean(points)
            ZI = np.reshape(res['f(x)'], XI.shape)
            # 0.8s/2,601pts = 0.3ms/pt
            # 6.3s/10,201pts = 0.6ms/pt
            # ERR/252,252pts ; uses >470 GiB RAM
            
        self.timing_end_msg('{} fill'.format(fill), threshold=4, indent=2)
                
                

        self.grid = grid
        self.xi, self.yi = xi, yi
        self.XI, self.YI = XI, YI
        self.ZI = ZI


    def load_hps(self, signal_name='model_1', source_dir='./', ext='.npy'):
        # TOCHANGE
        #source_dir = '../../../gpcamv4and5/data/current_data/'
        infile = Path(source_dir, 'hyperparameters_{}{}'.format(signal_name, ext))
        return np.load(infile, allow_pickle=True)
        
        
    def hps_initialize_default(self, params, signals, hps_guess=None):
        
        # Default kernel in 2D uses [[signal variance bounds], [length scale 1 bounds], [length scale 2 bounds]]
        # E.g. for 2D problem: hyperparameter_bounds = np.array([ [0.001,1e9], [0.001,100], [0.001,100] ])
        
        # For bounds, we take the spread in the signal/data and extend it
        #spread = np.square(np.max(signals[:,0])-np.min(signals[:,0])) # NB: Need to square the raw spread since variance = std^2
        spread = np.square(np.std(signals[:,0]))*20
        hps_bounds = np.array([ [spread*1e-3, spread] ])
        for i in range(params.shape[1]):
            spread = abs(np.max(params[:,i])-np.min(params[:,i]))
            hps_bounds = np.append( hps_bounds, [[spread*1e-4, spread*10]], axis=0 )
            
        if hps_guess is None:
            # For initial guess, we could simply use the midpoint of each range
            hps_guess = np.asarray([ 0.5*(lower+upper) for lower, upper in hps_bounds ])
            
            # A better estimate for signal variance, from standard deviation
            hps_guess[0] = np.square(np.std(signals[:,0]))*0.5 # The 0.5 is an empirical fudge-factor
            
            # A better estimate of lengthscales, from decorrelation lengths
            xlength, ylength = self.correlation_lengths_2D(params, signals)
            hps_guess[1] = xlength*0.75 # Empirical fudge-factor
            hps_guess[2] = ylength*0.75 # Empirical fudge-factor

            
        return None, hps_bounds, hps_guess


    def hps_initialize_periodic2Daniso(self, params, signals, hps_guess=None, code_PATH=None):
        
        # TOCHANGE: Make sure kernel definitions are available
        #code_PATH='../../../gpcamv4and5/scripts/'
        if code_PATH is not None:
            code_PATH in sys.path or sys.path.append(code_PATH)
        from kernel_definition import periodic_kernel_2d_isotropic, periodic_kernel_2d_anisotropic
        gp_kernel = periodic_kernel_2d_anisotropic
        
        # [[signal variance bounds], [length scale 1 bounds], [length scale 2 bounds], [intercept bounds], [slope bounds] ]
        # E.g. for 2D problem: hyperparameter_bounds = np.array([ [0.001,1e9], [0.001,100], [0.001,100] ])
        
        # For bounds, we take the spread in the signal/data and extend it
        #spread = np.square(np.max(signals[:,0])-np.min(signals[:,0])) # NB: Need to square the raw spread since variance = std^2
        spread = np.square(np.std(signals[:,0]))*20
        hps_bounds = np.array([ [spread*1e-3, spread] ])
        for i in range(params.shape[1]):
            spread = abs(np.max(params[:,i])-np.min(params[:,i]))
            hps_bounds = np.append( hps_bounds, [[spread*1e-4, spread*10]], axis=0 )

        #slope=1.59483070e-05, intercept=0.295558503        
        hps_bounds = np.append( hps_bounds, [[0.270, 0.33]], axis=0 ) # intercept
        hps_bounds = np.append( hps_bounds, [[-0.0004, +0.0004]], axis=0 ) # slope
            
        if hps_guess is None:
            # For initial guess, we could simply use the midpoint of each range
            hps_guess = np.asarray([ 0.5*(lower+upper) for lower, upper in hps_bounds ])
            
            # A better estimate for signal variance, from standard deviation
            hps_guess[0] = np.square(np.std(signals[:,0]))*0.5 # The 0.5 is an empirical fudge-factor
            
            # A better estimate of lengthscales, from decorrelation lengths
            xlength, ylength = self.correlation_lengths_2D(params, signals)
            hps_guess[1] = xlength*0.75 # Empirical fudge-factor
            hps_guess[2] = ylength*0.75 # Empirical fudge-factor

            
        return gp_kernel, hps_bounds, hps_guess


    def correlation_lengths_2D(self, params, signals, target_decay=0.5):
        '''Estimate the correlation length in the x and y directions.
        Conceptually, we are looking for the average length-scale in 
        the x-direction over which the selected signal de-correlates
        (bears no resemablance to itself).
        We compute this by creating a grid/matrix, and then comparing
        it to itself as we progressivley offset in this direction. By
        summing over this image, we naturally average over the different
        parts of the image.'''
        
        x, y, s = params[:,0], params[:,1], signals[:,0]
        
        # Create a grid version of the data
        grid, xi, yi, XI, YI = self.make_grid( n_grid=int(np.sqrt(len(s)))*2 ) # We densify by 2X2 = 4
        import scipy.interpolate
        ZI = scipy.interpolate.griddata(np.column_stack((x, y)), s, (XI, YI), rescale=True, method='linear')
        ZI = np.ma.masked_where( np.isnan(ZI), ZI)
        
        # Normalize
        ZI -= np.average(ZI)
        ZI /= np.std(ZI)
        
        # Calculate the baseline (maximum) correlation (zero offset)
        Zc = ZI*ZI
        corr_max = np.sum(Zc)/Zc.count() # Should be 1.0?
        
        h, w = ZI.shape
        for ix in range(w):
            # Compute correlation for offset ix by multiplying
            # the original and an offset matrix:
            Zc = ZI*np.roll(ZI, ix, axis=1)
            Zc = Zc[:,ix:] # Select just the overlap area
            corr_rel = ( np.sum(Zc)/Zc.count() )/corr_max
            if corr_rel<=target_decay:
                # We define a correlation length as the first crossing 
                # below the target value
                break

        for iy in range(h):
            Zc = ZI*np.roll(ZI, iy, axis=0)
            Zc = Zc[iy:,:] # Overlap area
            corr_rel = ( np.sum(Zc)/Zc.count() )/corr_max
            if corr_rel<=target_decay:
                break

        # Convert from indices to physical lengthscales
        xlength = (ix/w)*( np.max(xi)-np.min(xi) )
        ylength = (iy/h)*( np.max(yi)-np.min(yi) )
        
        return xlength, ylength


    def optimize_hps(self, hps_ranges, hps_guess, gp, opt_mode='siman', iterations=10, **kwargs):
        
        self.msg('Starting hps optimization', 4, 1)
        self.timing_start()
        
        from SimpleOptimizer import SimpleOptimizer
        SimpOpt = SimpleOptimizer(hps_ranges, initial_guess=hps_guess, evaluation_function=gp.gp.log_likelihood, verbosity=10)
        
        
        #hps, err = SimpOpt.optimize_random(iterations=iterations)
        #hps, err = SimpOpt.optimize_local_step(iterations=iterations)
        #hps, err = SimpOpt.optimize_siman()
        
        try:
            hps, err = getattr(SimpOpt, 'optimize_{}'.format(opt_mode))(iterations=iterations, **kwargs)
        except AttributeError as e:
            self.msg("ERROR: optimize_hps doesn't understand mode '{}'".format(opt_mode), 1, 0)
            raise e
        
        
        self.timing_end_msg('hps optimization', iterations=iterations, threshold=4, indent=1)
        
        return hps_guess


    # Plot
    ########################################
    def plot(self, outfile=None, title=None, faded=0.0, dpi=150, plot_buffers=[0.21,0.12,0.18,0.10], **kwargs):
        self.msg('Plotting {}'.format(self.z_name), 3, 0)
        
        d = Data2DSM()
        d.model = self
        d.data = self.ZI
        d.x_axis = self.xi
        d.y_axis = self.yi
        
        d.x_vals = self.x_vals
        d.y_vals = self.y_vals
        d.z_vals = self.z_vals
        
        d.set_z_display([None, None, 'linear', 1.0])
        
        d.x_rlabel = self.label(self.x_name)
        d.y_rlabel = self.label(self.y_name)
        d.z_rlabel = self.label(self.z_name)
        
        zmin, zmax, cmap = self.zscale(self.z_name)
        
        d.plot_args['rcParams'] = { 
                        'axes.labelsize': 50,
                        'xtick.labelsize': 40,
                        'ytick.labelsize': 40,    
                        }
        
        outfile = self.get_outfile(outfile)
            
        d.plot(save=outfile, show=False, cmap=cmap, zmin=zmin, zmax=zmax, title=title, plot_buffers=plot_buffers, plot_range=self.grid, plot_2D_type='pcolormesh', dpi=dpi, transparent=False, faded=faded, **kwargs)
        
        self.msg('Saved plot as: {}'.format(outfile), 4, 1)
        

    def plot_signal(self, x='x_pos', y='y_pos', signal=None, N_limit=None, **kwargs):
        
        self.select(x=x, y=y, signal=signal, N_limit=N_limit)
        self.interpolate(**kwargs)
        self.plot(**kwargs)
        
            
    def plot3D(self, outfile=None, title=None, faded=0.0, elev=30, azim=-60, dpi=150, plot_buffers=[0.05,0.10,0.05,0.05], **kwargs):
        self.msg('Plotting (3D) {}'.format(self.z_name), 3, 0)
        
        d = Data2DSM()
        d.model = self
        d.data = self.ZI
        d.X = self.XI
        d.Y = self.YI
        
        d.z_vals = self.z_vals
        
        d.x_rlabel = self.label(self.x_name)
        d.y_rlabel = self.label(self.y_name)
        d.z_rlabel = self.label(self.z_name)
        
        
        zmin, zmax, cmap = self.zscale(self.z_name)

        d.plot_args['rcParams'] = { 
                        'axes.labelsize': 40,
                        'xtick.labelsize': 20,
                        'ytick.labelsize': 20,    
                        }    
        
        outfile = self.get_outfile(outfile, extra_dir='3D')
        
        d.plot3D(save=outfile, show=False, cmap=cmap, zmin=zmin, zmax=zmax, title=title, plot_buffers=plot_buffers, plot_range=self.grid, elev=elev, azim=azim, dpi=dpi, transparent=False, **kwargs)        
    


    # Additional
    ########################################
    def get_outfile(self, outfile=None, extra='', show_N=True, subdir=True, extra_dir=None, ext='.png'):
        if outfile is None:
            protocol, signal = self.z_name.split('__')
            outdir = Path(self.output_dir)
            if subdir:
                outdir = outdir.joinpath(signal)
            if extra_dir is not None:
                outdir = outdir.joinpath(extra_dir)
            outdir.mkdir(parents=True, exist_ok=True)
            N = 'N{:04d}'.format(self.get_N()) if show_N else ''
            outfile = Path(outdir, '{}-{}-{}{}{}'.format(self.sample, signal, extra, N, ext))

        return outfile
        
        
    def copy_current(self, outfile=None, copy_to=None):
        '''Copies the most recently-created plot to another location.
        The intention is for this copied file to act as an updating
        status of the experiment.
        NOTE: This code is highly contingent. Various paths are hard-
        coded, and will need to be changed for a new setup.'''
        
        # NB: This code used to be called "status_online"
        
        outfile = self.get_outfile(outfile)
            
        if copy_to is None:
            # TOCHANGE: Add a valid path on the current machine
            search = [
                #'/home/kyager/Desktop/',
                '/home/kyager/software/statpage/',
                '/home/xf12id/software/statpage/',
                ]
            for s in search:
                copydir = Path(s)
                if copydir.is_dir():
                    copy_to = copydir.joinpath('current.png')
                    break
            if copy_to is None:
                self.msg("ERROR: copy_current couldn't find a valid path.",1,0)
                return
        
        import shutil
        shutil.copyfile(outfile, copy_to)
    
    
    
    ########################################
    # End: class SurrogateModel(Base)



# Animation(Base)
########################################
class Animation(Base):
    
    def __init__(self, model, name='anim', verbosity=3, **kwargs):
        
        self.kwargs = kwargs
        self.model = model
        self.name = name
        self.verbosity = verbosity


    def animation(self, anim_mode='gif', **kwargs):
        
        # Make the frames
        N_list = self.get_N_list(**kwargs)
        outfiles = self.plot_sequence(N_list, **kwargs)

        # Where to save?
        exts = {'gif':'.gif', 'mpeg':'.mp4'}
        outfile = self.model.get_outfile(extra='anim', show_N=False, subdir=False, ext=exts[anim_mode])

        try:
            getattr(self, 'animation_{}'.format(anim_mode))(outfile, outfiles, **kwargs)
        except AttributeError as e:
            self.msg("ERROR: animation doesn't understand mode '{}'".format(anim_mode), 1, 0)
            raise e
        
        
    def animation_gif(self, outfile, file_list, **kwargs):
        '''Generate an animated GIF based on the plotted frames.
        Users imagemagick on the commandline.'''
        # If you get a 'cache resources exhausted' error, you can increase the cache sizes:
        # sudo nano /etc/ImageMagick-6/policy.xml
        
        self.msg("Generating animated gif: {}".format(outfile), 3, 0)


        # Prepare command
        # (Animation is generated using imagemagick 'convert' bash command.)

        # -loop 0 makes the GIF loop forever; -loop 1 just plays once
        #cmd = "convert -delay 20 -resize 50% -fill white  -undercolor '#00000080'  -gravity NorthWest -annotate +0+5 ' Text ' "
        #cmd = "convert -delay 15 -loop 1 -resize 50% "
        #cmd = "convert -crop 450x450+60+220  +repage -delay 15 -loop 1 -resize 50% "
        cmd = "convert -dispose previous +repage -delay 30 -loop 1 -resize 30% "
                
        for infile in file_list:
            cmd += '{} '.format(infile)
        
        cmd += ' {}'.format(outfile)

        # Execute command
        os.system(cmd)
        
        # Add a white background
        #cmd = 'convert {} -coalesce -background white -alpha remove {}'.format(outfile, outfile[:-4]+'w.gif')
        #os.system(cmd)        


    def animation_mpeg(self, outfile, file_list, temp_file='aframe-Ahj5vaqu-', v_size=720, **kwargs):
        '''Generates a movie based on plotted frames.
        Uses ffmpeg on the commandline.'''
        
        self.msg("Generating mpeg movie: {}".format(outfile), 3, 0)
        
        # Prepare command
        # (Movie is generated using ffmpeg.)
        
        # We need a command like:
        #ffmpeg -r 10 -f image2 -i aframe%04d.png -filter:v scale=720:-1 -vcodec libx264 -crf 25  -pix_fmt yuv420p output.mp4
        # Where the frames are sequential. As a hack, we create symlinks:
        self.msg("Generating {} temporary frames".format(len(file_list)), 4, 1)
        for i, infile in enumerate(sorted(file_list)):
            linkfile = outfile.parent.joinpath( '{}{:04d}{}'.format(temp_file, i, infile.suffix))
            os.symlink(infile.resolve(), linkfile)
        
        file_pattern = '{}/{}%04d{}'.format(outfile.parent, temp_file, file_list[0].suffix)
        cmd = 'ffmpeg -r 10 -f image2 -i {} -filter:v scale={}:-1 -vcodec libx264 -crf 25  -pix_fmt yuv420p {}'.format(file_pattern, v_size, outfile)
        os.system(cmd)
        
        # Remove temporary files
        self.msg("Removing {} temporary frames".format(len(file_list)), 4, 1)
        for i, infile in enumerate(sorted(file_list)):
            linkfile = outfile.parent.joinpath( '{}{:04d}{}'.format(temp_file, i, infile.suffix))
            if linkfile.is_symlink():
                os.remove(linkfile)

    def plot_sequence(self, N_list, force=False, **kwargs):
        
        N_max = self.model.get_N_max()
        self.msg("Plotting sequence with {} frames ({:.1f}% of available)".format(len(N_list), 100.*len(N_list)/N_max), 3, 0)
        self.msg("N_list = {} ".format(N_list), 6, 1)
        
        outfiles = []
        for i, N in enumerate(N_list):
            self.model.change_N_limit(N)
            outfile = self.model.get_outfile()
            outfiles.append(outfile)
            if force or not outfile.exists():
                
                self.msgm("Plotting N {}/{} = {:.1f}% (frame {}/{} = {:.1f}%)".format(N, N_max, 100.*N/N_max, i+1, len(N_list), 100.*(i+1)/len(N_list)), 2, 1)
                
                self.model.interpolate(**kwargs)
                self.model.plot(**kwargs)
                
        return outfiles
        
        
    def get_N_list(self, spacing='power', N_total=None, N_spacing=None, exponent=5.0, **kwargs):
        
        N_length = len(self.model.z_vals)
        
        if spacing=='linear':
            if N_total is not None:
                N_spacing = int(N_length/N_total)
            elif N_spacing is None:
                N_spacing = 50
            N_list = np.arange(N_spacing, N_length, N_spacing)
        elif spacing=='power':
            if N_total is None:
                N_total = 140
            N_list = self.power_N_list(N_length, num=N_total, exponent=exponent)
        else:
            self.msg("ERROR: get_N_list doesn't recognize spacing='{}'".format(spacing), 1, 0)
            return None
        
        return N_list
    
            
    def power_N_list(self, N_max, N_min=3, num=40, exponent=3.0):
        '''Generates a list of integers that are spread out more logarithmically.
        That is, the list starts with small increments between integers, but finishes
        with large increments.'''
        
        #N_list = ( (np.exp( np.linspace(0, 1, num=40) ) - 1)/(np.exp(1)-1) )*len(z_vals)
        x = np.linspace(0, 1, num=num)
        N_list = np.power(x, exponent)*(N_max-N_min) + N_min
        N_list = np.unique(N_list.astype(int))
        #N_list = N_list[ (N_list>=N_min) & (N_list<=N_max) ] # Unnecessary
        
        return N_list
            
            

    ########################################
    # End: class Animation(Base)



# Tasker(Base)
########################################
class Tasker(Base):
    
    def __init__(self, name='task', distributed=False, processes=True, threads_per_worker=1, n_workers=20, verbosity=3, **kwargs):
        
        self.kwargs = kwargs
        self.name = name
        self.distributed = distributed
        self.verbosity = verbosity
        
        if self.distributed:
            from dask.distributed import Client, progress
            self.client = Client(processes=processes, threads_per_worker=threads_per_worker, n_workers=n_workers)
            

    def load_data(self, model, extractions, **kwargs):
        
        # Extract and save data
        
        #model.load_xml(extractions, use_cached=False)
        model.load_sql(extractions)
        
        model.sort_data()
        model.save_data()
        

    def prepare_data(self, model, **kwargs):
        
        # Prepare data
        model.load_data()
        #model.print_data()
        model.coord_transform_vy()
        #model.enforce_periodic()                     
        
        
    def plot_signals(self, model, signals, **kwargs):
        
        orig_name = model.name
        if self.distributed:
            futures = []
            
        for signal, zmin, zmax in signals:
            model.name = '{}-{}'.format(orig_name, signal)
            kwargs['signal'] = signal
            if self.distributed:
                self.msg('Launching dask job for signal {}'.format(signal), 3, 0)
                future = self.client.submit(model.plot_signal, **kwargs)
                futures.append(future)
            else:
                self.msg('Running plot_signal for {}'.format(signal), 3, 0)
                model.plot_signal(**kwargs)
                #if signal=='prefactor1': # TOCHANGE
                    #model.copy_current()
                
        
        if self.distributed:
            self.msg('Telling dask to gather results', 4, 1)
            results = self.client.gather(futures)
        
        model.name = orig_name


    def plot_sequence(self, model, force=False, x='x_pos', y='y_pos', signal=None, **kwargs):
        
        model.select(x=x, y=y, signal=signal)
        
        anim = Animation(model, verbosity=self.verbosity)
        N_list = anim.get_N_list(**kwargs)
        
        N_max = model.get_N_max()
        self.msg("Plotting sequence with {} frames ({:.1f}% of available)".format(len(N_list), 100.*len(N_list)/N_max), 3, 0)
        self.msg("N_list = {} ".format(N_list), 6, 1)
        
        
        orig_name = model.name
        if self.distributed:
            futures = []
            kwargs.update({'x':x, 'y':y, 'signal':signal})
        
        #outfiles = anim.plot_sequence(N_list, **kwargs)
        outfiles = []
        for i, N in enumerate(N_list):
            model.name = '{}-{}'.format(orig_name, N)
            model.change_N_limit(N)
            outfile = model.get_outfile()
            outfiles.append(outfile)
            if force or not outfile.exists():
                
                self.msgm("Plotting N {}/{} = {:.1f}% (frame {}/{} = {:.1f}%)".format(N, N_max, 100.*N/N_max, i+1, len(N_list), 100.*(i+1)/len(N_list)), 2, 1)

                if self.distributed:
                    self.msg('Launching dask job for N = {}'.format(N), 3, 0)
                    kwargs['N_limit'] = N
                    future = self.client.submit(model.plot_signal, **kwargs)
                    futures.append(future)
                else:
                    model.interpolate(**kwargs)
                    model.plot(**kwargs)

        if self.distributed:
            self.msg('Telling dask to gather results', 4, 1)
            results = self.client.gather(futures)

        model.name = orig_name
                
        return outfiles


    def animation(self, model, anim_mode='mpeg', force=False, **kwargs):
        
        # Create the frames
        self.plot_sequence(model, force=force, **kwargs)
        
        # Animation
        anim = Animation(model, verbosity=self.verbosity)
        anim.animation(anim_mode=anim_mode, force=False, **kwargs)
        

    ########################################
    # End: class Tasker(Base)









# RUN
########################################
# TOCHANGE: Update information below for the current experiment


verbosity = 5
model = SurrogateModel(sample='AE_PLA_acc25_run1_SM', sample_pattern='AE_PLA_acc25_run1_', source_dir='../', verbosity=verbosity)

# Signals to extract from SciAnalysis output
extractions = [ 
    ['metadata_extract', 
        [
        'x_position', 
        'y_position', 
        'sequence_ID',
        ]
    ],
    #['linecut_angle_fit',
        #[
        #'fit_eta_eta', 
        #'orientation_factor', 
        #'orientation_angle', 
        #'fit_eta_span_prefactor',
        #]
    ['circular_average_q2I_fit', 
        [
        'fit_peaks_prefactor1', 
        'fit_peaks_x_center1', 
        'fit_peaks_sigma1', 
        'fit_peaks_chi_squared', 
        'fit_peaks_d0', 
        'fit_peaks_grain_size',
        ]
    ],
    ]


    
    
    
# Define plotting behavior
signals = [
    ['prefactor1', 1e3, 6e3],
    ['grain_size', 32, 39],
    ['d0', 0.537, 0.539],
    ]
model.add_zscales(signals)

grid, n_grid = [None,None,None,None], 200 # Typical
#grid, n_grid = [0,65,-0.45,0.75], [1000,250]
figsize, plot_buffers = (10,10), [0.21,0.12,0.18,0.10] # Typical square
#figsize, plot_buffers = (20,10), [0.12,0.06,0.18,0.10] # Extended rectangle



# Define hyperparameters (for gpCAM)
hps_guess = None

# Saved by running gpCAM instance
#hps_guess = model.load_hps(signal_name='prefactor', source_dir = '../../../gpcamv4and5/data/current_data/')

# Manually defined
#hps_guess = [2764, 10, 0.1] # Example for a default anisotropic kernel
#hps_guess = [551116.4, 8.0, 0.094, 0.296, 0] # periodic_kernel_2d_anisotropic





if __name__ == '__main__':
    
    task = Tasker(distributed=False, verbosity=verbosity)
    
    #task.load_data(model, extractions)
    task.prepare_data(model)
    
    kwargs = {
        'x': 'v_print' ,
        'y': 'y_corrected' ,
        'interp_mode': 'griddata' , # 'griddata', 'gpcam'
        'grid': grid , 
        'n_grid': n_grid , 
        'hps_guess': hps_guess ,
        'gp_method': None , # 'global', 'local', None
        'pre_optimize': False ,
        'figsize': figsize ,
        'plot_buffers': plot_buffers ,
        }
    
    task.plot_signals(model, signals, **kwargs)
    
    kwargs['signal'] = 'd0'
    #task.plot_sequence(model, force=False, **kwargs)
    #task.animation(model, **kwargs)
    
    
