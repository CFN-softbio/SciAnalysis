#!/usr/bin/python3
# -*- coding: utf-8 -*-

# This generates a map (2D false-color plot or 3D height plot) for a set of
# experiments (that are presumptively defined in some (x,y) space). The code
# assumes you've already used SciAnalysis to process your data, such that you
# have XML files in your "results" sub-folder with the analysis results of
# interest. This code then compiles that data and generates the plot.

# The data can be interpolated in a variety of ways to yield a "surrogate
# model". For instance, a naive linear interpolation can be used, or 
# gpCAM (v7) can be invoked with a physics-aware kernel.

# The code can also be used to generate an animation of the sequence of
# measurements during the experiment.




# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
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
        #self.ax.scatter(self.x_vals, self.y_vals, s=80, c=self.z_vals, cmap=cmap, vmin=zmin, vmax=zmax, edgecolor='k', linewidth=0.4, zorder=100)
        
        
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


        self.ax.tick_params('both', which='major', length=10)
        
        #self.ax.axhline(0, linewidth=1, color='g', alpha=0.75)
        #self._plot_guides(**plot_args)
        #self._plot_guides_SM(**plot_args)
        #self._plot_guidebox(**plot_args)

        if 'aspect' in plot_args and plot_args['aspect'] is not None and plot_args['aspect'] is not False:
            self.ax.set_aspect('equal', 'box')
            # How to set? 'auto', 'equal', num (force ratio)
            # What to adjust? None, 'box', 'datalim'
        
        
    def _plot_guides(self, **plot_args):
        '''Example of overlaying some guide-lines'''
        
        xi, xf, yi, yf = self.ax.axis()

        intercept = 20
        slope = (130-intercept)/(6.88708-0)
        #lam = slope*x + intercept
        
        for lam in [50, 100]:
            x = (lam-intercept)/slope
            self.ax.axvline(x, color='k', linewidth=2, alpha=0.5)
            s = '$\Lambda = {:.0f} \, \mathrm{{nm}}$'.format(lam)
            self.ax.text(x, yf, s, size=30, horizontalalignment='center', verticalalignment='bottom')
        
        #self.ax2 = self.ax.twiny()
        #self.ax2.set_xlim(self.ax.get_xlim())
        #self.ax2.set_xlabel('$\Lambda \, (\mathrm{mm})$')
        
    def _plot_guides_SM(self, **plot_args):
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
        
    def _plot_guidebox(self, **plot_args):
        '''Example of overlaying some guide-box.'''
        origin, w, h = (-3.37, 0.57), 6.88708, 3
        
        rect = mpl.patches.Rectangle(origin, w, -h, linewidth=5, edgecolor='r', facecolor='none', alpha=0.5, zorder=20)
        rect.set_linestyle( (0, [2,4]) )
        self.ax.add_patch(rect)
        
        origin, w, n = -1.63, 0.5, 5
        for i in range(n):
            spacing = w/(n-1)
            x = origin - spacing*(n-1)/2 + i*spacing
            lw = 3*(1 - abs(x-origin)/w)
            self.ax.axvline(x, color='r', linewidth=lw, alpha=0.5, dashes=[2,4])
        
        
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



# DataLineSM(DataLine)
########################################
class DataLineSM(DataLine):
    
    def _plot_main(self, error=False, error_band=False, dashes=None, **plot_args):
        
        zmin, zmax = self.plot_range[2], self.plot_range[3]
        cmap = mpl.cm.get_cmap(self.cmap)
        
        for xc, yc in zip(self.x, self.y):
            yr = (yc-zmin)/(zmax-zmin)
            self.ax.scatter(xc, yc, s=60, color=cmap(yr), linewidths=0.5, edgecolors='k')
    
    
    def __plot_extra(self, **plot_args):
        '''Example of adding some guide-lines.'''
        
        xi, xf, yi, yf = self.ax.axis()
        x = np.linspace(xi, xf, endpoint=True, num=2000)
        y = x
        self.ax.plot(x, y, color='k', linewidth=2, dashes=[5,5], alpha=0.5)
        idx = np.where(y>=yf)[0][0]
        self.ax.text(x[idx], y[idx-2], '$1:1$', size=20, color='k', alpha=0.5, horizontalalignment='left', verticalalignment='top')
        
        y = x/2
        self.ax.plot(x, y, color='k', linewidth=2, dashes=[5,5], alpha=0.5)
        idx = np.where(y>=yf)[0][0]
        self.ax.text(x[idx], y[idx-2], '$1:2$', size=20, color='k', alpha=0.5, horizontalalignment='left', verticalalignment='top')
        
        #self.ax.set_aspect('equal', 'box')
    

    ########################################
    # End: class DataLineSM(DataLine)



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

    def msgm(self, txt=None, threshold=0, indent=0, indent_txt='  ', verbosity=None, mark='=', nmark=40, empty_lines=1, **kwargs):
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
        
        self.trim_list = None # Optional list of points to exclude
        
        # Selected values just before interpolating
        self.x_vals, self.y_vals, self.z_vals = None, None, None
        self.x_name, self.y_name, self.z_name = None, None, None
        self.z_errors = None
        
        
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
            [ 'x_position', '$x \, (\mathrm{mm})$', None, None, None ] ,
            [ 'y_position', '$y \, (\mathrm{mm})$', None, None, None ] ,
            [ 'sequence_ID', '$\mathrm{sID}$', 0, None, None ] ,
            [ 'prefactor', r'$p \, (\mathrm{a.u.})$', 0, None, cmap_vge ],
            [ 'chi_squared', '$\chi^2\, (\mathrm{a.u.})$', 0, None, 'plasma' ],            
            [ 'fit_peaks_x_center1', '$q_0 \, (\mathrm{\AA^{-1}})$', None, None, 'viridis' ] ,
            [ 'fit_peaks_d0', '$d_0 \, (\mathrm{nm})$', None, None, 'viridis' ] ,
            [ 'fit_peaks_sigma1', r'$\sigma_0 \, (\mathrm{\AA^{-1}})$', None, None, 'inferno' ] ,
            [ 'fit_peaks_grain_size', r'$\xi \, (\mathrm{nm})$', 0, None, 'inferno' ] ,
            [ 'fit_peaks_prefactor1', '$p \, (\mathrm{a.u.})$', 0, None, cmap_vge ] ,
            [ 'fit_eta_eta', r'$\eta$', 0, 1, 'inferno' ],
            [ 'fit_eta_span_eta', r'$\eta$', 0, 1, 'inferno' ],
            [ 'fit_MaierSaupe_m', r'$m_{\mathrm{MS}}$', 0, None, 'inferno' ],
            [ 'orientation_factor', r'$f_{\mathrm{ori}}$', -1, +1, 'gray' ],
            [ 'orientation_angle', r'$\chi \, (\mathrm{^{\circ}})$', -90, +90, cmap_cyclic_rb ],
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
        
        
    def zscale(self, name, ztrim=[0.02,0.02]):
        
        label, zmin, zmax, cmap = self.label(name, retall=True)
        
        if cmap is None: cmap = 'jet'
        
        for match, zmin_c, zmax_c in reversed(self.zscales):
            if match in name:
                if zmin_c is not None: zmin = zmin_c
                if zmax_c is not None: zmax = zmax_c
        
        if ztrim is not None:
            # Trim off the specified amount on either end of the histogram of values
            values = np.sort(self.z_vals)
            if zmin is None: zmin = values[ +int( len(values)*ztrim[0] ) ]
            if zmax is None:
                idx = -int( len(values)*ztrim[1] )
                if idx>=0:
                    idx = -1
                zmax = values[idx]
            if zmax<=zmin:
                zmax = max(values)

        else:
            # Set the limits to the min/max of the data
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

    def save_data(self, outfile=None, extra=''):
        '''Save the model data list to a npy file.'''
        if outfile is None:
            outfile = Path(self.output_dir, '{}-{}{}.npy'.format(self.name, self.sample, extra))
        np.save(outfile, self.data, allow_pickle=True)
        
        
    def load_data(self, infile=None, extra=''):
        '''Load data from npy file.
        Data should be formatted exactly as SurrogateModel expects it
        (e.g. from model.save_data()).'''
        if infile is None:
            infile = Path(self.output_dir, '{}-{}{}.npy'.format(self.name, self.sample, extra))
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
                if self.verbosity>=6:
                    # Output the keys that are being missed
                    for protocol, signals in extractions:
                        if protocol in result.keys():
                            for signal in signals:
                                if signal not in result[protocol]:
                                    self.msg("for protocol {}, missing signal: {}".format(protocol, signal), 6, 4)
                        else:
                            self.msg("missing protocol: {}".format(protocol), 6, 4)
                
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
        

    def select(self, x='x_pos', y='y_pos', signal=None, N_limit=None, e_signature='_error'):
        '''Decide what parameters to use for x, y and z/signal in subsequent plotting.'''
        
        self.N_limit = N_limit
        
        if signal is None:
            signal = list(self.data['signals'].keys())[0]

        self.msg("Selecting parameters/signals, based on x='{}', y='{}', signal='{}'".format(x, y, signal), 4, 0)
            
        # Get x values
        for name, data in self.data['parameters'].items():
            if x==name:
                self.msg("Using exact match '{}' ({} values) as x-coordinate".format(name, len(data)), 4, 1)
                self.x_name = name
                self.x_vals = self.data['parameters'][name]
                break
            if x in name:
                self.msg("Using '{}' ({} values) as x-coordinate".format(name, len(data)), 4, 1)
                self.x_name = name
                self.x_vals = self.data['parameters'][name]
            
        # Get y values
        for name, data in self.data['parameters'].items():
            if y==name:
                self.msg("Using exact match '{}' ({} values) as y-coordinate".format(name, len(data)), 4, 1)
                self.y_name = name
                self.y_vals = self.data['parameters'][name]
                break
            if y in name:
                self.msg("Using '{}' ({} values) as y-coordinate".format(name, len(data)), 4, 1)
                self.y_name = name
                self.y_vals = self.data['parameters'][name]
                
        # Get z values
        for name, data in self.data['signals'].items():
            if signal==name:
                self.msg("Using exact match '{}' ({} values) as z-values".format(name, len(data)), 4, 1)
                self.z_name = name
                self.z_vals = self.data['signals'][name]
                break
            if signal in name:
                self.msg("Using '{}' ({} values) as z-values".format(name, len(data)), 4, 1)
                self.z_name = name
                self.z_vals = self.data['signals'][name]
                
        
        # Avoid accidentally selecting an error
        if self.z_name.endswith(e_signature):
            name = self.z_name[:-len(e_signature)]
            if name in self.data['signals'].keys():
                self.msg("Using better match '{}' ({} values) as z-values".format(name, len(data)), 4, 1)
                self.z_name = name
                self.z_vals = self.data['signals'][name]
            
                
        # Get corresponding z errors
        if '{}{}'.format(self.z_name, e_signature) in self.data['signals'].keys():
            self.z_errors = self.data['signals']['{}{}'.format(self.z_name, e_signature)]
        else:
            self.z_errors = None


        if self.trim_list is not None:
            idx = self.trim_list
            self.x_vals, self.y_vals, self.z_vals = self.x_vals[idx], self.y_vals[idx], self.z_vals[idx]
            if self.z_errors is not None:
                self.z_errors = self.z_errors[idx]

        N = self.get_N()
        if self.N_limit is not None:
            self.x_vals, self.y_vals, self.z_vals = self.x_vals[:N], self.y_vals[:N], self.z_vals[:N]
            if self.z_errors is not None:
                self.z_errors = self.z_errors[:N]
        
        
        nans = self.count_nans(self.z_vals)
        if nans>0:
            self.msg("NOTE: selected signal includes NaN for {}/{} = {:.1f}% of the entries".format(nans, N, 100.*nans/N), 2, 2)
        nans = self.count_nans(self.z_errors)
        if nans>0:
            self.msg("NOTE: selected error includes NaN for {}/{} = {:.1f}% of the entries".format(nans, N, 100.*nans/N), 2, 2)


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
    
    def coord_transform_rotate(self, angle=0.0, origin=(0,0)):
        '''Rotate an (x,y) space by the specified angle.'''
        
        x_vals = self.select_single('x_position')
        y_vals = self.select_single('y_position')
        
        xo, yo = origin
        
        c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
        x_rotated = +c*(x_vals-xo) + -s*(y_vals-yo) + xo
        y_rotated = +s*(x_vals-xo) + +c*(y_vals-yo) + yo

        self.data['parameters']['x_rotated'] = x_rotated
        self.label_associations.append(['x_rotated', '$x_{\mathrm{r}} \, (\mathrm{mm})$', None, None, None])
        self.data['parameters']['y_rotated'] = y_rotated
        self.label_associations.append(['y_rotated', '$y_{\mathrm{r}} \, (\mathrm{mm})$', None, None, None])

    
    def coord_transform_translate(self, origin=(0, 0)):
        
        x_vals = self.select_single('x_position')
        y_vals = self.select_single('y_position')
        
        xo, yo = origin
        
        x_corrected = x_vals-xo
        y_corrected = y_vals-yo

        self.data['parameters']['x_corrected'] = x_corrected
        self.label_associations.append(['x_corrected', '$x_{\mathrm{c}} \, (\mathrm{mm})$', None, None, None])
        self.data['parameters']['y_corrected'] = y_corrected
        self.label_associations.append(['y_corrected', '$y_{\mathrm{c}} \, (\mathrm{mm})$', None, None, None])



    # Data transformations
    ########################################
    
    def add_trim(self, points_to_keep):
        '''Add an optional mask of points to keep or trim. The mask
        is a simple boolean array. True elements are retained, False
        are trimmed/removed.'''
        
        self.msg('trimming data from list:', 6, 4)
        
        
        if self.trim_list is None:
            # Initialize
            param = list(self.data['parameters'].values())[0] # Example
            self.trim_list = np.ones_like(param).astype(bool)
        
        count = np.count_nonzero(self.trim_list)
        self.msg('trim_list {}/{} = {:.1f}% True'.format(count, len(self.trim_list), 100.*count/len(self.trim_list)), 6, 5)
        
        count = np.count_nonzero(points_to_keep)
        self.msg('points_to_keep {}/{} = {:.1f}% True'.format(count, len(points_to_keep), 100.*count/len(points_to_keep)), 6, 5)
        
        self.trim_list = np.logical_and(self.trim_list, points_to_keep)

        count = np.count_nonzero(self.trim_list)
        self.msg('trim_list {}/{} = {:.1f}% True'.format(count, len(self.trim_list), 100.*count/len(self.trim_list)), 6, 5)
        
        

    def trim_to_inside_hull(self):
        # Exclude (x,y) points outside the convex hull where we have conversions defined
        
        self.msg('Trimming to restrict to convex hull', 4, 2)
        
        x_vals = self.select_single('x_corrected')
        y_vals = self.select_single('y_corrected')
        
        test_points = np.stack((x_vals,y_vals), axis=1)
        inside = self.limit_hull.find_simplex(test_points)>=0
        self.msg("    {:d}/{:d} = {:.1f}% points inside convex hull".format(np.sum(inside), len(x_vals), 100*np.sum(inside)/len(x_vals)), 4, 3)
        
        self.add_trim(inside)
        
        
    def trim_coord(self, **axes):
        '''Only include points within the given ranges of the gives axes.'''
        
        for axis_name, (a_min, a_max) in axes.items():
            self.msg('Trimming axis {}, restricting from {:.3g} to {:.3g}'.format(axis_name, a_min, a_max), 4, 2)
            vals = self.select_single(axis_name)
            points_to_keep = np.logical_and( vals>=a_min, vals<=a_max )
            self.add_trim(points_to_keep)


    def trim_signal(self, **axes):
        '''Only include points within the given signal range.'''
        
        for axis_name, (a_min, a_max) in axes.items():
            self.msg('Trimming signal {}, restricting from {:.3g} to {:.3g}'.format(axis_name, a_min, a_max), 4, 2)
            vals = self.select_single(axis_name, order=['signals'])
            points_to_keep = np.logical_and( vals>=a_min, vals<=a_max )
            self.add_trim(points_to_keep)


    def count_nans(self, vals):
        if vals is None:
            return 0
        return np.count_nonzero(np.isnan(vals))


    def fix_nans(self, nan_error=10, e_signature='_error'):
        
        for name, data in self.data['signals'].items():
            nans = self.count_nans(data)
            if nans>0:
                self.msg('Fixing NaNs: {}/{} = {:.1f}% in {}'.format(nans, len(data), 100*nans/len(data), name), 2, 2)
            data = np.nan_to_num(data)
            
            name_error = '{}{}'.format(name, e_signature)
            if name_error in self.data['signals'].keys():
                error = self.data['signals'][name_error]
                nans = self.count_nans(error)
                if nans>0:
                    self.msg('Fixing NaNs: {}/{} = {:.1f}% in {}'.format(nans, len(error), 100*nans/len(error), name_error), 2, 2)
                    max_error_rel = np.max(np.nan_to_num(error)/np.abs(data))
                    idx = np.isnan(error)
                    error[idx] = np.abs(data[idx])*max_error_rel
                    self.msg('NaNs set to {:.4g}×{:.3g} = {:.3g} (relative)'.format(nan_error, max_error_rel, nan_error*max_error_rel), 6, 3)                        
                
                    
        
    
    def compute_missing_errors(self, nan_error=10, e_signature='_error'):

        # We will be changing the signals dictionary during iteration,
        # so we use list() to unpack the current version here.
        for name, data in list(self.data['signals'].items()):
            
            if name[-len(e_signature):]!=e_signature:
            
                name_error = '{}{}'.format(name, e_signature)
                
                if name_error not in self.data['signals'].keys():
                    # Compute the error based on the related "i_name" signal
                    error = None
                    
                    if 'fit_peaks_d0' in name:
                        i_name = name.replace('fit_peaks_d0', 'fit_peaks_x_center1')
                        i_data = self.data['signals'][i_name]
                        i_error = self.data['signals']['{}{}'.format(i_name, e_signature)]
                        error = (data/i_data)*i_error
                        
                    elif 'fit_peaks_grain_size' in name:
                        i_name = name.replace('fit_peaks_grain_size', 'fit_peaks_sigma1')
                        i_data = self.data['signals'][i_name]
                        i_error = self.data['signals']['{}{}'.format(i_name, e_signature)]
                        error = (data/i_data)*i_error
                        
                    if error is not None:
                        nans = self.count_nans(error)
                        if nans>0:
                            self.msg('WARNING: {}/{} = {:.1f}% are NaN in computed errors for {}'.format(nans, len(error), 100*nans/len(error), name_error), 2, 2)
                            max_error_rel = np.max(np.nan_to_num(error)/np.abs(data))
                            idx = np.isnan(error)
                            error[idx] = np.abs(data[idx])*max_error_rel
                            self.msg('NaNs set to {:.4g}×{:.3g} = {:.3g} (relative)'.format(nan_error, max_error_rel, nan_error*max_error_rel), 6, 3)                        

                        self.data['signals'][name_error] = error
                        
            
    def find_point(self, **axes):
        
        elements = [ '{}={:.3g}'.format(axis_name, position) for axis_name, position in axes.items() ]
        position_str = ', '.join(elements)
        self.msg('Finding point close to: ({})'.format(position_str), 3, 1)
        
        distances = None
        for i, (axis_name, position) in enumerate(axes.items()):
            
            if axis_name not in self.data['parameters'].keys():
                self.msg('ERROR: Axis {} not recognized.'.format(axis_name), 2, 2)
                return
            
            data = self.data['parameters'][axis_name]
            
            if distances is None:
                # Initialize
                N = len(data)
                distances = np.zeros( (len(axes),N) )
                
            print(axis_name, position)
            distances[i] = np.square(data - position)
        
        # Compute overall distance
        distances = np.sum(distances, axis=0)
        distances = np.sqrt(distances)
        idx = np.argmin(distances)
        
        self.msg('Point #{}/{} (d = {:.4g})'.format(idx, len(data), distances[idx]), 3, 2)
        
        for axis_name, data  in self.data['parameters'].items():
            self.msg('{} = {:.5g}'.format(axis_name, data[idx]), 3, 3)
        


    # Interpolate
    ########################################
    
    def interpolate(self, interp_mode='griddata', **kwargs):
        
        self.msg('Interpolating using method: {}'.format(interp_mode), 3, 0)
        
        if self.verbosity>=4:
            tools.val_stats(self.z_vals, name='z_vals')
        nans = self.count_nans(self.z_vals)
        if nans>0:
            self.msg('WARNING: z_vals contains {}/{} = {:.1f}% NaNs'.format(nans, len(self.z_vals), 100*nans/len(self.z_vals)), 2, 1)
        nans = self.count_nans(self.z_errors)
        if nans>0:
            self.msg('WARNING: z_errors contains {}/{} = {:.1f}% NaNs'.format(nans, len(self.z_errors), 100*nans/len(self.z_errors)), 2, 1)
        
        getattr(self, 'interpolate_{}'.format(interp_mode))(**kwargs)

        if self.verbosity>=4:
            tools.val_stats(self.ZI, name='ZI')
    
    
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
        
        

    def interpolate_gpcam(self, hps_guess=None, gp_method='global', pre_optimize=False, fill='pixelwise', gpcam_PATH=None, error_relative=None, renormalize_signals=True, hps_lock=None, convex_clip=False, **kwargs):
        
        # gpCAM code:
        # https://bitbucket.org/MarcusMichaelNoack/gpcam.git
        
        # TOCHANGE: Make sure gpCAM is available
        #gpcam_PATH='/home/kyager/current/code/gpcam/main/'
        #gpcam_PATH='../../../gpcamv4and5/gpcam/'
        if gpcam_PATH is not None:
            gpcam_PATH in sys.path or sys.path.append(gpcam_PATH)
        from gpcam.gp_optimizer import GPOptimizer
        
        params = np.stack( (self.x_vals, self.y_vals), axis=1 )
        signals = np.stack( (self.z_vals, ), axis=1 )
        
        if self.z_errors is None:
            if error_relative is None:
                variances = None
            else:
                variances = np.stack( (np.square(np.abs(self.z_vals)*error_relative), ), axis=1 )
        else:
            variances = np.stack( (np.square(self.z_errors), ), axis=1 )
            
        if renormalize_signals:
            signal = signals[:,0]
            avg, std = np.average(signal), np.std(signal)
            signals[:,0] = (signal-avg)/std
            
            if variances is not None:
                variances[:,0] /= np.square(std)
        
        # input_space_bounds is the valid range for each parameter
        # We set this to the actual min/max range of the parameters.
        input_space_bounds = np.asarray([ [np.min(col), np.max(col)] for col in params.T ])

        # Hyperparameters depend on the kernel definition
        gp_kernel, hps_bounds, hps_guess = self.hps_initialize_default(params, signals, hps_guess=hps_guess)
        #gp_kernel, hps_bounds, hps_guess = self.hps_initialize_periodic2Daniso(params, signals, hps_guess=hps_guess)
        
        if hps_lock is not None:
            # Restrict tuning of some of the hyperparameters
            for i, lock in enumerate(hps_lock):
                if lock:
                    hps_bounds[i] = [hps_guess[i], hps_guess[i]]
        
        gp = GPOptimizer(
                input_space_dimension=params.shape[1], # Number of parameters (dimensionality of the search space)
                input_space_bounds=input_space_bounds, # Valid range for each parameter
                )
        

        gp.tell(params, signals, variances) # Load the experimental data
        gp.init_gp(hps_guess, compute_device='cpu', gp_kernel_function=gp_kernel)
        
        
        if pre_optimize:
            # Use a crude method to roughly guess the hyperparameters
            hps_guess = self.optimize_hps(hps_bounds, hps_guess, gp)
            
        if gp_method is not None:
            self.timing_start()
            gp.train_gp(hps_bounds,
                        method=gp_method, # 'global', 'local', 'hgdl'
                        pop_size=20,
                        tolerance=1e-4,
                        max_iter=100,
                        )
            self.timing_end_msg('gp.train_gp {} '.format(gp_method), threshold=4, indent=1)

        
        hps = gp.hyperparameters
        self.msg('gpCAM using hyperparameters: {}'.format(repr(hps)), 4, 1)
        self.msg('with log-likelihood: {:,.1f}'.format(gp.log_likelihood(hps)), 4, 2)
        
        
        grid, xi, yi, XI, YI = self.make_grid(**kwargs)
        self.msg("Interpolating {:,} points to {:,}×{:,} = {:,} points (densification = {:.0f})".format(len(self.z_vals), len(xi), len(yi), len(xi)*len(yi), len(xi)*len(yi)/len(self.z_vals)), 3, 1)

        
        ZI = np.zeros_like(XI)
        
        self.msg('Filling ZI grid ({})'.format(fill), 4, 1)
        self.timing_start()
        if fill=='pixelwise':
            for ix, x in enumerate(xi):
                self.timing_progress_msg(ix, len(xi), 4)
                for iy, y in enumerate(yi):
                    res = gp.posterior_mean( np.array([[x,y]]) )
                    ZI[iy,ix] = res['f(x)'][0]
            # 0.6s/2,601pts = 0.2ms/pt
            # 2.2s/10,201pts = 0.2ms/pt
            # 47.6s/252,252pts = 0.2ms/pt
        
        elif fill=='whole':
            points = np.column_stack((np.ravel(XI), np.ravel(YI)))
            res = gp.posterior_mean(points)
            ZI = np.reshape(res['f(x)'], XI.shape)
            # 0.8s/2,601pts = 0.3ms/pt
            # 6.3s/10,201pts = 0.6ms/pt
            # ERR/252,252pts ; uses >470 GiB RAM
            
        self.timing_end_msg('{} fill'.format(fill), threshold=4, indent=2)
        
        
        if renormalize_signals:
            ZI = (ZI*std)+avg
        
        if convex_clip:
            # Force the grid to only exist within the convex hull defined by the points
            import scipy.interpolate
            POINTS = np.column_stack((self.x_vals, self.y_vals))
            VALUES = self.z_vals
            grid, xi, yi, XI, YI = self.make_grid(**kwargs)
            clip = np.isnan(scipy.interpolate.griddata(POINTS, VALUES, (XI, YI), rescale=True, method='linear'))
            ZI = np.ma.masked_where( clip, ZI)

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
        spread = np.square(np.std(signals[:,0]))
        hps_bounds = np.array([ [spread*1e-2, spread*20] ])
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
        (bears no resemblance to itself).
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
        SimpOpt = SimpleOptimizer(hps_ranges, initial_guess=hps_guess, evaluation_function=gp.log_likelihood, verbosity=10)
        
        
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
    def plot_signal(self, x='x_pos', y='y_pos', signal=None, N_limit=None, **kwargs):
        
        self.select(x=x, y=y, signal=signal, N_limit=N_limit)
        self.interpolate(**kwargs)
        self.plot(**kwargs)
        #self.plot3D(**kwargs)

    
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
        
        self.msg('Saved plot3D as: {}'.format(outfile), 4, 1)
    
    
    def plot_signal1D(self, x='x_pos', signal=None, N_limit=None, **kwargs):

        self.select(x=x, signal=signal, N_limit=N_limit)
        self.plot1D(**kwargs)
        

    def plot1D(self, outfile=None, title=None, plot_buffers=[0.18,0.05,0.16,0.05], **kwargs):
        self.msg('Plotting (1D) {}'.format(self.z_name), 3, 0)
        
        d = DataLineSM(x=self.x_vals, y=self.z_vals, y_err=self.z_errors)
        d.model = self
        d.x_rlabel = self.label(self.x_name)
        d.y_rlabel = self.label(self.z_name)
        
        zmin, zmax, cmap = self.zscale(self.z_name)
        d.plot_range = kwargs['plot_range'] if 'plot_range' in kwargs else [None, None, zmin, zmax]
        kwargs['plot_range'] = d.plot_range
        d.cmap = cmap
        
        outfile = self.get_outfile(outfile, extra='1D-')
        d.plot(save=outfile, title=title, plot_buffers=plot_buffers, **kwargs)
        
        self.msg('Saved plot1D as: {}'.format(outfile), 4, 1)


    # Additional
    ########################################
    def get_outfile(self, outfile=None, extra='', show_N=True, subdir=True, extra_dir=None, ext='.png'):
        if outfile is None:
            try:
                protocol, signal = self.z_name.split('__')
            except:
                protocol, signal = '', self.z_name
                
            outdir = Path(self.output_dir)
            if subdir:
                #outdir = outdir.joinpath(signal)
                outdir = outdir.joinpath('{}__{}'.format(protocol,signal))
            if extra_dir is not None:
                outdir = outdir.joinpath(extra_dir)
            outdir.mkdir(parents=True, exist_ok=True)
            N = 'N{:04d}'.format(self.get_N()) if show_N else ''
            outfile = Path(outdir, '{}-{}-{}{}{}'.format(self.sample, signal, extra, N, ext))

        return outfile
        
        
    def copy_current(self, outfile=None, copy_to=None, online=True):
        '''Copies the most recently-created plot to another location.
        The intention is for this copied file to act as an updating
        status of the experiment.
        NOTE: This code is highly contingent. Various paths are hard-
        coded, and will need to be changed for a new setup.'''
        
        # NB: This code used to be called "status_online"
        outfile = self.get_outfile(outfile)
        if False:
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
        
        if online:
            # Put a copy into some kind of online storage
            code_PATH='../../../../' # TOCHANGE
            code_PATH in sys.path or sys.path.append(code_PATH)
            from CustomS3 import Queue_analyze as queue
            q = queue()
            # Publish for this experiment
            q.publish_status_file(outfile, 'current_map')
            # Publish for the generic current beamline status
            q.experiment = 'current'
            q.publish_status_file(outfile, 'current_map')

    
    ########################################
    # End: class SurrogateModel(Base)



# SurrogateModel variants
########################################
class SurrogateModelLithoArray(SurrogateModel):

    def coord_transform_ebeam(self, origin=(0, 0)):

        x_vals = self.select_single('x_corrected')
        y_vals = self.select_single('y_corrected')
        
        slope = (130-20)/(6.88708-0)
        pitch = (x_vals-origin[0])*slope + 20
        lw = y_vals
        
        
        interpolator_pitch, interpolator_lw = self.coord_transform_ebeam_helper()
        pitch = interpolator_pitch( np.stack((x_vals, y_vals), axis=1) )
        lw = interpolator_lw( np.stack((x_vals, y_vals), axis=1) )

        duty_cycle = 100*lw/pitch

        self.data['parameters']['pitch'] = pitch
        self.label_associations.append(['pitch', '$\Lambda \, (\mathrm{mm})$', None, None, None])
        self.data['parameters']['linewidth'] = lw
        self.label_associations.append(['linewidth', '$\mathrm{linewidth} \, (\mathrm{nm})$', None, None, None])
        self.data['parameters']['duty_cycle'] = duty_cycle
        self.label_associations.append(['duty_cycle', '$\mathrm{duty \, cycle} \, ( \% )$', None, None, None])
        

    def coord_transform_ebeam_helper(self, infile='design_lookup03.txt', origin=(0, 0)):
        # Patterned area: contiguous ~60x60um regions.
        # x: from 20 nm to 130 nm, 1 nm steps: i.e. 111 regions, total width 6.6 mm
        #    however, patterned areas are not exactly 60um wide, so patterned region is 6.88708 mm wide
        # y: 49 regions, total height 2.94 mm
        
        
        l_edge = 654 # um
        patch_sizes = np.asarray([59.97, 61.97, 63.97, 65.97, 67.97, 59.96, 61.67, 63.38, 65.10, 66.81, 59.96, 61.46, 62.96, 64.46, 65.96, 59.94, 61.27, 62.61, 63.94, 65.27, 59.95, 61.15, 62.35, 63.55, 64.75, 59.90, 60.99, 62.08, 63.16, 64.25, 59.94, 60.94, 61.94, 62.94, 63.94, 59.93, 60.85, 61.78, 62.70, 63.62, 59.92, 60.78, 61.63, 62.49, 63.35, 59.93, 60.73, 61.53, 62.32, 63.12, 59.92, 60.67, 61.42, 62.17, 62.92, 59.84, 60.55, 61.25, 61.95, 62.66, 59.85, 60.52, 61.18, 61.85, 62.51, 59.85, 60.48, 61.11, 61.74, 62.37, 59.9, 60.5, 61.1, 61.7, 62.3, 59.85, 60.42, 60.99, 61.56, 62.13, 59.84, 60.39, 60.93, 61.47, 62.02, 59.8, 60.32, 60.84, 60.84, 61.88, 59.88, 60.38, 60.88, 61.38, 61.88, 59.88, 60.36, 60.84, 61.31, 61.79, 59.8]) # Starting at 30 nm patch, which is 654um from left edge of patterned area        
        
        data = []
        
        with open(infile, 'r') as fin:
            for idose, line in enumerate(fin.readlines()):
                for ipitch, element in enumerate(line.split()):
                    if element!='NA':
                        
                        #x = origin[0] + (ipitch+10)*0.060
                        x = origin[0] + ( l_edge + np.sum(patch_sizes[:ipitch]) + patch_sizes[ipitch]*0.5 )/1000.0
                        y = origin[1] + idose*0.060
                        pitch = 30.0 + ipitch*1.0
                        linewidth = float(element)
                        
                        data.append( [x, y, pitch, linewidth] )
                        
        
        data = np.asarray(data)
        
        tools.val_stats(data[:,0], 'x')
        tools.val_stats(data[:,1], 'y')
        tools.val_stats(data[:,2], 'pitch')
        tools.val_stats(data[:,3], 'lw')
                        
        
        from scipy.interpolate import NearestNDInterpolator
        interpolator_pitch = NearestNDInterpolator( data[:,:2], data[:,2], rescale=True )
        interpolator_lw = NearestNDInterpolator( data[:,:2], data[:,3], rescale=True )
        
        if PLOT_LAMDC:
            from scipy.spatial import Delaunay
            self.limit_hull = Delaunay( data[:,:2] ) # Hull in the original (x,y) coordinate space
        
        return interpolator_pitch, interpolator_lw


    def signal_transform_angle(self):
        
        angle = self.select_single('orientation_angle') # -180 to +180
        
        idx = np.where(angle<0)
        angle[idx] = angle[idx] + 180 # 0 to +180 (taking advantage of SAXS 2-fold symmetry)
        
        idx = np.where(angle>90)
        angle[idx] = 180-angle[idx] # 0 to +90 (taking advantage of presumsed symmetry about qx)
        
        self.data['signals']['orientation_angle_corrected'] = angle
        self.label_associations.append(['orientation_angle_corrected', '$\chi_{\mathrm{c}} \, (^{\circ})$', 0, +90, 'Blues'])

    
class SurrogateModelCombi(SurrogateModel):
    
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
        

class SurrogateModel3DPrint(SurrogateModel):        
    
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
        
        # Apply desired transformations
        #model.coord_transform_rotate(angle=-2.25, origin=(-3.35+7/2, 0.57-3/2))
        #model.coord_transform_translate(origin=(-3.37, 0.57-3)) # Adjusted by eye
        #model.coord_transform_translate(origin=(-3.45, 0.57-3)) # Re-adjusted based on plot1D, matching pitch
        #model.coord_transform_ebeam()
        
        #model.signal_transform_angle()
        #model.compute_missing_errors()
        #model.fix_nans()
        
        #model.trim_coord(x_corrected=[0.28, 6.87], y_corrected=[0.1, 2.9]) # for plot1D
        
        
        
        
    def plot_signals(self, model, signals, **kwargs):
        
        
        orig_name = model.name
        if self.distributed:
            futures = []
            
        
        hps_guess = kwargs['hps_guess'] if 'hps_guess' in kwargs else None # Save the 'default'
            
        for signal, zmin, zmax in signals:
            
            if 'hps_list' in kwargs:
                if signal in kwargs['hps_list'].keys():
                    kwargs['hps_guess'] = kwargs['hps_list'][signal]
                else:
                    kwargs['hps_guess'] = hps_guess
            
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
model = SurrogateModel(sample='AE_DSA_SIS_sampleD_run2_SM', sample_pattern='AE_DSA_SIS_sampleD_run2_', source_dir='../', verbosity=verbosity)

# Signals to extract from SciAnalysis output
peak_extractions = [
        'fit_peaks_prefactor1',
        'fit_peaks_prefactor1_error', 
        #'fit_peaks_m',
        #'fit_peaks_b',
        'fit_peaks_x_center1',
        'fit_peaks_x_center1_error', 
        'fit_peaks_sigma1',
        'fit_peaks_sigma1_error', 
        'fit_peaks_d0',
        'fit_peaks_d0_error', 
        'fit_peaks_grain_size',
        'fit_peaks_grain_size_error',
        'fit_peaks_chi_squared', 
        ]
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
        #'max_position',
        #'max_height',
        #'fwhm',
        #'orientation_angle', 
        #'orientation_factor', 
        
        #'fit_eta_prefactor',
        #'fit_eta_prefactor_error',
        #'fit_eta_eta',
        #'fit_eta_eta_error', 
        
        #'fit_eta_span_prefactor',
        #'fit_eta_span_prefactor_error',
        #'fit_eta_span_eta',
        #'fit_eta_span_eta_error', 
        
        #'fit_MaierSaupe_m',
        #'fit_MaierSaupe_m_error',
        #]
    #]
    ['circular_average_q2I_fit', peak_extractions],
    #['line90', peak_extractions],
    #['line60', peak_extractions],
    #['line30', peak_extractions],
    #['line90thin', peak_extractions],
    #['line60thin', peak_extractions],
    #['line45thin', peak_extractions],
    #['line30thin', peak_extractions],
    ]


    
    
    
# Define plotting behavior
signals = [
    ['prefactor1', 1e3, 6e3],
    ['grain_size', 32, 39],
    ['d0', 0.537, 0.539],
    ]
model.add_zscales(signals)

grid, n_grid = [None,None,None,None], 200 # Typical
#grid, n_grid = [0, 65, -0.45, 0.75], [1000,250]
figsize, plot_buffers = (10,10), [0.21,0.12,0.18,0.10] # Typical square
#figsize, plot_buffers = (20,10), [0.12,0.06,0.18,0.10] # Extended rectangle


# Define hyperparameters (for gpCAM)
hps_guess = None # Will try to guess

# Saved by running gpCAM instance
#hps_guess = model.load_hps(signal_name='prefactor', source_dir = '../../../gpcamv4and5/data/current_data/')

# Manually defined
#hps_guess = [1.0, 10, 0.1] # Example for a default anisotropic kernel
#hps_guess = [1.0, 8.0, 0.094, 0.296, 0] # periodic_kernel_2d_anisotropic

# Accumulated list
#hps_list = {
    #'linecut_angle_fit__fit_eta_eta' : [1.0892445 , 2.88554136, 2.04346513],
    #'orientation_angle_corrected' : [1.01415859, 0.38701475, 0.29268221],
    #'orientation_factor' : [1.02035666, 0.39905006, 0.2893818],
    #}



if __name__ == '__main__':
    
    task = Tasker(distributed=False, verbosity=verbosity)
    
    task.load_data(model, extractions)
    task.prepare_data(model)
    
    kwargs = {
        'x': 'x_position' ,
        'y': 'y_position' ,
        'interp_mode': 'griddata' , # 'griddata', 'gpcam'
        'grid': grid , 
        'n_grid': n_grid , 
        'hps_guess': hps_guess ,
        #'hps_list' : hps_list ,
        #'hps_lock' : [False, True, True],
        'gp_method': None , # 'global', 'local', 'hgdl', None (don't optimize)
        'pre_optimize': False ,
        'figsize': figsize ,
        'plot_buffers': plot_buffers ,
        'gpcam_PATH': '/home/kyager/current/code/gpcam/main/',
        #'aspect': True,
        }
    
    
    task.plot_signals(model, signals, **kwargs)
    
    kwargs['signal'] = 'd0'
    #task.plot_sequence(model, force=False, **kwargs)
    #task.animation(model, **kwargs)
    
    
    #model.plot_signal1D(x='pitch', signal='orientation_angle_corrected', plot_range=[0, 140, 0, 90], yticks=[0, 30, 60, 90], transparent=True)
    #model.save_data(extra='-final')
    
    
    
    
