#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.CurveAnalysis.Data` - Base objects for analysing 1D curves
================================================
.. module:: SciAnalysis.CurveAnalysis
   :synopsis: Provides base classes for doing analysis of 1D (x,y) curves
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  This code defines some baseline objects for analysing 1D curves.
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################



#import sys
import re # Regular expressions

import numpy as np
import pylab as plt
import matplotlib as mpl
#from scipy.optimize import leastsq
#import scipy.special


from .. import tools
from ..Data import *




   
    
# DataLineStructured
################################################################################    
class DataLineStructured(DataLine):
    
    
    def analyze(self, outfile, plot=True, save=None, show=False, **run_args):
        
        results = {}
        
        if plot:
            self.plot(save=outfile, **run_args)
        
        return results
        
        
    def plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.15,0.05,0.15,0.05], **kwargs):
        '''Plots the scattering data.
        
        Parameters
        ----------
        save : str
            Set to 'None' to avoid saving to disk. Provide filename to save.
        show : bool
            Set to true to open an interactive window.
        plot_range : [float, float, float, float]
            Set the range of the plotting (None scales automatically instead).
        '''  
        
        self._plot(save=save, show=show, plot_range=plot_range, plot_buffers=plot_buffers, **kwargs)
        
        
    # End class DataLineStructured(DataLine)
    ########################################
    
    
# DataLineStructuredSort
################################################################################    
class DataLineStructuredSort(DataLineStructured):
    
    # Analysis
    ########################################
    
    def analyze(self, outfile, plot=True, save=None, show=False, **run_args):
        
        results = {}
        import scipy
        
        
        # Renormalize
        self.x -= np.min(self.x) # Rezero
        self.x *= 1.0/np.max(self.x) # Rescale
        #self.y -= np.min(self.y) # Rezero
        #self.y *= 1.0/np.max(self.y) # Rescale
        self.y -= np.average(self.y) # Set zero to average
        self.y /= np.std(self.y) # Set yscale to standard deviation
        
        
        # Fit primary data
        self.fit_lines = []
        color_list = ['b', 'purple', 'r', 'green', 'orange',]
        for i, fit_name in enumerate(['fit_sines', 'fit_eta', ]):
            
            #lm_result, fit_line, fit_line_e = self.fit_sines(self, **run_args)
            lm_result, fit_line, fit_line_e = getattr(self, fit_name)(self, **run_args)
            fit_line.plot_args['color'] = color_list[i%len(color_list)]
            self.fit_lines.append([fit_line, fit_line_e])
        
            prefactor_total = 0
            for param_name, param in lm_result.params.items():
                results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                if 'prefactor' in param_name:
                    prefactor_total += np.abs(param.value)
                
            results['{}_prefactor_total'.format(fit_name)] = prefactor_total
            results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree

        # Histogram
        hist, bin_edges = np.histogram(self.y, bins=20, range=[-4,+4])
        self.hist_x = bin_edges[:-1] + 0.5*(bin_edges[1]-bin_edges[0])
        self.hist_y = hist/np.max(hist)
        results['histogram'] = self.hist_y
        
        # Fit histogram
        self.fit_lines_hist = []
        color_list = ['darkblue', 'darkmagenta', 'darkred', 'darkgreen']
        for i, fit_name in enumerate(['fit_gauss', 'fit_lognormal']):
            
            lm_result, fit_line, fit_line_e = getattr(self, fit_name)(DataLine(x=self.hist_x, y=self.hist_y), **run_args)
            fit_line.plot_args['color'] = color_list[i%len(color_list)]
            self.fit_lines_hist.append([fit_line, fit_line_e])
        
            prefactor_total = 0
            for param_name, param in lm_result.params.items():
                results['histogram_{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                if 'prefactor' in param_name:
                    prefactor_total += np.abs(param.value)
                
            results['histogram_{}_prefactor_total'.format(fit_name)] = prefactor_total
            results['histogram_{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
            results['histogram_{}_residuals'.format(fit_name)] = self.hist_y-fit_line.y
        

        # Cumulative Distribution Function
        self.cdf_x = self.hist_x
        self.cdf_y = hist.cumsum()
        self.cdf_y = self.cdf_y/np.max(self.cdf_y)


        # Fit CDF
        self.fit_lines_cdf = []
        color_list = ['b', 'purple', 'r', 'green', 'orange',]
        for i, fit_name in enumerate(['fit_erf']):
            
            lm_result, fit_line, fit_line_e = getattr(self, fit_name)(DataLine(x=self.cdf_x, y=self.cdf_y), **run_args)
            fit_line.plot_args['color'] = color_list[i%len(color_list)]
            self.fit_lines_cdf.append([fit_line, fit_line_e])
        
            prefactor_total = 0
            for param_name, param in lm_result.params.items():
                results['cdf_{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                if 'prefactor' in param_name:
                    prefactor_total += np.abs(param.value)
                
            results['cdf_{}_prefactor_total'.format(fit_name)] = prefactor_total
            results['cdf_{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
            results['cdf_{}_residuals'.format(fit_name)] = self.cdf_y-fit_line.y


        

        # Create sorted_y data (monotonically decreasing curve)

        # TODO: Handle the case where the x-data is not uniformly spaced (e.g. gaps).
        #self.sort_x()
        
        idx = np.argsort(self.y)
        self.y_sort = self.y[idx][::-1]
        self.x_sort = self.x
        
        
        # Fit primary data
        self.fit_lines_sort = []
        color_list = ['b', 'purple', 'r', 'green', 'orange',]
        for i, fit_name in enumerate(['fit_erfinv', 'fit_nlog', 'fit_expdecay', 'fit_linear']):
            
            #lm_result, fit_line, fit_line_e = self.fit_sines(self, **run_args)
            lm_result, fit_line, fit_line_e = getattr(self, fit_name)(DataLine(x=self.x_sort, y=self.y_sort), **run_args)
            fit_line.plot_args['color'] = color_list[i%len(color_list)]
            self.fit_lines_sort.append([fit_line, fit_line_e])
        
            prefactor_total = 0
            for param_name, param in lm_result.params.items():
                results['sort_{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                if 'prefactor' in param_name:
                    prefactor_total += np.abs(param.value)
                
            results['sort_{}_prefactor_total'.format(fit_name)] = prefactor_total
            results['sort_{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
            
            f_interp = scipy.interpolate.interp1d(self.x_sort, self.y_sort-fit_line.y)
            results['sort_{}_residuals'.format(fit_name)] = f_interp(np.linspace(0, 1.0-1e-10, num=40, endpoint=True))
            # Note kludge of endpoint being "1.0-1e-10" instead of "1.0", because interp1d sometimes issues an error ("A value in x_new is below the interpolation") when the endpoints of the fit and prediction match too closely.
        
        
        if plot:
            self.plot(save=outfile, **run_args)
            
            
        return results
        
        
    def fit_peaks(self, line, num_curves=6, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Gaussian with linear background.'''
            m = v['b']
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('b', value=np.average(line.y))
        
        xspan = np.max(line.x) - np.min(line.x)
        for i in range(num_curves):
            xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            params.add('prefactor{:d}'.format(i+1), value=np.max(line.y))
            params.add('x_center{:d}'.format(i+1), value=xpos)
            params.add('sigma{:d}'.format(i+1), value=xspan/num_curves, min=0)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended
    
    def fit_sines(self, line, num_curves=3, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Gaussian with linear background.'''
            m = x*0.0
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.sin( 2.0*np.pi*(x-v['x_center{:d}'.format(i+1)]) /(v['period{:d}'.format(i+1)]) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('b', value=np.average(line.y))
        
        xspan = np.max(line.x) - np.min(line.x)
        for i in range(num_curves):
            xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            params.add('prefactor{:d}'.format(i+1), value=np.max(line.y)*0.5)
            params.add('x_center{:d}'.format(i+1), value=xpos)
            params.add('period{:d}'.format(i+1), value=0.5/(i+1), min=0.05, max=2)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended    
    
    def fit_eta(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['prefactor']*( 1 - (v['eta']**2) )/( ((1+v['eta'])**2) - 4*v['eta']*( np.square(np.cos(  (v['symmetry']/2.0)*(x-v['x_center'])*(2.0*np.pi)  )) ) ) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=1.0, min=0)
        params.add('x_center', value=0.5, min=np.min(line.x), max=np.max(line.x))
        params.add('eta', value=0.4, min=0, max=1)
        params.add('symmetry', value=4, min=0.5, max=20)
        params.add('baseline', value=np.min(line.y), vary=False)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})

        return lm_result, fit_line, fit_line_extended
    
    
    def fit_gauss(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Gaussian with linear background.'''
            m = v['prefactor']*np.exp( -np.square(x-v['x_center'])/(2*(v['sigma']**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=1.0, min=0)
        params.add('x_center', value=0.0, min=-4, max=+4)
        params.add('sigma', value=1.0, min=0)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended
    
    
    def fit_lognormal(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Gaussian with linear background.'''
            x = np.maximum( x-v['x_zero'], x*0.0+1e-20 )
            m = v['prefactor']*np.exp( -np.square( np.log(x) - v['mu'] )/(2*(v['sigma']**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=1.0, min=0)
        params.add('x_zero', value=-1.5, min=-4, max=+4)
        params.add('mu', value=0.0, min=0)
        params.add('sigma', value=0.5, min=0)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended
    
    
    def fit_erf(self, line, **run_args):
        
        import lmfit
        import scipy
        
        def model(v, x):
            '''Error function (integral of Gaussian).'''
            m = v['prefactor']*scipy.special.erf( (x-v['x_center'])/v['x_scale'] ) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=0.5, min=0, max=2.0)
        params.add('x_center', value=0.0, min=-4, max=+4)
        params.add('x_scale', value=1.0, min=0, max=8.0)
        params.add('baseline', value=0.5, min=0, max=1)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended    
    
    def fit_erfinv(self, line, **run_args):
        
        import lmfit
        import scipy
        
        def model(v, x):
            '''Inverse error function.'''
            x = np.clip(x, 0.001, 0.999) # Avoid infinities (asymptotes) at endpoints
            m = v['prefactor']*scipy.special.erfinv( -(x*2-1) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=1.4, min=0)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended    
    
    
    def fit_nlog(self, line, **run_args):
        
        import lmfit
        import scipy
        
        def model(v, x):
            x = np.clip(x, 0.001, 1.0)
            m = v['prefactor']*-1.0*np.log(x/v['x_scale']) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=1.0, min=0)
        params.add('x_scale', value=1.0, min=0)
        params.add('baseline', value=np.min(line.y))
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended
    
    
    def fit_expdecay(self, line, **run_args):
        
        import lmfit
        import scipy
        
        def model(v, x):
            x = np.clip(x, 0.001, 1.0)
            m = v['prefactor']*np.exp(-x/v['x_scale']) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=3.0-np.min(line.y), min=0)
        params.add('x_scale', value=0.5, min=0)
        params.add('baseline', value=np.min(line.y), max=np.max(line.y))
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended    
    
    
    def fit_linear(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            m = v['m']*x + v['b']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('m', value=-3.0)
        params.add('b', value=line.x[0])
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended        
    
    # Plotting
    ########################################
    
    def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.15,0.05,0.15,0.05], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, **kwargs):
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        
        # Layout
        ########################################
       
        self.fig = plt.figure( figsize=(10,10), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        
        subfig_gap = 0.02
        subfig_width = fig_width
        subfig_height = (fig_height-subfig_gap)*0.5
        
        vresid_size_r = 0.15
        hhist_size_r = 0.30
        subfig_fig_height = subfig_height*(1-vresid_size_r)
        subfig_fig_width = subfig_width*(1-hhist_size_r)
        subfig_resid_height = subfig_height*(vresid_size_r)
        subfig_resid_width = subfig_width*(1-hhist_size_r)
        
        subfig_hist_width = subfig_width*hhist_size_r
        
        
        


        # Upper plot (origina data)
        self.ax1f = self.fig.add_axes( [left_buf, bottom_buf+subfig_height+subfig_gap, subfig_fig_width, subfig_fig_height] )
        self.ax1r = self.fig.add_axes( [left_buf, bottom_buf+subfig_height+subfig_gap+subfig_fig_height, subfig_resid_width, subfig_resid_height] )
        
        # Upper plot histogram
        self.ax1fh = self.fig.add_axes( [left_buf+subfig_resid_width, bottom_buf+subfig_height+subfig_gap, subfig_hist_width*0.75, subfig_fig_height] )
        self.ax1fhr = self.fig.add_axes( [left_buf+subfig_resid_width+subfig_hist_width*0.75, bottom_buf+subfig_height+subfig_gap, subfig_hist_width*0.25, subfig_fig_height] )
        
        # Lower plot (sorted data)
        self.ax2f = self.fig.add_axes( [left_buf, bottom_buf, subfig_fig_width, subfig_fig_height] )
        self.ax2r = self.fig.add_axes( [left_buf, bottom_buf+subfig_fig_height, subfig_resid_width, subfig_resid_height] )

        # Lower plot histogram
        self.ax2fh = self.fig.add_axes( [left_buf+subfig_resid_width, bottom_buf, subfig_hist_width*0.75, subfig_fig_height] )
        self.ax2fhr = self.fig.add_axes( [left_buf+subfig_resid_width+subfig_hist_width*0.75, bottom_buf, subfig_hist_width*0.25, subfig_fig_height] )
        
        
        # Plotting
        ########################################
        
        # Upper plot (origina data)
        p_args = dict([(i, plot_args[i]) for i in self.plot_valid_keys if i in plot_args])
        l, = self.ax1f.plot(self.x, self.y, **p_args)
        if dashes is not None:
            l.set_dashes(dashes)

        self.ax1r.axhline(y=0, color='k', zorder=-10)
        for fit_line, fit_line_e in reversed(self.fit_lines):
            
            p_args_cur = dict([(i, fit_line.plot_args[i]) for i in fit_line.plot_valid_keys if i in fit_line.plot_args])
            self.ax1f.plot(fit_line_e.x, fit_line_e.y, zorder=-1, **p_args_cur)
            
            resid_y = fit_line.y-self.y
            self.ax1r.plot(self.x, resid_y, zorder=-1, **p_args_cur)


        # Upper plot histogram
        bin_spacing = self.hist_x[1]-self.hist_x[0]
        self.ax1fh.barh(self.hist_x-0.5*bin_spacing, self.hist_y, height=bin_spacing, color='0.5', alpha=0.75)
        self.ax1fhr.axvline(x=0, color='k', zorder=-10)
        for fit_line, fit_line_e in reversed(self.fit_lines_hist):
            
            p_args_cur = dict([(i, fit_line.plot_args[i]) for i in fit_line.plot_valid_keys if i in fit_line.plot_args])
            self.ax1fh.plot(fit_line_e.y, fit_line_e.x, zorder=-1, **p_args_cur)
            
            resid_y = fit_line.y-self.hist_y
            self.ax1fhr.plot(resid_y, fit_line.x, zorder=-1, **p_args_cur)


        # Lower plot (sorted data)
        l, = self.ax2f.plot(self.x_sort, self.y_sort, **p_args)
        
        self.ax2r.axhline(y=0, color='k', zorder=-10)
        for fit_line, fit_line_e in reversed(self.fit_lines_sort):
            
            p_args_cur = dict([(i, fit_line.plot_args[i]) for i in fit_line.plot_valid_keys if i in fit_line.plot_args])
            self.ax2f.plot(fit_line_e.x, fit_line_e.y, zorder=-1, **p_args_cur)
            
            resid_y = fit_line.y-self.y_sort
            self.ax2r.plot(self.x, resid_y, zorder=-1, **p_args_cur)
        
        
        # Lower plot histogram
        bin_spacing = self.cdf_x[1]-self.cdf_x[0]
        self.ax2fh.barh(self.cdf_x-0.5*bin_spacing, self.cdf_y, height=bin_spacing, color='0.5', alpha=0.75)
        self.ax2fhr.axvline(x=0, color='k', zorder=-10)
        for fit_line, fit_line_e in reversed(self.fit_lines_cdf):
            
            p_args_cur = dict([(i, fit_line.plot_args[i]) for i in fit_line.plot_valid_keys if i in fit_line.plot_args])
            self.ax2fh.plot(fit_line_e.y, fit_line_e.x, zorder=-1, **p_args_cur)
            
            resid_y = fit_line.y-self.cdf_y
            self.ax2fhr.plot(resid_y, fit_line.x, zorder=-1, **p_args_cur)

        
        self.ax1f.set_xlabel('')
        self.ax1f.set_xticks([])
        #self.ax1f.set_ylabel(self.y_rlabel)
        self.ax1f.set_ylabel(r'$\left ( y-\langle y \rangle \right ) / \sigma_y$')
        self.ax1r.set_xlabel('')
        self.ax1r.set_xticks([])
        self.ax1r.set_yticklabels([])
        
        self.ax1fh.set_xticklabels([])
        self.ax1fh.set_yticklabels([])
        self.ax1fhr.set_xticklabels([])
        self.ax1fhr.set_yticklabels([])
        
        self.ax2f.set_xlabel(self.x_rlabel)
        self.ax2f.set_ylabel(r'$y_{\mathrm{sorted}}$')
        self.ax2r.set_xlabel('')
        self.ax2r.set_xticks([])
        self.ax2r.set_yticklabels([])
        
        self.ax2fh.set_xticklabels([])
        self.ax2fh.set_yticklabels([])
        self.ax2fhr.set_xticklabels([])
        self.ax2fhr.set_yticklabels([])

        
        
        if xlog:
            plt.semilogx()
        if ylog:
            plt.semilogy()
        
        
        # Axis scaling
        xi, xf, yi, yf = self.ax1f.axis()

        xi = 0
        xf = 1
        yf = np.max( [np.max(self.y), np.abs(np.min(self.y))] )
        yi = -yf
        
        if plot_range[0] != None: xi = plot_range[0]
        if plot_range[1] != None: xf = plot_range[1]
        if plot_range[2] != None: yi = plot_range[2]
        if plot_range[3] != None: yf = plot_range[3]
        self.ax1f.axis( [xi, xf, yi, yf] )
        self.ax1r.axis( [xi, xf, yi, yf] )
        self.ax1fh.axis( [0, 1.1, yi, yf] )
        self.ax1fhr.axis( [-0.5, 0.5, yi, yf] )
        
        self.ax2f.axis( [xi, xf, yi, yf] )
        self.ax2r.axis( [xi, xf, yi, yf] )
        self.ax2fh.axis( [0, 1.1, yi, yf] )
        self.ax2fhr.axis( [-0.5, 0.5, yi, yf] )
        
        self._plot_extra(**plot_args)
        
        if save:
            plt.savefig(save)
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)
        
 
    

    # End class DataLineStructuredSort(DataLineStructured)
    ########################################
    
    
# DataLineStructuredStd
################################################################################    
class DataLineStructuredStd(DataLineStructured):
    
    def analyze(self, outfile, plot=True, save=None, show=False, **run_args):
        
        results = {}
        
        # Renormalize
        self.x -= np.min(self.x) # Rezero
        self.x *= 1.0/np.max(self.x) # Rescale
        self.y -= np.average(self.y) # Set zero to average
        self.y /= np.std(self.y) # Set yscale to standard deviation
        
        
        # Fit primary data
        self.fit_lines = []
        color_list = ['b', 'purple', 'r', 'green', 'orange',]
        for i, fit_name in enumerate(['fit_zones']):
            
            #lm_result, fit_line, fit_line_e = self.fit_sines(self, **run_args)
            lm_result, fit_line, fit_line_e = getattr(self, fit_name)(self, **run_args)
            fit_line.plot_args['color'] = color_list[i%len(color_list)]
            self.fit_lines.append([fit_line, fit_line_e])
        
            prefactor_total = 0
            for param_name, param in lm_result.params.items():
                results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                if 'prefactor' in param_name:
                    prefactor_total += np.abs(param.value)
                
            results['{}_prefactor_total'.format(fit_name)] = prefactor_total
            results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
        
        
        
        self.std_lines = []
        self.std_lines_x = []
        self.std_lines_y = []
        for half_points in range(2, int( 0.5*(len(self.x)-1) ), 2 ):
            num_points = half_points*2 + 1
            
            curve = np.zeros(len(self.x))
            for i, x in enumerate(self.x):
                subrange = np.roll(self.y, +half_points-i)[:num_points]
                curve[i] = np.std(subrange)
            
            self.std_lines.append(curve)
            self.std_lines_x.append(num_points/len(self.x))
            self.std_lines_y.append(np.average(curve))

        vals = np.clip(np.asarray(self.std_lines).ravel(), 0, 2)
        hist, bin_edges = np.histogram(vals, bins=20, range=[0,2])
        self.hist_s_x = bin_edges[:-1] + 0.5*(bin_edges[1]-bin_edges[0])
        self.hist_s_y = 1.*hist/np.max(hist)
        results['histogram_stds'] = self.hist_s_y

        vals = np.clip(np.asarray(self.std_lines[:min(3,len(self.std_lines))]).ravel(), 0, 2)
        #vals = np.clip(np.asarray(self.std_lines[0]).ravel(), 0, 2)
        hist, bin_edges = np.histogram(vals, bins=20, range=[0,2])
        self.hist_sl_x = bin_edges[:-1] + 0.5*(bin_edges[1]-bin_edges[0])
        self.hist_sl_y = 1.*hist/np.max(hist)
        results['histogram_stds_loc'] = self.hist_sl_y
        

        
        if plot:
            self.plot(save=outfile, **run_args)
        
        return results
    
    
    def fit_zones(self, line, num_curves=8, **run_args):
        
        import lmfit
        
        def model(v, x):
            xspan = 0.5*(np.max(line.x) - np.min(line.x))/num_curves
            m = x*0.0
            for i in range(num_curves):
                xpos = (i*2+1)*xspan
                m += np.where(np.abs(x-xpos)<=xspan, v['prefactor{:d}'.format(i+1)], 0.0)
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        
        
        for i in range(num_curves):
            params.add('prefactor{:d}'.format(i+1), value=0.0)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})        

        return lm_result, fit_line, fit_line_extended    
    
    
    # Plotting
    ########################################
        
    def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.15,0.05,0.15,0.05], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, **kwargs):
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        
        # Layout
        ########################################
       
        self.fig = plt.figure( figsize=(10,10), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        
        subfig_gap = 0.02
        subfig_width = fig_width
        subfig_height = (fig_height-subfig_gap)*0.5
        
        vresid_size_r = 0.15
        hhist_size_r = 0.30
        subfig_fig_height = subfig_height*(1-vresid_size_r)
        subfig_fig_width = subfig_width*(1-hhist_size_r)
        subfig_resid_height = subfig_height*(vresid_size_r)
        subfig_resid_width = subfig_width*(1-hhist_size_r)
        
        subfig_hist_width = subfig_width*hhist_size_r
        
        
        


        # Upper plot (origina data)
        self.ax1f = self.fig.add_axes( [left_buf, bottom_buf+subfig_height+subfig_gap, subfig_fig_width, subfig_fig_height] )
        self.ax1r = self.fig.add_axes( [left_buf, bottom_buf+subfig_height+subfig_gap+subfig_fig_height, subfig_resid_width, subfig_resid_height] )
        
        # Upper plot histogram
        self.ax1fh = self.fig.add_axes( [left_buf+subfig_resid_width, bottom_buf+subfig_height+subfig_gap, subfig_hist_width*0.75, subfig_fig_height] )
        self.ax1fhr = self.fig.add_axes( [left_buf+subfig_resid_width+subfig_hist_width*0.75, bottom_buf+subfig_height+subfig_gap, subfig_hist_width*0.25, subfig_fig_height] )
        
        # Lower plot
        self.ax2f = self.fig.add_axes( [left_buf, bottom_buf, subfig_fig_width, subfig_fig_height] )
        self.ax2r = self.fig.add_axes( [left_buf, bottom_buf+subfig_fig_height, subfig_resid_width, subfig_resid_height] )

        # Lower plot histogram
        self.ax2fh = self.fig.add_axes( [left_buf+subfig_resid_width, bottom_buf, subfig_hist_width*0.75, subfig_fig_height] )
        self.ax2fhr = self.fig.add_axes( [left_buf+subfig_resid_width+subfig_hist_width*0.75, bottom_buf, subfig_hist_width*0.25, subfig_fig_height] )
        
        self.ax1fhr.set_visible(False)
        self.ax2fhr.set_visible(False)
        
        
        # Plotting
        ########################################
        
        # Upper plot (origina data)
        p_args = dict([(i, plot_args[i]) for i in self.plot_valid_keys if i in plot_args])
        l, = self.ax1f.plot(self.x, self.y, **p_args)
        if dashes is not None:
            l.set_dashes(dashes)
            
        # Plot derivative
        self.ax1f.plot(self.x, np.gradient(self.y), color='green')
            
        ones = self.x*0+1.0
        self.ax1f.axhline(y=0, color='k', zorder=-10)
        self.ax1f.fill_between(self.x, -1*ones, +1*ones, color='0.5', alpha=0.1, zorder=-10)
        self.ax1f.fill_between(self.x, -2*ones, +2*ones, color='0.5', alpha=0.1, zorder=-10)
        self.ax1f.fill_between(self.x, -3*ones, +3*ones, color='0.5', alpha=0.1, zorder=-10)

        self.ax1r.axhline(y=0, color='k', zorder=-10)
        for fit_line, fit_line_e in reversed(self.fit_lines):
            
            p_args_cur = dict([(i, fit_line.plot_args[i]) for i in fit_line.plot_valid_keys if i in fit_line.plot_args])
            self.ax1f.plot(fit_line_e.x, fit_line_e.y, zorder=-1, **p_args_cur)
            
            resid_y = fit_line.y-self.y
            self.ax1r.plot(self.x, resid_y, zorder=-1, **p_args_cur)


        # Upper plot, histogram of derivative
        hist, bin_edges = np.histogram(np.gradient(self.y), bins=20, range=[-4,+4])
        bin_spacing = bin_edges[1]-bin_edges[0]
        hist = 1.*hist/np.max(hist)
        self.ax1fh.barh(bin_edges[:-1], hist, height=bin_spacing, color='green', alpha=0.75)
        self.ax1fhr.axvline(x=0, color='k', zorder=-10)
        


        # Lower plot
        avg = self.x*0.0
        for i, line in enumerate(self.std_lines):
            amt = 1.*i/len(self.std_lines)
            color = (.2, .2, 1-amt*0.8)
            self.ax2f.plot(self.x, line, color=color, alpha=0.5)
            avg += line
            
            xs = line*0 + self.std_lines_x[i]
            self.ax2r.plot(xs, line, 'o', color='0.5', alpha=0.05)
            
            
        # Average curve
        self.ax2f.plot(self.x, avg/len(self.std_lines), color='blue', alpha=0.75, linewidth=2.0)


        # Lower plot, upper part ("residuals") area used to plot how stdev varies with binning width
        self.ax2r.plot(self.std_lines_x, self.std_lines_y, '-o', color='k')

            
        # Lower plot histogram
        bin_spacing = self.hist_s_x[1]-self.hist_s_x[0]
        self.ax2fh.barh(self.hist_s_x-0.5*bin_spacing, self.hist_s_y, height=bin_spacing, color='0.5', alpha=0.75)

        bin_spacing = self.hist_sl_x[1]-self.hist_sl_x[0]
        self.ax2fh.barh(self.hist_sl_x-0.5*bin_spacing, self.hist_sl_y, height=bin_spacing, color='b', alpha=0.4)




        
        self.ax1f.set_xlabel('')
        self.ax1f.set_xticks([])
        #self.ax1f.set_ylabel(self.y_rlabel)
        self.ax1f.set_ylabel(r'$\left ( y-\langle y \rangle \right ) / \sigma_y$')
        self.ax1r.set_xlabel('')
        self.ax1r.set_xticks([])
        self.ax1r.set_yticklabels([])
        
        self.ax1fh.set_xticklabels([])
        self.ax1fh.set_yticklabels([])
        self.ax1fhr.set_xticklabels([])
        self.ax1fhr.set_yticklabels([])
        
        self.ax2f.set_xlabel(self.x_rlabel)
        self.ax2f.set_ylabel(r'$\sigma_y(x,\delta x)$')
        self.ax2r.set_xlabel('')
        self.ax2r.set_xticks([])
        self.ax2r.set_yticklabels([])
        
        self.ax2fh.set_xticklabels([])
        self.ax2fh.set_yticklabels([])
        self.ax2fhr.set_xticklabels([])
        self.ax2fhr.set_yticklabels([])

        
        
        if xlog:
            plt.semilogx()
        if ylog:
            plt.semilogy()
        
        
        # Axis scaling
        xi, xf, yi, yf = self.ax1f.axis()

        xi = 0
        xf = 1
        yf = np.max( [np.max(self.y), np.abs(np.min(self.y))] )
        yi = -yf
        
        if plot_range[0] != None: xi = plot_range[0]
        if plot_range[1] != None: xf = plot_range[1]
        if plot_range[2] != None: yi = plot_range[2]
        if plot_range[3] != None: yf = plot_range[3]
        self.ax1f.axis( [xi, xf, yi, yf] )
        self.ax1r.axis( [xi, xf, yi, yf] )
        self.ax1fh.axis( [0, 1.1, yi, yf] )
        self.ax1fhr.axis( [-0.5, 0.5, yi, yf] )
        
        self.ax2f.axis( [xi, xf, 0, 2.0] )
        self.ax2r.axis( [xi, xf, 0, 2.0] )
        self.ax2fh.axis( [0, 1.1, 0, 2.0] )
        self.ax2fhr.axis( [-0.5, 0.5, yi, yf] )
        
        self._plot_extra(**plot_args)
        
        if save:
            plt.savefig(save)
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)    
    

    # End class DataLineStructuredStd(DataLineStructured)
    ########################################



# DataLineStructuredFFT
################################################################################    
class DataLineStructuredFFT(DataLineStructured):
    
    
    def analyze(self, outfile, plot=True, save=None, show=False, **run_args):
        
        results = {}
        import scipy
        
        # Renormalize
        self.x -= np.min(self.x) # Rezero
        self.x *= 1.0/np.max(self.x) # Rescale
        self.y -= np.average(self.y) # Set zero to average
        self.y /= np.std(self.y) # Set yscale to standard deviation
        
        
        self.fft_y = np.fft.fft( self.y )
        self.fft_y /= np.max(np.abs(self.fft_y))
        self.fft_x = np.arange(len(self.fft_y))*2.0*np.pi/1.0
        
        f_interp = scipy.interpolate.interp1d(self.fft_x, np.abs(self.fft_y))
        results['histogram_abs'] = f_interp(np.linspace(0, np.max(self.fft_x), num=40, endpoint=True))
        
        results['spectral_spread'] = np.std(np.abs(self.fft_y))
        
        
        if plot:
            self.plot(save=outfile, **run_args)
        
        return results
        
        
    # Data analysis
    ########################################        
    def stats(self, prepend='stats_', threshold=0.25):
        
        results = {}
        
        y = np.abs(self.fft_y)
        
        results[prepend+'max'] = np.max(y)
        results[prepend+'min'] = np.min(y)
        results[prepend+'average'] = np.average(y)
        results[prepend+'std'] = np.std(y)
        results[prepend+'N'] = len(y)
        results[prepend+'total'] = np.sum(y)
        
        results[prepend+'skew'] = stats.skew(y)
        
        results[prepend+'spread'] = results[prepend+'max'] - results[prepend+'min']
        results[prepend+'std_rel'] = results[prepend+'std'] / results[prepend+'average']
        
        y -= threshold
        zero_crossings = np.where(np.diff(np.signbit(y)))[0]
        results[prepend+'threshold_crossings'] = len(zero_crossings)
        
        return results        
        
        
    # Plotting
    ########################################
    

    def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.15,0.05,0.15,0.05], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, **kwargs):
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        
        # Layout
        ########################################
       
        self.fig = plt.figure( figsize=(10,10), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        
        subfig_gap = 0.02
        subfig_width = fig_width
        subfig_height = (fig_height-subfig_gap)*0.5
        
        vresid_size_r = 0.15
        hhist_size_r = 0.30
        subfig_fig_height = subfig_height*(1-vresid_size_r)
        subfig_fig_width = subfig_width*(1-hhist_size_r)
        subfig_resid_height = subfig_height*(vresid_size_r)
        subfig_resid_width = subfig_width*(1-hhist_size_r)
        
        subfig_hist_width = subfig_width*hhist_size_r
        
        
        


        # Upper plot (origina data)
        self.ax1f = self.fig.add_axes( [left_buf, bottom_buf+subfig_height+subfig_gap, subfig_fig_width, subfig_fig_height] )
        self.ax1r = self.fig.add_axes( [left_buf, bottom_buf+subfig_height+subfig_gap+subfig_fig_height, subfig_resid_width, subfig_resid_height] )
        
        # Upper plot histogram
        self.ax1fh = self.fig.add_axes( [left_buf+subfig_resid_width, bottom_buf+subfig_height+subfig_gap, subfig_hist_width*0.75, subfig_fig_height] )
        self.ax1fhr = self.fig.add_axes( [left_buf+subfig_resid_width+subfig_hist_width*0.75, bottom_buf+subfig_height+subfig_gap, subfig_hist_width*0.25, subfig_fig_height] )
        
        # Lower plot
        self.ax2f = self.fig.add_axes( [left_buf, bottom_buf, subfig_fig_width, subfig_fig_height] )
        self.ax2r = self.fig.add_axes( [left_buf, bottom_buf+subfig_fig_height, subfig_resid_width, subfig_resid_height] )

        # Lower plot histogram
        self.ax2fh = self.fig.add_axes( [left_buf+subfig_resid_width, bottom_buf, subfig_hist_width*0.75, subfig_fig_height] )
        self.ax2fhr = self.fig.add_axes( [left_buf+subfig_resid_width+subfig_hist_width*0.75, bottom_buf, subfig_hist_width*0.25, subfig_fig_height] )
        
        
        self.ax1r.set_visible(False)
        self.ax1fh.set_visible(False)
        self.ax1fhr.set_visible(False)
        self.ax2r.set_visible(False)
        self.ax2fhr.set_visible(False)
        
        
        # Plotting
        ########################################
        
        # Upper plot (origina data)
        p_args = dict([(i, plot_args[i]) for i in self.plot_valid_keys if i in plot_args])
        l, = self.ax1f.plot(self.x, self.y, **p_args)
        if dashes is not None:
            l.set_dashes(dashes)
            
            
            
        # Lower plot (FFT)
        self.ax2f.plot(self.fft_x, np.abs(self.fft_y), 'k', linewidth=2.0)
        
        hist, bin_edges = np.histogram(np.abs(self.fft_y), bins=20, range=[-1,+1])
        bin_spacing = bin_edges[1]-bin_edges[0]
        hist = 1.*hist/np.max(hist)
        self.ax2fh.barh(bin_edges[:-1], hist, height=bin_spacing, color='0.5', alpha=0.75, zorder=1)
        
        
        self.ax2f.plot(self.fft_x, np.real(self.fft_y), 'b', linewidth=2.0)
        
        hist, bin_edges = np.histogram(np.real(self.fft_y), bins=20, range=[-1,+1])
        bin_spacing = bin_edges[1]-bin_edges[0]
        hist = 1.*hist/np.max(hist)
        self.ax2fh.barh(bin_edges[:-1], hist, height=bin_spacing, color='b', alpha=0.2, zorder=-1)

        
        self.ax2f.plot(self.fft_x, np.imag(self.fft_y), 'purple', linewidth=2.0)

        hist, bin_edges = np.histogram(np.imag(self.fft_y), bins=20, range=[-1,+1])
        bin_spacing = bin_edges[1]-bin_edges[0]
        hist = 1.*hist/np.max(hist)
        self.ax2fh.barh(bin_edges[:-1], hist, height=bin_spacing, color='purple', alpha=0.2, zorder=-2)
        


        
        self.ax1f.set_xlabel('')
        self.ax1f.set_xticks([])
        #self.ax1f.set_ylabel(self.y_rlabel)
        self.ax1f.set_ylabel(r'$\left ( y-\langle y \rangle \right ) / \sigma_y$')
        self.ax1r.set_xlabel('')
        self.ax1r.set_xticks([])
        self.ax1r.set_yticklabels([])
        
        self.ax1fh.set_xticklabels([])
        self.ax1fh.set_yticklabels([])
        self.ax1fhr.set_xticklabels([])
        self.ax1fhr.set_yticklabels([])
        
        self.ax2f.set_xlabel(r'$q$')
        self.ax2f.set_ylabel('FFT')
        self.ax2r.set_xlabel('')
        self.ax2r.set_xticks([])
        self.ax2r.set_yticklabels([])
        
        self.ax2fh.set_xticklabels([])
        self.ax2fh.set_yticklabels([])
        self.ax2fhr.set_xticklabels([])
        self.ax2fhr.set_yticklabels([])

        
        
        if xlog:
            plt.semilogx()
        if ylog:
            plt.semilogy()
        
        
        # Axis scaling
        xi, xf, yi, yf = self.ax1f.axis()

        xi = 0
        xf = 1
        yf = np.max( [np.max(self.y), np.abs(np.min(self.y))] )
        yi = -yf
        
        if plot_range[0] != None: xi = plot_range[0]
        if plot_range[1] != None: xf = plot_range[1]
        if plot_range[2] != None: yi = plot_range[2]
        if plot_range[3] != None: yf = plot_range[3]
        self.ax1f.axis( [xi, xf, yi, yf] )
        self.ax1r.axis( [xi, xf, yi, yf] )
        self.ax1fh.axis( [0, 1.1, yi, yf] )
        self.ax1fhr.axis( [-0.5, 0.5, yi, yf] )
        
        self.ax2f.axis( [0, 0.5*np.max(self.fft_x), -1, +1] )
        self.ax2r.axis( [xi, xf, 0, 2.0] )
        self.ax2fh.axis( [0, 1.1, -1, +1] )
        self.ax2fhr.axis( [-0.5, 0.5, yi, yf] )
        
        self._plot_extra(**plot_args)
        
        if save:
            plt.savefig(save)
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)  

        
    # End class DataLineStructuredFFT(DataLineStructured)
    ########################################
