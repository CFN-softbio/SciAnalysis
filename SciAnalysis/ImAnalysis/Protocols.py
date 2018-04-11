#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.ImAnalysis.Protocols` - Data analysis protocols
================================================
.. module:: SciAnalysis.ImAnalysis.Protocols
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


class ProcessorIm(Processor):

    
    def load(self, infile, **kwargs):

        data = Data2DImage(infile, **kwargs)
        data.infile = infile
        
        return data
        
        
        
class shell(Protocol):
    
    def __init__(self, name='shell_test', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : 0.5,
                        'blur' : 1.0,
                        'resize' : 0.5,
                        'cmap' : mpl.cm.bone,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'])
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize']) # Shrink
        
        data.set_z_display([None, None, 'gamma', 1.0])
        outfile = self.get_outfile(data.name, output_dir)
        data.plot_image(outfile, ztrim=[0.01, 0.005], **run_args)
        
        return results
    
    
    

class thumbnails(Protocol):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : 0.5,
                        'blur' : 1.0,
                        'resize' : 0.5,
                        'cmap' : mpl.cm.bone,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'])
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize']) # Shrink
        
        data.set_z_display([None, None, 'gamma', 1.0])
        outfile = self.get_outfile(data.name, output_dir)
        data.plot_image(outfile, ztrim=[0.01, 0.005], **run_args)
        
        return results
        
        
        
class fft(Protocol):
    
    def __init__(self, name='fft', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'blur' : None,
                        'Iqn_n' : 1.0,
                        'fourier_filter_invert' : False,
                        'fourier_filter_shift' : 0.3,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)
                
                
        # FFT
        data_fft = data.fft()
        results.update(data_fft.stats())
        
        
        # FFT plots
        if run_args['verbosity']>=1:
            data_fft.set_z_display( [None, None, 'gamma', 0.3] )
            outfile = self.get_outfile('fft', output_dir, ext='.png', ir=True)
            data_fft.plot(save=outfile, ztrim=[0.3,0.001], blur=run_args['blur'])
        if run_args['verbosity']>=2:
            data_fft.set_z_display( [None, None, 'gamma', 0.3] )
            outfile = self.get_outfile('fft_zoom', output_dir, ext='.png', ir=True)
            x_axis, y_axis = data_fft.xy_axes()
            q_max = np.max(x_axis)*0.25
            data_fft.plot(save=outfile, ztrim=[0.5,0.001], plot_range=[-q_max,+q_max,-q_max,+q_max], blur=run_args['blur'])
            if run_args['verbosity']>=4:
                outfile = self.get_outfile('fft_zoom_blur', output_dir, ext='.png', ir=True)
                data_fft.plot(save=outfile, ztrim=[0.5,0.001], plot_range=[-q_max,+q_max,-q_max,+q_max], blur=1.0,)
        if run_args['verbosity']>=3:
            data_fft.set_z_display( [None, None, 'linear', 0.3] )
            outfile = self.get_outfile('fft', output_dir, ext='.png', ir=True)
            x_axis, y_axis = data_fft.xy_axes()
            q_max = np.max(x_axis)*0.5
            data_fft.plot_components(save=outfile, plot_range=[-q_max,+q_max,-q_max,+q_max], blur=run_args['blur'])
            
            
        # 1D curve
        line = data_fft.circular_average()
        line.x_rlabel = '$q \, (\mathrm{nm}^{-1})$'
        outfile = self.get_outfile('fft_1d', output_dir, ext='.dat', ir=False)
        line.save_data(outfile)
        
        
        if 'q0' not in run_args:
            # Try to find a peak, in order to guess q0
            
            q0 = self.find_q0(line, output_dir, **run_args)
            run_args['q0'] = q0
            run_args['dq'] = q0*0.6
            
                
            
        # Peak fit near q0
        new_results, lines = self.analyze_q0(line, output_dir, **run_args)
        results.update(new_results)
        
        if run_args['verbosity']>=2:
            lines.x_label = 'q'
            lines.x_rlabel = '$q \, (\mathrm{nm}^{-1})$'
            lines.y_label = 'I'
            lines.y_rlabel = r'$\langle I \rangle \, (\mathrm{counts/pixel})$'
            
            outfile = self.get_outfile('fft_1d', output_dir, ext='.png', ir=True)
            plot_range = [0, np.max(line.x)*0.75, np.min(line.y[1:-1]), np.max(line.y[1:-1])]
            lines.plot(save=outfile, ylog=True, plot_range=plot_range)
        
        if run_args['verbosity']>=3:
            # I*q^2 vs. q (Kratky plot style)
            line_Iq2 = line.copy()
            line_Iq2.y = line_Iq2.y*np.square(line_Iq2.x)
            line_Iq2.y_label = 'Iq^2'
            line_Iq2.y_rlabel = r'$q^2 \langle I \rangle$'
            
            line_Iq2.trim(0, np.max(line_Iq2.x)*0.5)
            outfile = self.get_outfile('fft_1d_Iq2', output_dir, ext='.png', ir=True)
            plot_range = [0, np.max(line_Iq2.x), 0, None]
            line_Iq2.plot(save=outfile, plot_range=plot_range)        
            
            
        # Orientation analysis at q0
        q0 = results['q0']['value']
        q_spread = results['sigma_q0']['value']*3.0
        new_results = self.orientation_q0(data_fft, q0, q_spread, output_dir, **run_args)
        results.update(new_results)
        
        
        # Fourier-filtered image at q0
        q_spread = results['sigma_q0']['value']*5.0
        self.fourier_filter(data, q0, q_spread, output_dir, **run_args)            
            
            
        
        return results        
    
    
    def find_q0(self, line, output_dir, **run_args):
        '''Guess where the 'main' peak is.'''
        
        line_Iqn = line.copy()
        line_Iqn.y = line_Iqn.y*np.power(line_Iqn.x, run_args['Iqn_n'])
        line_Iqn.y_label = 'Iq^{:.2f}'.format(run_args['Iqn_n'])
        line_Iqn.y_rlabel = r'$q^{{ {} }} \langle I \rangle$'.format(run_args['Iqn_n'])

        # Focus on the 'good data' (avoid high-q, and very low-q)
        line_Iqn.x = line_Iqn.x[1:]
        line_Iqn.y = line_Iqn.y[1:]
        line_Iqn.trim(np.min(line_Iqn.x)*2, np.max(line_Iqn.x)*0.5)
        
        q0, I_max = line_Iqn.target_y(np.max(line_Iqn.y))
        
        if run_args['verbosity']>=4:
            
            class DataLines_current(DataLines):
                def _plot_extra(self, **plot_args):
                    xi, xf, yi, yf = self.ax.axis()
                    self.ax.text(q0, yf, r'$q_0$', color='b', size=30, verticalalignment='top', horizontalalignment='left')
                    self.ax.axvline(q0, color='b')
                    
            
            lines_current = DataLines_current([line_Iqn])
            lines_current.x_label = line_Iqn.x_label
            lines_current.x_rlabel = line_Iqn.x_rlabel
            lines_current.y_label = line_Iqn.y_label
            lines_current.y_rlabel = line_Iqn.y_rlabel
            
            
            outfile = self.get_outfile('fft_1d_Iqn', output_dir, ext='.png', ir=True)
            plot_range = [0, np.max(line_Iqn.x), 0, None]
            lines_current.plot(save=outfile, plot_range=plot_range, plot_buffers=[0.25,0.05,0.2,0.05])
            
            
        return q0
        
    
    def analyze_q0(self, line, output_dir, **run_args):
        '''Fit expected peak in 1D curve.'''
        
        results = {}
        
        q0 = run_args['q0']
        if 'dq' in run_args:
            dq = run_args['dq']
        else:
            dq = run_args['q0']*0.6
            
        sub_line = line.sub_range(q0-dq,q0+dq)
        
        lm_result, fit_line, fit_line_extended = self.fit_peak(sub_line, **run_args)
        
        q0 = lm_result.params['x_center'].value
        q0_err = lm_result.params['x_center'].stderr
        d0 = 2*np.pi/(q0)
        d0_err = (d0/q0)*q0_err

        if False:
            # Naive parameter copying
            params = lm_result.params.valuesdict()
            results.update(params)
        
        else:
            # 'Rich' results
            results['q0'] = { 'value': lm_result.params['x_center'].value, 
                                'units': 'nm^-1', 
                                'units_latex': r'\mathrm{nm}^{-1}', 
                                'unit_type': 'inverse distance', 
                                'error': lm_result.params['x_center'].stderr, 
                                'symbol': 'q0', 
                                'symbol_latex': 'q_0',
                                }
            results['d0'] = { 'value': d0,
                                'units': 'nm', 
                                'units_latex': r'\mathrm{nm}', 
                                'unit_type': 'distance', 
                                'error': d0_err, 
                                'symbol': 'd0', 
                                'symbol_latex': 'd_0',
                                }
            results['sigma_q0'] = { 'value': lm_result.params['sigma'].value,
                                'units': 'nm^-1', 
                                'units_latex': r'\mathrm{nm}^{-1}', 
                                'unit_type': 'inverse distance', 
                                'error': lm_result.params['sigma'].stderr, 
                                'symbol': 'sigma_q', 
                                'symbol_latex': '\sigma_q',
                                }
            
            results['peak0_prefactor'] = { 'value': lm_result.params['prefactor'].value, 'units': 'a.u.', 'units_latex': r'a.u.', 'unit_type': 'a.u.', 'error': lm_result.params['prefactor'].stderr, 'symbol': 'c', 'symbol_latex': 'c',}
            results['peak0_background_m'] = { 'value': lm_result.params['m'].value, 'units': 'a.u./(nm^-1)', 'units_latex': r'a.u./\mathrm{nm}^{-1}', 'unit_type': 'slope', 'error': lm_result.params['m'].stderr, 'symbol': 'm', 'symbol_latex': 'm',}
            results['peak0_background_b'] = { 'value': lm_result.params['b'].value, 'units': 'a.u.', 'units_latex': r'a.u.', 'unit_type': 'a.u.', 'error': lm_result.params['b'].stderr, 'symbol': 'b', 'symbol_latex': 'b',}

        

        
        
        text = r'$q_0 = {:.3g} \pm {:.2g} \, \mathrm{{nm}}^{{-1}}$'.format(q0, q0_err)
        text += '\n'
        text += r'$\sigma_q = {:.3g} \pm {:.2g} \, \mathrm{{nm}}^{{-1}}$'.format(lm_result.params['sigma'].value, lm_result.params['sigma'].stderr)
        text += '\n'
        text += r'$d_0 = {:.3g} \pm {:.2g} \, \mathrm{{nm}}$'.format(d0, d0_err)
        
        
        class DataLines_current(DataLines):
            def _plot_extra(self, **plot_args):
                xi, xf, yi, yf = self.ax.axis()
                plt.text(xf, yf, text, size=20, verticalalignment='top', horizontalalignment='right')

        if run_args['verbosity']>=4:
            outfile = self.get_outfile('fft_1d_peak', output_dir, ext='.png', ir=True)
            plot_range = [0, np.max(sub_line.x)*1.2, 0, None]
            lines = DataLines_current([sub_line, fit_line, fit_line_extended])
            lines.plot(save=outfile, ylog=False, plot_range=plot_range)                
        
        lines = DataLines_current([line, fit_line, fit_line_extended])        
        
        return results, lines
                   
                   
    
    def fit_peak(self, line, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peak(line, **run_args)
        
        import lmfit
        
        def model(v, x):
            '''Gaussian with linear background.'''
            m = v['prefactor']*np.exp( -np.square(x-v['x_center'])/(2*(v['sigma']**2)) ) + v['m']*x + v['b']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y), min=0)
        params.add('x_center', value=np.average(line.x))
        params.add('sigma', value=np.std(line.x), min=0)
        params.add('m', value=0)
        params.add('b', value=0)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        #fix_x = line.x
        #fit_y = line.y + lm_result.residual
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        

        x_span = abs(np.max(line.x) - np.min(line.x))
        xe = 0.5
        fit_x = np.linspace(np.min(line.x)-xe*x_span, np.max(line.x)+xe*x_span, num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':2.0})
            
        return lm_result, fit_line, fit_line_extended
        
        
    def orientation_q0(self, data, q_center, q_spread, output_dir, result_prepend='peak0_', **run_args):
        
        results = {}
        orig_data = data.data
        data.data = np.abs(data.data)
        
        line = data.linecut_angle(q_center, q_spread, **run_args)
        
        # Clean up the curve
        #line.remove_spurious(bins=5, tol=1.5)
        #line.smooth(1.0)
        for angle in [-180, -90, 0, +90, +180]:
            line.kill_x(angle, 1)

        outfile = self.get_outfile('orientation', output_dir, ext='.dat', ir=False)
        line.save_data(outfile)
        
        
        labels = []
        
        symmetry = 2
        lm_result, fit_line = self.angle_fit(line, symmetry_assumed=symmetry, color='b', **run_args)
        line_list = [line, fit_line]
        
        results['{}ori{}_{}'.format(result_prepend,symmetry,'eta')] = { 
                            'value': lm_result.params['eta'].value, 
                            'units': '', 
                            'units_latex': r'',
                            'unit_type': 'unitless', 
                            'error': lm_result.params['eta'].stderr, 
                            'symbol': 'eta', 
                            'symbol_latex': r'\eta',
                            }        
        results['{}ori{}_{}'.format(result_prepend,symmetry,'angle')] = { 
                            'value': lm_result.params['x_center'].value, 
                            'units': 'degrees', 
                            'units_latex': r'^{\circ}',
                            'unit_type': 'angle', 
                            'error': lm_result.params['x_center'].stderr, 
                            'symbol': 'chi0', 
                            'symbol_latex': r'\chi_0',
                            }   
        results['{}ori{}_{}'.format(result_prepend,symmetry,'prefactor')] = { 'value': lm_result.params['prefactor'].value, 'units': 'a.u.', 'units_latex': r'a.u.', 'unit_type': 'a.u.', 'error': lm_result.params['prefactor'].stderr, 'symbol': 'c', 'symbol_latex': 'c',}
        
        
        
        xp, yp = fit_line.target_x(lm_result.params['x_center'].value)
        text = '$\eta_{{ {} }} = {:.3g}$'.format(symmetry, lm_result.params['eta'].value)
        labels.append([xp, yp*1.2, text, 'b'])
        text = '${:.3g} ^{{ \circ }}$'.format(lm_result.params['x_center'].value)
        labels.append([xp, 0, text, 'b'])
        
        if 'symmetry' in run_args and run_args['symmetry']!=symmetry:
            # Run analysis again assuming a different symmetry
            
            symmetry = run_args['symmetry']
            lm_result, fit_line = self.angle_fit(line, symmetry_assumed=symmetry, color='purple', **run_args)
            line_list.append(fit_line)
            
            results['{}ori{}_{}'.format(result_prepend,symmetry,'eta')] = { 
                                'value': lm_result.params['eta'].value, 
                                'units': '', 
                                'units_latex': r'',
                                'unit_type': 'unitless', 
                                'error': lm_result.params['eta'].stderr, 
                                'symbol': 'eta', 
                                'symbol_latex': r'\eta',
                                }        
            results['{}ori{}_{}'.format(result_prepend,symmetry,'angle')] = { 
                                'value': lm_result.params['x_center'].value, 
                                'units': 'degrees', 
                                'units_latex': r'^{\circ}',
                                'unit_type': 'angle', 
                                'error': lm_result.params['x_center'].stderr, 
                                'symbol': 'chi0', 
                                'symbol_latex': r'\chi_0',
                                }   
            results['{}ori{}_{}'.format(result_prepend,symmetry,'prefactor')] = { 'value': lm_result.params['prefactor'].value, 'units': 'a.u.', 'units_latex': r'a.u.', 'unit_type': 'a.u.', 'error': lm_result.params['prefactor'].stderr, 'symbol': 'c', 'symbol_latex': 'c',}
            
            
            
            xp, yp = fit_line.target_x(lm_result.params['x_center'].value)
            text = '$\eta_{{ {} }} = {:.3g}$'.format(symmetry, lm_result.params['eta'].value)
            labels.append([xp, yp*1.2, text, 'purple'])
            text = '${:.3g} ^{{ \circ }}$'.format(lm_result.params['x_center'].value)
            labels.append([xp, yp*0.1, text, 'purple'])
        
        
        
        class DataLines_current(DataLines):
            def _plot_extra(self, **plot_args):
                xi, xf, yi, yf = self.ax.axis()
                for x, y, text, color in labels:
                    self.ax.text(x, y, text, size=20, color=color, verticalalignment='bottom', horizontalalignment='left')        
                    self.ax.axvline(x, color=color)
                    
        lines = DataLines_current(line_list)
        lines.x_label = 'angle'
        lines.x_rlabel = r'$\chi \, (^{\circ})$'
        lines.y_label = 'I'
        lines.y_rlabel = r'$I(\chi) \, (\mathrm{counts/pixel})$'
        
        
        if run_args['verbosity']>=2:
            outfile = self.get_outfile('orientation', output_dir, ext='.png', ir=True)
            plot_range = [-180, +180, 0, np.max(line.y)*1.2]
            lines.plot(save=outfile, plot_range=plot_range)

        if run_args['verbosity']>=3:
            outfile = self.get_outfile('orientation_polar', output_dir, ext='.png', ir=True)
            line.plot_polar(save=outfile)
            
            
        # TODO: Obtain order parameter from line
        # TODO: Add line.stats() to results
        
            
        data.data = orig_data
        
        return results
    
    
    def angle_fit(self, line, symmetry_assumed=2, color='b', **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['prefactor']*( 1 - (v['eta']**2) )/( ((1+v['eta'])**2) - 4*v['eta']*( np.square(np.cos((symmetry_assumed/2.0)*np.radians(x-v['x_center']))) ) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data        
        
        params = lmfit.Parameters()
        x_peak, y_peak = line.target_y(np.max(line.y))
        params.add('prefactor', value=y_peak, min=0)
        params.add('x_center', value=self.reduce_angle(x_peak, symmetry_assumed), min=-180/symmetry_assumed, max=+180/symmetry_assumed)
        params.add('eta', value=0.8, min=0, max=1)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = np.linspace(-180, +180, num=360, endpoint=True)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        
        fit_line = DataLineAngle(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':color, 'marker':None, 'linewidth':4.0})
        
        return lm_result, fit_line
    
    
    def reduce_angle(self, angle, symmetry):
        '''Reduce an angle to be minimal for the given symmetry.'''
        
        span = 180.0/symmetry
        # sym=1, repeats every 360deg, so span is -180 to +180
        # sym=2, repeats every 180deg, so span is -90 to +90
        # sym=4, repeats every 90deg, so span is -45 to +45
        # sym=6, repeats every 60deg, so span is -30 to +30
        
        while angle < span:
            angle += 2*span
        while angle > span:
            angle -= 2*span
            
        return angle
    
        
    def fourier_filter(self, data, q_center, q_spread, output_dir, **run_args):
        
        data.fourier_filter(q_center, q_spread)
        
        
        data.maximize_intensity_spread()
        if run_args['fourier_filter_invert']:
            data.invert()
        if run_args['fourier_filter_shift'] != 0.0:
            data.data = np.clip(data.data, 255*run_args['fourier_filter_shift'], 255)
            data.maximize_intensity_spread()
        
        
        if run_args['verbosity']>=2:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('filtered', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, zmin=0, zmax=255, cmap=mpl.cm.bone)


    # End class fft(Protocol)
    ########################################



class local_avg_realspace(Protocol):
    
    
    def __init__(self, name='local_avg_realspace', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'local_partition_image_size' : 75, # pixels
                        'local_partition_step' : 1.0, # relative to image_size
                        }
        self.run_args.update(kwargs)
        

        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)

        results = {}
        
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)
        


        # Pre-process
        data.equalize()
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        #data.lowkill(run_args['q0']*0.1)
        #data.blur(1.0)
        data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        data.maximize_intensity_spread()
        
        
        
        if run_args['verbosity']>=5:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('enhanced', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)
        
        
        exit()
        import imreg_dft as ird
        seed = Data2DImage(os.path.join(output_dir, '../seed.png'))
        local_image = Data2DImage()
        avg_image = np.zeros(seed.data.shape)
        n = 0
        
        height, width = data.data.shape
        par_size = run_args['local_partition_image_size']
        par_step = int(par_size*run_args['local_partition_step'])
        
        for ix in range(par_size, width-par_size, par_step):
            for iy in range(par_size, height-par_size, par_step):
                
                if run_args['verbosity']>=5:
                    ic = ix*height+iy
                    ip = 100.*ic/(width*height)
                    print('        ({},{}) (~{:.1f}% complete)'.format(ix,iy,ip))
                
                #local_image = data.data[iy-par_size:iy+par_size, ix-par_size:ix+par_size]
                local_image.data = data.data[iy-par_size:iy+par_size, ix-par_size:ix+par_size]
                #local_image.fourier_filter(0.1,0.01)
                #local_image.blur(2.0)
                #local_image.maximize_intensity_spread()
                #local_image.threshold_pixels(127)
                
                result = ird.similarity(seed.data, local_image.data, constraints={'scale':[1.0,0.01]})
                local_image.data = result['timg']
                
                avg_image += result['timg']
                n += 1
                
                if run_args['verbosity']>=3:
                    im = PIL.Image.fromarray( np.uint8(local_image.data) )
                    #outfile = self.get_outfile('im%03d_%03d'%(ix,iy), output_dir, ext='.png')
                    outfile = self.get_outfile('im{:0>4d}_{:0>4d}'.format(ix,iy), output_dir, ext='.png')
                    im.save(outfile)
                      
                      
                
        
        avg_image /= n
        if run_args['verbosity']>=2:
            im = PIL.Image.fromarray( np.uint8(avg_image) )
            outfile = self.get_outfile('average', output_dir, ext='.png', ir=True)
            im.save(outfile)
        

        results['test_result'] = { 'value': 1.0, 'units': 'A^-1', 'units_latex': r'\AA^{-1}', 'unit_type': 'inverse distance', 'error': 0.01, 'symbol': 'q', 'symbol_latex': 'q' }
        
        return results
    
    
    
    
    
    
class particles(Protocol):
    
    def __init__(self, name='particles', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        
        results = {}
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
        
        # Pre-process
        #data.equalize()
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        data.lowkill(run_args['q0']*0.1)
        data.blur(2.0) # lowpass
        #data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        #data.equalize()
        data.maximize_intensity_spread()

        if run_args['verbosity']>=4:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('processed', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
        data.threshold(run_args['threshold'], run_args['invert'])
        
        if run_args['verbosity']>=4:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('thresholded', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
            
        # Identify particles positions
        if 'diagonal_detection' in run_args and run_args['diagonal_detection']:
            #s = [[1,1,1],
            #    [1,1,1],
            #    [1,1,1]]
            s = ndimage.generate_binary_structure(2,2)
        else:
            s = [[0,1,0],
                [1,1,1],
                [0,1,0]]
        
        labeled_array, num_features = ndimage.measurements.label(data.data, structure=s)
        
        # TODO: Implement
        #if 'min_area' in run_args or 'max_area' in run_args:
        
        results['num_particles'] = num_features
            

        if run_args['verbosity']>=3 and False:
            # Colored image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            color_list = [ (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1),(0,1,1),(1,1,1),]
            color_list2 = [ (0.7*c[0], 0.7*c[1], 0.7*c[2]) for c in color_list ]
            color_list3 = [ (0.5*c[0], 1.0*c[1], 1.0*c[2]) for c in color_list ]
            color_list4 = [ (1.0*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list ]
            color_list5 = [ (1.0*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list ]
            color_list6 = [ (1.0*c[0], 0.7*c[1], 0.5*c[2]) for c in color_list ]
            color_list7 = [ (1.0*c[0], 0.5*c[1], 0.7*c[2]) for c in color_list ]
            color_list8 = [ (0.7*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list ]
            color_list9 = [ (0.5*c[0], 1.0*c[1], 0.7*c[2]) for c in color_list ]
            color_list10 = [ (0.7*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list ]
            color_list11 = [ (0.5*c[0], 0.7*c[1], 1.0*c[2]) for c in color_list ]
            color_list = color_list + color_list2 + color_list3 + color_list4 + color_list5 + color_list6 + color_list7 + color_list8 + color_list9 + color_list10 + color_list11
            
            h, w = data.data.shape
            for ix in range(w):
                for iy in range(h):
                    object_index = labeled_array[iy,ix]
                    if object_index==0:
                        c = (0, 0, 0)
                    else:
                        c = color_list[ object_index%len(color_list) ]
                        c = ( int(c[0]*255), int(c[1]*255), int(c[2]*255) )
                    pix[ix,iy] = c            
            
            
            outfile = self.get_outfile('colored', output_dir, ext='.png', ir=True)
            im.save(outfile)
            
        
        if run_args['verbosity']>=4 and False:
            # Boundary image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            c = ( 1*255, 0*255, 0*255 )
            h, w = data.data.shape
            for ix in range(w-1):
                for iy in range(h-1):
                    
                    num_zeros = np.bincount( labeled_array[iy:iy+2,ix:ix+2].flatten() )[0]
                    if not (num_zeros==0 or num_zeros==4):
                        pix[ix,iy] = c            
            

            outfile = self.get_outfile('boundaries', output_dir, ext='.png', ir=True)
            im.save(outfile)
            
        
        # Statistics on particles that have been found
        bins = np.bincount(labeled_array.flatten('C'))
        h, w = data.data.shape
        total_pixels = w*h
        background_pixels = bins[0]
        particle_pixels = total_pixels - background_pixels
        coverage = particle_pixels*1./(total_pixels*1.)
        results['coverage'] = coverage
        if run_args['verbosity']>=4:
            print('    Particle coverage: {:.1f}%'.format(coverage*100.))
        
        
        # Remove 'particles' of zero size
        idx = np.nonzero(bins)
        bins = bins[idx]
        # Remove the 'surrounding field' (index 0)
        bins = bins[1:]
        
        # Convert to physical sizes
        particle_sizes = bins*data.x_scale*data.y_scale # nm^2
        
        if 'area_min' in run_args:
            particle_sizes = particle_sizes[particle_sizes>run_args['area_min']]
        if 'area_max' in run_args:
            particle_sizes = particle_sizes[particle_sizes<run_args['area_max']]
        
        
        particle_radii = np.sqrt(particle_sizes/np.pi) # nm
        
        if 'radius_min' in run_args:
            particle_radii = particle_radii[particle_radii>run_args['radius_min']]
        if 'radius_max' in run_args:
            particle_radii = particle_radii[particle_radii<run_args['radius_max']]

        particle_sizes = np.pi*np.square(particle_radii)
        
        results['area_average'] = np.average(particle_sizes)
        results['area_std'] = np.std(particle_sizes)
        results['area_median'] = np.median(particle_sizes)
        
        results['radius_average'] = np.average(particle_radii)
        results['radius_std'] = np.std(particle_radii)
        results['radius_median'] = np.median(particle_radii)
        
        
        
        if run_args['verbosity']>=1:
            
            class DataHistogram_current(DataHistogram):
                def _plot_extra(self, **plot_args):
                    xi, xf, yi, yf = self.ax.axis()
                    yf = yf*1.2
                    
                    
                    self.ax.axvline(self.mean, color='b', linewidth=3.0, zorder=4)
                    self.ax.text(self.mean, yf, 'mean', color='b', rotation=90, verticalalignment='top', horizontalalignment='right', zorder=4)
                    
                    self.ax.axvspan( self.mean-self.std, self.mean+self.std, color='b', alpha=0.05, zorder=-10)
                    self.ax.axvspan( self.mean-2*self.std, self.mean+2*self.std, color='b', alpha=0.05, zorder=-10)

                    self.ax.axvline(self.median, color='purple', linewidth=2.0)
                    self.ax.text(self.median, yf, 'median', color='purple', rotation=90, verticalalignment='top', horizontalalignment='right')
                    
                    lm_result, fit_line, fit_line_e = self.fit_peak(self)
                    self.ax.axvline(lm_result.params['x_center'], color='r', linewidth=2.0)
                    self.ax.plot(fit_line_e.x, fit_line_e.y, 'r', linewidth=2.0)
                    
                    
                    els = self.x_rlabel.split('\,')
                    s = '{}= {:.1f} \pm {:.1f} \, {}'.format( els[0], self.mean, self.std, els[1].replace('(','').replace(')','') )
                    self.ax.text(xf, yf, s, size=30, verticalalignment='top', horizontalalignment='right')
                    
                    self.ax.axis( [xi, xf, yi, yf] )
                    
                    
                def fit_peak(self, line, **run_args):
                    import lmfit
                    
                    def model(v, x):
                        m = v['prefactor']*np.exp( -np.square(x-v['x_center'])/(2*(v['sigma']**2)) )
                        return m
                    
                    def func2minimize(params, x, data):
                        
                        v = params.valuesdict()
                        m = model(v, x)
                        
                        return m - data
                    
                    params = lmfit.Parameters()
                    params.add('prefactor', value=np.max(line.y), min=0, max=np.max(line.y)*2.0)
                    params.add('x_center', value=self.mean, min=np.min(line.x), max=np.max(line.x))
                    params.add('sigma', value=self.std, min=0)
                    
                    lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
                    
                    fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
                    fit_y = model(lm_result.params.valuesdict(), fit_x)
                    
                    fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
                    
                    x_span = abs(np.max(line.x) - np.min(line.x))
                    xe = 0.5
                    fit_x = np.linspace(np.min(line.x)-xe*x_span, np.max(line.x)+xe*x_span, num=2000)
                    fit_y = model(lm_result.params.valuesdict(), fit_x)
                    fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':2.0})
                        
                    return lm_result, fit_line, fit_line_extended
                    
            
            # Histogram of areas
            y, x = np.histogram(particle_sizes, bins=150, range=[0, max(particle_sizes)*1.05])
            
            # Instead of having x be ranges for each bar, center the x on the average of each range
            xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
            x = x[:-1]
            
            hist = DataHistogram_current(x=x, y=y, x_label='Area', x_rlabel='$A \, (\mathrm{nm}^{2})$', y_label='count')
            hist.mean = results['area_average']
            hist.std = results['area_std']
            hist.median = results['area_median']
            
            outfile = self.get_outfile('particle_areas', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
            
            
            # Histogram of radii
            y, x = np.histogram(particle_radii, bins=150, range=[0, max(particle_radii)*1.05])
            
            # Instead of having x be ranges for each bar, center the x on the average of each range
            xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
            x = x[:-1]
            
            hist = DataHistogram_current(x=x, y=y, x_label='Radius', x_rlabel='$r \, (\mathrm{nm})$', y_label='count')
            hist.mean = results['radius_average']
            hist.std = results['radius_std']
            hist.median = results['radius_median']
            
            outfile = self.get_outfile('particle_radii', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
            
            
        
        
        return results
        



            
            
class grain_size_hex(Protocol):
    
    def __init__(self, name='grain_size_hex', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        orig_data = data.data.copy()
        
        results = {}
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
        
        # Pre-process
        #data.equalize()
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        data.lowkill(run_args['q0']*0.1)
        data.blur(2.0) # lowpass
        #data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        #data.equalize()
        data.maximize_intensity_spread()

        if run_args['verbosity']>=4:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('processed', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
        data.threshold(run_args['threshold'], run_args['invert'])
        
        if run_args['verbosity']>=4:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('thresholded', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
            
        # Identify particles positions
        if 'diagonal_detection' in run_args and run_args['diagonal_detection']:
            #s = [[1,1,1],
            #    [1,1,1],
            #    [1,1,1]]
            s = ndimage.generate_binary_structure(2,2)
        else:
            s = [[0,1,0],
                [1,1,1],
                [0,1,0]]
        
        labeled_array, num_features = ndimage.measurements.label(data.data, structure=s)
        num_objects = num_features
        
        # Remove objects outside of size cutoffs
        bins = np.bincount(labeled_array.flatten('C'))
        bins = bins[1:] # Remove the 'surrounding field' (index 0)
        if 'area_min' in run_args:
            cutoff = run_args['area_min']/(data.x_scale*data.y_scale) # pix^2
            idx = np.where(bins<cutoff)[0]
            for object_index in idx:
                labeled_array[labeled_array==(object_index+1)] = 0
                
            #num_objects -= len(idx)
            
            # Update labeled_array
            labeled_array, num_features = ndimage.measurements.label(labeled_array, structure=s)
            num_objects = num_features
            

        # TODO:
        #if 'area_max' in run_args:
        
        
        
        results['num_particles'] = num_objects
            

        # Determine (x,y) position of each particle (center-of-mass of each particle)
        x_positions = np.zeros( (num_features) )
        y_positions = np.zeros( (num_features) )
        counts = np.zeros( (num_features) )
        for ix in range( len(labeled_array[0]) ):
            for iy in range( len(labeled_array) ):
                if labeled_array[iy,ix]!=0:
                    object_index = labeled_array[iy,ix] - 1
                    
                    x_positions[object_index] +=  1.0*ix
                    y_positions[object_index] +=  1.0*iy
                    counts[object_index] +=  1.0  
                    
        x_positions /= counts
        y_positions /= counts



        
        idx = np.where( np.nan_to_num(counts)>0 )[0]
        if(num_objects!=len(idx)):
            print('WARNING: disagreement in particle count.')
            print( '    num_features: {:d}, num_objects: {:d}, counts: {:d}'.format(num_features, num_objects, len(idx)) )
            
            
            
        if 'NN_cutoff_distance_pix' not in run_args:
            run_args['NN_cutoff_distance_pix'] = run_args['NN_cutoff_distance_nm']/( (data.x_scale + data.y_scale)*0.5 )
        NN_counts, angles, new_results = self.nearest_neighbor_count(x_positions, y_positions, counts, **run_args)
        results.update(new_results)



        
        
        
        
            
        if run_args['verbosity']>=5:
            # Colored image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            color_list = [ (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1),(0,1,1),(1,1,1),]
            color_list2 = [ (0.7*c[0], 0.7*c[1], 0.7*c[2]) for c in color_list ]
            color_list3 = [ (0.5*c[0], 1.0*c[1], 1.0*c[2]) for c in color_list ]
            color_list4 = [ (1.0*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list ]
            color_list5 = [ (1.0*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list ]
            color_list6 = [ (1.0*c[0], 0.7*c[1], 0.5*c[2]) for c in color_list ]
            color_list7 = [ (1.0*c[0], 0.5*c[1], 0.7*c[2]) for c in color_list ]
            color_list8 = [ (0.7*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list ]
            color_list9 = [ (0.5*c[0], 1.0*c[1], 0.7*c[2]) for c in color_list ]
            color_list10 = [ (0.7*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list ]
            color_list11 = [ (0.5*c[0], 0.7*c[1], 1.0*c[2]) for c in color_list ]
            color_list = color_list + color_list2 + color_list3 + color_list4 + color_list5 + color_list6 + color_list7 + color_list8 + color_list9 + color_list10 + color_list11
            
            h, w = data.data.shape
            for ix in range(w):
                for iy in range(h):
                    object_index = labeled_array[iy,ix]
                    if object_index==0:
                        c = (0, 0, 0)
                    else:
                        c = color_list[ object_index%len(color_list) ]
                        c = ( int(c[0]*255), int(c[1]*255), int(c[2]*255) )
                    pix[ix,iy] = c            
            
            
            outfile = self.get_outfile('colored', output_dir, ext='.png', ir=True)
            im.save(outfile)
            
        
        if run_args['verbosity']>=5:
            # Boundary image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            c = ( 1*255, 0*255, 0*255 )
            h, w = data.data.shape
            for ix in range(w-1):
                for iy in range(h-1):
                    
                    num_zeros = np.bincount( labeled_array[iy:iy+2,ix:ix+2].flatten() )[0]
                    if not (num_zeros==0 or num_zeros==4):
                        pix[ix,iy] = c            
            
            # Put a dot at center-of-mass (COM) of each object
            c = ( 0*255, 1*255, 0*255 )
            for object_index in range(len(counts)):
                #if counts[object_index]>0:
                ix = int(x_positions[object_index])
                iy = int(y_positions[object_index])
                pix[ix,iy] = c

            outfile = self.get_outfile('boundaries', output_dir, ext='.png', ir=True)
            im.save(outfile)

        if run_args['verbosity']>=3:
            # Color-coded image of nearest-neighbour (NN) count
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            for ix in range( len(labeled_array[0]) ):
                for iy in range( len(labeled_array) ):
                    object_index = labeled_array[iy,ix] - 1
                    
                    if object_index<0:
                        c = (0, 0, 0 )
                    else:
                        NN_count = NN_counts[object_index]
                        
                        if NN_count==0:
                            c = (20, 20, 20) # Dark grey
                        elif NN_count==1:
                            c = (40, 40, 40) # Dark grey
                        elif NN_count==2:
                            c = (60, 60, 60) # Dark grey
                        elif NN_count==3:
                            c = (0, 0, 150)
                        elif NN_count==4:
                            c = (0, 0, 200)
                        elif NN_count==5:
                            c = (0, 255, 255)
                        elif NN_count==6:
                            c = (0, 255, 0) # Greem
                        elif NN_count==7:
                            c = (255, 255, 0)
                        elif NN_count==8:
                            c = (100, 0, 0) # Dark red
                        elif NN_count==9:
                            c = (150, 0, 0) # Dark red
                        elif NN_count==10:
                            c = (200, 0, 0) # Dark red
                        else:
                            c = (255, 0, 0 ) # Red
                        
                    pix[ix,iy] = c
                    
            outfile = self.get_outfile('NN_count', output_dir, ext='.png', ir=True)
            im.save(outfile)


        angle_max_rad = 2.0*np.pi/run_args['symmetry']

        if run_args['verbosity']>=3:
            # False-color map of angles
            
            from scipy.interpolate import griddata
            
            positions = np.column_stack((x_positions,y_positions))
            grid_x, grid_y = np.mgrid[ 0:len(labeled_array) , 0:len(labeled_array[0]) ]
            angle_map = griddata(positions, angles, (grid_y, grid_x), method='nearest') # Avoids artifacts
            
            data_angles = Data2D()
            data_angles.data = angle_map
            data_angles.set_z_display([0, angle_max_rad, 'linear', 1.0])
            outfile = self.get_outfile('angles', output_dir, ext='.png', ir=True)
            data_angles.plot_image(outfile, ztrim=[0., 0.], cmap=cmap_cyclic_spectrum)
            
            
            # Overlay
            cur_data = orig_data - np.min(orig_data)
            cur_data = cur_data*(255.0/np.max(cur_data))
            img1 = PIL.Image.fromarray( np.uint8(cur_data) )
            img1 = img1.convert('RGBA')
            
            cmap = cmap_cyclic_spectrum
            Z = (data_angles.data-0.0)/(angle_max_rad-0.0)
            img2 = PIL.Image.fromarray(np.uint8(cmap(Z)*255))
            img2 = img2.convert('RGBA')
            
            img_blend = PIL.Image.blend(img1, img2, 0.5)
            outfile = self.get_outfile('angles', output_dir, ext='.png', ir=True)
            img_blend.save(outfile)

            
            
        # Angle histogram
        # TODO: Fix the angle histograms.
        hist, bin_edges = np.histogram(angles, bins=100, range=[0,angle_max_rad])
        bin_edges += bin_edges[1]-bin_edges[0]
        new_results = self.orientation_fit(np.degrees(bin_edges[:-1]), hist, output_dir, result_prepend='NN_', **run_args)
        results.update(new_results)            
        
        
        
        
        # Statistics on particles that have been found
        bins = np.bincount(labeled_array.flatten('C'))
        h, w = data.data.shape
        total_pixels = w*h
        background_pixels = bins[0]
        particle_pixels = total_pixels - background_pixels
        coverage = particle_pixels*1./(total_pixels*1.)
        results['coverage'] = coverage
        if run_args['verbosity']>=4:
            print('    Particle coverage: {:.1f}%'.format(coverage*100.))
        
        
        # Remove 'particles' of zero size
        idx = np.nonzero(bins)
        bins = bins[idx]
        # Remove the 'surrounding field' (index 0)
        bins = bins[1:]
        
        # Convert to physical sizes
        particle_sizes = bins*data.x_scale*data.y_scale # nm^2
        
        if 'area_min' in run_args:
            particle_sizes = particle_sizes[particle_sizes>run_args['area_min']]
        if 'area_max' in run_args:
            particle_sizes = particle_sizes[particle_sizes<run_args['area_max']]
        
        
        particle_radii = np.sqrt(particle_sizes/np.pi) # nm
        
        if 'radius_min' in run_args:
            particle_radii = particle_radii[particle_radii>run_args['radius_min']]
        if 'radius_max' in run_args:
            particle_radii = particle_radii[particle_radii<run_args['radius_max']]

        particle_sizes = np.pi*np.square(particle_radii)
        
        results['area_average'] = np.average(particle_sizes)
        results['area_std'] = np.std(particle_sizes)
        results['area_median'] = np.median(particle_sizes)
        
        results['radius_average'] = np.average(particle_radii)
        results['radius_std'] = np.std(particle_radii)
        results['radius_median'] = np.median(particle_radii)
        
        
        
        if run_args['verbosity']>=4:
            
            new_results = self.plot_particle_histograms(particle_radii, particle_sizes, output_dir, results, **run_args)
            results.update(new_results)
            
            
            
        # TODO: Compute correlation function.
      
            
        return results
            
            
            
    def plot_particle_histograms(self, particle_radii, particle_sizes, output_dir, prev_results, **run_args):
        
        results = {}
        
        class DataHistogram_current(DataHistogram):
            def _plot_extra(self, **plot_args):
                xi, xf, yi, yf = self.ax.axis()
                yf = yf*1.2
                
                
                self.ax.axvline(self.mean, color='b', linewidth=3.0, zorder=4)
                self.ax.text(self.mean, yf, 'mean', color='b', rotation=90, verticalalignment='top', horizontalalignment='right', zorder=4)
                
                self.ax.axvspan( self.mean-self.std, self.mean+self.std, color='b', alpha=0.05, zorder=-10)
                self.ax.axvspan( self.mean-2*self.std, self.mean+2*self.std, color='b', alpha=0.05, zorder=-10)

                self.ax.axvline(self.median, color='purple', linewidth=2.0)
                self.ax.text(self.median, yf, 'median', color='purple', rotation=90, verticalalignment='top', horizontalalignment='right')
                
                lm_result, fit_line, fit_line_e = self.fit_peak(self)
                self.ax.axvline(lm_result.params['x_center'], color='r', linewidth=2.0)
                self.ax.plot(fit_line_e.x, fit_line_e.y, 'r', linewidth=2.0)
                
                
                els = self.x_rlabel.split('\,')
                s = '{}= {:.1f} \pm {:.1f} \, {}'.format( els[0], self.mean, self.std, els[1].replace('(','').replace(')','') )
                self.ax.text(xf, yf, s, size=30, verticalalignment='top', horizontalalignment='right')
                
                self.ax.axis( [xi, xf, yi, yf] )
                
                
            def fit_peak(self, line, **run_args):
                import lmfit
                
                def model(v, x):
                    m = v['prefactor']*np.exp( -np.square(x-v['x_center'])/(2*(v['sigma']**2)) )
                    return m
                
                def func2minimize(params, x, data):
                    
                    v = params.valuesdict()
                    m = model(v, x)
                    
                    return m - data
                
                params = lmfit.Parameters()
                params.add('prefactor', value=np.max(line.y), min=0, max=np.max(line.y)*2.0)
                params.add('x_center', value=self.mean, min=np.min(line.x), max=np.max(line.x))
                params.add('sigma', value=self.std, min=0)
                
                lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
                
                fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
                fit_y = model(lm_result.params.valuesdict(), fit_x)
                
                fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
                
                x_span = abs(np.max(line.x) - np.min(line.x))
                xe = 0.5
                fit_x = np.linspace(np.min(line.x)-xe*x_span, np.max(line.x)+xe*x_span, num=2000)
                fit_y = model(lm_result.params.valuesdict(), fit_x)
                fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':2.0})
                    
                return lm_result, fit_line, fit_line_extended
                
        
        # Histogram of areas
        y, x = np.histogram(particle_sizes, bins=150, range=[0, max(particle_sizes)*1.05])
        
        # Instead of having x be ranges for each bar, center the x on the average of each range
        xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
        x = x[:-1]
        
        hist = DataHistogram_current(x=x, y=y, x_label='Area', x_rlabel='$A \, (\mathrm{nm}^{2})$', y_label='count')
        hist.mean = prev_results['area_average']
        hist.std = prev_results['area_std']
        hist.median = prev_results['area_median']
        
        outfile = self.get_outfile('particle_areas', output_dir, ext='.png')
        hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
        
        
        # Histogram of radii
        y, x = np.histogram(particle_radii, bins=150, range=[0, max(particle_radii)*1.05])
        
        # Instead of having x be ranges for each bar, center the x on the average of each range
        xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
        x = x[:-1]
        
        hist = DataHistogram_current(x=x, y=y, x_label='Radius', x_rlabel='$r \, (\mathrm{nm})$', y_label='count')
        hist.mean = prev_results['radius_average']
        hist.std = prev_results['radius_std']
        hist.median = prev_results['radius_median']
        
        outfile = self.get_outfile('particle_radii', output_dir, ext='.png')
        hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
        
        
        return results
        


    def nearest_neighbor_count(self, x_positions, y_positions, counts, **run_args):
        
        results = {}
        
        cutoff_pix = run_args['NN_cutoff_distance_pix']
        sym = run_args['symmetry']
        angle_max_deg = 360.0/sym # 60.0deg for 6-fold
        angle_max_rad = np.radians(angle_max_deg)
        
        num_features = len(counts)
        
        NN_counts = np.zeros(num_features)
        angles = np.zeros(num_features)
        
        for object_index in range(num_features):
            
            xpos = x_positions[object_index]
            ypos = y_positions[object_index]
            
            # Find nearby particles
            distances = np.sqrt( np.square( x_positions-xpos ) + np.square( y_positions-ypos ) )
            idx = np.where( distances<cutoff_pix )[0]
            idx = idx[ np.where( idx!=object_index )[0] ] # Excludes the particle itself
            
            NN_counts[object_index] = len(idx)
            
            cur_angles = ( np.arctan2( y_positions[idx]-ypos , x_positions[idx]-xpos ) )%angle_max_rad
            
            
            #angles[object_index] = np.sum(cur_angles)/NN_counts[object_index] # This naive average is wrong!
            angles[object_index] = ( self.average_angle( cur_angles*sym )/sym )%angle_max_rad
            
            if angles[object_index] <np.radians(0.0) or angles[object_index]>angle_max_rad:
                print( 'err', angles[object_index] )
                
                
                
        less = len(np.nonzero(NN_counts<sym)[0])
        equal = len(np.nonzero(NN_counts==sym)[0])
        more = len(np.nonzero(NN_counts>sym)[0])
        print("        <%d NN: %d (%.1f%%)" % ( sym, less, less*100.0/(less+equal+more) ) )
        print("        =%d NN: %d (%.1f%%)" % ( sym, equal, equal*100.0/(less+equal+more) ) )
        print("        >%d NN: %d (%.1f%%)" % ( sym, more, more*100.0/(less+equal+more) ) )
        
        results['NN_less%d'%(sym)] = less
        results['NN_less%d_fraction'%(sym)] = less*1.0/(less+equal+more)
        results['NN_equal%d'%(sym)] = equal
        results['NN_equal%d_fraction'%(sym)] = equal*1.0/(less+equal+more)
        results['NN_more%d'%(sym)] = more
        results['NN_more%d_fraction'%(sym)] = more*1.0/(less+equal+more)
        
        
        return NN_counts, angles, results


    def average_angle(self, list_of_angles):
        import cmath
        # from: http://stackoverflow.com/questions/491738/how-do-you-calculate-the-average-of-a-set-of-angles

        # make a new list of vectors
        vectors= [cmath.rect(1, angle) # length 1 for each vector
            for angle in list_of_angles]

        vector_sum = sum(vectors)

        # no need to average, we don't care for the modulus
        return cmath.phase(vector_sum)

            
    def orientation_fit(self, angle_vals, angle_ints, output_dir, result_prepend='NN_', **run_args):
        
        results = {}
        line = DataLineAngle(x=angle_vals, y=angle_ints)
        
        # Clean up the curve
        #line.remove_spurious(bins=5, tol=1.5)
        #line.smooth(1.0)
        for angle in [-180, -90, 0, +90, +180]:
            line.kill_x(angle, 1)

        outfile = self.get_outfile('orientation', output_dir, ext='.dat', ir=False)
        line.save_data(outfile)
        
        
        labels = []
        
        symmetry = 2
        lm_result, fit_line = self.angle_fit(line, symmetry_assumed=symmetry, color='b', **run_args)
        line_list = [line, fit_line]
        
        results['{}ori{}_{}'.format(result_prepend,symmetry,'eta')] = { 
                            'value': lm_result.params['eta'].value, 
                            'units': '', 
                            'units_latex': r'',
                            'unit_type': 'unitless', 
                            'error': lm_result.params['eta'].stderr, 
                            'symbol': 'eta', 
                            'symbol_latex': r'\eta',
                            }        
        results['{}ori{}_{}'.format(result_prepend,symmetry,'angle')] = { 
                            'value': lm_result.params['x_center'].value, 
                            'units': 'degrees', 
                            'units_latex': r'^{\circ}',
                            'unit_type': 'angle', 
                            'error': lm_result.params['x_center'].stderr, 
                            'symbol': 'chi0', 
                            'symbol_latex': r'\chi_0',
                            }   
        results['{}ori{}_{}'.format(result_prepend,symmetry,'prefactor')] = { 'value': lm_result.params['prefactor'].value, 'units': 'a.u.', 'units_latex': r'a.u.', 'unit_type': 'a.u.', 'error': lm_result.params['prefactor'].stderr, 'symbol': 'c', 'symbol_latex': 'c',}
        
        
        
        xp, yp = fit_line.target_x(lm_result.params['x_center'].value)
        text = '$\eta_{{ {} }} = {:.3g}$'.format(symmetry, lm_result.params['eta'].value)
        labels.append([xp, yp*1.2, text, 'b'])
        text = '${:.3g} ^{{ \circ }}$'.format(lm_result.params['x_center'].value)
        labels.append([xp, 0, text, 'b'])
        
        if 'symmetry' in run_args and run_args['symmetry']!=symmetry:
            # Run analysis again assuming a different symmetry
            
            symmetry = run_args['symmetry']
            lm_result, fit_line = self.angle_fit(line, symmetry_assumed=symmetry, color='purple', **run_args)
            line_list.append(fit_line)
            
            results['{}ori{}_{}'.format(result_prepend,symmetry,'eta')] = { 
                                'value': lm_result.params['eta'].value, 
                                'units': '', 
                                'units_latex': r'',
                                'unit_type': 'unitless', 
                                'error': lm_result.params['eta'].stderr, 
                                'symbol': 'eta', 
                                'symbol_latex': r'\eta',
                                }        
            results['{}ori{}_{}'.format(result_prepend,symmetry,'angle')] = { 
                                'value': lm_result.params['x_center'].value, 
                                'units': 'degrees', 
                                'units_latex': r'^{\circ}',
                                'unit_type': 'angle', 
                                'error': lm_result.params['x_center'].stderr, 
                                'symbol': 'chi0', 
                                'symbol_latex': r'\chi_0',
                                }   
            results['{}ori{}_{}'.format(result_prepend,symmetry,'prefactor')] = { 'value': lm_result.params['prefactor'].value, 'units': 'a.u.', 'units_latex': r'a.u.', 'unit_type': 'a.u.', 'error': lm_result.params['prefactor'].stderr, 'symbol': 'c', 'symbol_latex': 'c',}
            
            
            
            xp, yp = fit_line.target_x(lm_result.params['x_center'].value)
            text = '$\eta_{{ {} }} = {:.3g}$'.format(symmetry, lm_result.params['eta'].value)
            labels.append([xp, yp*1.2, text, 'purple'])
            text = '${:.3g} ^{{ \circ }}$'.format(lm_result.params['x_center'].value)
            labels.append([xp, yp*0.1, text, 'purple'])
        
        
        
        class DataLines_current(DataLines):
            def _plot_extra(self, **plot_args):
                xi, xf, yi, yf = self.ax.axis()
                for x, y, text, color in labels:
                    self.ax.text(x, y, text, size=20, color=color, verticalalignment='bottom', horizontalalignment='left')        
                    self.ax.axvline(x, color=color)
                    
        lines = DataLines_current(line_list)
        lines.x_label = 'angle'
        lines.x_rlabel = r'$\chi \, (^{\circ})$'
        lines.y_label = 'I'
        lines.y_rlabel = r'$I(\chi) \, (\mathrm{counts/pixel})$'
        
        
        if run_args['verbosity']>=3:
            outfile = self.get_outfile('orientation', output_dir, ext='.png', ir=True)
            plot_range = [0, 360.0/run_args['symmetry'], 0, np.max(line.y)*1.2]
            lines.plot(save=outfile, plot_range=plot_range)

        if run_args['verbosity']>=3:
            outfile = self.get_outfile('orientation_polar', output_dir, ext='.png', ir=True)
            line.plot_polar(save=outfile)
            
            
        # TODO: Obtain order parameter from line
        # TODO: Add line.stats() to results
        
            
        return results            
            
            
    def angle_fit(self, line, symmetry_assumed=6, color='b', **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['prefactor']*( 1 - (v['eta']**2) )/( ((1+v['eta'])**2) - 4*v['eta']*( np.square(np.cos((symmetry_assumed/2.0)*np.radians(x-v['x_center']))) ) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data        
        
        params = lmfit.Parameters()
        x_peak, y_peak = line.target_y(np.max(line.y))
        params.add('prefactor', value=y_peak, min=0)
        params.add('x_center', value=self.reduce_angle(x_peak, symmetry_assumed), min=-180/symmetry_assumed, max=+180/symmetry_assumed)
        params.add('eta', value=0.8, min=0, max=1)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = np.linspace(-180, +180, num=360, endpoint=True)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        
        fit_line = DataLineAngle(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':color, 'marker':None, 'linewidth':4.0})
        
        return lm_result, fit_line
    
                
    def reduce_angle(self, angle, symmetry):
        '''Reduce an angle to be minimal for the given symmetry.'''
        
        span = 180.0/symmetry
        # sym=1, repeats every 360deg, so span is -180 to +180
        # sym=2, repeats every 180deg, so span is -90 to +90
        # sym=4, repeats every 90deg, so span is -45 to +45
        # sym=6, repeats every 60deg, so span is -30 to +30
        
        while angle < span:
            angle += 2*span
        while angle > span:
            angle -= 2*span
            
        return angle
    
class grain_size(grain_size_hex):

    def __init__(self, name='grain_size', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'cmap' : mpl.cm.bone,
                        'blur_size_rel_d0' : 0.25,
                        'blur_orientation_image' : True,
                        'blur_orientation_image_num_passes' : 3,
                        'blur_orientation_image_size_rel' : 0.25,
                        }
        self.run_args.update(kwargs)
            
    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        
        #orig_data = data.data.copy()

        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
        
        # Pre-process
        #data.equalize()
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        #data.lowkill(run_args['q0']*0.1)
        data.blur(2.0) # lowpass
        data.enhance(contrast=1.3, contrast_passes=2, resharpen_passes=2)
        data.equalize()
        data.maximize_intensity_spread()

        if run_args['verbosity']>=4:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('processed', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
        orig_data = data.data.copy()
            
        angles = self.orientation_angle_map(data, output_dir, **run_args)
        angle_max_rad = +np.pi/2
        
        if run_args['verbosity']>=3:
            
            data_angles = Data2D()
            data_angles.data = angles
            data_angles.set_z_display( [None, None, 'linear', 1.0] )
            outfile = self.get_outfile('angles', output_dir, ext='.png', ir=True)
            data_angles.plot_image(save=outfile, zmin=-np.pi/2, zmax=+np.pi/2, cmap=cmap_cyclic_spectrum)
            
            # Overlay
            cur_data = orig_data - np.min(orig_data)
            cur_data = cur_data*(255.0/np.max(cur_data))
            img1 = PIL.Image.fromarray( np.uint8(cur_data) )
            img1 = img1.convert('RGBA')
            
            cmap = cmap_cyclic_spectrum
            Z = (data_angles.data+angle_max_rad)/(2*angle_max_rad)
            img2 = PIL.Image.fromarray(np.uint8(cmap(Z)*255))
            img2 = img2.convert('RGBA')
            
            img_blend = PIL.Image.blend(img1, img2, 0.5)
            outfile = self.get_outfile('angles', output_dir, ext='.png', ir=True)
            img_blend.save(outfile)


        # Angle histogram
        # TODO: Fix the angle histograms.
        hist, bin_edges = np.histogram(angles, bins=100, range=[-angle_max_rad,angle_max_rad])
        bin_edges += bin_edges[1]-bin_edges[0]
        new_results = self.orientation_fit(np.degrees(bin_edges[:-1]), hist, output_dir, result_prepend='ori_', **run_args)
        results.update(new_results)     

        
        return results



    def orientation_angle_map(self, data, output_dir, **run_args):
        
        if 'blur' in run_args and run_args['blur'] is not None:
            data.blur(run_args['blur'])
        elif 'q0' in run_args:
            blur_nm = 2*np.pi/run_args['q0']*run_args['blur_size_rel_d0']
            blur_pix = blur_nm/data.x_scale
            data.blur(blur_pix)

        if run_args['verbosity']>=5:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('blurred', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
        
        h, w = data.data.shape
        
        # Differentiate image in x and y directions
        dif_x = np.diff( data.data, axis=1 )
        line = np.zeros((h,1))
        dif_x = np.concatenate( (dif_x, line) , axis=1 )
        
        dif_y = np.diff( data.data, axis=0 )
        line = np.zeros((1,w))
        dif_y = np.concatenate( (dif_y, line) , axis=0 )
        
        numerator = 2.0*np.multiply( dif_x, dif_y )
        denominator = np.power(dif_x,2) - np.power(dif_y,2)
        

        if run_args['verbosity']>=4:
            display = Data2D()
            display.data = dif_x
            display.set_z_display( [None, None, 'linear', 1.0] )
            outfile = self.get_outfile('dif_x', output_dir, ext='.jpg', ir=True)
            dmax = max( np.abs(np.max(display.data)), np.abs(np.min(display.data)) )
            display.plot_image(save=outfile, zmin=-dmax, zmax=dmax, cmap=mpl.cm.seismic)

            display = Data2D()
            display.data = dif_y
            display.set_z_display( [None, None, 'linear', 1.0] )
            outfile = self.get_outfile('dif_y', output_dir, ext='.jpg', ir=True)
            dmax = max( np.abs(np.max(display.data)), np.abs(np.min(display.data)) )
            display.plot_image(save=outfile, zmin=-dmax, zmax=dmax, cmap=mpl.cm.seismic)

        
        
        if 'blur_orientation_image' in run_args and run_args['blur_orientation_image']:
            
            blur_nm = 2*np.pi/run_args['q0']*run_args['blur_orientation_image_size_rel']
            blur_pix = blur_nm/data.x_scale
            
            for i in range(run_args['blur_orientation_image_num_passes']):
                numerator = ndimage.filters.gaussian_filter( numerator, blur_pix )
                denominator = ndimage.filters.gaussian_filter( denominator, blur_pix )
            

        if run_args['verbosity']>=4:
            display = Data2D()
            display.data = numerator
            display.set_z_display( [None, None, 'linear', 1.0] )
            outfile = self.get_outfile('numerator', output_dir, ext='.jpg', ir=True)
            display.plot_image(save=outfile, cmap=mpl.cm.gnuplot2)

            display = Data2D()
            display.data = denominator
            display.set_z_display( [None, None, 'linear', 1.0] )
            outfile = self.get_outfile('denominator', output_dir, ext='.jpg', ir=True)
            display.plot_image(save=outfile, cmap=mpl.cm.gnuplot2)
            
            
            
        angles = 0.5*np.arctan2(numerator, denominator)
        # angles are in radians, from -pi/2 to +pi/2
        
        return angles

            
                            