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

import copy


class ProcessorIm(Processor):

    
    def load(self, infile, **kwargs):

        if 'load_param_file' in kwargs and kwargs['load_param_file']:
            # Load information about the image from the corresponding parameter file
            # For the Hitachi SEM, this is a txt file with the same name
            
            import re
            match_re = re.compile('^PixelSize=(\d+\.\d+)$')
            
            with open(infile[:-4]+'.txt') as fin:
                for line in fin.readlines():
                    m = match_re.match(line)
                    if m:
                        kwargs['scale'] = float(m.groups()[0])
            

        data = Data2DImage(infile, **kwargs)
        data.infile = infile
        
        if 'crop_edges' in kwargs:
            left, right, bottom, top = kwargs['crop_edges']
            data.crop_edges(left=left, right=right, bottom=bottom, top=top, relative=False)

        
        return data
        
        
        
        
class ProcessorImRGB(ProcessorIm):

    
    def load(self, infile, **kwargs):
        
        if 'defer_load' in kwargs and kwargs['defer_load']:
            data = Data2DImageRGB(**kwargs)
            data.name = tools.Filename(infile).get_filebase()
        else:
            data = Data2DImageRGB(infile, **kwargs)
        data.infile = infile
        
        if 'crop_edges' in kwargs:
            left, right, bottom, top = kwargs['crop_edges']
            data.crop_edges(left=left, right=right, bottom=bottom, top=top, relative=False)

        
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
        data = copy.deepcopy(data)
        
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
    
    
    
class preprocess(object):
    '''Base class to allow Protocols to run generic pre-processing
    routines on images.'''
    
    def preprocess(self, data, **run_args):
        # Usage: data = self.preprocess(data, **run_args)
        if 'preprocess' not in run_args:
            run_args['preprocess'] = 'default'
            
        if run_args['preprocess'] is None or run_args['preprocess']=='None':
            pass
        else:
            data = getattr(self, 'preprocess_{}'.format(run_args['preprocess']))(data, **run_args)
        
        return data


    def preprocess_default(self, data, **run_args):
        # Generic
        #data.equalize()
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        data.lowkill(run_args['q0']*0.1)
        data.blur(2.0) # lowpass
        #data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        #data.equalize()
        data.maximize_intensity_spread()
        
        return data
            

    def preprocess_AFM(self, data, **run_args):
        
        data.equalize()
        data.maximize_intensity_spread()
        
        return data

    def preprocess_SEM(self, data, **run_args):

        data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        data.blur(2.0) # lowpass
        data.blur(2.0) # lowpass
        data.equalize()
        data.maximize_intensity_spread()
        
        return data
        
    def preprocess_blur(self, data, **run_args):

        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(20):
            data.blur(2)
                            
        data.maximize_intensity_spread()
        
        return data

    def preprocess_blur_small(self, data, **run_args):

        for i in range(3):
            data.blur(0.5)
        data.maximize_intensity_spread()
        
        return data

    def preprocess_highloweq(self, data, **run_args):
        data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        
        data.blur(2.0) # lowpass
        data.blur(2.0) # lowpass
        data.equalize()
        data.maximize_intensity_spread()
        
        return data



class mask(object):
    '''Base shared class to allow Protocols to invoke a mask, which selects
    which regions of the image to analyze.'''

    def get_mask(self, data, output_dir='./', mask_threshold=127, mask_invert=False, **run_args):
        
        if 'mask' not in run_args or run_args['mask'] is None or run_args['mask'] in ['None', 'none',]:
            return None
            
        # Within a protocol, output_dir will be output_dir+self.name
        # We need to strip off the last part in order to get to the
        # output_dir being used for analysis results.
        output_dir = output_dir[:-len(self.name)]
            
        if run_args['mask']=='dots':
            infile = '{}/dots_vs_lines/{}/dots_vs_lines.png'.format(output_dir, data.name)
            mask_invert = not mask_invert
        
        elif run_args['mask']=='lines':
            infile = '{}/dots_vs_lines/{}/dots_vs_lines.png'.format(output_dir, data.name)
        
        else:
            infile = run_args['mask']
        
        
        mask_image = PIL.Image.open(infile).convert('I') # 'I' : 32-bit integer pixels
            
        mask = np.where( np.asarray(mask_image)>mask_threshold, 1, 0 ) # Binary array
        if mask_invert:
            mask = 1-mask
            
        return mask
    
    

class thumbnails(Protocol):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : 0.5,
                        'blur' : 1.0,
                        'resize' : 0.5,
                        'cmap' : mpl.cm.bone,
                        'preserve_data' : True,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['preserve_data']:
            # Avoid changing the data (which would disrupt downstream analysis of this data object)
            data = copy.deepcopy(data)
        
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
        
        
        
class thumb_view(thumbnails):
    '''Thumbnail variant intended to extract a square region of a 
    prescribed size (in nm).'''
    
    def __init__(self, name='thumb_view', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'size_nm' : 500,
                        'blur' : 1.0,
                        'resize' : 1.0,
                        'cmap' : mpl.cm.Greys_r,
                        'offset' : None,
                        'preserve_data' : True,
                        }
        self.run_args.update(kwargs)
        

    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['preserve_data']:
            # Avoid changing the data (which would disrupt downstream analysis of this data object)
            data = copy.deepcopy(data)
        
        # Crop out the correct region
        height, width = data.data.shape
        
        w_pix = run_args['size_nm']/data.x_scale
        h_pix = run_args['size_nm']/data.y_scale
        
        left = (width-w_pix)/2
        right = left
        bottom = (height-h_pix)/2
        top = bottom
        
        if run_args['offset'] is not None:
            left += run_args['offset'][0]
            right -= run_args['offset'][0]
            top += run_args['offset'][1]
            bottom -= run_args['offset'][1]
        
        data.crop_edges(left=left, right=right, bottom=bottom, top=top, relative=False)
        
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize']) # Shrink
        
        data.equalize()
        
        data.maximize_intensity_spread()
        #print(data.stats())
        
        data.set_z_display([None, None, 'linear', 1.0])
        outfile = self.get_outfile(data.name, output_dir)
        data.plot_image(outfile, zmin=0, zmax=255, **run_args)
        
        return results
    
    
    
class fft(Protocol, mask):
    
    def __init__(self, name='fft', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'blur' : None,
                        'Iqn_n' : 1.0,
                        'fourier_filter_invert' : False,
                        'fourier_filter_shift' : 0.3,
                        'q_max_multiple' : 0.25,
                        'q_max_q0_multiple' : 3.0,
                        }
        self.run_args.update(kwargs)


    def output_exists(self, name, output_dir):
        outfile = self.get_outfile('{}/01_fft'.format(name), output_dir, ext='.png')
        return os.path.isfile(outfile)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        run_args['mask'] = self.get_mask(data, output_dir=output_dir, **run_args)
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        
        
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
        
        
        
        if run_args['mask'] is not None:
            data = copy.deepcopy(data)
            avg = np.average(data.data)
            data.data = np.where(run_args['mask']==1, data.data, avg)
                
        
        
        if run_args['verbosity']>=5:
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
            
            if 'q0' in run_args:
                q_max = run_args['q0']*run_args['q_max_q0_multiple']
            else:
                q_max = np.max(x_axis)*run_args['q_max_multiple']
            
            data_fft.plot(save=outfile, ztrim=[0.5,0.0001], plot_range=[-q_max,+q_max,-q_max,+q_max], blur=run_args['blur'])
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
            lines.x_label = line.x_label
            lines.x_rlabel = line.x_rlabel
            lines.y_label = line.y_label
            lines.y_rlabel = line.y_rlabel
            
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
        params.add('x_center', value=np.average(line.x), min=0)
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



class local_d0(Protocol, preprocess):
    
    
    def __init__(self, name='local_d0', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'preprocess' : 'default',
                        'sub_region_size' : 75, # pixels
                        'sub_region_step_rel' : 0.5,
                        'blur' : 0.5,
                        }
        self.run_args.update(kwargs)

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)

        results = {}
        data = copy.deepcopy(data)
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)
        


        # Pre-process image
        data = self.preprocess(data, **run_args)
        
        
        if run_args['verbosity']>=5:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('enhanced', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)

        if 'd0' in run_args:
            q0 = 2*np.pi/run_args['d0']
        else:
            q0 = run_args['q0']
        dq = run_args['dq'] if 'dq' in run_args else q0*0.6
        s = 3
        plot_range = [-q0*s, +q0*s, -q0*s, +q0*s] if q0 is not None else [None,None,None,None]



        d0_x, d0_y, d0s = [], [], []

        sub_region_size = run_args['sub_region_size']
        sub_region_step = int(sub_region_size*run_args['sub_region_step_rel'])
        h, w = data.data.shape
        
        sub_image = Data2D()
        sub_image.x_scale, sub_image.y_scale = data.x_scale, data.y_scale
        
        
        
        for ix in range(sub_region_size, w-sub_region_size, sub_region_step):
            if run_args['verbosity']>=3:
                print('  ix {:d}/{:d} = {:.1f}%'.format(ix, w, 100.*ix/w))
            for iy in range(sub_region_size, h-sub_region_size, sub_region_step):
                
                sub_image.data = data.data[iy-sub_region_size:iy+sub_region_size , ix-sub_region_size:ix+sub_region_size]
                y_num, x_num = sub_image.data.shape
                
                if x_num>=sub_region_size*2 and y_num>=sub_region_size*2: # Full sub_images only
                    
                    if run_args['verbosity']>=5:
                        print('    sub_region ({:d}, {:d}) of size {:d}×{:d}'.format(ix, iy, x_num, y_num))


                    data_fft = sub_image.fft(update_origin=True)
                    
                    Re, Im = np.real(data_fft.data), np.imag(data_fft.data)
                    data_fft.data = np.abs(data_fft.data)
                    
                    if run_args['blur'] is not None:
                        data_fft.data = ndimage.filters.gaussian_filter(data_fft.data, run_args['blur'])
                        Re = ndimage.filters.gaussian_filter(Re, run_args['blur'])
                        Im = ndimage.filters.gaussian_filter(Im, run_args['blur'])
                        
                        
                    if run_args['verbosity']>=6:
                        outfile = self.get_outfile('sub_image_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        sub_image.plot_image(save=outfile, cmap=mpl.cm.bone, ztrim=[0,0])
                        
                        outfile = self.get_outfile('sub_FFT_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        data_fft.plot(save=outfile, plot_range=plot_range, ztrim=[0.25, 0.001], dpi=50)
                        
                    
                    outfile_extra = '_ix{:03d}iy{:03d}'.format(ix,iy)
                    q0local = self.analyze_q0(data_fft, q0i=q0, dqi=dq, outfile_extra=outfile_extra, output_dir=output_dir, **run_args)
                    d0local = 2*np.pi/q0local
                    
                    d0_x.append(ix)
                    d0_y.append(iy)
                    d0s.append(d0local)
                    
                    
        d0s = np.asarray(d0s)
        
        d0_values = np.column_stack([d0_x, d0_y, d0s])
        outfile = self.get_outfile('d0_local_values', output_dir, ext='.npy', ir=False)
        np.save(outfile, d0_values)
                    
        if run_args['verbosity']>=3:
            print_array(d0s, 'd0s')
            
            from scipy.interpolate import griddata
            grid = np.indices(data.data.shape)
            positions = np.column_stack((d0_y, d0_x))
            d0s_map = griddata(positions, d0s, tuple(grid), method='linear')
            
            d0min, d0max = np.min(d0s), np.max(d0s)
            d0s_map = (d0s_map-d0min)/(d0max-d0min)
            d0s_map = ndimage.filters.gaussian_filter(d0s_map, 1.0)
            print_array(d0s_map, 'd0s_map')
            outfile = self.get_outfile('d0s_interpolated', output_dir, ext='.npy', ir=False)
            np.save(outfile, d0s_map)
            
            d0s_interpolated = np.clip(d0s_map, 0, 1)
            cmap = mpl.cm.viridis
            img = PIL.Image.fromarray(np.uint8(cmap(d0s_interpolated)*255))
            outfile = self.get_outfile('d0s_interpolated', output_dir, ext='.png', ir=False)
            img.save(outfile)
                    
        
        
        return results



    def analyze_q0(self, data_fft, q0i, dqi, outfile_extra='', output_dir='./', **run_args):
        
        line = data_fft.circular_average()
        line.x_rlabel = '$q \, (\mathrm{nm}^{-1})$'
        
        sub_line = line.sub_range(q0i-dqi,q0i+dqi)
        lm_result, fit_line, fit_line_extended = self.fit_peak(sub_line, **run_args)
        q0 = lm_result.params['x_center'].value
        
        if run_args['verbosity']>=6:
            xpeak, ypeak = sub_line.target_y_max()
            lines = DataLines([line, fit_line, fit_line_extended])
            lines.copy_labels(line)
            outfile = self.get_outfile('fit{}'.format(outfile_extra), output_dir, ext='.png', ir=False)
            lines.plot(save=outfile, plot_range=[0, q0i*3, 0, ypeak*1.5], **run_args)
            
        
        return q0


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
        
        
        m = (line.y[-1]-line.y[0])/(line.x[-1]-line.x[0])
        b = line.y[0] - m*line.x[0]

        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y_max()
        
        params = lmfit.Parameters()
        params.add('prefactor', value=ypeak-(m*xpeak+b), min=0)
        params.add('x_center', value=xpeak, min=np.min(line.x), max=np.max(line.x))
        params.add('sigma', value=xspan/2, min=0)
        params.add('m', value=m)
        params.add('b', value=b)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=7:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        x_span = abs(np.max(line.x) - np.min(line.x))*0.5
        fit_x = np.linspace(np.min(line.x)-x_span, np.max(line.x)+x_span, num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':2.0})
            
        return lm_result, fit_line, fit_line_extended
                


    # End class fft(Protocol)
    ########################################



class local_avg_realspace(Protocol, preprocess):
    
    
    def __init__(self, name='local_avg_realspace', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        #'local_partition_image_size' : 75, # pixels
                        #'local_partition_step' : 1.0, # relative to image_size
                        'preprocess' : 'default',
                        'sub_region_size' : 75, # pixels
                        'sub_region_step_rel' : 0.5,
                        'blur' : 0.5,
                        }
        self.run_args.update(kwargs)


    def preprocess_default(self, data, **run_args):
        data.equalize()
        #data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        #data.lowkill(run_args['q0']*0.1)
        #data.blur(1.0)
        data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        data.maximize_intensity_spread()
        
        return data
            

        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)

        results = {}
        data = copy.deepcopy(data)
        
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)
        


        # Pre-process image
        data = self.preprocess(data, **run_args)
        
        
        if run_args['verbosity']>=5:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('enhanced', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)


        q0 = run_args['q0']
        dq = run_args['dq'] if 'dq' in run_args else q0*0.6
        s = 3
        plot_range = [-q0*s, +q0*s, -q0*s, +q0*s] if q0 is not None else [None,None,None,None]
                
        Z_accumluation = None
        Re_accumluation = None
        Im_accumluation = None
        
        realspace_accumulation = None
        num_images = 0
        
        sub_region_size = run_args['sub_region_size']
        sub_region_step = int(sub_region_size*run_args['sub_region_step_rel'])
        h, w = data.data.shape
        
        for ix in range(sub_region_size, w-sub_region_size, sub_region_step):
            for iy in range(sub_region_size, h-sub_region_size, sub_region_step):
                        
                sub_image = Data2D()
                sub_image.data = data.data[iy-sub_region_size:iy+sub_region_size , ix-sub_region_size:ix+sub_region_size]
                sub_image.x_scale, sub_image.y_scale = data.x_scale, data.y_scale
                y_num, x_num = sub_image.data.shape
                
                #if x_num>0 and y_num>0: # Images of size greater than 0
                if x_num>=sub_region_size*2 and y_num>=sub_region_size*2: # Full sub_images only
                    if run_args['verbosity']>=4:
                        print('    sub_region ({:d}, {:d}) of size {:d}×{:d}'.format(ix, iy, x_num, y_num))
                    
                    data_fft = sub_image.fft(update_origin=True)
                    
                    Re, Im = np.real(data_fft.data), np.imag(data_fft.data)
                    data_fft.data = np.abs(data_fft.data)
                    
                    if run_args['blur'] is not None:
                        data_fft.data = ndimage.filters.gaussian_filter(data_fft.data, run_args['blur'])
                        Re = ndimage.filters.gaussian_filter(Re, run_args['blur'])
                        Im = ndimage.filters.gaussian_filter(Im, run_args['blur'])
                        
                        
                    if run_args['verbosity']>=5:
                        outfile = self.get_outfile('sub_image_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        sub_image.plot_image(save=outfile, cmap=mpl.cm.bone, ztrim=[0,0])
                        
                        outfile = self.get_outfile('sub_FFT_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        data_fft.plot(save=outfile, plot_range=plot_range, ztrim=[0.25, 0.0005], dpi=50)
                    
                    line = data_fft.linecut_angle(d_center=q0, d_spread=dq, absolute_value=True)
                    
                    if run_args['verbosity']>=5:
                        outfile = self.get_outfile('sub_angle_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        line.plot(save=outfile, plot_range=[-180,+180,None,None])
                        outfile = self.get_outfile('sub_angle_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.dat', ir=False)
                        line.save_data(outfile)                        

                    if run_args['verbosity']>=10:
                        # Show where the linecut is being applied
                        data_fft.x_scale, data_fft.y_scale = 1, 1 # Hack because currently show_region assumes no coordinates have been applied
                        data_fft.origin = [0,0]
                        #plot_range = [None,None,None,None]
                        y_num, x_num = data_fft.data.shape
                        ws = 30
                        plot_range_c = [(x_num-ws)/2,(x_num+ws)/2,(y_num-ws)/2,(y_num+ws)/2]
                        data_fft.plot(save=False, plot_range=plot_range_c, dpi=50, show_region=True, show=True)
                        
                        
                    # Rotate the FFT to align the maximum
                    xm, ym = line.target_y_max()
                    Zabs = ndimage.interpolation.rotate(data_fft.data, -xm, reshape=False)
                    Re = ndimage.interpolation.rotate(Re, -xm, reshape=False)
                    Im = ndimage.interpolation.rotate(Im, -xm, reshape=False)
                    realspace = ndimage.interpolation.rotate(sub_image.data, -xm, reshape=False)
                    

                    if run_args['verbosity']>=5:
                        outfile = self.get_outfile('rot_FFT_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        Z = data_fft.data
                        data_fft.data = Zabs
                        data_fft.plot(save=outfile, plot_range=plot_range, ztrim=[0.25, 0.0005], dpi=50)
                        line2 = data_fft.linecut_angle(d_center=q0, d_spread=dq, absolute_value=True)
                        outfile = self.get_outfile('rot_angle_ix{:03d}iy{:03d}'.format(ix,iy), output_dir, ext='.jpg', ir=False)
                        line2.plot(save=outfile, plot_range=[-180,+180,None,None])
                        data_fft.data = Z


                    # Accumulate the rotated data
                    if Z_accumluation is None:
                        Z_accumluation = Zabs
                        Re_accumluation = Re
                        Im_accumluation = Im
                        realspace_accumulation = realspace
                    else:
                        Z_accumluation += Zabs
                        Re_accumluation += Re
                        Im_accumluation += Im
                        realspace_accumulation += realspace
                    num_images += 1
            
            
            
        Z_accumluation /= num_images
        Re_accumluation /= num_images
        Im_accumluation /= num_images
        realspace_accumulation /= num_images
            
        if run_args['verbosity']>=5:
            sub_image.data = realspace_accumulation
            outfile = self.get_outfile('avg_realspace', output_dir, ext='.jpg', ir=True)
            sub_image.plot_image(save=outfile, ztrim=[0,0], cmap=mpl.cm.bone)
            
        if run_args['verbosity']>=3:
            data_fft.data = Z_accumluation
            outfile = self.get_outfile('avg_FFT_abs', output_dir, ext='.jpg', ir=True)
            data_fft.plot(save=outfile, plot_range=plot_range, ztrim=[0.25,0.001])
            data_fft.data = Re_accumluation
            outfile = self.get_outfile('avg_FFT_Re', output_dir, ext='.jpg', ir=True)
            data_fft.plot(save=outfile, plot_range=plot_range, ztrim=[0.1,0.01])
            data_fft.data = Im_accumluation
            outfile = self.get_outfile('avg_FFT_Im', output_dir, ext='.jpg', ir=True)
            data_fft.plot(save=outfile, plot_range=plot_range, ztrim=[0.1,0.01])
            
        
        outfile = self.get_outfile('avg_realspace', output_dir, ext='.npy', ir=False)
        np.save(outfile, realspace_accumulation)
        outfile = self.get_outfile('avg_FFT_abs', output_dir, ext='.npy', ir=False)
        np.save(outfile, Z_accumluation)
        outfile = self.get_outfile('avg_FFT_Re', output_dir, ext='.npy', ir=False)
        np.save(outfile, Re_accumluation)
        outfile = self.get_outfile('avg_FFT_Im', output_dir, ext='.npy', ir=False)
        np.save(outfile, Im_accumluation)
        
        if run_args['verbosity']>=2:
            data_fft.data = Z_accumluation
            data_fft.recenter() # Put FFT origin back into corners
            realspace = Data2D()
            realspace.data = np.abs( np.fft.ifftn( data_fft.data ) )
            realspace.recenter(update_origin=True)
            y_num, x_num = data_fft.data.shape
            realspace.x_scale, realspace.y_scale = 2*np.pi/(data_fft.x_scale*x_num), 2*np.pi/(data_fft.x_scale*y_num)
            outfile = self.get_outfile('local_avg_realspace', output_dir, ext='.png', ir=False)
            realspace.x_rlabel, realspace.y_rlabel = r'$x \, (\mathrm{nm})$', r'$y \, (\mathrm{nm})$'
            realspace.plot(save=outfile, cmap='gist_heat', ztrim=[0.05, 0.03])
        
                    
                    
        return results

        
        
    def __deprecated_working_code(self):
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
    
    
    
    # End local_avg_realspace(Protocol, preprocess):
    ########################################


        
    
    
    
class particles(Protocol, preprocess):
    
    def __init__(self, name='particles', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        'preprocess' : 'default',
                        'hist_bins' : 100,
                        }
        self.run_args.update(kwargs)


    def preprocess_custom(self, data, **run_args):
        data.blur(0.6) # lowpass
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.blur(0.8) # lowpass
        data.equalize()
        data.maximize_intensity_spread()        
        
        return data



    def _find_objects(self, data, output_dir, results, **run_args):
        # results, labeled_array = self._find_objects(data, output_dir, results, **run_args)

        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
        
        
        # Pre-process image
        data = self.preprocess(data, **run_args)
        


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
        results['num_particles'] = num_features

        # Remove objects not meeting size criteria
        for i in range(np.max(labeled_array)+1):
            particle = (labeled_array==i).astype(np.uint8)
            area_pix = np.sum(particle)
            area_nm2 = area_pix*data.x_scale*data.y_scale # nm^2
            radius_nm = np.sqrt(area_nm2/np.pi) # nm

            # Select only desired particles
            exclude = False
            if 'area_min' in run_args and area_nm2<run_args['area_min']:
                exclude = True
            if 'area_max' in run_args and area_nm2>run_args['area_max']:
                exclude = True
            if 'radius_min' in run_args and radius_nm<run_args['radius_min']:
                exclude = True
            if 'radius_max' in run_args and radius_nm>run_args['radius_max']:
                exclude = True
            
            if exclude:
                idx = np.where(labeled_array==i)
                labeled_array[idx] = 0
            
        data.data = ( labeled_array>0 )*255

        if run_args['verbosity']>=4:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('excluded', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
            
        
        return results, labeled_array
    

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        data = copy.deepcopy(data)
        
        results, labeled_array = self._find_objects(data, output_dir, results, **run_args)


        if run_args['verbosity']>=3:
            # Colored image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
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
            
        
        if run_args['verbosity']>=4:
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
            print('    {} particles'.format(results['num_particles']))
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
        
        
        # Compute additional properties of of each particle
        import skimage.measure as measure
        analyzed_image = np.zeros((h,w,3))
        PrAs = []
        eccentricities = []
        for i in range(np.max(labeled_array)+1):
            
            if i>0: # Ignore the background object (index 0)
                
                # Extract just this particle
                #particle = np.where(labeled_array==i, 1, 0)
                particle = (labeled_array==i).astype(np.uint8)
                area_pix = np.sum(particle)
                area_nm2 = area_pix*data.x_scale*data.y_scale # nm^2
                radius_nm = np.sqrt(area_nm2/np.pi) # nm

                # Select only desired particles
                analyze = True
                if 'area_min' in run_args and area_nm2<run_args['area_min']:
                    analyze = False
                if 'area_max' in run_args and area_nm2>run_args['area_max']:
                    analyze = False
                if 'radius_min' in run_args and radius_nm<run_args['radius_min']:
                    analyze = False
                if 'radius_max' in run_args and radius_nm>run_args['radius_max']:
                    analyze = False


                h, w = labeled_array.shape
                

                if analyze:
                        
                    scale = (data.x_scale + data.y_scale)*0.5
                    perimeter_pix = measure.perimeter(particle)
                    perimeter_nm = perimeter_pix*scale
                    PrA = perimeter_nm*radius_nm/area_nm2
                    PrAs.append(PrA)
                
                    # Fit the particle to an ellipse
                    contour = measure.find_contours(particle, 0.5)[0]
                    ellipse = measure.EllipseModel()
                    ellipse.estimate(contour)
                    
                    if ellipse.params is None:
                        if run_args['verbosity']>=1:
                            print("WARNING: params is None for particle {:d}".format(i))
                    else:
                        xc, yc, a, b, theta = ellipse.params
                        if a>=b:
                            eccentricity = np.sqrt(1 - b**2/a**2)
                        else:
                            eccentricity = np.sqrt(1 - a**2/b**2)
                        eccentricities.append(eccentricity)

                        if run_args['verbosity']>=4:
                            print('    Particle {} ({} pixels)'.format(i, area_pix))
                            print('      A = {:.1f} nm^2; r = {:.1f} nm'.format(area_nm2, radius_nm))
                            print('      P = {:.1f} nm; P/A = {:.2g} 1/nm; Pr/A = {:.2f}'.format(perimeter_nm, perimeter_nm/area_nm2, PrA))
                            print('      e = {:.2f}'.format(eccentricity))

                        if run_args['verbosity']>=5:
                            analyzed_image += np.stack( (particle*255, particle*255, particle*255), axis=-1 )
                            xy = ellipse.predict_xy( np.linspace(0, 2*np.pi, 90) )
                            for y, x in xy:
                                if x>=0 and y>=0 and x<w and y<h:
                                    analyzed_image[int(y),int(x)] = [255, 0, 0]
                        
                        if run_args['verbosity']>=10:
                            # Output image of each particle separately (mostly for debugging)
                            outfile = self.get_outfile('particle{}'.format(i), output_dir, ext='.png', ir=False)
                            #import scipy.misc
                            #scipy.misc.toimage(particle*255).save(outfile) # Deprecated
                            particle_image = particle*255
                            Image.fromarray(particle_image.astype(np.uint8)).save(outfile)


        results['PrA_average'] = np.average(PrAs)
        results['PrA_std'] = np.std(PrAs)
        results['PrA_median'] = np.median(PrAs)

        results['eccentricity_average'] = np.average(eccentricities)
        results['eccentricity_std'] = np.std(eccentricities)
        results['eccentricity_median'] = np.median(eccentricities)


                    
        if run_args['verbosity']>=5:
            outfile = self.get_outfile('analyzed', output_dir, ext='.png', ir=True)
            #import scipy.misc
            #scipy.misc.toimage(analyzed_image).save(outfile) # Deprecated
            Image.fromarray(analyzed_image.astype(np.uint8)).save(outfile)
                
            
        
        
        if run_args['verbosity']>=1:
            
            class DataHistogram_current(DataHistogram):
                def _plot_extra(self, **plot_args):
                    
                    xi, xf, yi, yf = self.ax.axis()
                    yf *= 1.2
                    self.ax.axis([xi, xf, yi, yf])
                    
                    self._plot_stats(**plot_args)
                    
                    lm_result, fit_line, fit_line_e = self.fit_peak(self)
                    self.ax.axvline(lm_result.params['x_center'], color='r', linewidth=2.0)
                    self.ax.plot(fit_line_e.x, fit_line_e.y, 'r', linewidth=2.0)
                    
                    
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
            y, x = np.histogram(particle_sizes, bins=run_args['hist_bins'], range=[0, max(particle_sizes)*1.05])
            
            # Instead of having x be ranges for each bar, center the x on the average of each range
            xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
            x = x[:-1]
            
            hist = DataHistogram_current(x=x, y=y, x_label='Area', x_rlabel='$A \, (\mathrm{nm}^{2})$', y_label='count')
            hist.mean = results['area_average']
            hist.std = results['area_std']
            hist.median = results['area_median']
            
            outfile = self.get_outfile('particle_areas', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
            outfile = self.get_outfile('particle_areas', output_dir, ext='.dat')
            np.savetxt(outfile, np.transpose([x, y]))
            
            
            
            # Histogram of radii
            y, x = np.histogram(particle_radii, bins=run_args['hist_bins'], range=[0, max(particle_radii)*1.05])
            
            # Instead of having x be ranges for each bar, center the x on the average of each range
            xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
            x = x[:-1]
            
            hist = DataHistogram_current(x=x, y=y, x_label='Radius', x_rlabel='$r \, (\mathrm{nm})$', y_label='count')
            hist.mean = results['radius_average']
            hist.std = results['radius_std']
            hist.median = results['radius_median']
            
            outfile = self.get_outfile('particle_radii', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
            outfile = self.get_outfile('particle_radii', output_dir, ext='.dat')
            np.savetxt(outfile, np.transpose([x, y]))
            
           
            # Histogram of eccentricities
            y, x = np.histogram(eccentricities, bins=run_args['hist_bins'], range=[0, 1])
            
            # Instead of having x be ranges for each bar, center the x on the average of each range
            xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
            x = x[:-1]
            
            hist = DataHistogram_current(x=x, y=y, x_label='Eccentricity', x_rlabel='$\mathrm{eccentricity}$', y_label='count')
            hist.mean = results['eccentricity_average']
            hist.std = results['eccentricity_std']
            hist.median = results['eccentricity_median']
            
            outfile = self.get_outfile('eccentricities', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[0, 1, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
            outfile = self.get_outfile('eccentricities', output_dir, ext='.dat')
            np.savetxt(outfile, np.transpose([x, y]))
            
            
            # Histogram of PrAs
            y, x = np.histogram(PrAs, bins=run_args['hist_bins'], range=[2, max(PrAs)*1.05])
            
            # Instead of having x be ranges for each bar, center the x on the average of each range
            xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
            x = x[:-1]
            
            hist = DataHistogram_current(x=x, y=y, x_label='PrA', x_rlabel=r'$\frac{Pr}{A}$', y_label='count')
            hist.mean = results['PrA_average']
            hist.std = results['PrA_std']
            hist.median = results['PrA_median']
            
            outfile = self.get_outfile('PrA', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[2, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
            outfile = self.get_outfile('PrA', output_dir, ext='.dat')
            np.savetxt(outfile, np.transpose([x, y]))
                                    
        
        
        return results
        



class particles_annotated(particles):
    
    def __init__(self, name='particles', **kwargs):
        
        super().__init__(name=name, **kwargs)
        
        self.run_args['annotate_color'] = [1,0,0]
        self.run_args['annotate_threshold'] = 0.7
        self.run_args.update(kwargs)    
        
    
    def _find_objects(self, data, output_dir, results, **run_args):
        # results, labeled_array = self._find_objects(data, output_dir, results, **run_args)
        
        import scipy
        
        distance = np.linalg.norm(data.data_rgb - np.asarray(run_args['annotate_color'])*255.0, axis=-1)
        edges = ( distance<run_args['annotate_threshold']*255 )*255.0
        edges = edges.astype(np.uint8)
        
        filled = ( scipy.ndimage.morphology.binary_fill_holes(edges) )*255.0
        filled = filled.astype(np.uint8)

        
        if run_args['verbosity']>=6:
            outfile = self.get_outfile('distance', output_dir, ext='.png', ir=True)
            scipy.misc.toimage(distance).save(outfile)
        if run_args['verbosity']>=5:
            outfile = self.get_outfile('edges', output_dir, ext='.png', ir=True)
            scipy.misc.toimage(edges).save(outfile)
        if run_args['verbosity']>=4:
            outfile = self.get_outfile('filled', output_dir, ext='.png', ir=True)
            scipy.misc.toimage(filled).save(outfile)

        
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
        
        labeled_array, num_features = ndimage.measurements.label(filled, structure=s)
        results['num_particles'] = num_features
        
        
        return results, labeled_array
    
    
    
class skeletonize(particles):
    
    def __init__(self, name='skeletonize', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        'preprocess' : 'default',
                        'hist_bins' : 100,
                        }
        self.run_args.update(kwargs)


    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        orig_img = data.data
        data = copy.deepcopy(data)
        
        results, labeled_array = self._find_objects(data, output_dir, results, **run_args)
        
        results, skeleton = self._skeletonize(data, output_dir, results, labeled_array, **run_args)
        
        return results


    def _skeletonize(self, data, output_dir, results, labeled_array, **run_args):
        
        if run_args['verbosity']>=4:
            # Colored image
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
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

        # Statistics on particles that have been found
        bins = np.bincount(labeled_array.flatten('C'))
        h, w = data.data.shape
        total_pixels = w*h
        background_pixels = bins[0]
        particle_pixels = total_pixels - background_pixels
        coverage = particle_pixels*1./(total_pixels*1.)
        results['coverage'] = coverage
        if run_args['verbosity']>=4:
            print('    {} particles'.format(results['num_particles']))
            print('    Particle coverage: {:.1f}%'.format(coverage*100.))
                    
        # Remove 'particles' of zero size
        idx = np.nonzero(bins)
        bins = bins[idx]
        # Remove the 'surrounding field' (index 0)
        bins = bins[1:]
        
        # Exclude particles
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
        
        
        from skimage.morphology import skeletonize
        data.data = np.where(data.data>127, 1, 0)
        skeleton = skeletonize(data.data).astype(int)
        data.data = skeleton
        
        results['skeleton_coverage'] = np.sum(data.data)/total_pixels



        if run_args['verbosity']>=3:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('skeleton', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, zmin=0, zmax=1, cmap=mpl.cm.binary_r)
        
        
        #s = [[1,1,1],
        #    [1,1,1],
        #    [1,1,1]]
        s = ndimage.generate_binary_structure(2,2) # Detect along diagonals
        labeled_array, num_features = ndimage.measurements.label(data.data, structure=s)
        results['num_lines'] = num_features

        if run_args['verbosity']>=4:
            print('    {} skeleton lines'.format(num_features))
            print('    Skeleton coverage: {:.1f}%'.format(results['skeleton_coverage']*100.))

        if run_args['verbosity']>=5:
            # Colored image
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
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
            
        return results, skeleton


class skeleton_lines(skeletonize):
    
    def __init__(self, name='skeleton_lines', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        'preprocess' : 'default',
                        'hist_bins' : 100,
                        }
        self.run_args.update(kwargs)


    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        orig_img = data.data
        data = copy.deepcopy(data)
        
        results, labeled_array = self._find_objects(data, output_dir, results, **run_args)
        
        results, skeleton = self._skeletonize(data, output_dir, results, labeled_array, **run_args)
        
        results = self._skeletonize_lines(data, output_dir, results, skeleton, orig_img, **run_args)
        
        return results


    def _skeletonize_lines(self, data, output_dir, results, skeleton, orig_img, **run_args):


        # Neighbor image
        def filter_function(values):
            return values.sum()
        footprint = np.array([[1,1,1],
                            [1,0,1],
                            [1,1,1]])
        
        neighbors = ndimage.generic_filter(skeleton, filter_function, footprint=footprint, mode='constant')
        neighbors *= skeleton
        
        # Count certain kinds of defects
        results['line_ends_count'] = len(neighbors[neighbors==1])
        results['line_ends_density'] = results['line_ends_count']/data.image_area()
        
        junctions = np.where(neighbors>=3, 1, 0)
        s = ndimage.generate_binary_structure(2,2) # Detect along diagonals
        labeled_junctions, num_junctions = ndimage.measurements.label(junctions, structure=s)
        results['junctions_count'] = num_junctions
        results['junctions_density'] = results['junctions_count']/data.image_area()

        if run_args['verbosity']>=3:
            print('    {:d} line ends ({:.3g}/nm^2 = {:.2g}/μm^2)'.format(results['line_ends_count'], results['line_ends_density'], results['line_ends_density']*1e6))
            print('    {:d} junctions ({:.3g}/nm^2 = {:.2g}/μm^2))'.format(results['junctions_count'], results['junctions_density'], results['junctions_density']*1e6))

        if run_args['verbosity']>=6:
            #print_array(neighbors, 'neighbors')
            data.data = neighbors
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('neighbors', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, zmin=0, zmax=5, cmap=mpl.cm.jet)
            data.data = skeleton
            
        if run_args['verbosity']>=6:
            data.data = junctions
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('junctions', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, zmin=0, ztrim=[0,0], cmap=mpl.cm.binary_r)
            data.data = skeleton

        if run_args['verbosity']>=3:
            from PIL import ImageDraw
            # Two equivalent ways to extend the grayscale array into a RGB array
            #analyzed_image = np.repeat(orig_img[:, :, np.newaxis], 3, axis=2)
            analyzed_image = np.dstack([orig_img]*3)
            image = Image.fromarray(analyzed_image.astype(np.uint8))
            draw = ImageDraw.Draw(image)

            # Mark ends of lines
            line_end_idx = np.argwhere(neighbors==1)
            if 'd0' in run_args:
                r = 0.5*run_args['d0']/data.get_scale()
            elif 'q0' in run_args:
                r = 0.5*(2*np.pi/run_args['q0'])/data.get_scale()
            else:
                r = 10 # pixels
            for y, x in line_end_idx:
                draw.ellipse((x-r, y-r, x+r, y+r), outline=(255,0,0,0), fill=None)
                
            # Mark junctions
            import scipy
            com = scipy.ndimage.measurements.center_of_mass(junctions, labels=labeled_junctions, index=range(1,num_junctions+1))
            for y, x in com:
                draw.ellipse((x-r, y-r, x+r, y+r), outline=(0,0,255,0), fill=None)
            
            outfile = self.get_outfile('defects', output_dir, ext='.png', ir=True)
            image.save(outfile)

        if run_args['verbosity']>=6:
            # A black-and-white image that localizes areas with defects
            image = Image.fromarray(np.zeros_like(analyzed_image).astype(np.uint8))
            draw = ImageDraw.Draw(image)

            rc = r*3
            for y, x in line_end_idx:
                draw.ellipse((x-rc, y-rc, x+rc, y+rc), fill=(255,255,255,0))
            
            rc = r*4
            for y, x in com:
                draw.ellipse((x-rc, y-rc, x+rc, y+rc), fill=(255,255,255,0))
            
            outfile = self.get_outfile('defect_regions', output_dir, ext='.png', ir=True)
            image.save(outfile)            
            
        return results

        
class skeleton_angles(skeleton_lines):
    
    def __init__(self, name='skeleton_angles', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        'preprocess' : 'default',
                        'hist_bins' : 30,
                        'lookout_distance_pix' : 6,
                        }
        self.run_args.update(kwargs)


    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        orig_img = data.data
        data = copy.deepcopy(data)
        
        results, labeled_array = self._find_objects(data, output_dir, results, **run_args)
        
        results, skeleton = self._skeletonize(data, output_dir, results, labeled_array, **run_args)
        
        results = self._skeletonize_lines(data, output_dir, results, skeleton, orig_img, **run_args)
        
        results = self._skeletonize_angles(data, output_dir, results, labeled_array, skeleton, **run_args)
        
        return results        


    def _skeletonize_angles(self, data, output_dir, results, labeled_array, skeleton, **run_args):
        '''Use the computed skeleton, and calculate at each point along this skeleton the local
        "bending angle". To do this, for each pixel we look out to a local window. The bright
        pixels at the edge of this window are assumed to be extensions of our current line.
        (This is valid if the lines are thin compared to how close they are to each other.)
        Then we calculate the angles between these edge-pixels. Since we are centered on a
        particular line-pixel, these angle different between edge-pixels is the amount of
        bending.'''
        
        
        # Compute the local "bending" angle of skeleton lines
        # out to a "lookout" distance of l.
        l = run_args['lookout_distance_pix']
        s = 2*l+1 # Size of window
        
        # We define a "footprint" for the window that will search at any given point.
        # For instance, a footprint that searches 2 pixels away from the target point is a 5x5 matrix:
        #footprint = [ [1, 1, 1, 1, 1] ,
        #              [1, 0, 0, 0, 1] ,
        #              [1, 0, 0, 0, 1] ,
        #              [1, 0, 0, 0, 1] ,
        #              [1, 1, 1, 1, 1] ]
        footprint = np.ones((s,s)) # Set all values (including edges) to 1.0
        footprint[1:-1,1:-1] = 0 # Set inner (non-edge) values to 0.0

        # We pre-compute angles for the footprint, which has indices like:
        #  0  1  2  3  4
        #  5  -  -  -  6
        #  7  -  -  -  8
        #  9  -  -  - 10
        # 11 12 13 14 15
        #delta_x = [-2, -1, +0, +1, +2,    -2, +2, -2, +2, -2, +2,    -2, -1, +0, +1, +2]
        #delta_y = [+2, +2, +2, +2, +2,    +1, +1, +0, +0, -1, -1,    -2, -2, -2, -2, -2]
        grid = np.indices((s,s))
        delta_x = grid[1]-l
        delta_x = delta_x[footprint==1]
        delta_y = (grid[0]-l)*-1
        delta_y = delta_y[footprint==1]
        angle_lookup = np.degrees(np.arctan2(delta_y, delta_x))
        
        def filter_function(values):
            '''Calculates the local bending angle at this pixel.'''
            idx = np.where(values>0)[0]
            n = len(idx)
            if n<=1:
                return 0
            elif n==2:
                angle0 = angle_lookup[idx[0]]
                angle1 = angle_lookup[idx[1]]
                angle_diff = abs(angle1-angle0)
                angle_diff = min(angle_diff, 360-angle_diff)
                bend = 180-angle_diff
                return bend
            elif n==3:
                angle0 = angle_lookup[idx[0]]
                angle1 = angle_lookup[idx[1]]
                angle2 = angle_lookup[idx[2]]
                angle_diff01 = abs(angle1-angle0)
                angle_diff01 = min(angle_diff01, 360-angle_diff01)
                angle_diff12 = abs(angle2-angle1)
                angle_diff12 = min(angle_diff12, 360-angle_diff12)
                angle_diff20 = abs(angle0-angle2)
                angle_diff20 = min(angle_diff20, 360-angle_diff20)
                angles = sorted([angle_diff01, angle_diff12, angle_diff20])[:n-1]
                bends = 180-np.asarray(angles)
                return np.average(bends)
            
            else:
                # This block works for all n, but the above code might be slightly more efficient
                # (and is retained to document how this method works)
                angles = []
                for i in range(n):
                    angle_i = angle_lookup[idx[i]]
                    for j in range(i+1, n):
                        angle_j = angle_lookup[idx[j]]
                        diff = abs(angle_i-angle_j)
                        diff = min(diff, 360-diff)
                        angles.append(diff)
                angles = sorted(angles)[:n-1] # Exclude the largest outer angle
                bends = 180-np.asarray(angles)
                return np.average(bends)
                
        
        # Apply the local "angle calculator" to every pixel in the image.
        angles = ndimage.generic_filter(skeleton, filter_function, footprint=footprint, mode='constant', output=np.float)
        angles *= skeleton # Set to zero any non-skeleton pixels
        if run_args['verbosity']>=6:
            print_array(angles, 'angles')
            
        
        if run_args['verbosity']>=3:
            angles[skeleton==0] = -10 # We can set the non-skeleton values to some other number to help with visualizing.
            data.data = angles
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('angles', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, zmin=-10, zmax=90, _ztrim=[0,0], cmap=mpl.cm.jet)
            #data.plot(show=True, zmin=-90, zmax=90, _ztrim=[0,0], cmap=mpl.cm.jet)
            data.data = skeleton
            
            
        # Select only the valid pixels (along skeleton)
        angles = angles[skeleton==1]
            
            
        if run_args['verbosity']>=6:
            from scipy.interpolate import griddata
            grid = np.indices(skeleton.shape)
            positions = np.column_stack((grid[0][skeleton==1],grid[1][skeleton==1]))
            angle_map = griddata(positions, angles, tuple(grid), method='linear')
            
            angle_map = ndimage.filters.gaussian_filter(angle_map, l/4)
            outfile = self.get_outfile('angles_interpolated', output_dir, ext='.npy', ir=False)
            np.save(outfile, angle_map)
            
            angle_map = np.clip(angle_map/90, 0, 1)
            cmap = mpl.cm.jet
            img = PIL.Image.fromarray(np.uint8(cmap(angle_map)*255))
            outfile = self.get_outfile('angles_interpolated', output_dir, ext='.png', ir=True)
            img.save(outfile)
        
        
        
        # Calculate stats and histogram
        results['bending_angle_average'] = np.average(angles)
        results['bending_angle_std'] = np.std(angles)
        results['bending_angle_median'] = np.median(angles)
        
        y, x = np.histogram(angles, bins=run_args['hist_bins'], range=[0, 180])
        x = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ] # bin_edges to bin centers
        outfile = self.get_outfile('bending_angles', output_dir, ext='.dat')
        np.savetxt(outfile, np.stack([x, y], axis=1))
        

        if run_args['verbosity']>=2:
            class DataHistogram_current(DataHistogram):
                def _plot_extra(self, **plot_args):
                    xi, xf, yi, yf = self.ax.axis()
                    yf = np.max(self.y[1:])*0.25
                    self.ax.axis([xi, xf, yi, yf])
                    
                    self._plot_stats(**plot_args)
                    
                    self.ax.set_xticks(range(0, 120+30, 30))

            
            hist = DataHistogram_current(x=x, y=y, x_label='angle', x_rlabel=r'$\theta_{\mathrm{bend}} \, (^{\circ})$', y_label='count')
            hist.mean = results['bending_angle_average']
            hist.std = results['bending_angle_std']
            hist.median = results['bending_angle_median']
            
            outfile = self.get_outfile('bending_angles', output_dir, ext='.png')
            hist.plot(save=outfile, plot_range=[0, 120, 0, None], plot_buffers=[0.18,0.05,0.18,0.05],)

        
        return results


            
class grain_size_hex(Protocol, preprocess, mask):
    
    def __init__(self, name='grain_size_hex', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'symmetry': 6,
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        'correlation_edge_exclusion' : 10,
                        'correlation_step_size_points' : 5,
                        'trim_r_curve' : 0.8, # 1.0 doesn't trim anything; 0.8 trims the last 20% of the g(r)-curve
                        'preprocess' : 'default',
                        }
        self.run_args.update(kwargs)


    def preprocess_default(self, data, **run_args):
        #data.equalize()
        data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        #data.lowkill(run_args['q0']*0.1)
        
        data.blur(2.0) # lowpass
        data.blur(2.0) # lowpass
        #data.blur(1.0) # lowpass
        #data.blur(0.6) # lowpass
        #data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        data.equalize()
        data.maximize_intensity_spread()
        
        return data

    def preprocess_custom(self, data, **run_args):
        data.blur(0.6) # lowpass
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.blur(0.8) # lowpass
        data.equalize()
        data.maximize_intensity_spread()
        
        return data        
    
    def preprocess_SEM(self, data, **run_args):
        data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        data.blur(2.0) # lowpass
        data.blur(2.0) # lowpass
        data.equalize()
        data.maximize_intensity_spread()
        
        return data

    @run_default
    def run(self, data, output_dir, **run_args):
        
        run_args['mask'] = self.get_mask(data, output_dir=output_dir, **run_args)
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        
        
        orig_data = data.data.copy()
        
        results = {}
        data = copy.deepcopy(data)
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])

        
        # Pre-process image
        data = self.preprocess(data, **run_args)


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


        # Compute defect density
        h, w = data.data.shape
        im_area = (w*data.x_scale)*(h*data.y_scale)*1e-6 # um^2
        
        # Number of defects due to image edges (these are artifacts):
        L0 = 2*np.pi/run_args['q0']
        d = L0/(np.sqrt(3.0)/2.0)
        edge_defects = 2*w*data.x_scale/d # top/bottom edge
        edge_defects += 2*2*h*data.y_scale/L0 # left/right edge
        
        sym = run_args['symmetry']
        less = len(np.nonzero(NN_counts<sym)[0])
        equal = len(np.nonzero(NN_counts==sym)[0])
        more = len(np.nonzero(NN_counts>sym)[0])
        
        results['fraction_matching_symmetry'] = equal/(less+equal+more - edge_defects)
        results['number_defects'] = int(max( ( (less + more) - edge_defects ), 0 ))
        results['defect_density'] = { 'value': results['number_defects']/im_area, 'units': 'um^-2'}
        if run_args['verbosity']>=2:
            print("  {:d} (non-edge) defects / {:.2f} um^2 = ".format(results['number_defects'], im_area))
            print("  {:.4g} defects/um^2".format(results['defect_density']['value']))
            print("  f_{}NN = {:.3f}".format(run_args['symmetry'], results['fraction_matching_symmetry']))
        
        
        
        
        
            
        if run_args['verbosity']>=5:
            # Colored image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
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
                            c = (0, 255, 0) # Green
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


        # Compute angle map
        from scipy.interpolate import griddata
        positions = np.column_stack((x_positions,y_positions))
        grid_x, grid_y = np.mgrid[ 0:len(labeled_array) , 0:len(labeled_array[0]) ]
        angle_map = griddata(positions, angles, (grid_y, grid_x), method='nearest') # Avoids artifacts

        if run_args['verbosity']>=3:
            # False-color map of angles
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
            
            
        # Compute correlation function.
        if 'scale' not in run_args:
            run_args['scale'] = (data.x_scale+data.y_scale)*0.5
        new_results = self.correlation_function(angle_map, output_dir=output_dir, **run_args)
        results.update(new_results)

            
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
        
        outfile = self.get_outfile('particle_areas', output_dir, ext='.png', ir=True)
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
        
        outfile = self.get_outfile('particle_radii', output_dir, ext='.png', ir=True)
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
        if run_args['verbosity']>=2:
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
            r = 360.0/run_args['symmetry']
            if np.min(line.x)>=0:
                plot_range = [0, +r, 0, np.max(line.y)*1.2]
            else:
                plot_range = [-r/2, +r/2, 0, np.max(line.y)*1.2]
            lines.plot(save=outfile, plot_range=plot_range)

        if run_args['verbosity']>=3:
            outfile = self.get_outfile('orientation_polar', output_dir, ext='.png', ir=True)
            line.plot_polar(save=outfile, assumed_symmetry=run_args['symmetry'], symmetry_copy=True)
            
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



    def correlation_function(self, angles, output_dir='./', fit_curve=True, **run_args):
        '''Accumluate the pair-wise correlations into a 1D correlation curve.
        The decay of this curve reveals a characteristic correlation distance.
        This function computes the curve, fits it to an exponential decay, and
        plots the results.'''
        
        results = {}
        
        ex = run_args['correlation_edge_exclusion']
        step = run_args['correlation_step_size_points']
        symmetry = run_args['symmetry']
        scale = run_args['scale'] # nm/pixel
        
        h, w = angles.shape
        
        accumulator = np.zeros( (h*2,w*2) ) # Sum of correlation values
        count_accumulator = np.zeros( (h*2,w*2) ) # Keep track of counts (for normalization)
        one_field = np.ones( (h,w) )
        
        if 'mask' not in run_args or run_args['mask'] is None:
            # Normal analysis
            for ix in range(ex, w-ex, step):
                if run_args['verbosity']>=3:
                    print( '        Correlation analysis {:.1f}%'.format (100.*ix/w))
                for iy in range(ex, h-ex, step):
                    deltatheta = angles[iy,ix] - angles
                    orderparameter = np.cos(symmetry*deltatheta)
                    
                    accumulator[h-iy:-iy, w-ix:-ix] += orderparameter
                    count_accumulator[h-iy:-iy, w-ix:-ix] += one_field
                    
        else:
            # Analysis that ignores pixels excluded by the mask
            mask = run_args['mask']
            for ix in range(ex, w-ex, step):
                if run_args['verbosity']>=3:
                    print( '        Correlation analysis {:.1f}%'.format (100.*ix/w))
                for iy in range(ex, h-ex, step):
                    if mask[iy,ix]==1:
                        deltatheta = angles[iy,ix] - angles
                        orderparameter = np.cos(symmetry*deltatheta)
                        
                        accumulator[h-iy:-iy, w-ix:-ix] += orderparameter*mask
                        count_accumulator[h-iy:-iy, w-ix:-ix] += one_field*mask

        # Compute array of distances associated with the accumulator
        r_dist = np.zeros( (h*2,w*2) )
        ix_list = np.asarray( range(0, w*2) )
        for iy in range(0, h*2):
            r_dist[iy,:] = np.round( np.sqrt( (ix_list-w)**2 + (iy-h)**2 ) )

        
        
        # Create a g(r) curve by summing up the values in the accumulator
        # TODO: This section could be improved by using various numpy functions
        # (where, bincount, etc.)
        r_list = range( 0, int(np.max(r_dist)) )
        g_of_r = np.zeros( (len(r_list)) )
        count_list = np.zeros( (len(r_list)) )
        
        for iy in range(0, h*2):
            for ix in range(0, w*2):
                
                if count_accumulator[iy,ix]>0:
                    r = int( r_dist[iy,ix] )
                    g_of_r[r] += accumulator[iy,ix]
                    count_list[r] += count_accumulator[iy,ix]
        
        
        r_list_final = []
        r_nm_list_final = []
        g_of_r_final = []
        for r in range(len(r_list)):
            if count_list[r]>0 and not np.isnan(g_of_r[r]):
                
                # Convert from pixels to nm
                r_list_final.append( r_list[r] )
                r_nm_list_final.append( scale*r_list[r] )
                g_of_r_final.append( g_of_r[r]/count_list[r] )
                
                
        line = DataLine(x=r_nm_list_final, y=g_of_r_final)
        
        # Remove some of the line (the final values are not meaningful since they average over so few points)
        line.trim(0, np.max(line.x)*run_args['trim_r_curve'])

        if fit_curve:
            # Fit correlation curve
            
            # Get the 1/e distance
            x_1e, y_1e = line.target_y(1.0/np.e)
            results['xi_1_over_e_nm'] = x_1e
            
            lm_result, fit_line, fit_line_extended = self._fit_exp(line, **run_args)
            
            fit_name = 'fit_exp'
            for param_name, param in lm_result.params.items():
                results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            
            results['xi_nm'] = results['fit_exp_xi']['value']


            lm_result, fit_lineb, fit_line_extendedb = self._fit_expb(line, **run_args)
            fit_name = 'fit_expb'
            for param_name, param in lm_result.params.items():
                results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            results['xib_nm'] = results['fit_expb_xi']['value']
            
            
        if run_args['verbosity']>=3:
            # Plot curve
            
            class DataLines_current(DataLines):
                
                def _plot_extra(self, **plot_args):
                    self.ax.axhline(0, color='k')
                    
                    r, g, b = 0, 0, 1
                    self.ax.plot(x_1e, y_1e, 'o', color='b', markersize=20, markerfacecolor=(r,g,b,0.75), markeredgewidth=1.5, markeredgecolor=(r,g,b,1.0), zorder=10)
                    
                    if fit_curve:
                        self.ax.axvline(x_1e, color='0.5', dashes=[5,5], linewidth=1.0, zorder=-10)
                        s = r'$\xi_{{1/e}} = {:.1f} \, \mathrm{{nm}}$'.format(x_1e)
                        self.ax.text(x_1e, 1.0, s, size=20, color='0.5', verticalalignment='top', horizontalalignment='left')
                        
                        xi, xf, yi, yf = self.ax.axis()
                        s = r'$g(r) = e^{{-r/\xi}}$' + '\n' + r'$\xi = {:.1f} \pm {:.1f} \, \mathrm{{nm}}$'.format(results['fit_exp_xi']['value'], results['fit_exp_xi']['error'])
                        self.ax.text(xf, 1.0, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')


                        s = r'$g(r) = e^{{-r/\xi_b}} + b$' + '\n' + r'$\xi_b = {:.1f} \pm {:.1f} \, \mathrm{{nm}}$'.format(results['fit_expb_xi']['value'], results['fit_expb_xi']['error'])
                        s += '\n' + r'$b = {:.2f} \pm {:.2f}$'.format(results['fit_expb_b']['value'], results['fit_expb_b']['error'])
                        self.ax.text(xf, 0.75, s, size=20, color='purple', verticalalignment='center', horizontalalignment='right')
                        
                    
            lines = DataLines_current()
            
            
            lines.add_line(line)
            
            if fit_curve:
                lines.add_line(fit_line_extended)
                lines.add_line(fit_line_extendedb)
            
            
            lines.x_label = 'r'
            lines.x_rlabel = '$r \, \mathrm{(nm)}$'
            lines.y_label = 'g(r)'
            lines.y_rlabel = '$g(r)$'
            
            outfile = self.get_outfile('correlation', output_dir, ext='.png', ir=True)
            lines.plot(save=outfile, plot_range=[0, None, None, 1.0])
        
        
        
        return results
    
    def _fit_exp(self, line, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self._fit_exp(line, **run_args)
        
        x_1e, y_1e = line.target_y(1.0/np.e)

        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
            
        line.x = np.asarray(line.x)
        line.y = np.asarray(line.y)
        
        import lmfit
        
        def model(v, x):
            return np.exp( -x/v['xi'] )
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('xi', value=x_1e, min=0, max=np.max(line.x)*+10, vary=True)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        #fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=200)
        fit_x = np.linspace(0, np.max(line.x)*1.1, num=500)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'alpha':0.75, 'marker':None, 'linewidth':3.0})        

        return lm_result, fit_line, fit_line_extended             
    
    
    def _fit_expb(self, line, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self._fit_expb(line, **run_args)


        amt = int( len(line.y)*0.25 )
        b = np.average(line.y[-amt:])
        
        
        linet = line.copy()
        linet.y = (linet.y - b)/(1-b)
        x_1e, y_1e = linet.target_y(1.0/np.e)

        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
            
        line.x = np.asarray(line.x)
        line.y = np.asarray(line.y)
        
        import lmfit
        
        def model(v, x):
            return np.exp( -x/v['xi'] ) + v['b']
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('b', value=b, min=0, max=np.max(line.y)*0.8, vary=True)
        params.add('xi', value=x_1e, min=0, max=np.max(line.x)*+10, vary=True)
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'purple', 'marker':None, 'linewidth':4.0})
        
        #fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=200)
        fit_x = np.linspace(0, np.max(line.x)*1.1, num=500)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'purple', 'alpha':0.75, 'marker':None, 'linewidth':3.0})        

        return lm_result, fit_line, fit_line_extended    


    def orientation_angle_map(self, data, output_dir, **run_args):
        
        if 'blur' in run_args and run_args['blur'] is not None:
            data.blur(run_args['blur'])
        elif 'q0' in run_args:
            blur_nm = (2*np.pi/run_args['q0'])*run_args['blur_size_rel_d0']
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
            
            blur_nm = (2*np.pi/run_args['q0'])*run_args['blur_orientation_image_size_rel']
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


    
class grain_size(grain_size_hex):

    def __init__(self, name='grain_size', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'symmetry': 2,
                        'cmap' : mpl.cm.bone,
                        'blur_size_rel_d0' : 0.25,
                        'blur_orientation_image' : True,
                        'blur_orientation_image_num_passes' : 3,
                        'blur_orientation_image_size_rel' : 0.25,
                        'correlation_edge_exclusion' : 10,
                        'correlation_step_size_points' : 5,
                        'trim_r_curve' : 0.8, # 1.0 doesn't trim anything; 0.8 trims the last 20% of the g(r)-curve
                        'preprocess' : 'default',
                        }
        self.run_args.update(kwargs)
            
            
    def preprocess_SEM(self, data, **run_args):
        data.blur(2.0) # lowpass
        data.enhance(contrast=1.3, contrast_passes=2, resharpen_passes=2)
        data.equalize()
        data.maximize_intensity_spread()
        
        return data
    
            
    @run_default
    def run(self, data, output_dir, **run_args):
        
        run_args['mask'] = self.get_mask(data, output_dir=output_dir, **run_args)
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        
        
        results = {}
        data = copy.deepcopy(data)
        
        #orig_data = data.data.copy()

        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])
        
        
        
        # Pre-process image
        data = self.preprocess(data, **run_args)


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



        # Compute correlation function
        if 'scale' not in run_args:
            run_args['scale'] = data.get_scale()
        new_results = self.correlation_function(angles, output_dir=output_dir, **run_args)
        results.update(new_results)

        
        return results
    


            
class defects_lines(grain_size_hex):
    
    def __init__(self, name='defects_lines', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'symmetry': 2,
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        }
        self.run_args.update(kwargs)
        
        
    def preprocess_default(self, data, **run_args):
        #data.equalize()
        data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        #data.lowkill(run_args['q0']*0.1)
        
        data.blur(2.0) # lowpass
        data.blur(2.0) # lowpass
        #data.blur(1.0) # lowpass
        #data.blur(0.6) # lowpass
        #data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        data.equalize()
        data.maximize_intensity_spread()
        
        return data
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        orig_data = data.data.copy()
        
        results = {}
        data = copy.deepcopy(data)
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('initial', output_dir, ext='.jpg', ir=True)
            data.plot_image(save=outfile, ztrim=[0,0], cmap=run_args['cmap'])

        
        # Pre-process image
        data = self.preprocess(data, **run_args)


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
            
            
            

        # Compute defect density
        h, w = data.data.shape
        im_area = (w*data.x_scale)*(h*data.y_scale)*1e-6 # um^2
        
        # Number of particles expected in perfect case
        L0 = 2*np.pi/run_args['q0']
        num_perfect = h*data.y_scale/(L0)
        
        results['number_defects'] = int(max(num_objects - num_perfect, 0))
        results['defect_density'] = { 'value': results['number_defects']/im_area, 'units': 'um^-2'}
        if run_args['verbosity']>=2:
            print("  {:d} (non-edge) defects / {:.2f} um^2 = ".format(results['number_defects'], im_area))
            print("  {:.4g} defects/um^2".format(results['defect_density']['value']))
        
        
        
            
        if run_args['verbosity']>=5:
            # Colored image
            
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
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
            
            
            
            
        return results




class dots_vs_lines(particles):
    
    def __init__(self, name='dots_vs_lines', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'threshold' : 127,
                        'invert' : False,
                        'diagonal_detection' : False,
                        'cmap' : mpl.cm.bone,
                        'preprocess' : 'default',
                        'method' : 'nearest',
                        }
        self.run_args.update(kwargs)
        
    def preprocess_default(self, data, **run_args):
        #data.equalize()
        data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        for i in range(2):
            data.highpass(run_args['q0']*0.1, run_args['q0']*0.4)
        #data.lowkill(run_args['q0']*0.1)
        
        data.blur(2.0) # lowpass
        data.blur(2.0) # lowpass
        #data.blur(1.0) # lowpass
        #data.blur(0.6) # lowpass
        #data.enhance(contrast=1.3, contrast_passes=0, resharpen_passes=2)
        data.equalize()
        data.maximize_intensity_spread()
        
        return data        
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
        results = {}
        data = copy.deepcopy(data)
        
        results, labeled_array = self._find_objects(data, output_dir, results, **run_args)

        if run_args['verbosity']>=2:
            print("    Identifying regions using '{}' method.".format(run_args['method']))
        method = getattr(self, 'identify_regions_{}'.format(run_args['method']))
        new_results = method(data, output_dir, labeled_array, **run_args)
        results.update(new_results)



        return results


    def identify_regions_nearest(self, data, output_dir, labeled_array, grow_size=3, **run_args):
        
        results = {}
        
        scale_nm = (data.x_scale+data.y_scale)/2 # nm/pixel
        cutoff_pix = run_args['dot_size_cutoff_nm']/scale_nm
        cutoff_area = np.pi*np.square(cutoff_pix) # pixels
        
        if run_args['verbosity']>=5:
            print('        Cutoff {:.2f} nm ({:.1f} pixels); surface area {:.1f} pixels'.format(run_args['dot_size_cutoff_nm'], cutoff_pix, cutoff_area))
        
        particle_areas = np.bincount( labeled_array.flatten() )

        if run_args['verbosity']>=6:
            # Color-coded image of object types
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            
        # Sort labelled objects into the two size categories
        start_time = time.time()
        
        particles_sized = particle_areas[labeled_array] # Image where each pixel has the corresponding size of that object (we obtain this by returning the particle_areas as a lookup table and labeled_array for the object indices)
        
        labeled_array_categories = np.where( (particles_sized<cutoff_area) & (labeled_array>0), 1, 0 ) # Small particles
        labeled_array_categories += np.where( (particles_sized>=cutoff_area) & (labeled_array>0), 2, 0 ) # Big particles
        
        if run_args['verbosity']>=5:
            # Create a color-coded image of the objects
            color_background = [0, 0, 0] # Black
            color_small = [255, 0, 0 ] # Red
            color_big = [0, 255, 0] # Green
            
            red = np.where( (particles_sized<cutoff_area) & (labeled_array>0), color_small[0], color_background[0] )
            red = np.where( (particles_sized>=cutoff_area) & (labeled_array>0), color_big[0], red )
            
            green = np.where( (particles_sized<cutoff_area) & (labeled_array>0), color_small[1], color_background[0] )
            green = np.where( (particles_sized>=cutoff_area) & (labeled_array>0), color_big[1], green )
            
            blue = np.where( (particles_sized<cutoff_area) & (labeled_array>0), color_small[2], color_background[0] )
            blue = np.where( (particles_sized>=cutoff_area) & (labeled_array>0), color_big[2], blue )
            image_data = np.dstack( (red, green, blue) )
            
            image = PIL.Image.fromarray( np.uint8(image_data) )
            outfile = self.get_outfile('coded', output_dir, ext='.png', ir=True)
            image.save(outfile)
            
        
        if run_args['verbosity']>=4:
            print("    categorizing took {:.1f} s".format(time.time()-start_time))

        # Binarize image based on categories
        start_time = time.time()
        
        idx = np.where(labeled_array>0) # Only non-background objects
        values = labeled_array_categories[idx] - 1
        h, w = labeled_array.shape
        XI, YI = np.meshgrid(range(w), range(h))
        points = np.column_stack( (XI[idx], YI[idx]) )
        from scipy.interpolate import griddata
        labeled_array = griddata(points, values, (XI, YI), method='nearest')
        # labeled_array has 0 for regions of small objects; 1 for regions of big objects
        
        if run_args['verbosity']>=1:
            # Create black-and-white image
            im_data = np.uint8(labeled_array*255.0)
            im_data = np.dstack( (im_data, im_data, im_data) )
            im = PIL.Image.fromarray(im_data)
            outfile = self.get_outfile('dots_vs_lines', output_dir, ext='.png', ir=False)
            im.save(outfile)        
        
        if run_args['verbosity']>=3:
            print("    {} method took {:.1f} s".format(run_args['method'], time.time()-start_time))
        
        
        region_2_pixels = np.sum(labeled_array)
        region_1_pixels = labeled_array.size - region_2_pixels
        print( '    Black region: {:,d} pixels ({:.1f} %)'.format(region_1_pixels, region_1_pixels*100.0/labeled_array.size) )
        print( '    White region: {:,d} pixels ({:.1f} %)'.format(region_2_pixels, region_2_pixels*100.0/labeled_array.size) )
        
        results['dot_fractional_area'] = region_1_pixels*1.0/(labeled_array.size*1.0)
        results['line_fractional_area'] = region_2_pixels*1.0/(labeled_array.size*1.0)
        
        return results


    def identify_regions_grow(self, data, output_dir, labeled_array, grow_size=3, **run_args):
        
        
        results = {}
        
        scale_nm = (data.x_scale+data.y_scale)/2 # nm/pixel
        cutoff_pix = run_args['dot_size_cutoff_nm']/scale_nm
        cutoff_area = np.pi*np.square(cutoff_pix) # pixels
        
        if run_args['verbosity']>=5:
            print('        Cutoff {:.2f} nm ({:.1f} pixels); surface area {:.1f} pixels'.format(run_args['dot_size_cutoff_nm'], cutoff_pix, cutoff_area))
        
        particle_areas = np.bincount( labeled_array.flatten() )

        if run_args['verbosity']>=6:
            # Color-coded image of object types
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
        # Sort labelled objects into the two size categories
        start_time = time.time()
        for ix in range( len(labeled_array[0]) ):
            for iy in range( len(labeled_array) ):
                object_index = labeled_array[iy,ix]

                if object_index==0:
                    labeled_array[iy,ix] = 0
                elif particle_areas[object_index]<cutoff_area:
                    labeled_array[iy,ix] = 1
                else:
                    labeled_array[iy,ix] = 2

                if run_args['verbosity']>=6:
                    if object_index==0:
                        c = (0, 0, 0) # Black (no particles)
                    elif particle_areas[object_index]<cutoff_area:
                        c = (255, 0, 0) # Red (small particles)
                    else:
                        c = (0, 255, 0 ) # Green (big particles)
                    pix[ix,iy] = c
                    
        if run_args['verbosity']>=6:
            outfile = self.get_outfile('coded', output_dir, ext='.png', ir=True)
            im.save(outfile)

        if run_args['verbosity']>=4:
            print("    categorizing took {:.1f} s".format(time.time()-start_time))



        # Grow particle regions
        start_time = time.time()
        found_zeros = True
        igrow = 1
        h, w = data.data.shape
        while found_zeros:
            
            num_zero = labeled_array.size - np.count_nonzero(labeled_array)
            print( '        Growing regions: pass {:d} ({:,d} zeros left); {:.2f}% done...'.format(igrow, num_zero, 100.0*(1-num_zero/labeled_array.size) ) )
            igrow += 1
            found_zeros = False
            if num_zero>0:
                labeled_array_cur = np.copy(labeled_array)
                
                for ix in range(w-(grow_size-1)):
                    for iy in range(h-(grow_size-1)):
                        patch = labeled_array[iy:iy+grow_size,ix:ix+grow_size]
                        if 0 in patch:
                            found_zeros = True
                            if 1 in patch:
                                if 2 in patch:
                                    # 1's and 2's
                                    labeled_array_cur[iy:iy+grow_size,ix:ix+grow_size] = 1 # 1's are arbitrarily given priority...
                                else:
                                    # Just 1's and 0's
                                    labeled_array_cur[iy:iy+grow_size,ix:ix+grow_size] = 1
                            elif 2 in patch:
                                # Just 2's and 0's
                                labeled_array_cur[iy:iy+grow_size,ix:ix+grow_size] = 2
                            else:
                                pass # All 0's, do nothing
                                
                labeled_array = labeled_array_cur
                        
        if run_args['verbosity']>=3:
            print("    {} method took {:.1f} s".format(run_args['method'], time.time()-start_time))
                        
                        
        region_2_pixels = np.sum(labeled_array) - labeled_array.size
        region_1_pixels = labeled_array.size - region_2_pixels
        print( '    Black region: {:,d} pixels ({:.1f} %)'.format(region_1_pixels, region_1_pixels*100.0/labeled_array.size) )
        print( '    White region: {:,d} pixels ({:.1f} %)'.format(region_2_pixels, region_2_pixels*100.0/labeled_array.size) )
        
        results['dot_fractional_area'] = region_1_pixels*1.0/(labeled_array.size*1.0)
        results['line_fractional_area'] = region_2_pixels*1.0/(labeled_array.size*1.0)
        
        if run_args['verbosity']>=1:
            # Create black-and-white image
            im = PIL.Image.fromarray( np.uint8(data.data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            for ix in range(w):
                for iy in range(h):
                    object_index = labeled_array[iy,ix]
                    
                    if labeled_array[iy,ix]==0:
                        c = (255, 0, 0) # Red (error)
                    elif labeled_array[iy,ix]==1:
                        c = (0, 0, 0) # Black (small particles)
                    else:
                        c = (255, 255, 255 ) # White (big particles)
                    pix[ix,iy] = c
                    
            outfile = self.get_outfile('dots_vs_lines', output_dir, ext='.png', ir=False)
            im.save(outfile)
                        
        return results
    


color_list1 = [ (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1),(0,1,1),(1,1,1),]
color_list2 = [ (0.7*c[0], 0.7*c[1], 0.7*c[2]) for c in color_list1 ]
color_list3 = [ (0.5*c[0], 1.0*c[1], 1.0*c[2]) for c in color_list1 ]
color_list4 = [ (1.0*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list1 ]
color_list5 = [ (1.0*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list1 ]
color_list6 = [ (1.0*c[0], 0.7*c[1], 0.5*c[2]) for c in color_list1 ]
color_list7 = [ (1.0*c[0], 0.5*c[1], 0.7*c[2]) for c in color_list1 ]
color_list8 = [ (0.7*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list1 ]
color_list9 = [ (0.5*c[0], 1.0*c[1], 0.7*c[2]) for c in color_list1 ]
color_list10 = [ (0.7*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list1 ]
color_list11 = [ (0.5*c[0], 0.7*c[1], 1.0*c[2]) for c in color_list1 ]
color_list = color_list1 + color_list2 + color_list3 + color_list4 + color_list5 + color_list6 + color_list7 + color_list8 + color_list9 + color_list10 + color_list11
