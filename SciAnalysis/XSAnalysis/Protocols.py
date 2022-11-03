#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.Protocols` - Data analysis protocols
================================================
.. module:: SciAnalysis.XSAnalysis.Protocols
   :synopsis: Convenient protocols for data analysis.
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  Data analysis protocols.
################################################################################
# Known Bugs:
#  Using both 'flag_swaxs' and 'background' at the same time likely won't work.
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################

from SciAnalysis.XSAnalysis.Data import * #from .Data import *
from SciAnalysis.tools import * #from ..tools import *
#from SciAnalysis.IO_HDF import *

import copy, glob

class ProcessorXS(Processor):

    
    def load(self, infile, **kwargs):

        #calibration = kwargs['calibration'] if 'calibration' in kwargs else None
        #mask = kwargs['mask'] if 'mask' in kwargs else None
        #data = Data2DScattering(infile, calibration=calibration, mask=mask)
        
        #data = Data2DScattering(infile, **kwargs)
        data = self.handle_calibration(infile, **kwargs)
        data.infile = infile
        
        data.threshold_pixels(4294967295-1) # Eiger inter-module gaps
        
        if 'background' in kwargs:
            data.name = data.name+'_rmbkg'
            self.handle_background(data, **kwargs)
        
        if 'dezing' in kwargs and kwargs['dezing']:
            data.dezinger()
        
        if 'flip' in kwargs and kwargs['flip']:
            #if flip: self.im = self.im.transpose(Image.ROTATE_90).transpose(Image.FLIP_LEFT_RIGHT)
            data.data = np.rot90(data.data) # rotate CCW
            data.data = np.fliplr(data.data) # Flip left/right

        if 'rotCCW' in kwargs and kwargs['rotCCW']:
            data.data = np.rot90(data.data) # rotate CCW

        if 'rotCW' in kwargs and kwargs['rotCW']:
            data.data = np.rot90(data.data, k=3) # rotate CW

        if 'rot180' in kwargs and kwargs['rot180']:
            data.data = np.flipud(data.data) # Flip up/down
            data.data = np.fliplr(data.data) # Flip left/right


        if data.mask is not None:
            data.data *= data.mask.data


        return data
        
    
    def handle_background(self, data, scanid=None, **kwargs):
        
        verbosity = kwargs['verbosity'] if 'verbosity' in kwargs else 3
        
        if isinstance(kwargs['background'], (int, float)):
            # Constant background to be subtracted from whole image
            data.data -= kwargs['background']
            
        elif isinstance(kwargs['background'], str):
            # Subtract whole image as background
            infiles_background = glob.glob(kwargs['background'])
            if verbosity>=5:
                print('# {} Background Files: {}'.format(len(infiles_background), infiles_background))
            average_background_data = np.zeros(data.data.shape)
            for ii, infile_background in enumerate(infiles_background):
                data_background = Data2DScattering(infile_background, **kwargs)
                average_background_data += data_background.data
            average_background_data /= len(infiles_background)

            if isinstance(kwargs['transmission_int'], (str)):
                # Read from file
                import pandas as pd
                df = pd.read_csv(kwargs['transmission_int'])
                if verbosity>=6:
                    print(df)
                
                # Find the best matched name from CSV
                if scanid ==None:
                    emptyname, df0 = Filename(infile_background).get_best_match(df)
                else:
                    emptyname, df0 = Filename(infile_background).get_best_match(scanID)

                samplename, df1 = Filename(data.infile).get_best_match(df)
                print("# Found {} for background, {} for sample".format(df0['a_filename'].values, df1['a_filename'].values))
                print("# i.e. {} for background, {} for sample".format(df0['b_scanID'].values, df1['b_scanID'].values))
                # User the latest if more than one
                factor = df1.c_I0.to_numpy()[-1] / df0.c_I0.to_numpy()[-1]

            elif isinstance(kwargs['transmission_int'], (int, float)):
                # Specify value
                factor = kwargs['transmission_int']		  
            else:
                print('# WARNING: transmission_int invalid, use factor=1')
                factor = 1.0             

            if verbosity>=3: print("# factor = {:.3f}".format(factor))
            average_background_data[average_background_data>=0] *= factor
            if verbosity>=5:
                print("# Before: data MAX {:.3f}, MEAN {:.3f}".format(np.max(data.data), np.mean(data.data)))

            data.data -= average_background_data
            
            if verbosity>=5:
                print("# After: data MAX {:.3f}, MEAN {:.3f}".format(np.max(data.data),np.mean(data.data)))

            
        elif isinstance(kwargs['background'], (list, np.ndarray)):
            # Subtract supplied array as background
            data.data -= kwargs['background']
            
        else:
            print("ProcessorXS.load: Specified background type not recognized.")
        
        
    def handle_calibration(self, infile, **kwargs):
        # This is currently an ad-hoc definition tuned to a particular set of
        # assumptions about kwargs names.
        # TODO: This should probably be replaced with a general way to handle
        # multiple detectors (each with their own calibration). For instance, the
        # calibration object can internally have references to multipe calibration
        # objects for detectors. Then datasets can control which calibration they
        # are using.
        
        verbosity = kwargs['verbosity'] if 'verbosity' in kwargs else 3
        
        if 'flag_swaxs' in kwargs and kwargs['flag_swaxs']:
            if 'calibration2' in kwargs and 'mask2' in kwargs:
                if verbosity>=5:
                    print('    ProcessorXS.load using calibration2 and mask2')
                data = Data2DScattering(infile, calibration=kwargs['calibration2'], mask=kwargs['mask2'])
                # WARNING: This init doesn't pass other kwargs. This is probably fine, but means
                # if future kwargs behavior is added in the future, they won't be recongized here.
            else:
                print('ERROR: calibration2 and/or mask2 not provided.')
        else:
            if verbosity>=5:
                print('    ProcessorXS.load using default calibration')
            data = Data2DScattering(infile, **kwargs) # Default load
            
        return data


    def connect_databroker(self, beamline):
        '''Establish a connection to a Bluesky databroker database.
        This allows protocols to access databroker to get metadata
        that might be useful to their analysis.'''

        #from databroker import list_configs
        #print( list_configs() ) # List of available beamlines

        from databroker import Broker
        self.db = Broker.named(beamline)


    def get_db(self, recent_days=None, verbosity=3, **constraints):
        
        if recent_days is not None:
            constraints['since'] = time.time() - recent_days*24*60*60
        
        headers = self.db(**constraints)
        headers = [h for h in headers]
        
        if verbosity>=1 and len(headers)!=1:
            print('  WARNING: get_db got {} matches for constraints: {}'.format(len(headers), constraints))
        
        
        return headers[0]
            
        
class HDF5(Protocol):

    def __init__(self, name='HDF5', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.h5'
        self.run_args = {
            }
        self.run_args.update(kwargs)

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'hdf5' in run_args['save_results']:
            self.save_Data2D_HDF5(data, 'raw detector image', output_dir, results=results)
                                  
                                  
        return results


        

class thumbnails(Protocol):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
            'crop' : None,
            'shift_crop_up' : 0.0,
            'make_square' : False,
            'blur' : 2.0,
            'resize' : 0.2,
            'ztrim' : [0.05, 0.005],
            'preserve_data' : True,
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['preserve_data']:
            # Avoid changing the data (which would disrupt downstream analysis of this data object)
            data = copy.deepcopy(data)
            # TODO: This can raise errors (mpl or qt nodes may not copy cleanly).
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'], shift_crop_up=run_args['shift_crop_up'], make_square=run_args['make_square'])
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize']) # Shrink
        
        data.set_z_display([None, None, 'gamma', 0.3])
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)
        
        #print(data.stats())
        
        if 'plots' in run_args['save_results']:
            results['files_saved'] = [
                { 'filename': '{}'.format(outfile) ,
                'description' : 'quick view (thumbnail) image' ,
                'type' : 'plot' # 'data', 'plot'
                } ,
                ]
            data.plot_image(outfile, **run_args)
        
        return results
        
        
class circular_average(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'bins_relative' : 1.0,
            'markersize' : 0,
            'linewidth' : 1.5,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)
        
        
        line = data.circular_average_q_bin(error=True, bins_relative=run_args['bins_relative'])
        #line.smooth(2.0, bins=10)
        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])

        if 'twotheta' in run_args and run_args['twotheta']:
            line.x = data.calibration.q_to_angle(line.x)
            line.x_label = '2theta (deg)'
            line.x_rlabel = '$2 \theta \, (^{\circ})$'
        
        # TODO: Add options for a second x-axis scaled as twotheta or d=2*pi/q

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)

        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, **run_args)

        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)

       
        return results



class circular_average_sum(circular_average):

    def __init__(self, name='circular_average_sum', **kwargs):

        self.name = self.__class__.__name__ if name is None else name

        self.default_ext = '.png'
        self.run_args = {
            'bins_relative' : 1.0,
            'markersize' : 0,
            'linewidth' : 1.5,
            }
        self.run_args.update(kwargs)


    @run_default
    def run(self, data, output_dir, **run_args):

        results = {}

        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)


        line = data.circular_average_q_bin(error=True, bins_relative=run_args['bins_relative'])
        #line.smooth(2.0, bins=10)

        line_sub = line.sub_range(run_args['sum_range'][0], run_args['sum_range'][1])
        results['values_sum'] = np.sum(line_sub.y)

        


        # Plot and save data
        class DataLines_current(DataLines):

            def _plot_extra(self, **plot_args):

                xi, xf, yi, yf = self.ax.axis()
                v_spacing = (yf-yi)*0.10

                yp = yf
                s = '$S = \, {:.1f} $'.format(self.results['values_sum'])
                self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')
                
                self.ax.axvspan(run_args['sum_range'][0], run_args['sum_range'][1], color='b', alpha=0.1)


        lines = DataLines_current([line])
        lines.copy_labels(line)
        lines.results = results

        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            lines.plot(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)

        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)


        return results
                
                
                
class circular_average_q2I(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'bins_relative' : 1.0,
            'markersize' : 0,
            'linewidth' : 1.5,
            'qn_power' : 2.0,
            'num_curves' : 1, # For (optional) fitting
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.circular_average_q_bin(error=True, bins_relative=run_args['bins_relative'])
        
        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
        
        
        if run_args['qn_power']==2.0:
            line.y *= np.square(line.x)
            line.y_label = 'q^2*I(q)'
            line.y_rlabel = '$q^2 I(q) \, (\AA^{-2} \mathrm{counts/pixel})$'
            
        else:
            line.y *= np.power(line.x, run_args['qn_power'])
            line.y_label = 'q^n*I(q)'
            line.y_rlabel = '$q^n I(q) \, (\AA^{-n} \mathrm{counts/pixel})$'
            
        results['qn_power'] = run_args['qn_power']
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
            line.plot(save=outfile, **run_args)
        
        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='_q2I.dat')
            line.save_data(outfile)

        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)
        
        
        return results
                       

    def output_exists(self, name, output_dir):

        if 'file_extension' in self.run_args:
            ext = '_q2I{}'.format(self.run_args['file_extension'])
        else:
            ext = '_q2I{}'.format(self.default_ext)

        outfile = self.get_outfile(name, output_dir, ext=ext)
        return os.path.isfile(outfile)


class fit_peaks(Protocol):
    
    def _fit(self, line, results, **run_args):
        
        # Fit
        lm_result, fit_line, fit_line_extended, fit_line_curves = self._fit_peaks(line, **run_args)
        
        fit_name = 'fit_peaks'
        prefactor_total = 0
        for param_name, param in lm_result.params.items():
            results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            if 'prefactor' in param_name:
                prefactor_total += np.abs(param.value)
            
        results['{}_prefactor_total'.format(fit_name)] = prefactor_total
        results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
        
        # Calculate some additional things
        for i in range(run_args['num_curves']):
            q = results['{}_x_center{}'.format(fit_name, i+1)]['value']
            d = 0.1*2.*np.pi/q
            err = results['{}_x_center{}'.format(fit_name, i+1)]['error']
            if err is None:
                err = 0
            d_err = err*(d/q)
            #results['{}_d0{}'.format(fit_name, i+1)] = d
            results['{}_d0{}'.format(fit_name, i+1)] = { 'value': d, 'error': d_err }
            
            sigma = results['{}_sigma{}'.format(fit_name, i+1)]['value']
            if 'instrumental_resolution' in run_args:
                sigma = np.sqrt( np.square(sigma) - np.square(run_args['instrumental_resolution']) )
            xi = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/sigma
            err = results['{}_sigma{}'.format(fit_name, i+1)]['error']
            if err is None:
                err = 0
            xi_err = err*(xi/sigma)            
            #results['{}_grain_size{}'.format(fit_name, i+1)] = xi
            results['{}_grain_size{}'.format(fit_name, i+1)] = { 'value': xi, 'error': xi_err }
            
        results['{}_d0'.format(fit_name)] = results['{}_d01'.format(fit_name)]
        results['{}_grain_size'.format(fit_name)] = results['{}_grain_size1'.format(fit_name)]
        
        
        # Plot and save data
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                
                if 'fit_range' in self._run_args:
                    xstart, xend = self._run_args['fit_range']
                    line = self.lines[0].sub_range(xstart, xend)
                else:
                    line = self.lines[0]
                
                yf = np.max(line.y)*1.5
                self.ax.axis([xi, xf, yi, yf])

                color = 'b'
                font_size = self._run_args['font_size'] if 'font_size' in self._run_args else 18
                v_spacing = (yf-yi)*0.065*(font_size/20)

                s = '$\chi^2 = \, {:.4g}$'.format(self.results['fit_peaks_chi_squared'])
                self.ax.text(xi, yi, s, size=font_size, color=color, verticalalignment='bottom', horizontalalignment='left')


                
                for i in range(self._run_args['num_curves']):
                    
                    self.ax.axvline(self.results['fit_peaks_x_center{}'.format(i+1)]['value'], linewidth=1, color=color, alpha=0.5)
                    
                    if i<=1:
                        yp = yf
                    else:
                        yp -= v_spacing*1.5
                    if i==0:
                        ha, xp = 'right', xf
                    else:
                        ha, xp = 'left', xi
        
                    s = '$p_{{ {:d} }} = \, {:.3g}$'.format(i+1, self.results['fit_peaks_prefactor{}'.format(i+1)]['value'])
                    self.ax.text(xp, yp, s, size=font_size, color=color, verticalalignment='top', horizontalalignment=ha)

                    yp -= v_spacing
                    s = '$q = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_x_center{}'.format(i+1)]['value'])
                    self.ax.text(xp, yp, s, size=font_size, color=color, verticalalignment='top', horizontalalignment=ha)

                    yp -= v_spacing
                    s = r'$d \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_d0{}'.format(i+1)]['value'])
                    self.ax.text(xp, yp, s, size=font_size, color=color, verticalalignment='top', horizontalalignment=ha)

                    yp -= v_spacing
                    s = '$\sigma = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_sigma{}'.format(i+1)]['value'])
                    self.ax.text(xp, yp, s, size=font_size, color=color, verticalalignment='top', horizontalalignment=ha)
                    
                    yp -= v_spacing
                    s = r'$\xi \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_grain_size{}'.format(i+1)]['value'])
                    self.ax.text(xp, yp, s, size=font_size, color=color, verticalalignment='top', horizontalalignment=ha)        
                    
        
        
        lines = DataLines_current([line, fit_line, fit_line_extended])
        if 'num_curves' in run_args and run_args['num_curves']>1 and 'show_curves' in run_args and run_args['show_curves']:
            for curve in fit_line_curves:
                lines.add_line(curve)
        lines.results = results
        lines._run_args = run_args
        lines.copy_labels(line)
        
        # Note that the results dictionary is modified within this function.
        # Thus although it is not returned, it is part of the set of returned
        # information.
        return lines
    
    
    def _fit_peaks(self, line, q0=None, num_curves=1, **run_args):
        # Usage: lm_result, fit_line, fit_line_extended = self._fit_peaks(line, **run_args)

        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            # Linear background
            m = v['m']*x + v['b']
            # Power-law background
            m += v['qp']*np.power( np.abs(x), v['qalpha'] )
            
            # Gaussian peaks
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()


        m = (line.y[-1]-line.y[0])/(line.x[-1]-line.x[0])
        b = line.y[0] - m*line.x[0]

        xs = np.abs(line.x)
        ys = line.y
        qalpha = (np.log(ys[0])-np.log(ys[-1]))/(np.log(xs[0])-np.log(xs[-1]))
        qp = np.exp( np.log(ys[0]) - qalpha*np.log(xs[0]) )

        if True:
            # Linear background
            params.add('m', value=m, min=abs(m)*-10, max=abs(m)*+10+1e-12, vary=False)
            params.add('b', value=b, vary=False)
            
            params.add('qp', value=0, vary=False)
            params.add('qalpha', value=1.0, vary=False)
            
        else:
            # Power-law background
            params.add('m', value=0, vary=False)
            params.add('b', value=0, vary=False)
            
            params.add('qp', value=qp, vary=False)
            params.add('qalpha', value=qalpha, vary=False)
            
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y_max()



        # Best guess for peak position
        xs = np.asarray(line.x)
        ys = np.asarray(line.y)
        
        if isinstance(q0, (list, tuple, np.ndarray)):
            # q0 may be a list of q0 positions for a set of peaks
            q0s = q0
            q0 = q0s[0]
        else:
            q0s = None
        
        if q0 is not None:
            # Sort
            indices = np.argsort(xs)
            x_sorted = xs[indices]
            y_sorted = ys[indices]
            
            idx = np.where( x_sorted>=q0 )[0][0]

            xpeak = x_sorted[idx]
            ypeak = y_sorted[idx]
        
        else:
            # Sort
            indices = np.argsort(ys)
            x_sorted = xs[indices]
            y_sorted = ys[indices]
            
            target = np.max(ys)
            idx = np.where( y_sorted>=target )[0][0]
            
            xpeak = x_sorted[idx]
            ypeak = y_sorted[idx]
            
        xpeak, ypeak = line.target_x(xpeak)
                                 

        prefactor = ypeak - ( m*xpeak + b )
        if 'sigma' in run_args:
            sigma = run_args['sigma']
        else:
            sigma = 0.1*xspan
        
        for i in range(num_curves):
            
            params.add('prefactor{:d}'.format(i+1), value=prefactor, min=0, max=max(np.max(line.y)*1.5,0)+1e-12, vary=False)
            if i==0:
                # 1st peak should be at max location
                params.add('x_center{:d}'.format(i+1), value=xpeak, min=np.min(line.x), max=np.max(line.x), vary=False)
            else:
                # Additional peaks can be spread out
                # (or use q0s value if available)
                if q0s is not None and len(q0s)>i:
                    xpos = q0s[i]
                else:
                    xpos = np.min(line.x) + (xspan/num_curves)*i
                params.add('x_center{:d}'.format(i+1), value=xpos, min=np.min(line.x), max=np.max(line.x), vary=False)
                
            params.add('sigma{:d}'.format(i+1), value=sigma, min=0.00001, max=xspan*0.5, vary=False)
        
        
        # Fit only the peak width
        params['sigma1'].vary = True
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if True:
            # Tweak peak position
            lm_result.params['sigma1'].vary = False
            lm_result.params['x_center1'].vary = True
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
        
        if True:
            # Relax entire fit
            lm_result.params['m'].vary = True
            lm_result.params['b'].vary = True
            #lm_result.params['qp'].vary = True
            #lm_result.params['qalpha'].vary = True
            
            for i in range(num_curves):
                lm_result.params['prefactor{:d}'.format(i+1)].vary = True
                lm_result.params['sigma{:d}'.format(i+1)].vary = True
                lm_result.params['x_center{:d}'.format(i+1)].vary = True
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
            
            #lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y), method='nelder')
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        x_span = abs(np.max(line.x)-np.min(line.x))
        fit_x = np.linspace(np.min(line.x)-x_span, np.max(line.x)+x_span, num=2000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'alpha':0.5, 'marker':None, 'linewidth':2.0})
        
        # Generate component curves
        prefactors = [lm_result.params['prefactor{:d}'.format(i+1)].value for i in range(num_curves) ]
        fit_line_curves = []
        for i in range(num_curves):
            # Set all but one prefactors to zero
            for j, prefactor in enumerate(prefactors):
                if j!=i:
                    prefactor = 0
                lm_result.params['prefactor{:d}'.format(j+1)].value = prefactor

            fit_y = model(lm_result.params.valuesdict(), fit_x)
            fit_line_curve = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'purple', 'alpha':0.5, 'marker':None, 'linewidth':1.0})
            fit_line_curves.append(fit_line_curve)
        
        # Return the model to the correct state
        for i, prefactor in enumerate(prefactors):
            lm_result.params['prefactor{:d}'.format(i+1)].value = prefactor
            

        return lm_result, fit_line, fit_line_extended, fit_line_curves


class fit_FormFactor_Sphere(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'load_from' : 'circular_average',
            'markersize' : 0,
            'linewidth' : 1.5,
            'ylog' : True,
            }
        self.run_args.update(kwargs)


    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['load_from'] is not None:
            # Load data from corresponding circular_average
            load_dir = Path(output_dir).parent.joinpath(run_args['load_from'])
            infile = load_dir.joinpath(data.name + '.dat')
            data_1d = np.loadtxt(infile) # q qerr I(q) I(q)err
            
            line = DataLine(x=data_1d[:,0], y=data_1d[:,2], x_err=data_1d[:,1], y_err=data_1d[:,3], x_label='q', y_label='P(q)', x_rlabel='$q \, (\mathrm{\AA})$', y_rlabel='P(q)')
            
        if False:
            # Testing using CSV file
            load_dir = Path(output_dir).parent.joinpath('csv')
            infile = load_dir.joinpath('SP38nm_Free' + '.csv')
            data_1d = np.loadtxt(infile, delimiter=',', skiprows=1) # q P(q)
            line = DataLine(x=data_1d[:,1], y=data_1d[:,2], x_label='q', y_label='P(q)', x_rlabel='$q \, (\mathrm{\AA})$', y_rlabel='$P(q)$')


        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
            
            
        # Do the fit
        lines = self._fit(line, results, **run_args)


        # Save results
        if 'save_fit' in run_args and run_args['save_fit']:
            outfile = self.get_outfile(data.name, output_dir, ext='_fit.dat')
            lines.lines[1].save_data(outfile)            
        
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name, output_dir, ext=self.default_ext)
            lines.plot(save=outfile, **run_args)
        
        
        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)
        
        return results

    
    def _fit(self, line, results, **run_args):
        
        # Do the fit
        lm_result, fit_line, fit_line_extended = self._fit_FF(line, **run_args)


        # Store results
        fit_name = 'fit_FF'
        prefactor_total = 0
        for param_name, param in lm_result.params.items():
            results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            
        results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree



        # Plot and save data
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                
                if 'fit_range' in self._run_args:
                    xstart, xend = self._run_args['fit_range']
                    line = self.lines[0].sub_range(xstart, xend)
                else:
                    line = self.lines[0]
                
                yi = np.min(line.y)*0.5
                yf = np.max(line.y)*2.0
                self.ax.axis([xi, xf, yi, yf])

                font_size = self._run_args['font_size'] if 'font_size' in self._run_args else 18
                v_spacing = (np.log(yf)-np.log(yi))*(font_size/20)*0.2

                s = '$\chi^2 = \, {:.4g}$'.format(self.results['fit_FF_chi_squared'])
                self.ax.text(xi, yi, s, size=font_size, color='b', verticalalignment='bottom', horizontalalignment='left')
                
                xp, yp = xf, yf
                for k, v in self._run_args['initial_guess'].items():
                    r = self.results['fit_FF_{}'.format(k)]
                    if r['error'] is None:
                        e = 'NA'
                    else:
                        e = '{:.2g}'.format(r['error'])
                    s = '{} = {:.3g} ± {}'.format(k, r['value'], e)
                    self.ax.text(xp, yp, s, size=font_size, color='b', verticalalignment='top', horizontalalignment='right')
                    
                    yp /= v_spacing
                    


                    
        
        
        lines = DataLines_current([line, fit_line, fit_line_extended])
        lines.results = results
        lines._run_args = run_args
        lines.copy_labels(line)
        
        # Note that the results dictionary is modified within this function.
        # Thus although it is not returned, it is part of the set of returned
        # information.
        return lines        
        
        
    def _fit_FF(self, line, **run_args):
        
        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            
            background = v['qpower_scale']*np.power(x, v['qpower']) + v['background']
            
            #P = self.sphere_form_factor_intensity(q=x, R=v['radius'], delta_rho=v['delta_rho'] )
            P = self.poly_sphere_form_factor_intensity( q=x, R=v['radius'], sigma=v['sigma'], delta_rho=v['delta_rho'] )
            
            return P + background
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            #return m - data
            return np.log(data/m)
        
        params = lmfit.Parameters()
        
        for k, v in run_args['initial_guess'].items():
            v_min, v_max = run_args['limits'][k]
            params.add(k, value=v, min=v_min, max=v_max, vary=True)

        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))

        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        x_span = abs(np.max(line.x)-np.min(line.x))
        fit_x = np.linspace(0, np.max(line.x)+x_span*0.2, num=2000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'alpha':0.5, 'marker':None, 'linewidth':2.0})

        
        return lm_result, fit_line, fit_line_extended


    def sphere_form_factor(self, q, R, delta_rho=1):
        """Calculate the form factor (fq) of a sphere."""
        from scipy.special import spherical_jn
        volume = (4.0/3.0)*np.pi*(R**3)
        qR = R*q
        F =  delta_rho  * volume * 3 * spherical_jn(1, qR) / qR      
        return F

    def sphere_form_factor_intensity(self, q, R, delta_rho=1):
        """Calculate the form factor intensity (Pq=fq**2) of a sphere."""  
        F =  self.sphere_form_factor(q=q, R=R, delta_rho=delta_rho)   
        return np.abs(F)**2 

    def poly_sphere_form_factor_intensity(self, q, R, sigma, delta_rho=1):
        """Calculate the form factor intensity of a polydispersed distribution
        of spheres.
        Pq = sum( wi * fqi**2 )"""    
        disx, disw, disdx = self.distribution_gaussian( R, sigma * R )
        f1_ = np.zeros_like(q)
        cts=0
        for i, ri in enumerate(disx):
            f1_ += self.sphere_form_factor_intensity(q, R=ri, delta_rho= delta_rho ) * disw[i] * disdx
            cts += disw[i] *disdx  
        return f1_/cts
    
    def distribution_gaussian(self, radius=1.0, sigma=0.01, num_points=30, spread=3, only_positive=True): 
        ''' Create a gaussian distribution'''    
        #distribution_list = []  
        x, step= np.linspace( radius - spread* sigma, radius + spread*sigma, num_points,retstep=True)
        #print(x,step)
        prefactor = 1/( sigma*np.sqrt(2*np.pi) ) 
        delta = radius - x
        wt = prefactor*np.exp( - (delta**2)/(2*( sigma**2 ) ) )         
        if only_positive:
            wP = x>0
            return x[wP], wt[wP], step
        else:
            return x, wt, step      
        

class circular_average_q2I_fit(circular_average_q2I, fit_peaks):
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.circular_average_q_bin(error=True, bins_relative=run_args['bins_relative'])
        
        if run_args['qn_power']==2.0:
            line.y *= np.square(line.x)
            line.y_label = 'q^2*I(q)'
            line.y_rlabel = '$q^2 I(q) \, (\AA^{-2} \mathrm{counts/pixel})$'
            
        else:
            line.y *= np.power(line.x, run_args['qn_power'])
            line.y_label = 'q^n*I(q)'
            line.y_rlabel = '$q^n I(q) \, (\AA^{-n} \mathrm{counts/pixel})$'
            
        results['qn_power'] = run_args['qn_power']
            

        if 'txt' in run_args['save_results']:        
            outfile = self.get_outfile(data.name, output_dir, ext='_q2I.dat')
            line.save_data(outfile)        

        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
        
        lines = self._fit(line, results, **run_args)
        #lines = DataLines([line])
        
        if 'save_fit' in run_args and run_args['save_fit']:
            # lines.lines contains: line, fit_line, fit_line_extended, fit_curve1, fit_curve2, ...
            if 'show_curves' in run_args and run_args['show_curves']:
                for i, line in enumerate(lines.lines[2:]):
                    outfile = self.get_outfile(data.name, output_dir, ext='_q2I_fit{}.dat'.format(i))
                    line.save_data(outfile)            
            else:
                outfile = self.get_outfile(data.name, output_dir, ext='_q2I_fit.dat')
                lines.lines[1].save_data(outfile)            
        
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
            lines.plot(save=outfile, **run_args)
        
        
        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)
        
        return results
         
        #End class circular_average_q2I_fit(Protocol, fit_peaks)
        ########################################
         
         
class sector_average(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'bins_relative' : 1.0,
            'markersize' : 0,
            'linewidth' : 1.5,
            'error' : True, 
            'show_region' : False,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)
        
        
        line = data.sector_average_q_bin(**run_args)
        #line.smooth(2.0, bins=10)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
        
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)

        if 'hdf5' in run_args['save_results']:          
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results) 

        return results
    

class sector_average_fit(sector_average, fit_peaks):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'bins_relative' : 1.0,
            'markersize' : 0,
            'linewidth' : 1.5,
            'error' : False, 
            #'elinewidth' : 0.25,
            #'ecolor': '0.5',
            #'capsize': 1.5,
            'show_region' : False,
            'qn_power' : None,
            'num_curves' : 1,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)
        
        
        line = data.sector_average_q_bin(**run_args)
        #line.smooth(2.0, bins=10)
        
        if run_args['qn_power'] is not None:
            if run_args['qn_power']==2.0:
                line.y *= np.square(line.x)
                line.y_label = 'q^2*I(q)'
                line.y_rlabel = '$q^2 I(q) \, (\AA^{-2} \mathrm{counts/pixel})$'
            else:
                line.y *= np.power(line.x, run_args['qn_power'])
                line.y_label = 'q^n*I(q)'
                line.y_rlabel = '$q^n I(q) \, (\AA^{-n} \mathrm{counts/pixel})$'
        
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
        
        
        #if 'plots' in run_args['save_results']:
            #self.label_filename(data, line, **run_args)
            #outfile = self.get_outfile(data.name, output_dir)
            #line.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)

        if 'hdf5' in run_args['save_results']:          
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results) 


        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])


        # Do fit
        lines = self._fit(line, results, **run_args)
        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            lines.plot(save=outfile, **run_args)        
        

        return results
    

                
class roi(Protocol):

    def __init__(self, name='roi', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.txt'
        self.run_args = {
            'extra' : '',
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        results.update( data.roi_q(**run_args) )
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        if run_args['verbosity']>=3:
            print('ROI stats:')
            print(results)
            
        outfile = self.get_outfile(data.name, output_dir)
        with open(outfile, 'w') as fout:
            for k, v in results.items():
                fout.write('{} : {}\n'.format(k, v))
        
        return results                
                
                
class linecut_angle(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'plot_range' : [-180, 180, 0, None],
            'markersize' : 0,
            'linewidth' : 1.5,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_angle(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile, ztrim=run_args['ztrim'])
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        if 'smooth' in run_args:
            line.smooth(run_args['smooth'])
        
        if 'polarization_correction' in run_args and run_args['polarization_correction']:
            chi_deg = line.x
            two_theta_rad = 2.0*np.arcsin(run_args['q0']/(2.*data.calibration.get_k()))
            P_h = 1 - np.square(np.sin(two_theta_rad))*np.square(np.sin(np.radians(chi_deg)))
            line.y *= P_h


        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, **run_args)
            
            #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
            #line.plot_polar(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        return results
                                


class linecut_angle_graininess(linecut_angle):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'ztrim' : [0.1, 0.01],
            'plot_range' : [-180, 180, 0, None],
            'markersize' : 3,
            'linewidth' : 1.5,
            'smooth' : 0.25, # Data-smoothing, as a fraction of sigma_chi
            'multiplicity' : 6,
            'scattering_volume' : None,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        # Graininess analysis is based on:
        # https://doi.org/10.1107/S1600576714020822
        
        results = {}
        
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_angle(**run_args)
        
        q0 = run_args['q0']
        results['f_chi'] = line.f_chi

        if 'sigma_q' not in run_args or run_args['sigma_q'] is None:
            # Determine sigma_q based on circular_average_q2I_fit
            
            head, tail = os.path.split(output_dir)
            result_file = os.path.join(head, 'results', '{}.xml'.format(data.name))
            prev_results = tools.get_result_xml(result_file, 'circular_average_q2I_fit')
            run_args['sigma_q'] = prev_results['fit_peaks_sigma1']
        
        results['sigma_q'] = run_args['sigma_q']
        sigma_chi_rad = np.arctan2(run_args['sigma_q'], q0)
        #sigma_chi_rad = run_args['sigma_q']/q0 # Small angle approximation
        results['sigma_chi'] = np.degrees(np.abs(sigma_chi_rad))

        # Subtract background
        if 'background' not in run_args or run_args['background'] is None:
            # Determine background based on circular_average_q2I_fit
            
            head, tail = os.path.split(output_dir)
            result_file = os.path.join(head, 'results', '{}.xml'.format(data.name))
            prev_results = tools.get_result_xml(result_file, 'circular_average_q2I_fit')
            m, b = prev_results['fit_peaks_m'], prev_results['fit_peaks_b']
            qp, qalpha = prev_results['fit_peaks_qp'], prev_results['fit_peaks_qalpha']
            qn = prev_results['qn_power']
            
            background = (m*q0 + b) + qp*np.power(np.abs(q0), qalpha)
            run_args['background'] = background/np.power(q0, qn)
            
        line.y -= run_args['background']*0.9 # Slightly 'under-estimate' the background (since there is a q-dependency that we are averaging over in the binning process)

        # First point often spurious
        line.x = line.x[1:]
        line.y = line.y[1:]
        line.mask_fractions = line.mask_fractions[1:]
        
        # Smooth curve (but not by more than sigma_chi)
        if run_args['smooth'] is not None and run_args['smooth']>0:
            dbins = run_args['smooth']*results['sigma_chi']/line.dchi
            line.smooth(dbins)
        
        # Compute the standard deviation of the intensity variation along the ring
        sigma_R = np.std(line.y)/np.average(line.y)
        
        # N_g = c_sigmaR * (sigma_R)^beta
        # c.f. Eq. (30) (page 5) of https://doi.org/10.1107/S1600576714020822  
        beta = -2
        c_sigmaR_numerator = 2*( np.sqrt(np.pi) - results['sigma_chi'])*q0
        c_sigmaR_denominator = run_args['multiplicity']*np.sqrt(2*np.pi)*run_args['sigma_q']*sigma_chi_rad*results['f_chi']
        c_sigmaR = c_sigmaR_numerator/c_sigmaR_denominator
        
        N_g = c_sigmaR*np.power( sigma_R, beta )
        
        results['average'] = np.average(line.y)
        results['sigma'] = np.std(line.y)
        results['sigma_R'] = sigma_R
        results['N_g'] = N_g
        
        if run_args['scattering_volume'] is not None:
            results['V'] = run_args['scattering_volume'] # mm
            results['N_g_per_V'] = N_g/run_args['scattering_volume'] # mm^-3
            results['grain_size_V'] = np.power(3*run_args['scattering_volume']/(4*np.pi*N_g), 1/3)*1e3 # um


        
        
        
        
        
        # TODO: I_sort and related analysis
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile, ztrim=run_args['ztrim'])
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        if 'smooth' in run_args:
            line.smooth(run_args['smooth'])
        
        if 'polarization_correction' in run_args and run_args['polarization_correction']:
            chi_deg = line.x
            two_theta_rad = 2.0*np.arcsin(q0/(2.*data.calibration.get_k()))
            P_h = 1 - np.square(np.sin(two_theta_rad))*np.square(np.sin(np.radians(chi_deg)))
            line.y *= P_h


        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            #line.plot(save=outfile, **run_args)
            
            outfile = self.get_outfile(data.name, output_dir) #, ext='_graininess.png')
            line.results = results
            line.plot_graininess(save=outfile, **run_args)
            
            #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
            #line.plot_polar(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        return results                
                


class linecut_qr(Protocol):

    def __init__(self, name='linecut_qr', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'plot_range' : [None, None, 0, None]
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qr(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        #line.smooth(2.0, bins=10)

        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        return results
                                
                                
class linecut_qz(Protocol):

    def __init__(self, name='linecut_qz', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name

        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'plot_range' : [None, None, 0, None]
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qz(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        return results                


class linecut_q(Protocol):

    def __init__(self, name='linecut_q', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'plot_range' : [None, None, 0, None]
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_q(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)

        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)

        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)          
        
        return results




class linecut_qr_fit(linecut_qr, fit_peaks):
    '''Takes a linecut along qr, and fits the data to a simple model
    (Gaussian peak with background).'''
    

    def __init__(self, name='linecut_qr_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False ,
            'plot_range' : [None, None, 0, None] ,
            'num_curves' : 1 ,
            }
        self.run_args.update(kwargs)    
    
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qr(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        #line.smooth(2.0, bins=10)
        
        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
        
        lines = self._fit(line, results, **run_args)
        #lines = DataLines([line])

        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            lines.plot(save=outfile, **run_args)        
        
        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)          
        
        
        return results

        #End class linecut_qr_fit(linecut_qr)
        ########################################



# TODO: Remove
class _old_linecut_qr_fit(linecut_qr):
    '''Takes a linecut along qr, and fits the data to a simple model
    (Gaussian peak with background).'''
    

    def __init__(self, name='linecut_qr_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False ,
                         'plot_range' : [None, None, 0, None] ,
                         'auto_plot_range_fit' : True ,
                         }
        self.run_args.update(kwargs)    
    
    
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qr(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        line.plot(save=outfile, **run_args)
        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        
        # Fit data
        #if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
            #line.trim(run_args['fit_range'][0], run_args['fit_range'][1])
        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
        
        lm_result, fit_line, fit_line_extended = self._fit_peaks(line, **run_args)
        
        
        # Save fit results
        fit_name = 'fit_peaks'
        prefactor_total = 0
        for param_name, param in lm_result.params.items():
            results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            if 'prefactor' in param_name:
                prefactor_total += np.abs(param.value)
            
        results['{}_prefactor_total'.format(fit_name)] = prefactor_total
        results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
        
        # Calculate some additional things
        d = 0.1*2.*np.pi/results['{}_x_center1'.format(fit_name)]['value']
        results['{}_d0'.format(fit_name)] = d
        xi = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/results['{}_sigma1'.format(fit_name)]['value']
        results['{}_grain_size'.format(fit_name)] = xi       
        
        
        # Plot and save data
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                v_spacing = (yf-yi)*0.10
                
                yp = yf
                s = '$q_0 = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_x_center1']['value'])
                self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')

                yp -= v_spacing
                s = r'$d_0 \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_d0'])
                self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')

                yp -= v_spacing
                s = '$\sigma = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_sigma1']['value'])
                self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')
                
                yp -= v_spacing
                s = r'$\xi \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_grain_size'])
                self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')

        
        lines = DataLines_current([line, fit_line, fit_line_extended])
        lines.copy_labels(line)
        lines.results = results

        outfile = self.get_outfile(data.name+'-fit', output_dir, ext='.png')
        
        # Tweak the plotting range for the fit-plot
        run_args_cur = run_args.copy()
        if run_args['auto_plot_range_fit']:
            run_args_cur['plot_range'] = [ run_args['plot_range'][0] , run_args['plot_range'][1] , run_args['plot_range'][2] , run_args['plot_range'][3] ]
            if 'fit_range' in run_args_cur:
                span = abs(run_args['fit_range'][1]-run_args_cur['fit_range'][0])
                run_args_cur['plot_range'][0] = run_args['fit_range'][0]-span*0.25
                run_args_cur['plot_range'][1] = run_args_cur['fit_range'][1]+span*0.25
            
            run_args_cur['plot_range'][2] = 0
            run_args_cur['plot_range'][3] = max(fit_line.y)*1.3
        
        
        try:
            #lines.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
            lines.plot(save=outfile, **run_args_cur)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        #print(results)
        return results        
        
        
    def _fit_peaks(self, line, num_curves=1, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self._fit_peaks(line, **run_args)

        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            
            # Linear background
            m = v['m']*x + v['b']
            
            # Power-law background
            m += v['qp']*np.power( np.abs(x), v['qalpha'] )
            
            # Gaussian peaks
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()


        m = (line.y[-1]-line.y[0])/(line.x[-1]-line.x[0])
        b = line.y[0] - m*line.x[0]

        xs = np.abs(line.x)
        ys = line.y
        qalpha = (np.log(ys[0])-np.log(ys[-1]))/(np.log(xs[0])-np.log(xs[-1]))
        qp = np.exp( np.log(ys[0]) - qalpha*np.log(xs[0]) )

        if True:
            # Linear background
            params.add('m', value=m, min=abs(m)*-4, max=abs(m)*+4+1e-12, vary=False)
            params.add('b', value=b, min=0, max=np.max(line.y)*100+1e-12, vary=False)
            
            params.add('qp', value=0, vary=False)
            params.add('qalpha', value=1.0, vary=False)
            
        else:
            # Power-law background
            params.add('m', value=0, vary=False)
            params.add('b', value=0, vary=False)
            
            params.add('qp', value=qp, vary=False)
            params.add('qalpha', value=qalpha, vary=False)
            
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        
        # Best guess for peak position
        if True:
            # Account for power-law scaling (Kratky-like)
            xs = np.asarray(line.x)
            ys = np.asarray(line.y)
            
            ys = ys*np.power( np.abs(xs), np.abs(qalpha) ) # Kratky-like
            
            # Sort
            indices = np.argsort(ys)
            x_sorted = xs[indices]
            y_sorted = ys[indices]
            
            target = np.max(ys)

            # Search through y for the target
            idx = np.where( y_sorted>=target )[0][0]
            xpeak = x_sorted[idx]
            ypeak = y_sorted[idx]
            
            xpeak, ypeak = line.target_x(xpeak)
                                 

        prefactor = ypeak - ( m*xpeak + b )
        sigma = 0.05*xspan
        
        for i in range(num_curves):
            
            params.add('prefactor{:d}'.format(i+1), value=prefactor, min=0, max=np.max(line.y)*1.5, vary=False)
            params.add('x_center{:d}'.format(i+1), value=xpeak, min=np.min(line.x), max=np.max(line.x), vary=False)
            params.add('sigma{:d}'.format(i+1), value=sigma, min=0, max=xspan*0.75, vary=False)
        
        
        # Fit only the peak width
        params['sigma1'].vary = True
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if True:
            # Tweak peak position
            lm_result.params['sigma1'].vary = False
            lm_result.params['x_center1'].vary = True
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
        
        if True:
            # Relax entire fit
            lm_result.params['m'].vary = True
            lm_result.params['b'].vary = True
            #lm_result.params['qp'].vary = True
            #lm_result.params['qalpha'].vary = True
            
            lm_result.params['prefactor1'].vary = True
            lm_result.params['sigma1'].vary = True
            lm_result.params['x_center1'].vary = True
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        #fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=200)
        fit_x = np.linspace(np.average( [np.min(line_full.x), np.min(line.x)] ), np.average( [0, np.max(line.x)] ), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'alpha':0.5, 'marker':None, 'linewidth':2.0})        

        return lm_result, fit_line, fit_line_extended         
    

        #End class linecut_qr_fit(linecut_qr)
        ########################################
        



class linecut_qz_fit(linecut_qz): # TODO: Use class fit_peaks

    def __init__(self, name='linecut_qz_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'plot_range' : [None, None, 0, None],
            'show_guides' : True,
            'markersize' : 0,
            'linewidth' : 1.5,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        # Usage example:
        # linecut_qz_fit(show_region=False, show=False, qr=0.009, dq=0.0025, q_mode='qr', fit_range=fit_range, q0_expected=q0_expected, plot_range=[0, 0.08, 0, None]) ,

        
        results = {}
        
        line = data.linecut_qz(**run_args)
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
        
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, line, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            line.plot(save=outfile, **run_args)

        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)

        
        if 'incident_angle' not in run_args:
            run_args['incident_angle'] = data.calibration.incident_angle
            
            import re
            filename_re = re.compile('^.+_th(-?\d+\.\d+)_.+$')
            m = filename_re.match(data.name)
            if m:
                run_args['incident_angle'] = float(m.groups()[0])
                
        if 'verbosity' in run_args and run_args['verbosity']>=4:
            print('    Using incident_angle = {:.3f} degrees'.format(run_args['incident_angle']))
                
        #if 'critical_angle_film' not in run_args:
            #run_args['critical_angle_film'] = 0
        #if 'critical_angle_substrate' not in run_args:
            #run_args['critical_angle_substrate'] = 0
                
        
        
        
        # Fit data
        lm_result, fit_line, fit_line_extended = self._fit_peaks(line, **run_args)
        
        
        # Save fit results
        fit_name = 'fit_peaks'
        prefactor_total = 0
        for param_name, param in lm_result.params.items():
            results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            if 'prefactor' in param_name:
                prefactor_total += np.abs(param.value)
            
        results['{}_prefactor_total'.format(fit_name)] = prefactor_total
        results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
        
        # Calculate some additional things
        q0 = results['{}_x_center1'.format(fit_name)]['value']
        d = 0.1*2.*np.pi/q0
        results['{}_d0'.format(fit_name)] = d
        xi = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/results['{}_sigma1'.format(fit_name)]['value']
        results['{}_grain_size'.format(fit_name)] = xi       
        

        def angle_to_q(two_theta_s_rad):
            k = data.calibration.get_k()
            qz = 2*k*np.sin(two_theta_s_rad/2.0)
            return qz      
        def q_to_angle(q):
            k = data.calibration.get_k()
            two_theta_s_rad = 2.0*np.arcsin(q/(2.0*k))
            return two_theta_s_rad
        
        if 'critical_angle_film' in run_args:
            # Account for refraction distortion
            
            theta_incident_rad = np.radians(run_args['incident_angle'])
            theta_c_f_rad = np.radians(run_args['critical_angle_film'])
            #theta_c_s_rad = np.radians(run_args['critical_angle_substrate'])
            
            alpha_i_rad = np.arccos( np.cos(theta_incident_rad)/np.cos(theta_c_f_rad) )
            
            
            # Scattering from incident (refracted) beam
            two_theta_s_rad = q_to_angle(q0)
            theta_f_rad = two_theta_s_rad - theta_incident_rad
            
            alpha_f_rad = np.arccos( np.cos(theta_f_rad)/np.cos(theta_c_f_rad) )
            
            two_alpha_s_rad = alpha_i_rad + alpha_f_rad
            qT = angle_to_q(two_alpha_s_rad)
            results['{}_qT'.format(fit_name)] = qT
            results['{}_dT'.format(fit_name)] = 0.1*2.*np.pi/qT
            
            
            # Scattering from reflected beam
            two_alpha_s_rad = abs( alpha_f_rad-alpha_i_rad )
            qR = angle_to_q(two_alpha_s_rad)
            results['{}_qR'.format(fit_name)] = qR
            results['{}_dR'.format(fit_name)] = 0.1*2.*np.pi/qR

            
        
        
        
        # Plot and save data
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                y_fit_max = np.max(self.lines[1].y)
                yf = y_fit_max*2.0
                v_spacing = (yf-yi)*0.06
                
                q0 = self.results['fit_peaks_x_center1']['value']
                color = 'purple'
                
                yp = yf
                s = '$q_0 = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(q0)
                self.ax.text(xf, yp, s, size=15, color=color, verticalalignment='top', horizontalalignment='right')

                yp -= v_spacing
                s = r'$d_0 \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_d0'])
                self.ax.text(xf, yp, s, size=15, color=color, verticalalignment='top', horizontalalignment='right')

                yp -= v_spacing
                s = '$\sigma = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_sigma1']['value'])
                self.ax.text(xf, yp, s, size=15, color=color, verticalalignment='top', horizontalalignment='right')
                
                yp -= v_spacing
                s = r'$\xi \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_grain_size'])
                self.ax.text(xf, yp, s, size=15, color=color, verticalalignment='top', horizontalalignment='right')
                
                
                self.ax.axvline(q0, color=color, linewidth=0.5)
                self.ax.text(q0, yf, '$q_0$', size=20, color=color, horizontalalignment='center', verticalalignment='bottom')
                
                if 'critical_angle_film' in self.run_args:
                    s = '$q_T = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$ \n $d_T = \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_qT'], self.results['fit_peaks_dT'])
                    self.ax.text(q0, y_fit_max, s, size=15, color='b', horizontalalignment='left', verticalalignment='bottom')
                    self.ax.plot( [self.results['fit_peaks_qT'], q0], [y_fit_max, y_fit_max], '-', color='b' )

                    s = '$q_R = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$ \n $d_R = \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_qR'], self.results['fit_peaks_dR'])
                    self.ax.text(q0, 0, s, size=15, color='r', horizontalalignment='left', verticalalignment='bottom')
                    self.ax.plot( [self.results['fit_peaks_qR'], q0], [yi, yi], '-', color='r' )
                
                
                
                if self.run_args['show_guides']:
                    # Show various guides of scattering features
                    theta_incident_rad = np.radians(self.run_args['incident_angle'])
                    
                    # Direct
                    qz = 0
                    self.ax.axvline( qz, color='0.25' )
                    self.ax.text(qz, yf, '$\mathrm{D}$', size=20, color='0.25', horizontalalignment='center', verticalalignment='bottom')

                    # Horizon
                    qz = angle_to_q(theta_incident_rad)
                    l = self.ax.axvline( qz, color='0.5' )
                    l.set_dashes([10,6])
                    self.ax.text(qz, yf, '$\mathrm{H}$', size=20, color='0.5', horizontalalignment='center', verticalalignment='bottom')

                    # Specular beam
                    qz = angle_to_q(2*theta_incident_rad)
                    self.ax.axvline( qz, color='r' )
                    self.ax.text(qz, yf, '$\mathrm{R}$', size=20, color='r', horizontalalignment='center', verticalalignment='bottom')


                    if 'critical_angle_film' in self.run_args:
                        theta_c_f_rad = np.radians(self.run_args['critical_angle_film'])
                    
                        # Transmitted (direct beam refracted by film)
                        if theta_incident_rad<=theta_c_f_rad:
                            qz = angle_to_q(theta_incident_rad) # Horizon
                        else:
                            alpha_i_rad = np.arccos( np.cos(theta_incident_rad)/np.cos(theta_c_f_rad) )
                            two_theta_s_rad = theta_incident_rad - alpha_i_rad
                            qz = angle_to_q(two_theta_s_rad)
                        l = self.ax.axvline( qz, color='b' )
                        l.set_dashes([4,4])

                        # Yoneda
                        qz = angle_to_q(theta_incident_rad+theta_c_f_rad)
                        self.ax.axvline( qz, color='gold' )
                        self.ax.text(qz, yf, '$\mathrm{Y}_f$', size=20, color='gold', horizontalalignment='center', verticalalignment='bottom')
                    
                    if 'critical_angle_substrate' in self.run_args:
                        theta_c_s_rad = np.radians(self.run_args['critical_angle_substrate'])

                        # Transmitted (direct beam refracted by substrate)
                        if theta_incident_rad<=theta_c_s_rad:
                            qz = angle_to_q(theta_incident_rad) # Horizon
                        else:
                            alpha_i_rad = np.arccos( np.cos(theta_incident_rad)/np.cos(theta_c_s_rad) )
                            two_theta_s_rad = theta_incident_rad - alpha_i_rad
                            qz = angle_to_q(two_theta_s_rad)
                        self.ax.axvline( qz, color='b' )
                        self.ax.text(qz, yf, '$\mathrm{T}$', size=20, color='b', horizontalalignment='center', verticalalignment='bottom')

                        # Yoneda
                        qz = angle_to_q(theta_incident_rad+theta_c_s_rad)
                        self.ax.axvline( qz, color='gold' )
                        self.ax.text(qz, yf, '$\mathrm{Y}_s$', size=20, color='gold', horizontalalignment='center', verticalalignment='bottom')
                        

                
                self.ax.axis([xi, xf, yi, yf])
                

        lines = DataLines_current([line, fit_line, fit_line_extended])
        lines.copy_labels(line)
        lines.results = results
        lines.run_args = run_args

        
        
        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name+'-fit', output_dir, ext='.png')
            #lines.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
            lines.plot(save=outfile, **run_args)

        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)  
        
        #print(results)
        
        return results        
        


                
    def _fit_peaks(self, line, num_curves=1, **run_args):
        # TODO: Use class fit_peaks
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)
        
        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
            
        
        import lmfit
        
        def model(v, x):
            '''Gaussians with constant background.'''
            m = v['m']*x + v['b']
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        
        m = (line.y[-1]-line.y[0])/(line.x[-1]-line.x[0])
        b = line.y[0] - m*line.x[0]
        
        #params.add('m', value=0)
        #params.add('b', value=np.min(line.y), min=0, max=np.max(line.y))
        params.add('m', value=m, min=abs(m)*-10, max=abs(m)*+10)
        #params.add('m', value=m, min=abs(m)*-5, max=1e-12) # Slope must be negative
        params.add('b', value=b, min=0)
        
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        
        if 'q0_expected' in run_args and run_args['q0_expected'] is not None:
            xpeak, ypeak = line.target_x(run_args['q0_expected'])
            
        prefactor = ypeak - (m*xpeak+b)
        
        for i in range(num_curves):
            #xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            #xpos, ypos = line.target_x(xpeak*(i+1))
            xpos, ypos = xpeak, ypeak
            
            params.add('prefactor{:d}'.format(i+1), value=prefactor, min=0, max=np.max(line.y)*4)
            params.add('x_center{:d}'.format(i+1), value=xpos, min=np.min(line.x), max=np.max(line.x))
            #params.add('x_center{:d}'.format(i+1), value=-0.009, min=np.min(line.x), max=np.max(line.x))
            params.add('sigma{:d}'.format(i+1), value=0.003, min=0, max=xspan*0.5)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        # https://lmfit.github.io/lmfit-py/fitting.html
        #lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y), method='nelder')
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'purple', 'marker':None, 'linewidth':4.0})
        
        span = abs( np.max(line.x) - np.min(line.x) )
        fit_x = np.linspace(np.min(line.x)-0.5*span, np.max(line.x)+0.5*span, num=500)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'purple', 'alpha':0.5, 'marker':None, 'linewidth':2.0})        

        return lm_result, fit_line, fit_line_extended       


        #End class linecut_qz_fit(linecut_qz)
        ########################################





class linecut_angle_fit(linecut_angle):

    def __init__(self, name='linecut_angle_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'show_region' : False,
            'plot_range' : [-90, 90, None, None],
            'do_max' : True,
            'do_fits' : True,
            'markersize' : 0,
            'linewidth' : 1.5,
            'between_convention' : False,
            }
        self.run_args.update(kwargs)
        
        
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'q0' not in run_args or run_args['q0'] is None:
            # Determine q based on circular_average_q2I_fit
            
            head, tail = os.path.split(output_dir)
            result_file = os.path.join(head, 'results', '{}.xml'.format(data.name))
            prev_results = tools.get_result_xml(result_file, 'circular_average_q2I')
            
            run_args['q0'] = prev_results['fit_peaks_x_center1']
        
        # Extract circular average
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_angle(**run_args)
        
        
        # Fine-tune data
        #line.kill_x(0, 3) # Remove the data near chi=0 since this is artificially high (specular ridge)
        line.x = line.x[1:]
        line.y = line.y[1:]
        line.smooth(2.0)
        
        
        if 'show_region' in run_args:
            if run_args['show_region']=='save':
                outfile = self.get_outfile(data.name, output_dir, ext='_region.png')
                data.plot(save=outfile)
            elif run_args['show_region']:
                data.plot(show=True)
            
        
        if 'txt' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        
        
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                #v_spacing = 0.1
                xi, xf, yi, yf = self.ax.axis()
                #yi = min( self.lines[-1].y )*0.5
                #yf = max( self.lines[-1].y )*1.5
                #self.ax.axis( [xi, xf, yi, yf] )

                if run_args['do_fits']:
                    yp = yf
                    s = '$\eta = {:.2f}$'.format(self.results['fit_eta_eta']['value'])
                    self.ax.text(xi, yp, s, size=30, color='b', verticalalignment='top', horizontalalignment='left')

                    yp = yf
                    s = '$m = {:.2f}$'.format(self.results['fit_MaierSaupe_m']['value'])
                    self.ax.text(xf, yp, s, size=30, color='purple', verticalalignment='top', horizontalalignment='right')
                    
                    
                if run_args['do_max']:
                    xt, yt = self.results['max_position'], self.results['max_height']
                    fwhm = self.results['fwhm']
                    self.ax.plot([xt], [yt], 'o', color='b', markersize=8)
                    s = r'$\chi_{{ \mathrm{{max}} }} = {:.1f}^{{\circ}}$'.format(xt)
                    self.ax.text(xt, yt, s, color='b', size=20, verticalalignment='bottom', horizontalalignment='left')
                    
                    self.ax.plot([xt-fwhm/2, xt+fwhm/2], [yt/2, yt/2], 'o-', color='b', markersize=8, linewidth=1.0)
                    s = r'$\mathrm{{FWHM}} = {:.1f}^{{\circ}}$'.format(fwhm)
                    self.ax.text(xt+fwhm/2, yt/2, s, color='b', size=20, verticalalignment='center', horizontalalignment='left')
                    
                
                angle = self.results['orientation_angle']
                self.ax.axvline(angle, color='r', linewidth=1.0)
                
                s = '$\chi = {:.1f} ^{{\circ}} $'.format(angle)
                self.ax.text(angle, yf, s, size=20, color='r', verticalalignment='top', horizontalalignment='right')
                s = '$f = {:.3f}$'.format(self.results['orientation_factor'])
                self.ax.text(angle, yi, s, size=20, color='r', verticalalignment='bottom', horizontalalignment='right')
                
                #self.ax.set_xticks([-180, -90, 0, +90, +180])
                
                
                
        lines = DataLines_current()
        lines.add_line(line)
        lines.copy_labels(line)
        angle = None
            
            
        xt, yt = line.target_y(max(line.y))
        angle = xt
            
        if run_args['do_max']:
            # FWHM
            hm = yt*0.5

            # Right (+) side
            line_temp = line.sub_range(xt, xt+60)
            if len(line_temp.y)>0:
                xtr, ytr = line_temp.target_y(hm)
            else:
                xtr = None
            
            # Left (-) side
            line_temp = line.sub_range(xt-60, xt)
            if len(line_temp.y)>0:
                xtl, ytl = line_temp.target_y(hm)
                if xtr is None:
                    xtr, ytr = 2*xt-xtl, ytl
            else:
                xtl, ytl = 2*xt-xtr, ytr
                
            xtl, ytl = 2*xt-xtr, ytr
            
            results['max_position'] = xt
            results['max_height'] = yt
            results['fwhm'] = abs(xtr-xtl)
            
            
            
            
        if run_args['do_fits']:
            color_list = ['b', 'purple', 'r', 'green', 'orange',]
            for i, fit_name in enumerate(['fit_eta', 'fit_MaierSaupe', 'fit_eta_span']):
                
                lm_result, fit_line, fit_line_e = getattr(self, fit_name)(line, **run_args)
                fit_line_e.plot_args['color'] = color_list[i%len(color_list)]
                lines.add_line(fit_line_e)
                if fit_name=='fit_eta_span':
                    lines.add_line(fit_line)
                
            
                #prefactor_total = 0
                for param_name, param in lm_result.params.items():
                    results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                    #if 'prefactor' in param_name:
                        #prefactor_total += np.abs(param.value)
                    
                #results['{}_prefactor_total'.format(fit_name)] = prefactor_total
                results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
                
                #if fit_name=='fit_eta':
                if fit_name=='fit_eta' or fit_name=='fit_eta_span':
                    angle = lm_result.params['x_center'].value
                
        
        if angle<0:
            angle += 180
        # angle is now 0 to +180
        
        if run_args['between_convention']:
            # Convert to the 'in between the peaks' convention
            angle -= 90 
            # angle is now -90 to +90        
            
        qz_unit = np.cos(np.radians(angle))
        qx_unit = np.sin(np.radians(angle))
        orientation_factor = 2*np.square(qz_unit) - 1
        results['orientation_angle'] = angle
        results['orientation_factor'] = orientation_factor
            
            
        # Output
        if 'plots' in run_args['save_results']:
            self.label_filename(data, lines, **run_args)
            outfile = self.get_outfile(data.name, output_dir)
            lines.lines.reverse() # Puts the data on top
            lines.results = results
            lines.plot(save=outfile, **run_args)

            if run_args['verbosity']>=4:
                outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
                line.plot_polar(save=outfile, **run_args)

        if 'hdf5' in run_args['save_results']:
            self.save_DataLine_HDF5(line, data.name, output_dir, results=results)  

        return results
    
    
    def fit_eta(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['prefactor']*( 1 - (v['eta']**2) )/( ((1+v['eta'])**2) - 4*v['eta']*( np.square(np.cos(  (v['symmetry']/2.0)*np.radians(x-v['x_center'])  )) ) ) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data

        
        xpeak, ypeak = line.target_y(np.max(line.y))
        
        params = lmfit.Parameters()
        params.add('prefactor', value=ypeak*0.5, min=0)
        params.add('x_center', value=xpeak, min=np.min(line.x), max=np.max(line.x), vary=True)
        params.add('eta', value=0.4, min=0, max=1)
        params.add('symmetry', value=2, min=0.5, max=20, vary=False)
        params.add('baseline', value=0, min=0, max=np.max(line.y)+1e-10, vary=False)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=1000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})

        return lm_result, fit_line, fit_line_extended
                    
                    
    def fit_eta_span(self, line, span=30, **run_args):
        '''Fit the data with an "eta orientation" function, but over a limited span of angle.
        This guards against spurious fits for distributions that do not strictly match the
        eta function.'''
        
        import lmfit
        
        xpeak, ypeak = line.target_y(np.max(line.y))
        if xpeak<0:
            xpeak += 180
        line_full = line
        
        # Expand the line so that we don't have to worry about edges
        x = np.asarray(line.x)
        x = np.concatenate( [x-360, x, x+360] )
        y = np.asarray(line.y)
        y = np.concatenate( [y, y, y] )
        line_extended = DataLine(x=x, y=y)
        
        region_right = line_extended.sub_range( xpeak-span/2, xpeak+span/2 )
        region_left = line_extended.sub_range( (xpeak-180)-span/2, (xpeak-180)+span/2 )
        
        x = np.concatenate( [region_left.x, region_right.x] )
        y = np.concatenate( [region_left.y, region_right.y] )
        line = DataLine( x=x, y=y )
        
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['prefactor']*( 1 - (v['eta']**2) )/( ((1+v['eta'])**2) - 4*v['eta']*( np.square(np.cos(  (v['symmetry']/2.0)*np.radians(x-v['x_center'])  )) ) ) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        
        params = lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y)*0.5, min=0)
        params.add('x_center', value=xpeak, min=np.min(line.x), max=np.max(line.x), vary=True)
        params.add('eta', value=0.4, min=0, max=1)
        params.add('symmetry', value=2, min=0.5, max=20, vary=False)
        params.add('baseline', value=0, min=0, max=np.max(line.y)+1e-10, vary=False)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':'o', 'markersize':8.0, 'linewidth':0.0})
        
        fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=1000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':1.5})

        return lm_result, fit_line, fit_line_extended
    
           
    def fit_MaierSaupe(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''orientation function.'''
            m = v['prefactor']*np.exp(v['m']*np.square(np.cos(np.radians(x-v['x_center'])))) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        
        xpeak, ypeak = line.target_y(np.max(line.y))
        
        params = lmfit.Parameters()
        params.add('prefactor', value=ypeak*0.1, min=0)
        params.add('x_center', value=xpeak, min=np.min(line.x), max=np.max(line.x))
        params.add('m', value=2.0, min=0)
        params.add('baseline', value=0, min=0, max=np.max(line.y)+1e-10, vary=False)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line.x), np.max(line.x), num=1000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})

        return lm_result, fit_line, fit_line_extended

        #End class linecut_angle_fit(linecut_angle)
        ########################################

                
                
class calibration_check(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'dq': 0.01,
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'resize' in run_args and (run_args['resize'] is not None):
            data.resize(run_args['resize']) # Shrink
        
        
        outfile = self.get_outfile('{}_full'.format(data.name), output_dir)
        data.plot_image(outfile, **run_args)

        
        #data.blur(2.0)
        
        dq = run_args['dq']

        if 'AgBH' in run_args and run_args['AgBH']:
            q0 = 0.1076 # A^-1
            
            for i in range(11):
                data.overlay_ring(q0*(i+1), q0*(i+1)*dq)
                
        if 'CeO2' in run_args and run_args['CeO2']:
            # FCC structure, a = 5.411 Angstroms
            # Ref: https://community.dur.ac.uk/john.evans/topas_workshop/data/ceo2.cif
            q0 = (2*np.pi/5.411)*np.sqrt(3) # A^-1, (111)

            qlist = ( q0/np.sqrt(3) )*np.array((np.sqrt(3), 2, np.sqrt(8),np.sqrt(11),np.sqrt(12),np.sqrt(16),np.sqrt(19),  np.sqrt(20)))
            data.overlay_ring(q0, q0*dq)
            for q in qlist:
                data.overlay_ring(q, q*dq)                
                
        if 'q0'  in run_args:
            q0 = run_args['q0']
            if 'num_rings' in run_args:
                num_rings = run_args['num_rings']
            else:
                num_rings = 5

            for i in range(num_rings):
                data.overlay_ring(q0*(i+1), q0*(i+1)*dq)

        
        outfile = self.get_outfile(data.name, output_dir)

        data.plot(save=outfile, **run_args)
        
        return results
                                       
                
                
                
                
                
                
                
                
class fit_calibration(Protocol):
    # Work in progress
    ################################################################################                    

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'material' : 'AgBH01',
            }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        # WARNING: This procedure doesn't work very well.
        
        results = {}
        
        if run_args['material'] is 'AgBH01':
            
            import lmfit
            # https://lmfit.github.io/lmfit-py/parameters.html#simple-example
            def fcn2min(params, x, data):
                '''Gaussian with linear background.'''
                v = params.valuesdict()
                model = v['prefactor']*np.exp( -np.square(x-v['x_center'])/(2*(v['sigma']**2)) ) + v['m']*x + v['b']
                
                return model - data
            
            self.lmfit = lmfit
            self._fcn2min = fcn2min
            
            q_peak = 0.1076 # A^-1
            dq = q_peak*0.3
            
            self._find_local_q_match(data, 'distance_m', 0.02, 0.002, q_peak, dq)
            self._find_local_q_match(data, 'distance_m', 0.001, 0.0005, q_peak, dq)
            
            #self._find_local_minimum_width(data, 'x0', 40, 2, q_peak, dq)
            #self._find_local_minimum_width(data, 'y0', 40, 2, q_peak, dq)
            #self._find_local_minimum_width(data, 'x0', 5, 0.5, q_peak, dq)
            #self._find_local_minimum_width(data, 'y0', 5, 0.5, q_peak, dq)
            
            self._find_local_q_match(data, 'distance_m', 0.005, 0.0002, q_peak, dq)
            
            
            print('Final values:\n    x0 = %.1f\n    y0 = %.1f\n    dist = %g'%(data.calibration.x0, data.calibration.y0, data.calibration.distance_m))
            
            
            line = data.circular_average_q_range(q_peak, dq, error=False)
            result = self._fit_peak(line)
            self.lmfit.report_fit(result.params)
            fit_data = line.y + result.residual
            fit_line = DataLine(x=line.x, y=fit_data, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
            
            outfile = self.get_outfile(data.name, output_dir)
            lines = DataLines( [line, fit_line] )
            lines.plot(save=outfile, show=False)
            outfile = self.get_outfile(data.name, output_dir, ext='.dat')
            line.save_data(outfile)
        
        
        return results
    
    
    def _fit_peak(self, line):
        
        params = self.lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y), min=0)
        params.add('x_center', value=np.average(line.x))
        params.add('sigma', value=np.std(line.x), min=0)
        params.add('m', value=0)
        params.add('b', value=0)
        
        result = self.lmfit.minimize(self._fcn2min, params, args=(line.x, line.y))
        
        return result
        
    
    def _find_local_minimum_width(self, data, attr, spread, step, q_peak, dq):
        
        v_start = getattr(data.calibration, attr)
        v_min = None
        v_min_value = None
        
        for v_displacement in np.arange(-spread, +spread, step):
            
            setattr(data.calibration, attr, v_start + v_displacement)
            data.calibration.clear_maps()
            
            line = data.circular_average_q_range(q_peak, dq, error=False)
            result = self._fit_peak(line)
            width = result.params['sigma'].value
            
            print('  %s: %.1f, width: %g' % (attr, v_start+v_displacement, width))
            
            if v_min is None:
                v_min = v_start+v_displacement
                v_min_value = width
            elif width<v_min_value:
                v_min = v_start+v_displacement
                v_min_value = width
                
        # Go to the best position
        setattr(data.calibration, attr, v_min)
        data.calibration.clear_maps()
               
               
    def _find_local_q_match(self, data, attr, spread, step, q_peak, dq):
                
        v_start = getattr(data.calibration, attr)
        v_min = None
        v_min_value = None
        
        for v_displacement in np.arange(-spread, +spread, step):
            
            setattr(data.calibration, attr, v_start + v_displacement)
            data.calibration.clear_maps()
            
            line = data.circular_average_q_range(q_peak, dq, error=False)
            result = self._fit_peak(line)
            err = abs(result.params['x_center'].value - q_peak)
            
            print('  %s: %g, pos: %g, err: %g' % (attr, v_start+v_displacement, result.params['x_center'], err))
            
            if v_min is None:
                v_min = v_start+v_displacement
                v_min_value = err
            elif err<v_min_value:
                v_min = v_start+v_displacement
                v_min_value = err
                
        # Go to the best position
        setattr(data.calibration, attr, v_min)
        data.calibration.clear_maps()                
       
       
       
       

class q_image(Protocol):
    
    def __init__(self, name='q_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'blur' : None,
            'ztrim' : [0.05, 0.005],
            'method' : 'nearest',
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        
        q_data = data.remesh_q_bin(**run_args)
        
        if run_args['verbosity']>=10:
            # Diagnostic
            
            # WARNING: These outputs are not to be trusted.
            # The maps are oriented relative to data.data (not q_data.data)
            # These are for diagnostic purposes only.
            data_temp = Data2DReciprocal()
            
            data_temp.data = data.calibration.qx_map()
            outfile = self.get_outfile('qx-{}'.format(data.name), output_dir, ext='.png', ir=True)
            r = np.max( np.abs(data_temp.data) )
            data_temp.set_z_display([-r, +r, 'linear', 0.3])
            data_temp.plot(outfile, cmap='bwr', **run_args)

            data_temp.data = data.calibration.qy_map()
            outfile = self.get_outfile('qy-{}'.format(data.name), output_dir, ext='.png', ir=True)
            r = np.max( np.abs(data_temp.data) )
            data_temp.set_z_display([-r, +r, 'linear', 0.3])
            data_temp.plot(outfile, cmap='bwr', **run_args)
            
            data_temp.data = data.calibration.qz_map()
            outfile = self.get_outfile('qz-{}'.format(data.name), output_dir, ext='.png', ir=True)
            r = np.max( np.abs(data_temp.data) )
            data_temp.set_z_display([-r, +r, 'linear', 0.3])
            data_temp.plot(outfile, cmap='bwr', **run_args)
            
        
        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)

        if 'q_max' in run_args and run_args['q_max'] is not None:
            q_max = run_args['q_max']
            run_args['plot_range'] = [-q_max, +q_max, -q_max, +q_max]
        
        q_data.set_z_display([None, None, 'gamma', 0.3])
        q_data.plot_args = { 'rcParams': {'axes.labelsize': 55,
                                    'xtick.labelsize': 40,
                                    'ytick.labelsize': 40,
                                    'xtick.major.pad': 10,
                                    'ytick.major.pad': 10,
                                    },
                            } 


        if 'plot_buffers' not in run_args:
            run_args['plot_buffers'] = [0.30,0.05,0.25,0.05]
        self.label_filename(data, q_data, **run_args)
        q_data.plot(outfile, **run_args)
        
        #if 'save_data' in run_args and run_args['save_data']: # Deprecated
        if 'npz' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)

        if 'hdf5' in run_args['save_results']:
            q_data.name = data.name
            self.save_Data2D_HDF5(q_data, 'q_image', output_dir, results=results)
        
        
        return results
    
    
class qr_image(Protocol):
    
    def __init__(self, name='qr_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'blur' : None,
            'ztrim' : [0.05, 0.005],
            'method' : 'nearest',
            #'save_data' : False, # Deprecated
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'dezing_fill' in run_args and run_args['dezing_fill']:
            data.dezinger(sigma=3, tol=5, mode='median', mask=False, fill=True)
        
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        
        q_data = data.remesh_qr_bin(**run_args)
        
        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)

        if 'q_max' in run_args and run_args['q_max'] is not None:
            q_max = run_args['q_max']
            run_args['plot_range'] = [-q_max, +q_max, -q_max, +q_max]
        
        q_data.set_z_display([None, None, 'gamma', 0.3])
        q_data.plot_args = { 'rcParams': {'axes.labelsize': 55,
                                    'xtick.labelsize': 40,
                                    'ytick.labelsize': 40,
                                    'xtick.major.pad': 10,
                                    'ytick.major.pad': 10,
                                    },
                            } 
        q_data.x_label = 'qr'
        q_data.x_rlabel = '$q_r \, (\mathrm{\AA^{-1}})$'

        if 'plot_buffers' not in run_args:
            run_args['plot_buffers'] = [0.30,0.05,0.25,0.05]
        self.label_filename(data, q_data, **run_args)
        q_data.plot(outfile, **run_args)
        
        #if 'save_data' in run_args and run_args['save_data']: # Deprecated
        if 'npz' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)
        
        if 'hdf5' in run_args['save_results']:
            q_data.name = data.name
            self.save_Data2D_HDF5(q_data, 'q_image', output_dir, results=results)
        
        
        return results
        
        
class q_image_special(q_image):
    
    def __init__(self, name='q_image_special', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'blur' : None,
            'ztrim' : [0.05, 0.005],
            'method' : 'nearest',
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        
        q_data = data.remesh_q_bin(**run_args)
        
        data_hold = data.data
        data.data = data.mask.data
        q_mask = data.remesh_q_bin(**run_args)
        data.data = data_hold
        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)

        if 'q_max' in run_args and run_args['q_max'] is not None:
            q_max = run_args['q_max']
            run_args['plot_range'] = [-q_max, +q_max, -q_max, +q_max]
            
        
        # Determine incident angle
        if 'incident_angle' not in run_args:
            import re
            filename_re = re.compile('^.+_th(-?\d+)_.+')
            m = filename_re.match(data.name)
            if m:
                run_args['incident_angle'] = float(m.groups()[0])/100.0
            else:
                print("ERROR: Couldn't identify theta from filename: {}".format(data.name))
                run_args['incident_angle'] = 0
        
        # Tweak the plotting methods of our q_data object/instance
        q_data.incident_angle = run_args['incident_angle']
        q_data.critical_angle = run_args['critical_angle']
        
        
        def _plot(self, save=None, show=False, ztrim=[0.01, 0.01], size=10.0, plot_buffers=[0.1,0.1,0.1,0.1], **kwargs):
            
            # Data2D._plot()
            
            plot_args = self.plot_args.copy()
            plot_args.update(kwargs)
            self.process_plot_args(**plot_args)
            
            
            self.fig = plt.figure( figsize=(size,size), facecolor='white' )
            left_buf, right_buf, bottom_buf, top_buf = plot_buffers
            fig_width = 1.0-right_buf-left_buf
            fig_height = 1.0-top_buf-bottom_buf
            self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
            
            
            
            # Set zmin and zmax. Top priority is given to a kwarg to this plot function.
            # If that is not set, the value set for this object is used. If neither are
            # specified, a value is auto-selected using ztrim.
            
            values = np.sort( self.data.flatten() )
            if 'zmin' in plot_args and plot_args['zmin'] is not None:
                zmin = plot_args['zmin']
            elif self.z_display[0] is not None:
                zmin = self.z_display[0]
            else:
                zmin = values[ +int( len(values)*ztrim[0] ) ]
                
            if 'zmax' in plot_args and plot_args['zmax'] is not None:
                zmax = plot_args['zmax']
            elif self.z_display[1] is not None:
                zmax = self.z_display[1]
            else:
                idx = -int( len(values)*ztrim[1] )
                if idx>=0:
                    idx = -1
                zmax = values[idx]
                
            if zmax==zmin:
                zmax = max(values)
                
            print( '        data: %.1f to %.1f\n        z-scaling: %.1f to %.1f\n' % (np.min(self.data), np.max(self.data), zmin, zmax) )
            
            self.z_display[0] = zmin
            self.z_display[1] = zmax
            self._plot_z_transform()
                
            
            shading = 'flat'
            #shading = 'gouraud'
            
            if 'cmap' in plot_args:
                cmap = plot_args['cmap']
                
            else:
                # http://matplotlib.org/examples/color/colormaps_reference.html
                #cmap = mpl.cm.RdBu
                #cmap = mpl.cm.RdBu_r
                #cmap = mpl.cm.hot
                #cmap = mpl.cm.gist_heat
                cmap = mpl.cm.jet
            
            x_axis, y_axis = self.xy_axes()
            extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
            
            Zm = np.ma.masked_where(q_mask.data < 0.5, self.Z)
            self.im = plt.imshow(Zm, vmin=0, vmax=1, cmap=cmap, interpolation='nearest', extent=extent, origin='lower')
            #plt.pcolormesh( self.x_axis, self.y_axis, self.Z, cmap=cmap, vmin=zmin, vmax=zmax, shading=shading )
            
            if self.regions is not None:
                for region in self.regions:
                    plt.imshow(region, cmap=mpl.cm.spring, interpolation='nearest', alpha=0.75)
                    #plt.imshow(np.flipud(region), cmap=mpl.cm.spring, interpolation='nearest', alpha=0.75, origin='lower')

            x_label = self.x_rlabel if self.x_rlabel is not None else self.x_label
            y_label = self.y_rlabel if self.y_rlabel is not None else self.y_label
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            
            if 'xticks' in kwargs and kwargs['xticks'] is not None:
                self.ax.set_xticks(kwargs['xticks'])
            if 'yticks' in kwargs and kwargs['yticks'] is not None:
                self.ax.set_yticks(kwargs['yticks'])
            
            
            if 'plot_range' in plot_args:
                plot_range = plot_args['plot_range']
                # Axis scaling
                xi, xf, yi, yf = self.ax.axis()
                if plot_range[0] != None: xi = plot_range[0]
                if plot_range[1] != None: xf = plot_range[1]
                if plot_range[2] != None: yi = plot_range[2]
                if plot_range[3] != None: yf = plot_range[3]
                self.ax.axis( [xi, xf, yi, yf] )
            
            if 'title' in plot_args:
                #size = plot_args['rcParams']['axes.labelsize']
                size = plot_args['rcParams']['xtick.labelsize']
                plt.figtext(0, 1, plot_args['title'], size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
            
            self._plot_extra(**plot_args)
            
            if save:
                if 'transparent' not in plot_args:
                    plot_args['transparent'] = True
                if 'dpi' in plot_args:
                    plt.savefig(save, dpi=plot_args['dpi'], transparent=plot_args['transparent'])
                else:
                    plt.savefig(save, transparent=plot_args['transparent'])
            
            if show:
                self._plot_interact()
                plt.show()
                
            plt.close(self.fig.number)        
            
        
        def _plot_extra(self, **plot_args):
            '''This internal function can be over-ridden in order to force additional
            plotting behavior.'''
            
            self.ax.get_yaxis().set_tick_params(which='both', direction='out')
            self.ax.get_xaxis().set_tick_params(which='both', direction='out')       
            
            xi, xf, yi, yf = self.ax.axis()
            
            xmin, xmax = 0.5, 1.0
            
            # Horizon
            s = '$\mathrm{H}$'
            qz = data.calibration.angle_to_q(self.incident_angle)
            self.ax.axhline(qz, xmin=xmin, xmax=xmax, color='0.5', linewidth=4.0, dashes=[15,15])
            self.ax.text(xf, qz, s, size=30, color='0.5', horizontalalignment='left', verticalalignment='center')
            # Specular
            s = '$\mathrm{R}$'
            qz = data.calibration.angle_to_q(2.0*self.incident_angle)
            self.ax.axhline(qz, xmin=xmin, xmax=xmax, color='r', linewidth=4.0)
            self.ax.text(xf, qz, s, size=30, color='r', horizontalalignment='left', verticalalignment='center')
            # Yoneda
            s = '$\mathrm{Y}$'
            qz = data.calibration.angle_to_q(self.incident_angle+self.critical_angle)
            self.ax.axhline(qz, xmin=xmin, xmax=xmax, color='#bfbf00', linewidth=4.0)
            self.ax.text(xf, qz, s, size=30, color='#bfbf00', horizontalalignment='left', verticalalignment='center')
            # Transmitted beam
            s = '$\mathrm{T}$'
            if self.incident_angle < self.critical_angle:
                qz = data.calibration.angle_to_q(self.incident_angle)
            else:
                numerator = np.cos(np.radians(self.incident_angle))
                denominator = np.cos(np.radians(self.critical_angle))
                alpha_incident = np.degrees(np.arccos(numerator/denominator))
                qz = data.calibration.angle_to_q(self.incident_angle - alpha_incident)
            self.ax.axhline(qz, xmin=xmin, xmax=xmax, color='#5555ff', linewidth=4.0)
            self.ax.text(xf, qz, s, size=30, color='#0000ff', horizontalalignment='left', verticalalignment='center')
                
            
                                 
            
        import types
        q_data._plot = types.MethodType(_plot, q_data)
        q_data._plot_extra = types.MethodType(_plot_extra, q_data)
    
    
        q_data.set_z_display([None, None, 'gamma', 0.3])
        q_data.plot_args = { 'rcParams': {'axes.labelsize': 55,
                                    'xtick.labelsize': 40,
                                    'ytick.labelsize': 40,
                                    'xtick.major.pad': 10,
                                    'ytick.major.pad': 10,
                                    },
                            } 

        if 'plot_buffers' not in run_args:
            run_args['plot_buffers'] = [0.30,0.08,0.25,0.05]
        self.label_filename(data, q_data, **run_args)
        q_data.plot(outfile, **run_args)
        
        
        #if 'save_data' in run_args and run_args['save_data']: # Deprecated
        if 'npz' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)
        
        if 'hdf5' in run_args['save_results']:
            q_data.name = data.name
            self.save_Data2D_HDF5(q_data, 'q_image', output_dir, results=results)
        
        
        return results
            
    
    
class q_phi_image(Protocol):
    
    def __init__(self, name='q_phi_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
            'blur' : None,
            'bins_relative' : 0.5,
            'bins_phi' : 360.0/1.0,
            'ztrim' : [0.05, 0.005],
            'method' : 'nearest',
            'yticks' : [-180, -90, 0, 90, 180],
            #'save_data' : False, # Deprecated
            'save_data_pickle' : True,
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        
        q_data = data.remesh_q_phi(**run_args)
        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)

        if 'q_max' in run_args and run_args['q_max'] is not None:
            run_args['plot_range'] = [0, +run_args['q_max'], -180, +180]
        
        q_data.set_z_display([None, None, 'gamma', 0.3])
        q_data.plot_args = { 'rcParams': {'axes.labelsize': 55,
                                    'xtick.labelsize': 40,
                                    'ytick.labelsize': 40,
                                    },
                            } 

        if 'plot_buffers' not in run_args:
            run_args['plot_buffers'] = [0.20,0.05,0.20,0.05]
        self.label_filename(data, q_data, **run_args)
        q_data.plot(outfile, **run_args)
        
        #if 'save_data' in run_args and run_args['save_data']: # Deprecated
        if 'npz' in run_args['save_results']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)         

        if 'hdf5' in run_args['save_results']:
            q_data.name = data.name
            self.save_Data2D_HDF5(q_data, 'q_image', output_dir, results=results)


        # TODO: Deprecate in favor of 'save_data' .npz file
        #if 'save_data_pickle' in run_args and run_args['save_data_pickle']:
        if 'pkl' in run_args['save_results']:
            # Save Data2DQPhi() object
            import pickle
            outfile = self.get_outfile(data.name, output_dir, ext='.pkl')
            with open(outfile, 'wb') as fout:
                out_data = q_data.data, q_data.x_axis, q_data.y_axis
                pickle.dump(out_data, fout)

        
        return results





class export_STL(Protocol):
    
    def __init__(self, name='export_STL', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.stl'
        self.run_args = {
            'crop' : None ,
            'shift_crop_up' : 0.0 ,
            'blur' : 0.5 ,
            'resize' : 1.0 ,
            'ztrim' : [0.05, 0.005] ,
            'stl_pedestal' : 75.0 ,
            'stl_zscale' : 200.0 ,
            'crop_zone' : None ,
            'crop_beam' : None , 
            'crop_GI' : None ,
            'logo_file' : None ,
            'logo_resize' : 1.0 ,
            }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        # Usage:
        #Protocols.export_STL( stl_zscale=150, stl_pedestal=8, blur=1.5, resize=1.0, crop_GI=700)
        #Protocols.export_STL( stl_zscale=150, stl_pedestal=8, blur=1.5, resize=1.0, crop_beam=[400,400], ztrim=[0.05, 0.006], logo_file=SciAnalysis_PATH+'./examples/STL/logo.png', logo_resize=0.3 )
        
        results = {}
        
        data.data = data.data.astype(float)
        calibration = data.calibration
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'], shift_crop_up=run_args['shift_crop_up'])
        if run_args['crop_zone'] is not None:
            xi, xf, yi, yf = run_args['crop_zone']
        elif run_args['crop_beam'] is not None:
            xw, yw = run_args['crop_beam']
            xi = calibration.x0 - xw/2
            xf = calibration.x0 + xw/2
            yi = calibration.y0 - yw/2
            yf = calibration.y0 + yw/2
            data.data = data.data[ int(yi):int(yf), int(xi):int(xf) ]
        elif run_args['crop_GI'] is not None:
            w = run_args['crop_GI']
            xi = calibration.x0 - w/2
            xf = calibration.x0 + w/2
            yi = calibration.y0 - w
            yf = calibration.y0
            data.data = data.data[ int(yi):int(yf), int(xi):int(xf) ]
            
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize']) # Shrink
        
        
        
        
        # Adjust scaling of data
        data.set_z_display([None, None, 'gamma', 0.2])
        
        if run_args['verbosity']>=3:
        
            outfile = self.get_outfile(data.name, output_dir, ext='.jpg')
            results['files_saved'] = [
                { 'filename': '{}'.format(outfile) ,
                'description' : 'quick view (thumbnail) image' ,
                'type' : 'plot' # 'data', 'plot'
                } ,
                ]
            data.plot_image(outfile, cmap=cmap_vge_hdr, **run_args)
                
        

        if run_args['verbosity']>=4:
            print('        data.data from {:.2f} to {:.2f} ({:.2f} ± {:.2f})'.format(np.min(data.data), np.max(data.data), np.average(data.data), np.std(data.data)))

        if 'clip' in run_args:
            data.data = np.clip(data.data, run_args['clip'][0], run_args['clip'][1])
                                                                                      
        data._plot_z_transform()
        
        if run_args['verbosity']>=4:
            print('        data.Z from {:.2f} to {:.2f} ({:.2f} ± {:.2f})'.format(np.min(data.Z), np.max(data.Z), np.average(data.Z), np.std(data.Z)))
                                                                                      
                                                                                      
        height_map = (data.Z - np.min(data.Z))/(np.max(data.Z)-np.min(data.Z)) # Data now from 0...1
        
        
        if run_args['verbosity']>=4:
            print('        height_map (normed) from {:.2f} to {:.2f} ({:.2f} ± {:.2f})'.format(np.min(height_map), np.max(height_map), np.average(height_map), np.std(height_map)))
            
        
        # Readjust height map output
        height_map = ( height_map*run_args['stl_zscale'] ) + run_args['stl_pedestal']
        height_map *= run_args['resize']
        if run_args['verbosity']>=4:
            print('        height_map (STL scale) from {:.2f} to {:.2f} ({:.2f} ± {:.2f})'.format(np.min(height_map), np.max(height_map), np.average(height_map), np.std(height_map)))
        
        
        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)
        
        #print(data.stats())
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'STL version of data (for 3D printing)' ,
             'type' : 'print' # 'data', 'plot'
            } ,
            ]


        if run_args['logo_file'] is not None:
            
            # Load a logo file
            import scipy.ndimage as ndimage
            logo_im = ndimage.imread(run_args['logo_file'])[:,:,0] # 1st channel
            logo_im = logo_im.astype(float)/255.0
            logo_im = ndimage.interpolation.zoom(logo_im, zoom=run_args['logo_resize'])
            
            # Expand size to match the height map
            hp, wp = logo_im.shape
            hi, wi = height_map.shape
            logo_im = np.pad(logo_im, ((0, hi-hp), (0, wi-wp)), mode='edge')
            
            # Adjust height map
            height_map = height_map*logo_im + (1 - logo_im)*0.5*run_args['stl_pedestal']
            
        
        # Kludge: We clip a corner to enforce a definition of the 'bottom' (minimum height)
        height_map[0,0] = 0
        
        if run_args['verbosity']>=4:
            print('        height_map (STL scale final) from {:.2f} to {:.2f} ({:.2f} ± {:.2f})'.format(np.min(height_map), np.max(height_map), np.average(height_map), np.std(height_map)))
        
        
        from stl_tools import numpy2stl
        numpy2stl( height_map, outfile, scale=1.0, mask_val=0, solid=True)

        
        
        return results            
        






class metadata_extract(Protocol):
    
    def __init__(self, name='metadata_extract', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        
        patterns = [
                    ['theta', '.+_th(\d+\.\d+)_.+'] ,
                    ['x_position', '.+_x(-?\d+\.\d+)_.+'] ,
                    ['y_position', '.+_y(-?\d+\.\d+)_.+'] ,
                    ['annealing_temperature', '.+_T(\d+\.\d\d\d)C_.+'] ,
                    ['annealing_time', '.+_(\d+\.\d)s_T.+'] ,
                    ['exposure_time', '.+_(\d+\.\d+)s_\d+_[a-zA-Z]axs.+'] ,
                    ['sequence_ID', '.+_(\d+)_[a-zA-Z]axs.+'] ,
                    ]            
            
        self.run_args = {
            'patterns' : patterns,
            }
        self.run_args.update(kwargs)

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        infile = data.infile
        f = tools.Filename(infile)
        filepath, filename, filebase, ext = f.split()
        
        results['infile'] = infile
        results['filepath'] = filepath
        results['filename'] = filename
        results['filebase'] = filebase
        results['fileext'] = ext
        
        results['sample_name'] = data.name
        results['file_ctime'] = os.path.getctime(infile)
        results['file_modification_time'] = os.path.getmtime(infile)
        results['file_access_time'] = os.path.getatime(infile)
        results['file_size'] = os.path.getsize(infile)
        
        patterns = run_args['patterns']
        
        for pattern_name, pattern_string in patterns:
            pattern = re.compile(pattern_string)
            m = pattern.match(filename)
            if m:
                if run_args['verbosity']>=5:
                    print('  matched: {} = {}'.format(pattern_name, m.groups()[0]))
                results[pattern_name] = float(m.groups()[0])
        
        
        outfile = self.get_outfile(data.name, output_dir, ext='.npy')
        np.save(outfile, results)
        
        
        return results



class databroker_extract(Protocol):

    def __init__(self, name='databroker_extract', **kwargs):

        self.name = self.__class__.__name__ if name is None else name

        self.default_ext = '.npy'
        self.run_args = {
            'constraints': None, # Set of md constraints to help find the files of interest
            'timestamp': True, # Restrict db searching based on file timestamp
            'timestamp_window': 3600, # We provide a finite window since the db entry is at the start of the measure, while the file timestamp is at the end
            'recent_days': None, # Restrict db searching to recent history
            'section': 'start', # Put results from a particular document section into the top-level
            }
        self.run_args.update(kwargs)


    @run_default
    def run(self, data, output_dir, **run_args):

        results = {}


        db = self._processor.db # databroker

        results['filename_resolved'] = Path(data.infile).resolve()
        results['filename'] = Path(data.infile).stem

        # This makes assumptions about what the filename looks like:
        filename = results['filename'][:-5] # The [:-5] removes the trailing "_saxs"
        scan_id = int(filename.split('_')[-1])


        # We will do a search like:
        start_time = time.time()
        #headers = db(filename=filename, scan_id=scan_id)

        # There are various ways to define the search
        constraints = run_args['constraints'] if run_args['constraints'] is not None else {}
        constraints['filename'] = filename
        constraints['scan_id'] = scan_id
        #constraints['measure_type'] = 'measure'
        if run_args['timestamp']:
            constraints['since'] = timestamp(results['filename_resolved'])-run_args['timestamp_window']
            constraints['until'] = timestamp(results['filename_resolved'])+run_args['timestamp_window']
        if run_args['recent_days'] is not None:
            constraints['since'] = time.time() - run_args['recent_days']*24*60*60


        # Search for this filename in databroker
        headers = db(**constraints)


        headers = [h for h in headers]
        if run_args['verbosity']>=1 and len(headers)!=1:
            print('  WARNING: databroker_extract got {} matches for filename: {}'.format(len(headers), filename))
        header = headers[0]

        if run_args['verbosity']>=4:
            print('        databroker lookup took {:.2f}s'.format(time.time()-start_time))
        results['uid'] = header['start']['uid']
        if run_args['section'] is not None and run_args['section'] in header:
            # Put the results of this section into the top-level
            results.update(header[run_args['section']])

        # Add all db metadata to the results
        for section_name in header.keys():
            results[section_name] = header[section_name]

        if run_args['verbosity']>=10:
            print_n(results)

        return results










# DEPRECATED
# Protocols that operate on multiple files
# These methods are being moved to a separate file (Multiple.py)
class _deprecated_sum_images(Protocol):
    
    def __init__(self, name='sum_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : None,
                        'blur' : None,
                        'resize' : None,
                        'ztrim' : [0.05, 0.005],
                        'pattern_re' : '^.+\/([a-zA-Z0-9_]+_)(\d+)(\.+)$',
                        'file_extension' : 'sum.npy',
                        'processor' : None
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['processor'] is None:
            print('ERROR: {} needs a reference to the Processor (so that it can load files).'.format(self.name))
            return {}
        
        data = self.transform(data, **run_args)
        
        infiles = run_args['infiles']
        import re
        pattern_re = re.compile(run_args['pattern_re'])
        
        # Find data among the infiles
        for infile in infiles:
            if data.name in infile:
                mainfile = infile
                break
        m = pattern_re.match(mainfile)
        if m:
            search_base = m.groups()[0]
            search_ext = m.groups()[1]
            if run_args['verbosity']>=5:
                print('    base = "{}"; ext = "{}"'.format(search_base, search_ext))
        else:
            print('ERROR: No match for {}'.format(infile))
            return results
        
        
        outfile = self.get_outfile(search_base, output_dir, ext=run_args['file_extension'])
        
        if (not run_args['force']) and os.path.isfile(outfile):
            print(' Skipping (internal check) {} for {}'.format(self.name, data.name))
            return results
        
        # Find all files that match
        for infile in infiles:
            if infile==mainfile:
                pass
            else:
                m = pattern_re.match(infile)
                if m:
                    add_base = m.groups()[0]
                    add_ext = m.groups()[1]
                    if run_args['verbosity']>5:
                        print('    base = "{}"; ext = "{}"'.format(add_base, add_ext))
                        
                    # Add this new image to the data...
                    newdata = run_args['processor'].load(infile, **run_args['load_args'])
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
    
class _deprecated_merge_images_tiling(Protocol):
    
    def __init__(self, name='merge_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'pattern_re' : '^.+\/([a-zA-Z0-9_]+_)(\d+)(\.+)$',
                        'file_extension' : 'merged.npy',
                        'processor' : None,
                        'normalizations' : None,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['processor'] is None:
            print('ERROR: {} needs a reference to the Processor (so that it can load files).'.format(self.name))
            return {}
        
        #data = self.transform(data, **run_args)
        
        infiles = run_args['infiles']
        import re
        pattern_re = re.compile(run_args['pattern_re'])
        
        # Find data among the infiles
        for infile in infiles:
            if data.name in infile:
                mainfile = infile
                break
        m = pattern_re.match(mainfile)
        if m:
            search_base = m.groups()[0]
            search_ext = m.groups()[1]
            if run_args['verbosity']>=5:
                print('    base = "{}"; ext = "{}"'.format(search_base, search_ext))
        else:
            print('ERROR: No match for {}'.format(infile))
            return results
        
        
        outfile = self.get_outfile(search_base, output_dir, ext=run_args['file_extension'])
        
        if (not run_args['force']) and os.path.isfile(outfile):
            print(' Skipping (internal check) {} for {}'.format(self.name, data.name))
            return results
        
        
        # Prepare a large-area q matrix
        q_range = run_args['q_range']
        dq = run_args['dq']
        qxs = np.arange(q_range[0], q_range[1], dq)
        qzs = np.arange(q_range[2], q_range[3], dq)
        QXs, QZs = np.meshgrid(qxs, qzs)
        
        Intensity_map = np.zeros( (len(qzs), len(qxs)) )
        count_map = np.zeros( (len(qzs), len(qxs)) )
        
        if run_args['verbosity']>=5:
            print('      Expecting array size num_qx = {}, num_qz = {}'.format(len(qxs), len(qzs)))
            print('      Coordinate matrices sized {}'.format(QXs.shape))
            #print('      Coordinate matrices sized {}'.format(QZs.shape))
            print('      Data matrices sized {}'.format(Intensity_map.shape))
            #print('      Data matrices sized {}'.format(count_map.shape))

        # Find all files that match
        for i, (infile, (x0,y0)) in enumerate(zip(infiles, run_args['beam_positions'])):
            
            m = pattern_re.match(infile)
            if m:
                add_base = m.groups()[0]
                add_ext = m.groups()[1]
                if run_args['verbosity']>5:
                    print('    base = "{}"; ext = "{}"'.format(add_base, add_ext))
                    
                # Add this new image to the data...
                newdata = run_args['processor'].load(infile, **run_args['load_args'])

                if run_args['verbosity']>=4:
                    print('    Updating calibration to (x0, y0) = ({:.3f}, {:.3f}) for {}'.format(x0, y0, infile))
                
                newdata.calibration.clear_maps()
                newdata.calibration.x0 = x0
                newdata.calibration.y0 = y0
                
                remesh_data, num_per_pixel = newdata.remesh_q_bin_explicit(qx_min=q_range[0], qx_max=q_range[1], num_qx=len(qxs), qz_min=q_range[2], qz_max=q_range[3], num_qz=len(qzs), **run_args)
                #newdata = self.transform(newdata, **run_args)
                #data.data += newdata.data

                if run_args['verbosity']>=5:
                    print('      remesh_data matrix sized {}'.format(remesh_data.shape))
                
                if run_args['normalizations'] is not None:
                    remesh_data *= run_args['normalizations'][i]
                
                Intensity_map += remesh_data
                count_map += num_per_pixel
                
            
            
        Intensity_map = np.nan_to_num( Intensity_map/count_map )
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'sum of multiple images (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        np.save(outfile, Intensity_map)    
            
            
        if True:
            import pickle
            outfile = self.get_outfile(search_base, output_dir, ext='merged.pkl')
            
            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'qx and qz axes of output data (Python pickle format)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                'comment' : "reload using: qxs, qzs = pickle.load(open('{}merged.pkl','rb'))".format(search_base)
                } , )
            with open(outfile, 'wb') as fout:
                pickle.dump([qxs, qzs], fout)
                
            # Reload using:
            #qxs, qzs = pickle.load(open('infile.pkl','rb'))
            
            
        
        return results
    
class _deprecated_merge_images_gonio_phi(Protocol):
    
    def __init__(self, name='merge_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'pattern_re' : '^.+\/([a-zA-Z0-9_]+_)(\d+)(\.+)$',
                        'file_extension' : 'merged.npy',
                        'processor' : None
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['processor'] is None:
            print('ERROR: {} needs a reference to the Processor (so that it can load files).'.format(self.name))
            return {}
        
        #data = self.transform(data, **run_args)
        
        infiles = run_args['infiles']
        import re
        pattern_re = re.compile(run_args['pattern_re'])
        
        # Find data among the infiles
        for infile in infiles:
            if data.name in infile:
                mainfile = infile
                break
        m = pattern_re.match(mainfile)
        if m:
            search_base = m.groups()[0]
            search_ext = m.groups()[1]
            if run_args['verbosity']>=5:
                print('    base = "{}"; ext = "{}"'.format(search_base, search_ext))
        else:
            print('ERROR: No match for {}'.format(infile))
            return results
        
        
        outfile = self.get_outfile(search_base, output_dir, ext=run_args['file_extension'])
        
        if (not run_args['force']) and os.path.isfile(outfile):
            print(' Skipping (internal check) {} for {}'.format(self.name, data.name))
            return results
        
        
        # Prepare a large-area q matrix
        q_range = run_args['q_range']
        dq = run_args['dq']
        qxs = np.arange(q_range[0], q_range[1], dq)
        qzs = np.arange(q_range[2], q_range[3], dq)
        QXs, QZs = np.meshgrid(qxs, qzs)
        
        Intensity_map = np.zeros( (len(qzs), len(qxs)) )
        count_map = np.zeros( (len(qzs), len(qxs)) )
        
        if run_args['verbosity']>=5:
            print('      Expecting array size num_qx = {}, num_qz = {}'.format(len(qxs), len(qzs)))
            print('      Coordinate matrices sized {}'.format(QXs.shape))
            #print('      Coordinate matrices sized {}'.format(QZs.shape))
            print('      Data matrices sized {}'.format(Intensity_map.shape))
            #print('      Data matrices sized {}'.format(count_map.shape))

        # Find all files that match
        for infile, phi in zip(infiles, run_args['phis']):
            m = pattern_re.match(infile)
            if m:
                add_base = m.groups()[0]
                add_ext = m.groups()[1]
                if run_args['verbosity']>5:
                    print('    base = "{}"; ext = "{}"'.format(add_base, add_ext))
                    
                # Add this new image to the data...
                newdata = run_args['processor'].load(infile, **run_args['load_args'])

                if run_args['verbosity']>=4:
                    print('    Updating calibration to phi = {} for {}'.format(phi, infile))
                
                newdata.calibration.clear_maps()
                newdata.calibration.set_angles(det_phi_g=phi, det_theta_g=0.)
                
                remesh_data, num_per_pixel = newdata.remesh_q_bin_explicit(qx_min=q_range[0], qx_max=q_range[1], num_qx=len(qxs), qz_min=q_range[2], qz_max=q_range[3], num_qz=len(qzs), **run_args)
                #newdata = self.transform(newdata, **run_args)
                #data.data += newdata.data

                if run_args['verbosity']>=5:
                    print('      remesh_data matrix sized {}'.format(remesh_data.shape))
                
                Intensity_map += remesh_data
                count_map += num_per_pixel
                
            
            
        Intensity_map = np.nan_to_num( Intensity_map/count_map )
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'sum of multiple images (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        np.save(outfile, Intensity_map)    
            
            
        if True:
            import pickle
            outfile = self.get_outfile(search_base, output_dir, ext='merged.pkl')
            
            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'qx and qz axes of output data (Python pickle format)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                'comment' : "reload using: qxs, qzs = pickle.load(open('{}merged.pkl','rb'))".format(search_base)
                } , )
            with open(outfile, 'wb') as fout:
                pickle.dump([qxs, qzs], fout)
                
            # Reload using:
            #qxs, qzs = pickle.load(open('infile.pkl','rb'))
            
            
        
        return results
# DEPRECATED
