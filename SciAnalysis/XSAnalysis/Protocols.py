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
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################

from .Data import *
from ..tools import *


class ProcessorXS(Processor):

    
    def load(self, infile, **kwargs):

        #calibration = kwargs['calibration'] if 'calibration' in kwargs else None
        #mask = kwargs['mask'] if 'mask' in kwargs else None
        #data = Data2DScattering(infile, calibration=calibration, mask=mask)
        
        data = Data2DScattering(infile, **kwargs)
        data.infile = infile
        
        data.threshold_pixels(4294967295-1) # Eiger inter-module gaps
        
        if 'dezing' in kwargs and kwargs['dezing']:
            data.dezinger()
        
        if 'flip' in kwargs and kwargs['flip']:
            #if flip: self.im = self.im.transpose(Image.ROTATE_90).transpose(Image.FLIP_LEFT_RIGHT)
            data.data = np.rot90(data.data) # rotate CCW
            data.data = np.fliplr(data.data) # Flip left/right

        if 'rotCCW' in kwargs and kwargs['rotCCW']:
            data.data = np.rot90(data.data) # rotate CCW

        if 'rot180' in kwargs and kwargs['rot180']:
            data.data = np.flipud(data.data) # Flip up/down
            data.data = np.fliplr(data.data) # Flip left/right


        if data.mask is not None:
            data.data *= data.mask.data


        return data
        
        

class thumbnails(Protocol):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : None,
                        'shift_crop_up' : 0.0,
                        'blur' : 2.0,
                        'resize' : 0.2,
                        'ztrim' : [0.05, 0.005]
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'], shift_crop_up=run_args['shift_crop_up'])
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
        self.run_args = {}
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)
        
        
        line = data.circular_average_q_bin(error=True)
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        
        try:
            line.plot(save=outfile, **run_args)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        # TODO: Fit 1D data
        
        return results
                
                
                
class circular_average_q2I(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {}
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.circular_average_q_bin(error=True)
        
        line.y *= np.square(line.x)
        line.y_label = 'q^2*I(q)'
        line.y_rlabel = '$q^2 I(q) \, (\AA^{-2} \mathrm{counts/pixel})$'
        
        
        outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
        line.plot(save=outfile, show=False, **run_args)
        
        outfile = self.get_outfile(data.name, output_dir, ext='_q2I.dat')
        line.save_data(outfile)        
        
        # TODO: Fit 1D data
        
        return results
                       

    def output_exists(self, name, output_dir):

        if 'file_extension' in self.run_args:
            ext = '_q2I{}'.format(self.run_args['file_extension'])
        else:
            ext = '_q2I{}'.format(self.default_ext)

        outfile = self.get_outfile(name, output_dir, ext=ext)
        return os.path.isfile(outfile)



         
class sector_average(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
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
        
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        outfile = self.get_outfile(data.name, output_dir)
        
        try:
            line.plot(save=outfile, show=False, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        return results
                                
                
                
                
                
class linecut_angle(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [-180, 180, 0, None]
                         }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir,ztrim=[0.02, 0.01], **run_args):
        
        results = {}
        
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_angle(**run_args)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True, ztrim=ztrim)
        
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        line.plot(save=outfile, **run_args)

        #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
        #line.plot_polar(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        return results
                                
                
                


class linecut_qr(Protocol):

    def __init__(self, name='linecut_qr', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [None, None, 0, None]
                         }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qr(**run_args)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        line.plot(save=outfile, **run_args)

        #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
        #line.plot_polar(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        return results
                                
                                
class linecut_qz(Protocol):

    def __init__(self, name='linecut_qz', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [None, None, 0, None]
                         }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qz(**run_args)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        line.plot(save=outfile, **run_args)

        #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
        #line.plot_polar(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        return results                


class linecut_q(Protocol):

    def __init__(self, name='linecut_q', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [None, None, 0, None]
                         }
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_q(**run_args)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        line.plot(save=outfile, **run_args)

        #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
        #line.plot_polar(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        return results






class linecut_qr_fit(linecut_qr):
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
        
        if 'show_region' in run_args and run_args['show_region']:
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
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)

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
            params.add('m', value=m, min=0, max=abs(m)*+4, vary=False) # Slope must be positive
            params.add('b', value=b, min=0, max=np.max(line.y)*100, vary=False)
            
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
    





                
                
class calibration_check(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'dq': 0.01}
        self.run_args.update(kwargs)
    
        
    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        outfile = self.get_outfile('{}_full'.format(data.name), output_dir)
        data.plot_image(outfile, **run_args)

        
        #data.blur(2.0)
        
        dq = run_args['dq']

        if 'AgBH' in run_args and run_args['AgBH']:
            q0 = 0.1076 # A^-1
            
            for i in range(11):
                data.overlay_ring(q0*(i+1), q0*(i+1)*dq)
                
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
                                       
                
                
                
                
                
                
                
# Work in progress
################################################################################                    
                
class fit_calibration(Protocol):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = { 'material' : 'AgBH01' }
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
        q_data.plot(outfile, plot_buffers=[0.30,0.05,0.25,0.05], **run_args)


        if 'plot_buffers' not in run_args:
            run_args['plot_buffers'] = [0.30,0.05,0.25,0.05]
        q_data.plot(outfile, **run_args)
        
        if 'save_data' in run_args and run_args['save_data']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)
        
        
        return results
    
    
class qr_image(Protocol):
    
    def __init__(self, name='qr_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'blur' : None,
                        'ztrim' : [0.05, 0.005],
                        'method' : 'nearest',
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, image_output=False, **run_args):
        
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

        image_outfile = self.get_outfile(data.name, output_dir)
        #matrix_outfile = '/GPFS/xf11bm/data/2018_3/BOcko3/waxs/analysis/qr_image/test'

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
        q_data.x_rlabel = '$q_r \, (\AA^{-1})$'

        if 'plot_buffers' not in run_args:
            run_args['plot_buffers'] = [0.30,0.05,0.25,0.05]
        q_data.plot(outfile, **run_args)
        
        if 'save_data' in run_args and run_args['save_data']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)
        
        if image_output==True:
            q_data.save_image(image_outfile)
        
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
        q_data.plot(outfile, **run_args)
        
        
        if 'save_data' in run_args and run_args['save_data']:
            outfile = self.get_outfile(data.name, output_dir, ext='.npz')
            q_data.save_data(outfile)
        
        
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
                        'save_data_pickle': True,
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
        q_data.plot(outfile, plot_buffers=[0.20,0.05,0.20,0.05], **run_args)
        
        if 'save_data_pickle' in run_args and run_args['save_data_pickle']:
            # Save Data2DQPhi() object
            import pickle
            outfile = self.get_outfile(data.name, output_dir, ext='.pkl')
            with open(outfile, 'wb') as fout:
                out_data = q_data.data, q_data.x_axis, q_data.y_axis
                pickle.dump(out_data, fout)
            
        
        return results











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
    
    
    
class merge_images_tiling(Protocol):
    
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
            
    
class merge_images_gonio_phi(Protocol):
    
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
        
class test_tiling(Protocol):
    
    def __init__(self, name='merge_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.tiff'
        self.run_args = {
                        #'pattern_re' : '^.+\/([a-zA-Z0-9_]+_)(\d+)(\.+)$',
                        'file_extension' : 'merged.npy',
                        'processor' : None,
                        'normalizations' : None,
                        'mask1':None,
                        'mask2':None
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, infiles, output_dir,  mask1=None, mask2=None, **run_args):
        
        results = {}
        print('1st check\n')
       
        #outfile = self.preliminary(infiles, output_dir, **run_args)
        #if outfile is None:
            #return {}

        print('2st check\n')
        
        processor = run_args['processor']
        load_args = processor.load_args
        
        #print('processer = {}\n'.format(processer.name))
        print(infiles.name)
        print('=====')
        # Load first image
        #data = processor.load(infiles[0], **load_args, mask=mask1)
        #data = processor.load(infiles)
        #data = self.transform(data, **run_args)
        #data = Data2DScattering(infiles, mask=mask1)
        #data = Data2DScattering(infiles)
        
        # Iterate through remaining images
        #for infile in infiles[1:]:
            ## Add this new image to the data
            #newdata = run_args['processor'].load(infile, **load_args, mask=mask_pos2)
            #newdata = self.transform(newdata, **run_args)
            #data.data += newdata.data
        ## Iterate through remaining images
        #for infile in infiles[1:]:
            ## Add this new image to the data
            #newdata = run_args['processor'].load(infile, **load_args)
            #newdata = self.transform(newdata, **run_args)
            #data.data += newdata.data
                
            
        
        #results['files_saved'] = [
            #{ 'filename': '{}'.format(outfile) ,
             #'description' : 'sum of multiple images (npy format)' ,
             #'type' : 'data' # 'data', 'plot'
            #} ,
            #]
            
        #np.save(outfile, data.data)
        
        
        return results
        