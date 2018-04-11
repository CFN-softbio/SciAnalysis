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


        return data
        
        

class thumbnails(Protocol):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : None,
                        'blur' : 2.0,
                        'resize' : 0.2,
                        'ztrim' : [0.05, 0.005]
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
        
        data.set_z_display([None, None, 'gamma', 0.3])
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)
        
        
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
            line.plot(save=outfile, show=False, **run_args)
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
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        #line = data.linecut_angle(q0=run_args['q0'], dq=run_args['dq'])
        line = data.linecut_angle(**run_args)
        
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
        q_data.x_rlabel = '$q_r \, (\AA^{-1})$'

        q_data.plot(outfile, plot_buffers=[0.30,0.05,0.25,0.05], **run_args)
        
        
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
        
        if True:
            # Save Data2DQPhi() object
            import pickle
            outfile = self.get_outfile(data.name, output_dir, ext='.pkl')
            with open(outfile, 'wb') as fout:
                out_data = q_data.data, q_data.x_axis, q_data.y_axis
                pickle.dump(out_data, fout)
            
        
        
        
        return results
    
    
class sum_images(Protocol):
    
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
            print('ERROR: {} needs a reference to the Processor (so that it can load files.'.format(self.name))
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