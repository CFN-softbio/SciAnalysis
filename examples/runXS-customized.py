#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols




# Define some custom analysis routines
########################################

class circular_average_q2I_fit(Protocols.circular_average_q2I):

    def __init__(self, name='circular_average_q2I', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {}
        self.run_args.update(kwargs)
    
        
    @Protocols.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.circular_average_q_bin(error=True)
        
        # Convert to q2I
        line.y *= np.square(line.x)
        line.y_label = 'q^2*I(q)'
        line.y_rlabel = '$q^2 I(q) \, (\AA^{-2} \mathrm{counts/pixel})$'

        outfile = self.get_outfile(data.name, output_dir, ext='_q2I.dat')
        line.save_data(outfile)        

        # Fit
        lm_result, fit_line, fit_line_extended = self.fit_peaks(line, fit_range=[0, 0.05], **run_args)
        
        fit_name = 'fit_peaks'
        prefactor_total = 0
        for param_name, param in lm_result.params.items():
            results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            if 'prefactor' in param_name:
                prefactor_total += np.abs(param.value)
            
        results['{}_prefactor_total'.format(fit_name)] = prefactor_total
        results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
        
        
        lines = DataLines([line, fit_line, fit_line_extended])
        lines.copy_labels(line)
        
        outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
        lines.plot(save=outfile, show=False, **run_args)
        
        
        return results
    
    
        
    def fit_peaks(self, line, num_curves=2, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)
        
        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            '''Gaussians with constant background.'''
            m = v['b']
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('b', value=np.min(line.y), min=0)
        
        #xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        for i in range(num_curves):
            #xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            xpos, ypos = line.target_x(xpeak*(i+1))
            params.add('prefactor{:d}'.format(i+1), value=ypos, min=0)
            params.add('x_center{:d}'.format(i+1), value=xpos, min=np.min(line.x), max=np.max(line.x))
            params.add('sigma{:d}'.format(i+1), value=xpos*0.25, min=0)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':3.0})        

        return lm_result, fit_line, fit_line_extended



class linecut_angle_fit(Protocols.linecut_angle):

    def __init__(self, name='linecut_angle', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [-180, 180, 0, None]
                         }
        self.run_args.update(kwargs)
        
        self.load_fit_libraries()
        
        
        
    def load_fit_libraries(self):
        
        load_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/order_parameter/simulation_1D/eta_library/'
        
        filename_re = re.compile('^.+I(-?\d+\.\d+)_S(-?\d+\.\d+)\.dat$')
        
        infiles = glob.glob('{}/*.dat'.format(load_dir))
        infiles.sort()
        
        self.fit_library_eta = []
        for infile in infiles:
            m = filename_re.match(infile)
            if m:
                value = float(m.groups()[0])
                S = float(m.groups()[1])
                
                data = np.loadtxt(infile, skiprows=1, comments='#')
                self.fit_library_eta.append( [value, data, S] )
                
                
                
        load_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/order_parameter/simulation_1D/MS_library/'
        
        filename_re = re.compile('^.+I(-?\d+\.\d+)_S(-?\d+\.\d+)\.dat$')
        
        infiles = glob.glob('{}/*.dat'.format(load_dir))
        infiles.sort()
        
        self.fit_library_MaierSaupe = []
        for infile in infiles:
            m = filename_re.match(infile)
            if m:
                value = float(m.groups()[0])
                S = float(m.groups()[1])
                
                data = np.loadtxt(infile, skiprows=1, comments='#')
                self.fit_library_MaierSaupe.append( [value, data, S] )                
    
        
    @Protocols.run_default
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
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        # Basic output
        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        
        # Fine-tune data
        line.x = line.x[1:]
        line.y = line.y[1:]
        line.smooth(2.0)
        

        
        
        # Fit
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                v_spacing = 0.1
                
                xi, xf, yi, yf = self.ax.axis()

                yp = yf
                s = '$\eta = {:.2f}$'.format(self.results['fit_eta_eta'])
                self.ax.text(xi, yp, s, size=30, color='b', verticalalignment='top', horizontalalignment='left')

                yp -= (yf-yi)*v_spacing
                s = '$S_{{\eta}} = {:.3f}$'.format(self.results['fit_library_eta_S'])
                self.ax.text(xi, yp, s, size=30, color='b', verticalalignment='top', horizontalalignment='left')


                yp = yf
                s = '$m = {:.2f}$'.format(self.results['fit_MaierSaupe_m'])
                self.ax.text(xf, yp, s, size=30, color='purple', verticalalignment='top', horizontalalignment='right')

                yp -= (yf-yi)*v_spacing
                s = '$S_{{m}} = {:.3f}$'.format(self.results['fit_library_MaierSaupe_S'])
                self.ax.text(xf, yp, s, size=30, color='purple', verticalalignment='top', horizontalalignment='right')
                
                
                self.ax.set_xticks([-180, -90, 0, +90, +180])
                
                
        lines = DataLines_current([]) # Not sure why the lines=[] argument needs to be specified here...
        lines.add_line(line)
        lines.copy_labels(line)
        lines.results = {}
            
            
        color_list = ['b', 'purple', 'r', 'green', 'orange',]
        for i, fit_name in enumerate(['fit_eta', 'fit_MaierSaupe']):
            
            lm_result, fit_line, fit_line_e = getattr(self, fit_name)(line, **run_args)
            fit_line_e.plot_args['color'] = color_list[i%len(color_list)]
            lines.add_line(fit_line_e)
        
            #prefactor_total = 0
            for param_name, param in lm_result.params.items():
                results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                lines.results['{}_{}'.format(fit_name, param_name)] = param.value
                #if 'prefactor' in param_name:
                    #prefactor_total += np.abs(param.value)
                
            #results['{}_prefactor_total'.format(fit_name)] = prefactor_total
            results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
            
            
            
        # Special fit
        fit_result, fit_line = self.fit_library(line, library='eta', **run_args)
        fit_line.plot_args['color'] = 'blue'
        fit_line.plot_args['linewidth'] = 1.0
        lines.results['fit_library_eta_S'] = fit_result['S']
        lines.add_line(fit_line)
        results.update(self.prepend_keys(fit_result, 'fit_library_eta_'))
        
        
        fit_result, fit_line = self.fit_library(line, library='MaierSaupe', **run_args)
        fit_line.plot_args['color'] = 'purple'
        fit_line.plot_args['linewidth'] = 1.0
        lines.results['fit_library_MaierSaupe_S'] = fit_result['S']
        results.update(self.prepend_keys(fit_result, 'fit_library_MaierSaupe_'))
        lines.add_line(fit_line)
        
        


        # Output
        if run_args['verbosity']>=1:
            outfile = self.get_outfile(data.name, output_dir)
            lines.lines.reverse()
            lines.plot(save=outfile, **run_args)

        if run_args['verbosity']>=4:
            outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
            line.plot_polar(save=outfile, **run_args)

        return results
    
    
    def fit_library(self, line, library, **run_args):
        
        import lmfit
        from scipy.interpolate import interp1d
        
        fit_library_use = getattr(self, 'fit_library_{}'.format(library))

        
        num_curves = len(fit_library_use)
        
        best_ic = -1
        best_chi_squared = 1e20
        best_y = []
        
        for ic in range(num_curves):
        
            # 'Expand' curve
            value, data, S = fit_library_use[ic]
            x = data[:,0]
            y = data[:,1]
            
            x = np.concatenate( (x, 180-x[::-1]) )
            y = np.concatenate( (y, y[::-1]) )
            
            x = np.concatenate( (-1*x[::-1], x) )
            y = np.concatenate( (y[::-1], y) )
            
            y /= np.max(y)
            
            # Register with experimental data
            f = interp1d(x, y)
            xm = line.x
            ym = f(xm)
            to_fit = DataLine(x=xm, y=ym)
            
            
            
            def model(v, x):
                m = v['prefactor']*ym
                return m
            
            def func2minimize(params, x, data):
                
                v = params.valuesdict()
                m = model(v, x)
                
                return m - data        

            params = lmfit.Parameters()
            params.add('prefactor', value=np.max(line.y), min=0)

            lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
            
            chi_squared = lm_result.chisqr/lm_result.nfree
            
            if chi_squared<best_chi_squared:
                best_ic = ic
                best_chi_squared = chi_squared
                best_y = model(lm_result.params.valuesdict(), line.x)
                
            #print('        Curve #{:d} (S = {:.3f}), chi = {:f}'.format(ic, fit_library_use[ic][2], chi_squared))
        
        
        
        fit_result = {}
        fit_line = DataLine(x=line.x, y=best_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        fit_result['S'] = fit_library_use[best_ic][2]
        fit_result['chi_squared'] = best_chi_squared
    
        return fit_result, fit_line
    
    
                        
                        

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
        
        params = lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y)*0.5, min=0)
        params.add('x_center', value=-90, min=np.min(line.x), max=np.max(line.x))
        params.add('eta', value=0.2, min=0, max=1)
        params.add('symmetry', value=2, min=0.5, max=20, vary=False)
        params.add('baseline', value=0, min=0, max=np.max(line.y), vary=False)
        
        
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
                    
           
    def fit_MaierSaupe(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['prefactor']*np.exp(v['m']*np.square(np.cos(np.radians(x-v['x_center'])))) + v['baseline']
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y)*0.1, min=0)
        params.add('x_center', value=-90, min=np.min(line.x), max=np.max(line.x))
        params.add('m', value=2.0, min=0)
        params.add('baseline', value=0, min=0, max=np.max(line.y), vary=False)
        
        
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








# qm is "q_medium", where q0 = 0.05, dq = 0.025
# This is diffuse (form factor) scattering (probably), so we need a different 'library' fit
class linecut_angle_fit_qm(linecut_angle_fit):

    def __init__(self, name='linecut_angle_qm', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [-180, 180, 0, None]
                         }
        self.run_args.update(kwargs)
        
        self.load_fit_libraries()
        

    def load_fit_libraries(self):
        
        load_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/order_parameter/simulation_1D_qm/eta_library/'
        
        filename_re = re.compile('^.+I(-?\d+\.\d+)_S(-?\d+\.\d+)\.dat$')
        
        infiles = glob.glob('{}/*.dat'.format(load_dir))
        infiles.sort()
        
        self.fit_library_eta = []
        for infile in infiles:
            m = filename_re.match(infile)
            if m:
                value = float(m.groups()[0])
                S = float(m.groups()[1])
                
                data = np.loadtxt(infile, skiprows=1, comments='#')
                self.fit_library_eta.append( [value, data, S] )
                
                
                
        load_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/order_parameter/simulation_1D_qm/MS_library/'
        
        filename_re = re.compile('^.+I(-?\d+\.\d+)_S(-?\d+\.\d+)\.dat$')
        
        infiles = glob.glob('{}/*.dat'.format(load_dir))
        infiles.sort()
        
        self.fit_library_MaierSaupe = []
        for infile in infiles:
            m = filename_re.match(infile)
            if m:
                value = float(m.groups()[0])
                S = float(m.groups()[1])
                
                data = np.loadtxt(infile, skiprows=1, comments='#')
                self.fit_library_MaierSaupe.append( [value, data, S] )                
    


class sector_average_fit(Protocols.sector_average):

    def __init__(self, name='sector_average_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'error' : True, 
                        'show_region' : False,
                        }
        self.run_args.update(kwargs)
    
        
    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)
        
        
        line = data.sector_average_q_bin(**run_args)
        #line.smooth(2.0, bins=10)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)


        # Fit data
        lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)

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
        xi = (2.*np.pi/np.sqrt(2.*np.pi))/results['{}_sigma1'.format(fit_name)]['value']
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

        outfile = self.get_outfile(data.name, output_dir)
        
        try:
            lines.plot(save=outfile, show=False, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        print(results)
        
        return results


    def fit_peaks(self, line, num_curves=1, **run_args):
        
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
        params.add('m', value=0)
        params.add('b', value=np.min(line.y), min=0, max=np.max(line.y))
        
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        for i in range(num_curves):
            #xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            #xpos, ypos = line.target_x(xpeak*(i+1))
            xpos, ypos = xpeak, ypeak
            
            params.add('prefactor{:d}'.format(i+1), value=ypos-np.min(line.y), min=0, max=np.max(line.y)*10)
            params.add('x_center{:d}'.format(i+1), value=xpos, min=np.min(line.x), max=np.max(line.x))
            params.add('sigma{:d}'.format(i+1), value=xspan*0.20, min=0, max=xspan*0.50)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'alpha':0.5, 'marker':None, 'linewidth':2.0})        

        return lm_result, fit_line, fit_line_extended
                                
         

# Experimental parameters
########################################

calibration = Calibration(wavelength_A=1.5418) # 8.04 keV (Cu K-alpha e.g. Nanostar)
calibration.set_image_size(2048)
calibration.set_pixel_size(width_mm=140.0)
calibration.set_distance(1.150)
calibration.set_beam_position(1016.0, 1021.5)

mask_dir = '/home/kyager/current/code/SciAnalysis/main/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Bruker_Nanostar_SAXS-2015Oct-mask.png')



# Files to analyze
########################################

root_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/2015_10Oct_18-SAXS_all/'

source_dir = os.path.join(root_dir, 'CNC_825_1p2T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_825_0p56T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_825_0T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_413_1p2T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_165_1p2T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_165_0T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_1p65_0p56T_1_kin/')


#source_dir = os.path.join(root_dir, 'Hydrogels/', 'Static tests & cals/')

output_dir = os.path.join(source_dir, 'SciAnalysis/')

infiles = glob.glob(os.path.join(source_dir, '*.dat'))
#infiles = glob.glob(os.path.join(source_dir, 'AgBH_kineticCal.dat'))
infiles.sort()
#infiles = [infiles[0]]
#infiles = [infiles[-1]]


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)

protocols = [ Protocols.calibration_check(show=False, AgBH=True, q0=0.015, num_rings=2) ]

#protocols = [ Protocols.thumbnails(crop=0.8)] # mini
#protocols = [ Protocols.thumbnails(crop=None, resize=1.0, blur=None, file_extension='.png') ] # full-size
#protocols = [ Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge, ztrim=[0, 0.001], file_extension='.png') ] # custom

#protocols = [ Protocols.q_image(q_max=0.14, blur=2.0, bins_relative=0.25, xticks=[-.1, 0, .1], ztrim=[0.01, 0.001])]

#protocols = [ Protocols.circular_average(ylog=True, plot_range=[0, 0.2, None, None]) ]
#protocols = [ Protocols.circular_average_q2I(plot_range=[0, 0.2, 0, None]) ]
#protocols = [ circular_average_q2I_fit(plot_range=[0, 0.2, 0, None]) ] # local protocol

#protocols = [ Protocols.linecut_angle(q0=0.01687, dq=0.00455*1.5, show_region=False) ]
#protocols = [ linecut_angle_fit(dq=0.00455*1.5) ]


protocols = [
    #Protocols.thumbnails(crop=0.8) ,
    #circular_average_q2I_fit(plot_range=[0, 0.2, 0, None]) ,
    #linecut_angle_fit(dq=0.00455*1.5) , # for q0
    linecut_angle_fit_qm(q0=0.05, dq=0.025) , # for qm
    ]
    



# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=True)



