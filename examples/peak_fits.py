#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis2/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis2018C1/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols



class linecut_qr_fit(Protocols.linecut_qr):
    

    def __init__(self, name='linecut_qr_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [None, None, 0, None]
                         }
        self.run_args.update(kwargs)    
    
    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.linecut_qr(**run_args)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        #line.plot(save=outfile, **run_args)

        #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
        #line.plot_polar(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        
        
        # Fit data
        if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
            line.trim(run_args['fit_range'][0], run_args['fit_range'][1])
        
        lm_result, fit_line, fit_line_extended = self._fit_peaks_special(line, **run_args)
        
        
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
        
        try:
            #lines.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
            lines.plot(save=outfile, **run_args)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        
        return results


    

    def _fit_peaks(self, line, num_curves=1, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)
        
        line_full = line
        #if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
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
                
                
    def _fit_peaks_special(self, line, num_curves=1, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)
        
        line_full = line
        #if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
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
        params.add('m', value=m, min=abs(m)*-10, max=0) # Slope must be negative
        params.add('b', value=b, min=0)
        
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        for i in range(num_curves):
            #xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            #xpos, ypos = line.target_x(xpeak*(i+1))
            xpos, ypos = xpeak, ypeak
            
            params.add('prefactor{:d}'.format(i+1), value=np.max(line.y), min=0, max=np.max(line.y)*10)
            params.add('x_center{:d}'.format(i+1), value=xpos, min=np.min(line.x), max=np.max(line.x))
            #params.add('x_center{:d}'.format(i+1), value=-0.009, min=np.min(line.x), max=np.max(line.x))
            params.add('sigma{:d}'.format(i+1), value=0.005, min=0, max=xspan*0.5)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        # https://lmfit.github.io/lmfit-py/fitting.html
        #lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y), method='nelder')
        
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



class linecut_angle_fit(Protocols.linecut_angle):

    def __init__(self, name='linecut_angle_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [-90, 90, None, None]
                         }
        self.run_args.update(kwargs)
        
        
        
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
        
        
        
        # Remove the data near chi=0 since this is artificially high (specular ridge)
        line.kill_x(0, 3)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        # Basic output
        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        #line.save_data(outfile)
        
        
        # Fine-tune data
        line.x = line.x[1:]
        line.y = line.y[1:]
        line.smooth(2.0)
        
        do_FWHM = True
        do_fits = False
        
        
        # Fit
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                v_spacing = 0.1
                
                xi, xf, yi, yf = self.ax.axis()

                #yi = min( self.lines[-1].y )*0.5
                #yf = max( self.lines[-1].y )*1.5
                #self.ax.axis( [xi, xf, yi, yf] )

                if do_fits:
                    yp = yf
                    s = '$\eta = {:.2f}$'.format(self.results['fit_eta_eta'])
                    self.ax.text(xi, yp, s, size=30, color='b', verticalalignment='top', horizontalalignment='left')

                    yp = yf
                    s = '$m = {:.2f}$'.format(self.results['fit_MaierSaupe_m'])
                    self.ax.text(xf, yp, s, size=30, color='purple', verticalalignment='top', horizontalalignment='right')
                    
                    
                if do_FWHM:
                    # FWHM
                    line = self.lines[0]
                    xt, yt = line.target_y(max(line.y))
                    hm = yt*0.5
                    
                    # Right (+) side
                    line_temp = line.sub_range(xt, +60)
                    xtr, ytr = line_temp.target_y(hm)
                    
                    
                    # Left (-) side
                    line_temp = line.sub_range(-60, xt)
                    xtl, ytl = line_temp.target_y(hm)
                    
                    
                    self.ax.plot([xtl, xtr], [ytl, ytr], 'o-', color='b', markersize=8, linewidth=1.0)
                    s = r'$\mathrm{{FWHM}} = {:.1f}^{{\circ}}$'.format( abs(xtr-xtl) )
                    self.ax.text(xtr, ytr, s, color='b', size=20, verticalalignment='center', horizontalalignment='left')

                
                #self.ax.set_xticks([-180, -90, 0, +90, +180])
                
                
                
        lines = DataLines_current([]) # Not sure why the lines=[] argument needs to be specified here...
        lines.add_line(line)
        lines.copy_labels(line)
        lines.results = {}
            
            
        if do_fits:
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
            
            
        # Output
        if run_args['verbosity']>=1:
            outfile = self.get_outfile(data.name, output_dir)
            lines.lines.reverse()
            lines.plot(save=outfile, **run_args)

        if run_args['verbosity']>=4:
            outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
            line.plot_polar(save=outfile, **run_args)

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
        
        params = lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y)*0.5, min=0)
        params.add('x_center', value=0, min=np.min(line.x), max=np.max(line.x), vary=True)
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
        params.add('x_center', value=0, min=np.min(line.x), max=np.max(line.x))
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



class linecut_qz_fit(Protocols.linecut_qz):

    def __init__(self, name='linecut_qz_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [None, None, 0, None],
                         'show_guides' : True,
                         'markersize' : 0,
                         'linewidth' : 1.5,
                         }
        self.run_args.update(kwargs)
    
        
    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        # Usage example:
        # linecut_qz_fit(show_region=False, show=False, qr=0.009, dq=0.0025, q_mode='qr', fit_range=fit_range, q0_expected=q0_expected, plot_range=[0, 0.08, 0, None]) ,

        
        results = {}
        
        line = data.linecut_qz(**run_args)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        
        #line.smooth(2.0, bins=10)
        
        outfile = self.get_outfile(data.name, output_dir)
        line.plot(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)

        
        if 'incident_angle' not in run_args:
            run_args['incident_angle'] = 0
            
            import re
            filename_re = re.compile('^.+_th(-?\d+\.\d+)_.+$')
            m = filename_re.match(data.name)
            if m:
                run_args['incident_angle'] = float(m.groups()[0])
                
                
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
        

        
        if 'critical_angle_film' in run_args:
            # Account for refraction distortion
            
            theta_incident_rad = np.radians(run_args['incident_angle'])
            theta_c_f_rad = np.radians(run_args['critical_angle_film'])
            #theta_c_s_rad = np.radians(run_args['critical_angle_substrate'])
            
            alpha_i_rad = np.arccos( np.cos(theta_incident_rad)/np.cos(theta_c_f_rad) )
            
            def angle_to_q(two_theta_s_rad):
                k = data.calibration.get_k()
                qz = 2*k*np.sin(two_theta_s_rad/2.0)
                return qz      
            def q_to_angle(q):
                k = data.calibration.get_k()
                two_theta_s_rad = 2.0*np.arcsin(q/(2.0*k))
                return two_theta_s_rad
            
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
                
                
                s = '$q_T = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$ \n $d_T = \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_qT'], self.results['fit_peaks_dT'])
                self.ax.text(q0, y_fit_max, s, size=15, color='b', horizontalalignment='left', verticalalignment='bottom')
                self.ax.plot( [q0-self.results['fit_peaks_qT'], q0], [y_fit_max, y_fit_max], '-', color='b' )

                s = '$q_R = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$ \n $d_R = \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_qR'], self.results['fit_peaks_dR'])
                self.ax.text(q0, 0, s, size=15, color='r', horizontalalignment='left', verticalalignment='bottom')
                self.ax.plot( [q0-self.results['fit_peaks_qR'], q0], [yi, yi], '-', color='r' )
                
                
                
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

        outfile = self.get_outfile(data.name+'-fit', output_dir, ext='.png')
        
        try:
            #lines.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
            lines.plot(save=outfile, **run_args)
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
        #params.add('m', value=m, min=abs(m)*-10, max=abs(m)*+10)
        params.add('m', value=m, min=abs(m)*-5, max=0) # Slope must be negative
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







# Experimental parameters
########################################
#cms.SAXS.setCalibration([464.0, 552.0], 5.038, [35.00, 35.00]) # 2017-06-18, 13.5 keV
# Data collected at SAXSx = 10, SAXSy = 19
# i.e. 

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'

if True:
    # SAXS detector on CMS
    calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
    calibration.set_pixel_size(pixel_size_um=172.0)
    calibration.set_image_size(1475, height=1679) # Pilatus2M
    #calibration.set_beam_position(752, 1080) # SAXSx = -65, SAXSy = -73
    calibration.set_beam_position(764, 1279) # SAXSx = -60, SAXSy = -42
    
    calibration.set_distance(3.015)
    
    #mask = Mask(mask_dir+'Pilatus300k_main_gaps-mask.png')
    #mask.load('./Pilatus300k_current-mask.png')

    mask = Mask(mask_dir+'/Dectris/Pilatus2M_gaps-mask.png')
    mask.load('./Pilatus2M_CMS_3m-mask.png')
    

else:
    # WAXS detector on CMS
    from SciAnalysis.XSAnalysis.DataRQconv import *
    calibration = CalibrationRQconv(wavelength_A=0.9184) # 13.5 keV
    calibration.set_image_size(1042) # psccd Photonic Sciences CCD
    calibration.set_pixel_size(pixel_size_um=101.7)
    calibration.set_distance(0.232) # Bigger number moves theory rings outwards (larger spacing)
    calibration.set_beam_position(22.0, 1042-22.0)
    calibration.set_angles(det_orient=45, det_tilt=-21, det_phi=0, incident_angle=0., sample_normal=0.)
    print('ratio Dw = {:.3f}'.format(calibration.get_ratioDw()))
    
    mask = Mask(mask_dir+'psccd_generic-mask.png')



#match_re = re.compile('^.+_th0.110_\d+\.\ds_T(\d\d\d\.\d\d\d)C_5\.00s_\d+_saxs.jpg$')



# Files to analyze
########################################

#root_dir = '/GPFS/xf11bm/Pilatus300/'
#root_dir = '/GPFS/xf11bm/Pilatus300/2016-3Pilatus2M_current-mask.png/CFN_aligned-BCP/'
#source_dir = os.path.join(root_dir, '')

source_dir = '../'


#output_dir = os.path.join(source_dir, 'analysis/')
output_dir = './'


#pattern = 'AgBH*10.00s*'
pattern = '*'

#pattern = '*_th0.150_*'

pattern = '*_th0.150_*'


#pattern = '20170830_1x_60GdF3_PMMA_II_th0.150_5101.3s_x-0.000_15.00s_423535_saxs'


pattern, fit_range = '*_th0.150_*', [0.009, 0.033]
pattern, fit_range = '*_th0.150_*422099_saxs', [0.010, 0.020]



# qz_fit
pattern, fit_range = '20170823_48hBCP_I_th0.200_510.9s_x-1.500_15.00s_421982_saxs', [0.045, 0.070]

pattern, fit_range, q0_expected = '20170823_48hBCP_I_*', [0.045, 0.070], None
pattern, fit_range, q0_expected = '20170823_48hBCP_I_th0.100_*', [0.05, 0.080], 0.634
pattern, fit_range, q0_expected = '20170823_48hBCP_I_th0.150_*', [0.046, 0.065], 0.0567
#pattern, fit_range, q0_expected = '20170823_48hBCP_I_th0.200_*', [0.046, 0.065], 0.0560
#pattern, fit_range, q0_expected = '20170823_48hBCP_I_th0.250_*', [0.05, 0.063], 0.0560


pattern, fit_range, q0_expected = '20170823_48hBCP_I_th0.150_*421981*', [0.046, 0.065], 0.0567


infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))

infiles.sort()



# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             'rot180' : False,
             }
run_args = { 'verbosity' : 3,
            'critical_angle_film' : 0.098, # PMMA at 13.5 keV
            'critical_angle_substrate' : 0.132, # Si at 13.5 keV
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)



protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.01, num_rings=4, ztrim=[0.0, 0.01], zmin=0 ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    #Protocols.sector_average(angle=0, dangle=20, ylog=True, plot_range=[0, 0.3, None, None], show_region=True) ,
    #Protocols.linecut_angle(q0=0.094, dq=0.015, show_region=False) ,
    #linecut_qr_fit(show_region=False, show=False, qz=0.024, dq=0.006, fit_range=[0.01, 0.02], plot_range=[0, 0.03, 0, None]) ,
    #linecut_qr_fit(show_region=False, show=False, qz=0.024, dq=0.006, fit_range=[0.007, 0.015], plot_range=[0, 0.03, 0, None]) ,
    #linecut_qr_fit(show_region=False, show=False, qz=0.037, dq=0.01, fit_range=[-0.14, -0.07], plot_range=[-0.25, 0, 0, None]) , # C67
    #Protocols.q_image(show=True, blur=None, bins_relative=0.5, plot_range=[0.0, 0.05, 0, 0.05], _xticks=[-0.04, -0.02, 0, 0.02,], cmap=cmap_vge_hdr, ztrim=[0.05, 0.001]) ,
    #Protocols.qr_image(blur=None, bins_relative=0.8, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.38, 0.002], dezing_fill=True) ,
    #Protocols.q_phi_image(bins_relative=0.25, plot_range=[0, 3.0, 0, +90]) ,
    #Protocols.linecut_qr(show_region=True, show=True, qz=0.024, dq=0.006, plot_range=[0, -0.3, 0, None]) ,
    
    #Protocols.q_image(show=False, blur=None, bins_relative=0.5, plot_range=[-0.20, 0.20, 0, 0.22], xticks=[-0.15, 0.0, 0.15], cmap=cmap_vge_hdr, ztrim=[0.04, 0.0003], zmin=0) ,
    
    
    #Protocols.linecut_qr(qz=0.031, dq=0.008, ylog=True, plot_range=[0, 0.15, None, None]) ,
    #Protocols.linecut_qz(show_region=False, qr=0, dq=0.01, q_mode='qr', plot_range=[0, 0.25, 0, None]) ,
    
    
    #linecut_qr_fit(show_region=False, show=False, qz=0.031, dq=0.008, fit_range=fit_range, plot_range=[0, 0.05, 0, None]) ,
    linecut_qz_fit(show_region=False, show=False, qr=0.009, dq=0.0025, q_mode='qr', fit_range=fit_range, q0_expected=q0_expected, plot_range=[0, 0.08, 0, None]) ,
    #linecut_angle_fit(show_region=True, show=False, ylog=True, q0=q0, dq=0.02) ,
    
    
    #Protocols.thumbnails(crop=None, resize=0.5, blur=None, cmap=cmap_vge_hdr, ztrim=[0.01, 0.0001]) , # For GISAXS
    #Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge_hdr, ztrim=[0.001, 0.0001]) , # For SAXS
    ]
    


print('\n\n\n')
# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=True)


# Loop
########################################
# This code is typically only used at the beamline (it loops forever, watching for new files).
import time
donefiles = []
while False:

    #infiles = glob.glob(os.path.join(source_dir, '*.tiff'))
    infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))

    for infile in infiles:
        if infile in donefiles:
            pass

        else:
            process.run([infile], protocols, output_dir=output_dir, force=False)

        donefiles.append(infile)

    time.sleep(4)




print('\n\n\n')

