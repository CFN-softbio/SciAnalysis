#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH='/GPFS/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import time

import glob
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols



# Define some custom analysis routines
########################################

from SciAnalysis.Result import * # The Results() class allows one to extract data from XML files.


        

class update_autonomous_data(Protocols.Protocol):
    
    def __init__(self, name='autonomous', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'crop' : None,
                        'shift_crop_up' : 0.0,
                        'blur' : 2.0,
                        'resize' : 0.2,
                        'ztrim' : [0.05, 0.005]
                        }
        self.run_args.update(kwargs)
        

    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        import time
        time.sleep(1) # Kludge to avoid corrupting XML file?
        
        
        # Compile results
        results_xml_file = self.get_outfile(data.name, './results/', ext='.xml')
        
        if run_args['verbosity']>=5:
            print('    Extracting results from: {}'.format(results_xml_file))
        
        
        
        infiles = [results_xml_file]
        #extractions = [ [ 'metadata_extract', ['x_position','y_position'] ] ,
                    #['circular_average_q2I_fit', ['fit_peaks_d0', 'fit_peaks_grain_size', 'fit_peaks_prefactor1'] ],
                    #]
        extractions = [ [ 'metadata_extract', ['x_position', 'y_position', 'cost'] ] ,
                    ['circular_average_q2I', ['fit_peaks_prefactor1', 'fit_peaks_x_center1', 'fit_peaks_sigma1', 'fit_peaks_chi_squared' ] ],
                    ]

        
        extracted_results = Results().extract_multi(infiles, extractions)
        
        data_vector = {}
        for ir, result in enumerate(extracted_results):
            data_vector[ir] = {
                    'x_position' : float(result[1]),
                    'y_position' : float(result[2]),
                    'cost' : float(result[3]),
                    'fit_peaks_prefactor1' : float(result[4]),
                    'fit_peaks_x_center1' : float(result[5]),
                    'fit_peaks_sigma1' : float(result[6]),
                    'fit_peaks_chi_squared' : float(result[7]),
                    
                }
        
        
        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)
        
        #print(data.stats())
        
        np.save(outfile, data_vector)
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'data vector (npy)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        
        
        
        # Update the file for SMART to analyze
        outfile = SMART_DIRECTORY+'new_experiment_result/experiment_result.npy'
        
        
        print('    Saving {}'.format(outfile))
        np.save(outfile, data_vector)
        if os.path.isfile(outfile):
            print('    Saved (file exists).')
        else:
            print('    Saved (file missing).')
        
        #if os.path.isfile(outfile):
            ## File exists; append data instead of over-writing
            #old_data_vector = np.load(outfile)
            #data_vector = np.append( old_data_vector, [data_vector], axis=0 )
            #np.save(outfile, data_vector)
            
        #else:
            #np.save(outfile, [data_vector])
        
        
        return results    
        


    
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
        
        #line.y *= np.abs(line.x)
        #line.y *= np.square(line.x)
        line.y *= np.power(np.abs(line.x), 3.0)
        
        line.y_label = 'q^n*I(q)'
        line.y_rlabel = '$q^n I(q) \, (\AA^{-n} \mathrm{counts/pixel})$'

        outfile = self.get_outfile(data.name, output_dir, ext='_q2I.dat')
        line.save_data(outfile)        

        if True:
            # Fit
            lm_result, fit_line, fit_line_extended = self._fit_peaks(line, **run_args)
            #print(lm_result.params)
            
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
                    
                    xstart, xend = self._run_args['fit_range']
                    line = self.lines[0].sub_range(xstart, xend)
                    #print(line.y)
                    yf = np.max(line.y)*1.5
                    self.ax.axis([xi, xf, yi, yf])
                    
                    
                    v_spacing = (yf-yi)*0.10

                    yp = yf
                    s = '$\chi^2 = \, {:.4g}$'.format(self.results['fit_peaks_chi_squared'])
                    self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')
                    
                    yp -= v_spacing
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
            lines.results = results
            lines._run_args = run_args
            lines.copy_labels(line)
            
        else:
            lines = DataLines([line])
            
            
           
        
        outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
        lines.plot(save=outfile, **run_args)
        
        return results
    
    
        
    def _fit_peaks(self, line, q0, vary=True, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)
        
        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            '''Gaussians with constant background.'''
            m = v['b']
            
            i = 0
            q0 = v['x_center1']
            m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-q0)/(2*(v['sigma{:d}'.format(i+1)]**2)) )

            #i = 1
            #q0 = v['x_center1']*np.sqrt(3.0)
            #m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-q0)/(2*(v['sigma{:d}'.format(i+1)]**2)) )

            #i = 2
            #q0 = v['x_center1']*2.0
            #m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-q0)/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        b = np.min(line.y)
        params.add('b', value=b, min=b*0.98, max=b*1.02+1e-10, vary=vary)
        
        
        if 'sigma' in run_args:
            sigma = run_args['sigma']
        else:
            sigma = 0.001
            
        xspan = np.max(line.x) - np.min(line.x)
        if q0 is None:
            xpeak, ypeak = line.target_y(np.max(line.y))
            q0 = xpeak
        else:
            xpeak, ypeak = line.target_x(q0)
            
        
        i = 0
        params.add('prefactor{:d}'.format(i+1), value=(ypeak-b), min=0, max=np.max(line.y), vary=vary)
        params.add('x_center{:d}'.format(i+1), value=q0, min=q0*0.5, max=q0*1.5, vary=vary)
        params.add('sigma{:d}'.format(i+1), value=sigma, min=0, vary=vary)

        if False:
            i = 1
            xpeak, ypeak = line.target_x(q0*np.sqrt(3))
            params.add('prefactor{:d}'.format(i+1), value=ypeak-b, min=0, max=np.max(line.y), vary=vary)
            params.add('sigma{:d}'.format(i+1), value=sigma, min=0, vary=vary)

            i = 2
            xpeak, ypeak = line.target_x(q0*2)
            params.add('prefactor{:d}'.format(i+1), value=ypeak-b, min=0, max=np.max(line.y), vary=vary)
            params.add('sigma{:d}'.format(i+1), value=sigma, min=0, vary=vary)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=2000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':1.0})        

        return lm_result, fit_line, fit_line_extended



class linecut_angle_fit_custom(Protocols.linecut_angle):

    def __init__(self, name='linecut_angle_fit', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False,
                         'plot_range' : [-180, 180, 0, None]
                         }
        self.run_args.update(kwargs)
            
            
        
    @tools.run_default
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
        #line.kill_x(0, 3)
        
        if 'show_region' in run_args and run_args['show_region']:
            data.plot(show=True)
        
        # Basic output
        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        #line.save_data(outfile)
        
        
        # Fine-tune data
        line.x = line.x[1:]
        line.y = line.y[1:]
        line.smooth(2.0)
        
        
        do_max = True
        do_FWHM = False
        do_fits = True
        
        
        # Fit
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                v_spacing = 0.1
                
                xi, xf, yi, yf = self.ax.axis()

                #yi = min( self.lines[-1].y )*0.5
                #yf = max( self.lines[-1].y )*1.5
                #self.ax.axis( [xi, xf, yi, yf] )

                self.ax.axvline(self.angle, color='r', linewidth=1.0)
                
                s = '$\chi = {:.1f} ^{{\circ}} $'.format(self.angle)
                self.ax.text(self.angle, yf, s, size=30, color='r', verticalalignment='top', horizontalalignment='left')
                s = '$f = {:.3f}$'.format(self.orientation_factor)
                self.ax.text(self.angle, yi, s, size=30, color='r', verticalalignment='bottom', horizontalalignment='left')
                
                
                
                if do_max:
                    line = self.lines[-1]
                    xpeak, ypeak = line.target_y(np.max(line.y))
                    self.ax.plot([xpeak], [ypeak], 'o-', color='b', markersize=12, linewidth=1.0)
                    self.ax.axvline(xpeak, color='b', linewidth=1.0)
                    
                    s = '$\chi = {:.1f} ^{{\circ}} $'.format(xpeak)
                    self.ax.text(xpeak, ypeak, s, size=30, color='b', verticalalignment='bottom', horizontalalignment='left')
                    

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
            for i, fit_name in enumerate(['fit_eta', 'fit_MaierSaupe', 'fit_eta_span']):
                
                lm_result, fit_line, fit_line_e = getattr(self, fit_name)(line, **run_args)
                fit_line.plot_args['color'] = color_list[i%len(color_list)]
                fit_line_e.plot_args['color'] = color_list[i%len(color_list)]
                
                lines.add_line(fit_line_e)
                if fit_name=='fit_eta_span':
                    lines.add_line(fit_line)
                    
            
                #prefactor_total = 0
                for param_name, param in lm_result.params.items():
                    results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
                    lines.results['{}_{}'.format(fit_name, param_name)] = param.value
                    #if 'prefactor' in param_name:
                        #prefactor_total += np.abs(param.value)
                    
                #results['{}_prefactor_total'.format(fit_name)] = prefactor_total
                results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
            
            
        angle = results['fit_eta_span_x_center']['value']
        if angle<0:
            angle += 180
        # angle is now 0 to +190
        # Convert to the 'in between the peaks' convention
        angle -= 90 
        # angle is now -90 to +90
        lines.angle = angle
        
        qz_unit = np.cos(np.radians(angle))
        qx_unit = np.sin(np.radians(angle))
        orientation_factor = 2*np.square(qz_unit) - 1
        lines.orientation_factor = orientation_factor
        
        
        results['orientation_angle'] = angle
        results['orientation_factor'] = orientation_factor
            
        # Output
        if run_args['verbosity']>=1:
            outfile = self.get_outfile(data.name, output_dir)
            lines.lines.reverse()
            lines.plot(save=outfile, **run_args)

        if run_args['verbosity']>=4:
            outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
            line.plot_polar(save=outfile, **run_args)

        return results
    
    
    def fit_eta_span(self, line, span=30, **run_args):
        
        import lmfit
        
        xpeak, ypeak = line.target_y(np.max(line.y))
        if xpeak<0:
            xpeak += 180
        line_full = line
        
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
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':'o', 'linewidth':0.0})
        
        fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=1000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':2.0})

        return lm_result, fit_line, fit_line_extended
                    


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
        
        xpeak, ypeak = line.target_y(np.max(line.y))
        
        params = lmfit.Parameters()
        params.add('prefactor', value=np.max(line.y)*0.1, min=0)
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




def autonomous_result(xml_file):
    
    time.sleep(1) # Kludge to avoid corrupting XML file?
    
    extractions = [ [ 'metadata_extract', ['x_position', 'y_position'] ] ,
                ['circular_average_q2I', ['fit_peaks_prefactor1', 'fit_peaks_x_center1', 'fit_peaks_sigma1', 'fit_peaks_chi_squared' ] ],
                ]
    #extractions = [ [ 'metadata_extract', ['x_position', 'y_position'] ] , # 1, 2
                #['circular_average_q2I_fit', ['fit_peaks_prefactor1', 'fit_peaks_x_center1', 'fit_peaks_sigma1', 'fit_peaks_chi_squared', 'fit_peaks_prefactor1_error', 'fit_peaks_x_center1_error', 'fit_peaks_sigma1_error' ] ], # 3, ... 9
                #['linecut_angle_fit', ['fit_eta_eta', 'orientation_factor', 'orientation_angle', 'fit_eta_span_prefactor'] ], # 10, ... 13
                #]


    infiles = [xml_file]
    extracted_results = Results().extract_multi(infiles, extractions)
    
    result = extracted_results[0]
    data_vector = {
            #'x_position' : float(result[1]),
            #'y_position' : float(result[2]),
            'prefactor' : float(result[3]),
            'q' : float(result[4]),
            'sigma' : float(result[5]),
            'chi_squared' : float(result[6]),
            
        }
    #data_vector = {
            ##'x_position' : float(result[1]),
            ##'y_position' : float(result[2]),
            #'prefactor' : float(result[3]),
            #'prefactor variance' : float(result[7]),
            #'q' : float(result[4]),
            #'q variance' : float(result[8]),
            #'sigma' : float(result[5]),
            #'sigma variance' : float(result[9]),
            #'chi_squared' : float(result[6]),
            #'chi_squared variance' : 0,
            #'eta' : float(result[10]),
            #'orientation_factor' : float(result[11]),
            #'orientation_angle' : float(result[12]),
            #'eta_prefactor' : float(result[13]),
            
        #}    
    
    return data_vector
    
    

# Tweaking data
########################################
# This code can be modified to tweak the data being passed around to control SMART

def analyze_more(communication_directory='../../Code/smart/data/', source_dir='../', pattern='*', save=False, verbosity=3):
    
    
    #input_file = communication_directory+'new_experiment_command/analysis_command.npy'
    input_file = communication_directory+'new_experiment_command/analysis_command_bak.npy'
    output_file = communication_directory+'new_experiment_result/experiment_result.npy'
    output_file = './experiment_result.npy'


    # Determine files to analyze
    infiles = glob.glob(source_dir+pattern+'.tiff')
    
    if verbosity>=3:
        print('Found {} infiles...'.format(len(infiles)))
        
        
    # Select only the 'right' files
    filename_re = re.compile('.+_(\d+)_saxs.+') # sequence_ID
    
    ##AE_LGX_C67blend_P045_run1_clock40304.16_localT492.0_th0.120_5.00s_2532754_saxs.tiff
    #filename_re = re.compile('.+_clock(\d+\.\d\d)_localT(\d+\.\d)_.+_(\d+)_saxs.+') # clock, T, sID
    
    to_analyze = []
    for infile in infiles:
        m = filename_re.match(infile)
        if m:
            sequence_ID = int(m.groups()[0])
            #clock, temperature, sequence_ID = float(m.groups()[0]), float(m.groups()[1]), int(m.groups()[2])
            if verbosity>=5:
                print('considering file: {}'.format(infile))

                print('    sID {:d}'.format(sequence_ID))
                #print('    clock: {:.2f}, T {:.1f}, sID {:d}'.format(clock, temperature, sequence_ID))
            #if sequence_ID==2533156:
            to_analyze.append(infile)
                
    if verbosity>=3:
        print('RE match {} files...'.format(len(to_analyze)))


    # Load existing analysis data
    #results = np.load(input_file).item()
    #index = max(results.keys()) + 1
    
    results = {}
    index = 0
    
    for infile in to_analyze:
        if False:
            if verbosity>=3:
                print('        Analysis for index {}, file: {}'.format(index, infile))
                print('\n')
            
            process.run([infile], protocols, output_dir=output_dir, force=False)
            
            if verbosity>=3:
                print('\n')
        
        #xml_file = '{}{}{}'.format( './results/', infile[len(source_dir):-5], '.xml' )
        xml_file = '{}{}{}{}'.format( source_dir, '../analysis/results/', infile[len(source_dir):-5], '.xml' )
        result = autonomous_result(xml_file)
        
        if verbosity>=5:
            print(result)
        
        m = filename_re.match(infile)
        if m:
            sequence_ID = int(m.groups()[0])
            #clock, temperature, sequence_ID = float(m.groups()[0]), float(m.groups()[1]), int(m.groups()[2])
        
            result['filename'] = infile
            #result['cost'] = 1
            
            #result['laser_power'] = 0
            #result['sweep_index'] = 0
            #result['l_position'] = lpos
            #result['o_position'] = opos
            
            result['measured'] = True
            result['analyzed'] = True
            results[index] = result
            index += 1
            
        else:
            if verbosity>=1:
                print("WARNING: RE didn't match for: {}".format(infile))
                
    if verbosity>=4:
        print( 'Results contains {} results, from index {} to {}'.format( len(results), min(results.keys()), max(results.keys()) ) )
    if verbosity>=5:
        for index in results:
            print(index)
            for key in results[index]:
                print ('    ', key, ':', results[index][key])

    if save:
        # Save analysis results for SMART to consider
        if os.path.isfile(output_file):
            time.sleep(5)
        if os.path.isfile(output_file):
            if verbosity>=3:
                print('    Saving experiment_result.npy (file exists).')
        else:
            if verbosity>=3:
                print('    Saving experiment_result.npy (file missing).')
        np.save(output_file, results)
            


def fix_results(infile, save=False, verbosity=3):
    
    
    results = np.load(infile).item()
    index = max(results.keys()) + 1
    
    filename_re = re.compile('.+_opos(-?\d+\.\d+)_lpos(-?\d+\.\d+)_.+_(\d+)_saxs.+')
    
    for index, result in results.items():
        
        if verbosity>=4:
            print('doing index {}, infile {}'.format(index, result['filename']))
        
        m = filename_re.match(result['filename'])
        if m:
            result['l_position'] = m.groups()[0]
            result['o_position'] = m.groups()[1]
            
        else:
            if verbosity>=1:
                print("WARNING: RE didn't match for: {}".format(result['filename']))

    if save:
        # Save analysis results for SMART to consider
        if os.path.isfile(infile):
            time.sleep(5)
        if os.path.isfile(infile):
            if verbosity>=3:
                print('    Saving experiment_result.npy (file exists).')
        else:
            if verbosity>=3:
                print('    Saving experiment_result.npy (file missing).')
        np.save(infile, results)



# Special data instantiation
if False:
    
    pattern = 'AMJ269_S20_e_AM_*'
    #pattern = 'AMJ269_S20_e_AM_x7.876_yy-19.400_1.00s_2218885_saxs'
    infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
    infiles.sort()
    
    print('    SciAnalysis processing {} files...'.format(len(infiles)))
    result_vector = {}
    for index, infile in enumerate(infiles):
        
        print('        Analysis for index {}, file: {}'.format(index, infile))
        #print('\n')
        #process.run([infile], protocols, output_dir=output_dir, force=True)
        #print('\n')
        
        xml_file = '{}{}{}'.format( './results/', infile[len(source_dir):-5], '.xml' )
        result = autonomous_result(xml_file)
        result['filename'] = infile
        
        result_vector[index] = result
        
    
    # Update the file for SMART to analyze
    outfile = SMART_DIRECTORY+'new_experiment_result/initialize_experiment_result.npy'
    np.save(outfile, result_vector)





# Autonomous execution
########################################
def load_retry(infile, max_retries=10, wait_time=0.1):
    attempts = 0
    while attempts < max_retries:
        try:
            results = np.load(infile, allow_pickle=True).item()
            break
        except Exception as ex:
            attempts += 1
            print("    Attempt {}: Exception {} occurred while loading {}".format(attempts+1, ex.__class__.__name__, infile))
            print(ex)
            
            time.sleep(wait_time)
            wait_time *= 2
            
            if attempts>=max_retries:
                raise ex
            
    return results

def run_autonomous_loop(protocols, communication_directory='../../Code/smart/data/', verbosity=3):

    if verbosity>=3:
        print('\n\n\n')
        print('=============================')
        print('==  Autonomous Experiment  ==')
        print('=============================')
        print('\n')
        
    
    input_file = communication_directory+'new_experiment_command/analysis_command.npy'
    output_file = communication_directory+'new_experiment_result/experiment_result.npy'

    
    #filename_re = re.compile('.+_(\d+)_saxs.+') # In case we need the sequence_ID
    #filename_re = re.compile('.+_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)_saxs.+')
    
    
    # Loop forever
    while True:
    
        # Wait for analysis command to appear
        start_time = time.time()
        while not os.path.isfile(input_file):
            if verbosity>=3:
                print(' Waiting for command file (analysis_command.npy) [{:.1f} minutes]'.format( (time.time()-start_time)/60 ))
            time.sleep(10)
            
        if verbosity>=2:
            print('New command file found (waited {:.1f} minutes)'.format((time.time()-start_time)/60) )
            
        # Open command file
        time.sleep(0.2)
        results = load_retry(input_file)
  
        
        if verbosity>=3:
            num_items = sum(1 for v in results.values() if 'analyzed' in v and v['analyzed'] is False)
            print('    SciAnalysis processing {} files...'.format(num_items))
        if verbosity>=10:
            print_dictionary(results)
        
        for index, result in results.items():
            
            # Ad hoc code to control analysis
            analyze = True
            #m = filename_re.match(result['filename'])
            #if m:
                #sequence_ID = int(m.groups()[0])
                #if sequence_ID>=863 and sequence_ID<=871:
                    #analyze = False        
            
            
            #if 'filename' in result: # Update all analysis
            if 'analyzed' in result and result['analyzed'] is False and 'filename' in result and analyze:

                infile = source_dir + result['filename']
                if verbosity>=3:
                    print('        Analysis for index {}, file: {}'.format(index, infile))
                    print('\n')
                
                process.run([infile], protocols, output_dir=output_dir, force=False)
                
                if verbosity>=3:
                    print('\n')
                
                xml_file = '{}{}{}'.format( './results/', infile[len(source_dir):-4], '.xml' )
                new_result = autonomous_result(xml_file)
                result.update(new_result)
                
                if verbosity>=5:
                    print(result)
                
                results[index] = result
                results[index]['analyzed'] = True
            
            
        
        # Update the file for SMART to analyze
        if verbosity>=3:
            print('    Removing command file (analysis_command.npy)')
        
        #os.remove(input_file)
        os.replace(input_file, communication_directory+'/new_experiment_command/analysis_command_bak.npy')
        
        
        
        if os.path.isfile(output_file):
            time.sleep(10)
            
        if os.path.isfile(output_file):
            if verbosity>=3:
                print('    Saving experiment_result.npy (file exists).')
        else:
            if verbosity>=3:
                print('    Saving experiment_result.npy (file missing).')
        np.save(output_file, results)
            
        time.sleep(1)
    


# Debug code
########################################
# This code is convenient for inspecting the inter-process files being passed around
# during a SMART autonomous experiment.

def print_dictionary(results):
    for index in results:
        print(index)
        for key in results[index]:
            print ('    ', key, ':', results[index][key])
            
# Example structure of a SMART dictionary
results = { 
    0 : {'x': 1.0, 'y': 1.0, 'prefactor': 10.0, 'prefactor variance': 0.1, 'cost': 1.0, 'measured': True, 'analyzed': True} ,
    1 : {'x': 2.0, 'y': 1.0, 'prefactor': 11.0, 'prefactor variance': 0.1, 'cost': 1.0, 'measured': True, 'analyzed': True} ,
    3 : {'x': 3.0, 'y': 2.0, 'prefactor': 10.0, 'prefactor variance': 0.1, 'cost': 1.0, 'measured': True, 'analyzed': False} ,
    }        
            
if False:
    # Examples of commands that can be run in an ipython session to inspect SMART dictionaries    
    import numpy as np
    results = np.load('experiment_command.npy.bak', allow_pickle=True).item() ; print_dictionary(results)
    results = np.load('analysis_command_bak.npy', allow_pickle=True).item() ; print_dictionary(results)
    max_cost = max( [result['cost'] for result in results.values()] )
    num_results = len( [ key for key in results.keys() ] )

    # Print dictionary
    def dict_str(d): return ''.join([ '{}\n{}'.format(i, sv(v)) for i, v in d.items() ])
    def sv(v): return ''.join([ se(k, e) for k, e in v.items() ])
    def se(k,e): return ''.join(['    ', str(k), ': ', str(e), '\n'])
    print(dict_str(results))    
    





# Experimental parameters
########################################

calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
calibration.set_image_size(1475, height=1679) # Pilatus2M
calibration.set_pixel_size(pixel_size_um=172.0)
#calibration.set_beam_position(765.0, 1680-579) # SAXSx -60, SAXSy -71
calibration.set_beam_position(757.0-1, 1680-530) # SAXSx -65, SAXSy -60

#calibration.set_distance(5.038) # 5m
#calibration.set_distance(2.001) # 2m
calibration.set_distance(5.09) # 5m

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'
mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
mask.load('./Pilatus2M_CMS_5m-circ.png')



# Files to analyze
########################################
source_dir = '../raw/'
output_dir = './'

pattern = '*'

infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))
infiles.sort()


# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)




patterns = [
            ['theta', '.+_th(\d+\.\d+)_.+'] ,
            ['x_position', '.+_x(-?\d+\.\d+)_.+'] ,
            ['y_position', '.+_yy(-?\d+\.\d+)_.+'] ,
            ['cost', '.+_Cost(\d+\.\d+)_.+'] ,
            ['annealing_temperature', '.+_T(\d+\.\d\d\d)C_.+'] ,
            ['annealing_time', '.+_(\d+\.\d)s_T.+'] ,
            ['exposure_time', '.+_(\d+\.\d+)c_\d+_saxs.+'] ,
            ['sequence_ID', '.+_(\d+)_saxs.+'] ,
            ]


protocols = [
    #Protocols.thumbnails(name='thumbnails2', crop=None, resize=1.0, blur=None, cmap=cmap_vge, ztrim=[0.0, 0.01])
    #circular_average_q2I_fit(show=False, q0=None, plot_range=[0, 0.025, 0, None], fit_range=[0.010, 0.022], vary=True) ,
    #circular_average_q2I_fit(show=False, q0=None, plot_range=[0, 0.025, 0, None], fit_range=[q0-0.006, q0+0.004], vary=True) ,
    linecut_qr_fit(show_region=False, show=False, qz=0.0285, dq=0.006, q0=0.0030, sigma=0.0010, fit_range=[0.00, 0.010], plot_range=[0, 0.050, 0, None], trim_range=[0, 0.05]) ,
    linecut_angle_fit_custom(q0=0.0135, dq=0.0035, show=False, show_region=False) ,
    Protocols.metadata_extract(patterns=patterns) ,
    ]
    

# Normal SciAnalysis loop
if False:

    # Run
    ########################################
    print('Processing {} infiles...'.format(len(infiles)))
    process.run(infiles, protocols, output_dir=output_dir, force=True)


    # Loop
    ########################################
    # This code is typically only used at the beamline (it loops forever, watching for new files).
    # This is typically only used at the beamline (it loops forever, watching for new files).
    #process.monitor_loop(source_dir=source_dir, pattern='*.tiff', protocols=protocols, output_dir=output_dir, force=False)


if False:
    run_autonomous_loop(protocols)





