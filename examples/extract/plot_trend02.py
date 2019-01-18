#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
import numpy as np
from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols


# Select file/conditions
root_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/2015_10Oct_18-SAXS_all/'
source_dir = os.path.join(root_dir, 'CNC_825_1p2T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_825_0p56T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_825_0T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_413_1p2T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_165_1p2T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_165_0T_kin/')
#source_dir = os.path.join(root_dir, 'CNC_1p65_0p56T_1_kin/')

results_dir = os.path.join(source_dir, 'SciAnalysis/', 'results/')


output_dir = os.path.join(source_dir, 'SciAnalysis', 'trend/')

t_initial = 3 # min
t_step = 1 # min


#show = ['I', 'q', 'd', 'sigma', 'xi', 'eta', 'm', 'S'] # all
show = ['I', 'd', 'xi', 'eta', 'S'] # most useful
#show = ['q', 'sigma', 'm'] # aux


# Custom plot

class DataLinesStacked_current(DataLinesStacked):
    
    def analyze(self, **run_args):
        
        self.fit_lines = []
        self.fit_results = []

        for i, line in enumerate(self.lines):
            
            lm_result, fit_line, fit_line_extended = self.fit_sigmoid(line, **run_args)
            self.fit_lines.append(fit_line_extended)
            
            result = [lm_result.params['x_center'].value,
                      lm_result.params['x_center'].stderr,
                      lm_result.params['x_scale'].value,
                      lm_result.params['x_scale'].stderr,
                      lm_result.params['b'].value,
                      lm_result.params['b'].stderr,
                      lm_result.params['prefactor'].value,
                      lm_result.params['prefactor'].stderr,
                      ]
            self.fit_results.append(result)
            
            print( line.y_label )
            tau1 = lm_result.params['x_center'].value
            tau2 = np.abs(lm_result.params['x_scale'].value)
            tau = lm_result.params['x_center'].value+np.abs(lm_result.params['x_scale'].value)
            tauerr = np.sqrt( np.square(lm_result.params['x_center'].stderr) + np.square(lm_result.params['x_scale'].stderr) )
            
            print( '  tau1 = {:.2f} min tau2 = {:.2f} min'.format(tau1, tau2) )
            print( '  tau = {:.2f} +/- {:.2f} min'.format(tau, tauerr) )
            If = lm_result.params['b'].value+lm_result.params['prefactor'].value
            Iferr = np.sqrt( np.square(lm_result.params['b'].stderr) + np.square(lm_result.params['prefactor'].stderr) )
            print( '  I_i, I,f = {:.3f}+/-{:.3f}, {:.3f}+/-{:.3f}'.format(lm_result.params['b'].value, lm_result.params['b'].stderr, If, Iferr ) )
            print( '  <I>+/- = {} +/- {}'.format(np.average(line.y), np.std(line.y)) )
            
            


    def fit_sigmoid(self, line, **run_args):
        
        import lmfit
        
        def model(v, x):
            '''Eta orientation function.'''
            m = v['b']
            m += v['prefactor']*( 1./(1 + np.exp(-(x-v['x_center'])/v['x_scale'])) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        x_span = np.max(line.x)-np.min(line.x)
        # Determine 'polarity'
        line.sort_x()
        y_start = np.average(line.y[:5])
        y_end = np.average(line.y[-5:])
        if y_start>y_end:
            polarity = -1
        else:
            polarity = +1
       
        params = lmfit.Parameters()
        params.add('b', value=np.min(line.y), min=0, max=np.max(line.y))
        params.add('prefactor', value=np.max(line.y)-np.min(line.y), min=0)
        params.add('x_center', value=x_span, min=np.min(line.x), max=np.max(line.x))
        params.add('x_scale', value=x_span*0.25*polarity)
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        fit_x = np.linspace(min(np.min(line.x), 0), np.max(line.x)*2.0, num=1000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})

        return lm_result, fit_line, fit_line_extended

    
    def _plot_extra(self, **plot_args):
        
        
        if hasattr(self, 'fit_lines') and self.fit_lines is not None:
            for i, fline in enumerate(self.fit_lines):
                
                ax = getattr(self, 'ax{}'.format(i+1))
                xi, xf, yi, yf = ax.axis()
                
                p_args = dict([(i, fline.plot_args[i]) for i in self.plot_valid_keys if i in fline.plot_args])
                p_args['color'] = 'b'
                p_args['linewidth'] = 3.0
                ax.plot(fline.x, fline.y, **p_args)
                
                result = self.fit_results[i]
                
                xpos = result[0]
                xerr = result[1]
                ax.axvline(xpos, color='b', linewidth=1.0)
                ax.text(xpos, yf, '${:.0f}\pm{:.0f} \, \mathrm{{min}}$'.format(xpos, xerr), verticalalignment='top', horizontalalignment='left', color='b')
            
                xpos = result[2] + xpos
                xerr = result[3]
                ax.axvline(xpos, color='purple', linewidth=1.0)
                ax.text(xpos, yi, '${:.0f}\pm{:.0f} \, \mathrm{{min}}$'.format(xpos, xerr), verticalalignment='bottom', horizontalalignment='left', color='purple')

                xpos = result[0] + result[2]*0.5
                ax.text(xpos, yi+0.25*(yf-yi), '$+{:.0f} \, \mathrm{{min}}$'.format(result[2]), verticalalignment='center', horizontalalignment='center', color='purple')


        for i, line in enumerate(self.lines):
            
            ax = getattr(self, 'ax{}'.format(i+1))
            xi, xf, yi, yf = ax.axis()
            
            if line.y_label=='m':
                ax.axis( [xi, xf, 0, yf] )
            elif line.y_label=='eta':
                ax.axis( [xi, xf, 0, yf] )
            elif line.y_label=='S':
                ax.axis( [xi, xf, 0, 1] )


            #yticks = ax.get_yticks() # Positions
            yticks = ax.yaxis.get_major_ticks() # Objects
            yticks[0].label1.set_visible(False) # First
            #yticks[-1].label1.set_visible(False) # Last


def load_data(infile):
    
    # Extract data
    import re
    filename_re = re.compile('^.+_(\d+)\.xml$')

    data = []
    with open(infile) as fin:
        
        for i, line in enumerate(fin.readlines()):
            
            if i==0 and line[0]=='#':
                headers = line[1:].split()
            
            els = line.split()
            
            m = filename_re.match(els[0].strip())
            if m:
                exposure_id = int(m.groups()[0])
                time = t_initial + (exposure_id-1)*t_step
                
                els[0] = time
                els = [float(el) for el in els]
                data.append(els)
                
    return headers, np.asarray(data)
                



infile = os.path.join(source_dir, 'SciAnalysis', 'trend', 'q2I.txt')
headers, data = load_data(infile)
lines = DataLinesStacked_current()

times = data[:,0]
time_errs = np.ones(len(data))*t_step

#headers = ['filename', 'fit_peaks_prefactor1', 'fit_peaks_prefactor1_error', 'fit_peaks_x_center1', 'fit_peaks_x_center1_error', 'fit_peaks_sigma1', 'fit_peaks_sigma1_error']


# Create lines
if 'I' in show:
    # Intensity
    idx = headers.index('fit_peaks_prefactor1')
    idxe = headers.index('fit_peaks_prefactor1_error')
    line = DataLine(x=times, y=data[:,idx], x_err=time_errs, y_err=data[:,idxe], y_label='intensity', y_rlabel='$I(q_0)$')
    
    #line.remove_spurious(bins=5, tol=0.25)
    
    lines.add_line(line)

if 'q' in show:
    # q position
    idx = headers.index('fit_peaks_x_center1')
    idxe = headers.index('fit_peaks_x_center1_error')
    line = DataLine(x=times, y=data[:,idx], x_err=time_errs, y_err=data[:,idxe], y_label='pos', y_rlabel='$q_{0} \, (\AA^{-1})$')
    lines.add_line(line)

if 'd' in show:
    # d spacing
    idx = headers.index('fit_peaks_x_center1')
    idxe = headers.index('fit_peaks_x_center1_error')
    y = 0.1*2.*np.pi/data[:,idx]
    yerr = 0.1*2.*np.pi*data[:,idxe]/np.square(data[:,idx])
    
    y_rlabel = '$d_{0} = 2 \pi / q_0 \, (\mathrm{nm})$'
    y_rlabel = '$d_{0} \, (\mathrm{nm})$'
    line = DataLine(x=times, y=y, x_err=time_errs, y_err=yerr, y_label='pos', y_rlabel=y_rlabel)
    lines.add_line(line)

if 'sigma' in show:
    # sigma_q
    idx = headers.index('fit_peaks_sigma1')
    idxe = headers.index('fit_peaks_sigma1_error')
    line = DataLine(x=times, y=data[:,idx], x_err=time_errs, y_err=data[:,idxe], y_label='width', y_rlabel='$\sigma_{q_{0}}  \, (\AA^{-1})$')
    lines.add_line(line)

if 'xi' in show:
    # grain size
    idx = headers.index('fit_peaks_sigma1')
    idxe = headers.index('fit_peaks_sigma1_error')
    sigma_total = data[:,idx]
    sigma_total_err = data[:,idxe]
    
    sigma_instrumental = 0.0010 # A^-1
    sigma_instrumental_err = sigma_instrumental*0.25
    
    sigma = np.sqrt( np.square(sigma_total) - np.square(sigma_instrumental) )
    sigma_err = 0.5*np.sqrt(np.square(2*sigma_total*sigma_total_err) + np.square(2*sigma_instrumental*sigma_instrumental_err))/sigma
    
    sigma_to_fwhm = 2.0*np.sqrt(2.0*np.log(2.0))
    K = 0.9394
    
    fwhm = sigma*sigma_to_fwhm
    fwhm_err = sigma_err*sigma_to_fwhm
    
    xi = 0.1*2.*np.pi*K/(fwhm)
    xi_err = 0.1*2.*np.pi*K*sigma_err/np.square(sigma)
    
    line = DataLine(x=times, y=xi, x_err=time_errs, y_err=xi_err, y_label='grain size', y_rlabel=r'$\xi_{0}  \, (\mathrm{nm})$')
    lines.add_line(line)






infile = os.path.join(source_dir, 'SciAnalysis', 'trend', 'angle_fit.txt')
headers, data = load_data(infile)

if 'eta' in show:
    # eta
    idx = headers.index('fit_eta_eta')
    idxe = headers.index('fit_eta_eta_error')
    line = DataLine(x=times, y=data[:,idx], x_err=time_errs, y_err=data[:,idxe], y_label='eta', y_rlabel='$\eta$')
    lines.add_line(line)

if 'm' in show:
    idx = headers.index('fit_MaierSaupe_m')
    idxe = headers.index('fit_MaierSaupe_m_error')
    line = DataLine(x=times, y=data[:,idx], x_err=time_errs, y_err=data[:,idxe], y_label='m', y_rlabel='$m$')
    lines.add_line(line)


if 'S' in show:
    idx = headers.index('fit_library_eta_S')
    y1 = data[:,idx]
    idx = headers.index('fit_library_MaierSaupe_S')
    y2 = data[:,idx]
    
    y = (y1+y2)/2.0
    y_err = np.maximum( np.std( (y1, y2), axis=0 ), np.abs(y1-y2) )
    
    y /= -0.5
    y_err *= 0.5
    
    line = DataLine(x=times, y=y, x_err=time_errs, y_err=y_err, y_label='S', y_rlabel='$S/(-1/2)$')
    lines.add_line(line)



# Plot everything
lines.analyze(verbosity=2)

lines.x_label = 'time (min)'
lines.x_rlabel = 'time (min)'
plot_args = { 'a' : 2.0,
             'rcParams': {'axes.labelsize': 25,
                            'xtick.labelsize': 25,
                            'ytick.labelsize': 15,
                            },
            } 



els = source_dir.strip('/').split('/')
outfile = '{}/trend_q0-{}.png'.format(output_dir, els[-1])
lines.plot(save=outfile, show=False, xticks=range(0, 401, 50), error_band=True, plot_range=[0, 400, None, None], plot_buffers=[0.18,0.05,0.12,0.05], **plot_args)


                      
outfile = '{}/trend_q0-{}.txt'.format(output_dir, els[-1])
with open(outfile, 'w') as fout:
    
    values = ['x_center', 'x_center_error', 'x_scale', 'x_scale_error', 'b', 'b_error', 'prefactor', 'prefactor_error']
    fout.write( '#property\t{}\n'.format( '\t'.join(values) ) )
    for i, line in enumerate(lines.lines):
        
        results = lines.fit_results[i]
        results_str = [str(result) for result in results]
        fout.write( '{}\t{}\n'.format( line.y_label, '\t'.join(results_str) ) )
        
        
        