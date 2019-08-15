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
    

    def __init__(self, name='linecut_qr', **kwargs):
        
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
        line.plot(save=outfile, **run_args)

        #outfile = self.get_outfile(data.name, output_dir, ext='_polar.png')
        #line.plot_polar(save=outfile, **run_args)

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)
        
        
        
        # Fit data
        if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
            line.trim(run_args['fit_range'][0], run_args['fit_range'][1])
        
        lm_result, fit_line, fit_line_extended = self.fit_peaks_special(line, **run_args)
        
        
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
        
        #print(results)
        
        return results        
        
        
        
        return results


    

    def fit_peaks(self, line, num_curves=1, **run_args):
        
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
                
                
    def fit_peaks_special(self, line, num_curves=1, **run_args):
        
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
        params.add('m', value=m, min=abs(m)*+10, max=0) # Slope must be positive
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
#cms.SAXS.setCalibration([464.0, 552.0], 5.038, [35.00, 35.00]) # 2017-06-18, 13.5 keV
# Data collected at SAXSx = 10, SAXSy = 19
# i.e. 

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'

if True:
    # SAXS detector on CMS
    calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
    calibration.set_pixel_size(pixel_size_um=172.0)
    calibration.set_image_size(1475, height=1679) # Pilatus2M
    #calibration.set_beam_position(761.0, 1680-579)
    #calibration.set_beam_position(760.0, 1680-591)
    calibration.set_beam_position(765.0, 1680-579) # 2019C1
    
    calibration.set_distance(2.001)
    
    #mask = Mask(mask_dir+'Pilatus300k_main_gaps-mask.png')
    #mask.load('./Pilatus300k_current-mask.png')

    mask = Mask(mask_dir+'Dectris/Pilatus2M_gaps-mask.png')
    #mask.load('./Pilatus2M_CMS_2m-mask_2019C1.png')
    mask.load('./Pilatus2M_CMS_2m-circ-mask.png')
    
    

else:
    # WAXS detector on CMS
    from SciAnalysis.XSAnalysis.DataRQconv import *
    calibration = CalibrationRQconv(wavelength_A=0.9184) # 13.5 keV
    calibration.set_image_size(1042) # psccd Photonic Sciences CCcalibration.set_beam_position(769.0, 1680-589)D
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



pattern = 'AgBH_calibration_standard_2m_10.00s_798020_saxs*'
pattern = '*'
#pattern = 'run10*'
pattern = '*_2003752_saxs'
pattern = '*_2001081_saxs'
pattern = '*_2000674_saxs'







infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))

infiles.sort()



# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             #'rot180' : False,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)




protocols = [
    #Protocols.export_STL( stl_zscale=150, stl_pedestal=8, blur=1.0, resize=1.0, crop_GI=700, logo_file='logo.png', logo_resize=0.5),
    Protocols.export_STL( stl_zscale=150, stl_pedestal=8, blur=1.0, resize=1.0, crop_beam=[900,900], logo_file='logo.png', logo_resize=0.5),
    ]
    



# Run
########################################
print('Processing {} infiles...'.format(len(infiles)))
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






