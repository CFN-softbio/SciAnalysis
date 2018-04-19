#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
# Update this to point to the directory where you copied the SciAnalysis base code
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis/'
#SciAnalysis_PATH='/home/xf11bm/software/SciAnalysis2/'
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
                
                
    def fit_peaks_special(self, line, num_curves=1, alpha=False, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)
        
        line_full = line
        #if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            
            if alpha:
                m = v['prealpha']*np.power(x, -v['alpha'])
            else:
                m = v['m']*x + v['b'] # linear background
                
            
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()



        if alpha:
            params.add('alpha', value=3, min=0, max=5.0)
            p = (line.y[0]+line.y[1])*0.5/( np.power(np.average(line.x), 3) )
            params.add('prealpha', value=p, min=0, max=p*100)
        else:
            m = (line.y[-1]-line.y[0])/(line.x[-1]-line.x[0])
            b = line.y[0] - m*line.x[0]
            
            params.add('m', value=m, min=abs(m)*-10, max=0) # Slope must be negative
            params.add('b', value=b, min=0)

        
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        for i in range(num_curves):
            #xpos = np.min(line.x) + xspan*(1.*i/num_curves)
            #xpos, ypos = line.target_x(xpeak*(i+1))
            xpos, ypos = xpeak, ypeak
            
            params.add('prefactor{:d}'.format(i+1), value=np.max(line.y), min=0, max=np.max(line.y)*10)
            params.add('x_center{:d}'.format(i+1), value=0.016, min=np.min(line.x), max=np.max(line.x)) # L74
            #params.add('x_center{:d}'.format(i+1), value=0.018, min=np.min(line.x), max=np.max(line.x)) #C67
            params.add('sigma{:d}'.format(i+1), value=0.001, min=0.0002, max=xspan*0.5)


            #params.add('prefactor{:d}'.format(i+1), value=1, min=0, max=np.max(line.y)*10, vary=True)
            #params.add('x_center{:d}'.format(i+1), value=0.022, min=np.min(line.x), max=np.max(line.x), vary=True) #C67
            #params.add('sigma{:d}'.format(i+1), value=0.0026, min=0.0002, max=xspan*0.5, vary=True)
        
        
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



class metadata_extract(tools.Protocol):
    
    def __init__(self, name='metadata_extract', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.dat'
        self.run_args = {}
        self.run_args.update(kwargs)
    
        
    @Protocols.run_default
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
        
        
        patterns = [
                    ['theta', '.+_th(\d+\.\d+)_.+'] ,
                    ['annealing_temperature', '.+_T(\d+\.\d\d\d)C_.+'] ,
                    ['annealing_time', '.+_(\d+\.\d)s_T.+'] ,
                    ['exposure_time', '.+_(\d+\.\d+)c_\d+_saxs.+'] ,
                    ['sequence_ID', '.+_(\d+)_saxs.+'] ,
                    ]
                    
        
        for pattern_name, pattern_string in patterns:
            pattern = re.compile(pattern_string)
            m = pattern.match(filename)
            if m:
                if run_args['verbosity']>=5:
                    print('  matched: {} = {}'.format(pattern_name, m.groups()[0]))
                results[pattern_name] = float(m.groups()[0])
        
        
        #outfile = self.get_outfile(data.name, output_dir)
        
        
        return results
    



# Experimental parameters
########################################
#cms.SAXS.setCalibration([464.0, 552.0], 5.038, [35.00, 35.00]) # 2017-06-18, 13.5 keV
# Data collected at SAXSx = 10, SAXSy = 19
# i.e. 

mask_dir = SciAnalysis_PATH + '/SciAnalysis/XSAnalysis/masks/'

if True:
    # SAXS detector on CMS
    calibration = Calibration(wavelength_A=0.9184) # 13.5 keV
    calibration.set_image_size(487, height=619) # Pilatus300k
    calibration.set_pixel_size(pixel_size_um=172.0)
    #calibration.set_beam_position(464.0, 552.0) # SAXSx = 35, SAXSy = 35
    calibration.set_beam_position(464.0-145.35, 552.0-93.02) # SAXSx = 10, SAXSy = 19
    calibration.set_distance(5.068)
    
    mask = Mask(mask_dir+'Pilatus300k_main_gaps-mask.png')
    mask.load('./Pilatus300k_current-mask.png')
    

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



match_re = re.compile('^.+_th0.110_\d+\.\ds_T(\d\d\d\.\d\d\d)C_5\.00s_\d+_saxs.jpg$')



# Files to analyze
########################################

#root_dir = '/GPFS/xf11bm/Pilatus300/'
#root_dir = '/GPFS/xf11bm/Pilatus300/2016-3/CFN_aligned-BCP/'
#source_dir = os.path.join(root_dir, '')

source_dir = '../'


#output_dir = os.path.join(source_dir, 'analysis/')
output_dir = './'


pattern = '*'
#pattern = 'AgBH_*'
#pattern = 'L74S3M3hp50_th0.080_2475.9s_T150.007C_5.00s_56900_saxs'
#pattern = 'L74S3M3hp30_th0.110_1897.2s_T199.994C_5.00s_57093_saxs'
#pattern = '*_th0.110_*'
#pattern = '*NoPB_th0.110_*'

#pattern = 'L74S6M6hp70_NoPB_th0.110_*'
#pattern = 'C67*_th0.110*'



# run01 L74
#pattern = 'L74*_th0.110_*_5.00s_*_saxs'
#pattern = 'L74S3M3hp70_th0.110_*_5.00s_*_saxs'

# run02 L74_NoPB
pattern = 'L74*NoPB_th0.110_*_5.00s_*_saxs'
pattern = 'L74S6M6hp70_NoPB_th0.110_*_5.00s_*_saxs'

# run03 L74_NoPB_S (S8-S12)
pattern = 'L74*NoPB_S8_th0.110_*_5.00s_*_saxs'
#pattern = 'L74*NoPB_S*_th0.110_*_5.00s_58318_saxs'

# run04 C67
#pattern = 'C67_S13_th0.110_*_5.00s_*_saxs'
#pattern = 'C67_S13_th0.110_*_5.00s_58581_saxs'
#pattern = 'C67_S13_th0.110_*_5.00s_58719_saxs'

#pattern = 'C67S3M3hp30_S14_th0.110_*_5.00s_*_saxs'
#pattern = 'C67S3M3hp30_S14_th0.110_*_5.00s_58582_saxs'
#pattern = 'C67S3M3hp30_S14_th0.110_*_5.00s_58603_saxs'
#pattern = 'C67S3M3hp30_S14_th0.110_*_5.00s_58606_saxs'
#pattern = 'C67S3M3hp30_S14_th0.110_2040.1s_T200.001C_5.00s_58720_saxs'

#pattern = 'C67S3M3hp30_PGMEA_S15_th0.110_*_5.00s_*_saxs'
#pattern = 'C67S3M3hp30_PGMEA_S15_th0.110_36.1s_T90.000C_5.00s_58583_saxs'
#pattern = 'C67S3M3hp30_PGMEA_S15_th0.110_72.0s_T109.123C_5.00s_58586_saxs'
#pattern = 'C67S3M3hp30_PGMEA_S15_th0.110_138.0s_T143.741C_5.00s_58592_saxs'
#pattern = 'C67S3M3hp30_PGMEA_S15_th0.110_2053.1s_T200.003C_5.00s_58721_saxs'



# run05 C67speed
#pattern = 'C67*speed*_th0.110_*_5.00s_*_saxs'
#pattern = 'C67*speed*_th0.110_*_5.00s_59062_saxs'


# run06 L74speed
pattern = 'L74*speed*_th0.110_*_5.00s_59986_saxs'




pattern = '*'
infiles = glob.glob(os.path.join(source_dir, pattern+'.tiff'))

infiles.sort()

True
# Analysis to perform
########################################

load_args = { 'calibration' : calibration, 
             'mask' : mask,
             }
run_args = { 'verbosity' : 3,
            }

process = Protocols.ProcessorXS(load_args=load_args, run_args=run_args)



protocols = [
    #Protocols.calibration_check(show=False, AgBH=True, q0=0.010, cmap=cmap_vge_hdr, num_rings=4, ztrim=[0.0, 0.01], ) ,
    #Protocols.circular_average(ylog=True, plot_range=[0, 0.12, None, None]) ,
    #Protocols.sector_average(angle=0, dangle=20, ylog=True, plot_range=[0, 0.3, None, None], show_region=True) ,
    #Protocols.linecut_angle(q0=0.094, dq=0.015, show_region=False) ,
    #linecut_qr_fit(show_region=False, show=False, qz=0.029, dq=0.010, fit_range=[0.009, 0.024], plot_range=[0, 0.035, 0, None], xticks=[0, 0.01, 0.02, 0.03]) , # L74
    #linecut_qr_fit(show_region=False, show=False, qz=0.027, dq=0.010, fit_range=[0.008, 0.028], plot_range=[0, 0.035, 0, None], xticks=[0, 0.01, 0.02, 0.03]) , # L74 tweak
    #linecut_qr_fit(show_region=False, show=False, qz=0.029, dq=0.010, fit_range=[0.0125, 0.035], plot_range=[0, 0.035, 0, None], xticks=[0, 0.01, 0.02, 0.03]) , # C67
    #Protocols.q_image(blur=1.0, bins_relative=0.5, plot_range=[-0.05, 0.05, 0, 0.10], xticks=[-0.04, -0.02, 0, 0.02,], cmap=cmap_vge, ztrim=[0.2, 0.01]) ,
    #Protocols.qr_image(blur=None, bins_relative=0.8, plot_range=[-0.1, 3.0, 0, 3.0], _xticks=[0, 1.0, 2.0, 3.0], ztrim=[0.38, 0.002], dezing_fill=True) ,
    #Protocols.q_phi_image(bins_relative=0.25, plot_range=[0, 3.0, 0, +90]) ,
    #Protocols.thumbnails(crop=None, resize=1.0, blur=None, cmap=cmap_vge_hdr, ztrim=[0.0, 0.01]) ,
    metadata_extract(),
    ]
    



# Run
########################################
process.run(infiles, protocols, output_dir=output_dir, force=True)


# Loop
########################################
# This code is typically only used at the beamline (it loops forever, watching for new files).
import time
donefiles = []
while False:

    infiles = glob.glob(os.path.join(source_dir, '*.tiff'))

    for infile in infiles:
        if infile in donefiles:
            pass

        else:
            process.run([infile], protocols, output_dir=output_dir, force=False)

        donefiles.append(infile)

    time.sleep(4)






