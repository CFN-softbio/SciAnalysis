#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys # To get commandline arguments
import glob

SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)
from SciAnalysis import tools
#from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.ImAnalysis.Data import *
from SciAnalysis.ImAnalysis import Protocols


class thumb_view(Protocols.thumbnails):
    
    def __init__(self, name='thumb_view', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'size_nm' : 500,
                        'blur' : 1.0,
                        'resize' : 1.0,
                        'cmap' : mpl.cm.Greys_r,
                        'offset' : None,
                        }
        self.run_args.update(kwargs)
        

    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
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
    



class fft_custom(Protocols.fft):
    
    def __init__(self, name='fft_custom', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'blur' : None,
                        'Iqn_n' : 1.0,
                        'fourier_filter_invert' : False,
                        'fourier_filter_shift' : 0.3,
                        'offset' : None,
                        }
        self.run_args.update(kwargs)
        
        
    def output_exists(self, name, output_dir):
        
        outfile = self.get_outfile('{}/01_fft'.format(name), output_dir, ext='.png')
        return os.path.isfile(outfile)
        

    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        tools.make_dir(output_dir)
        
        results = {}
        
        
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(data.data) )
            outfile = self.get_outfile('original', output_dir, ext='.png', ir=True)
            im.save(outfile)
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
                q_max = run_args['q0']*2.0
            else:
                q_max = np.max(x_axis)*0.25
                
            data_fft.plot(save=outfile, ztrim=[0.3,0.00004], plot_range=[-q_max,+q_max,-q_max,+q_max], blur=run_args['blur'], plot_buffers=[0,0,0,0])
            if run_args['verbosity']>=4:
                outfile = self.get_outfile('fft_zoom_blur', output_dir, ext='.png', ir=True)
                data_fft.plot(save=outfile, ztrim=[0.5,0.001], plot_range=[-q_max,+q_max,-q_max,+q_max], blur=1.0,)
        if run_args['verbosity']>=4:
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
        
        #results['q0']['value'] = 2*np.pi/36
        
        
        if run_args['verbosity']>=2:
            lines.x_label = 'q'
            lines.x_rlabel = '$q \, (\mathrm{nm}^{-1})$'
            lines.y_label = 'I'
            lines.y_rlabel = r'$\langle I \rangle \, (\mathrm{counts/pixel})$'
            
            outfile = self.get_outfile('fft_1d', output_dir, ext='.png', ir=True)
            plot_range = [0, np.max(line.x)*0.75, np.min(line.y[1:-1]), np.max(line.y[1:-1])]
            lines.plot(save=outfile, ylog=True, plot_range=plot_range)
        
        if run_args['verbosity']>=3:
            
            class DataLine_current(DataLine):
                def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.2,0.05,0.2,0.05], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, **kwargs):
                    
                    # DataLine._plot()
                    
                    plot_args = self.plot_args.copy()
                    plot_args.update(kwargs)
                    self.process_plot_args(**plot_args)
                    
                    self.fig = plt.figure( figsize=(7,7), facecolor='white' )
                    left_buf, right_buf, bottom_buf, top_buf = plot_buffers
                    fig_width = 1.0-right_buf-left_buf
                    fig_height = 1.0-top_buf-bottom_buf
                    self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
                    
                    p_args = dict([(i, plot_args[i]) for i in self.plot_valid_keys if i in plot_args])
                    self._plot_main(error=error, error_band=error_band, dashes=dashes, **p_args)
                    
                    
                    plt.xlabel(self.x_rlabel)
                    plt.ylabel(self.y_rlabel)
                    
                    if xlog:
                        plt.semilogx()
                    if ylog:
                        plt.semilogy()
                    if xticks is not None:
                        self.ax.set_xticks(xticks)
                    if yticks is not None:
                        self.ax.set_yticks(yticks)
                    
                    
                    # Axis scaling
                    xi, xf, yi, yf = self.ax.axis()
                    if plot_range[0] != None: xi = plot_range[0]
                    if plot_range[1] != None: xf = plot_range[1]
                    if plot_range[2] != None: yi = plot_range[2]
                    if plot_range[3] != None: yf = plot_range[3]
                    self.ax.axis( [xi, xf, yi, yf] )
                    
                    self._plot_extra(**plot_args)
                    
                    s = '${:.1f} \, \mathrm{{nm}}$'.format( 2.*np.pi/(self.q0) )
                    self.ax.text(xi, yf, s, size=50, color='b', verticalalignment='top', horizontalalignment='left')
                    
                    self.ax.axvline(self.q0, color='b', linewidth=4)
                    self.ax.text(self.q0, yf, '$1$', size=18, color='b', verticalalignment='top', horizontalalignment='left')
                    
                    self.ax.axvline(np.sqrt(3.0)*self.q0, color='b', linewidth=2, dashes=[5,5])
                    self.ax.text(np.sqrt(3.0)*self.q0, yf, '$\sqrt{3}$', size=18, color='b', verticalalignment='top', horizontalalignment='left')
                    
                    self.ax.axvline(2.0*self.q0, color='b', linewidth=2)
                    self.ax.text(2*self.q0, yf, '$2$', size=18, color='b', verticalalignment='top', horizontalalignment='left')

                    self.ax.axvline(np.sqrt(7.0)*self.q0, color='b', linewidth=2, dashes=[5,5])
                    self.ax.text(np.sqrt(7.0)*self.q0, yf, '$\sqrt{7}$', size=18, color='b', verticalalignment='top', horizontalalignment='left')
                    
                    
                    
                    if save:
                        if 'dpi' in plot_args:
                            plt.savefig(save, dpi=plot_args['dpi'])
                        else:
                            plt.savefig(save)
                    
                    if show:
                        self._plot_interact()
                        plt.show()
                        
                    plt.close(self.fig.number)
                                    
            
            # I*q^2 vs. q (Kratky plot style)
            line_Iq2 = DataLine_current()
            line_Iq2.x = line.x
            line_Iq2.y = line.y*np.square(line_Iq2.x)
            line_Iq2.y_label = 'Iq^2'
            line_Iq2.y_rlabel = r'$q^2 \langle I \rangle$'
            line_Iq2.plot_args['linestyle'] = '-'
            line_Iq2.plot_args['linewidth'] = 4.0
            line_Iq2.plot_args['marker'] = None


            #line_Iq2.q0 = run_args['q0']
            line_Iq2.q0 = results['q0']['value']
            if 'q0' in run_args:
                q_max = np.sqrt(2.0)*run_args['q0']*2.0
            else:
                q_max = np.max(line_Iq2.x)*0.5

            q_max = np.sqrt(2.0)*( 2*np.pi/36 )*2.0

            
            line_Iq2.trim(0, q_max)
            outfile = self.get_outfile('fft_1d_Iq2', output_dir, ext='.png', ir=True)
            plot_range = [0, np.max(line_Iq2.x), 0, None]
            line_Iq2.plot(save=outfile, plot_range=plot_range, plot_buffers=[0,0,0,0])        
            
        if run_args['verbosity']>=3:
            # Orientation analysis at q0
            q0 = results['q0']['value']
            q_spread = results['sigma_q0']['value']*3.0
            new_results = self.orientation_q0(data_fft, q0, q_spread, output_dir, **run_args)
            results.update(new_results)
        
        if run_args['verbosity']>=3:
            # Fourier-filtered image at q0
            q_spread = results['sigma_q0']['value']*5.0
            self.fourier_filter(data, q0, q_spread, output_dir, **run_args)            
            
            
        
        return results        


    def fourier_filter(self, data, q_center, q_spread, output_dir, **run_args):
        
        data.fourier_filter(q_center, q_spread)
        
        
        data.maximize_intensity_spread()
        if run_args['fourier_filter_invert']:
            data.invert()
        if run_args['fourier_filter_shift'] != 0.0:
            data.data = np.clip(data.data, 255*run_args['fourier_filter_shift'], 255)
            data.maximize_intensity_spread()
        
        
        # Custom crop
        height, width = data.data.shape
        
        w_pix = run_args['size_nm']/data.x_scale
        h_pix = run_args['size_nm']/data.y_scale
        
        left = (width-w_pix)/2
        right = left
        bottom = (height-h_pix)/2
        top = bottom
        
        if run_args['offset'] is not None:
            left += offset[0]
            right -= offset[0]
            top += offset[0]
            bottom -= offset[0]
        
        data.crop_edges(left=left, right=right, bottom=bottom, top=top, relative=False)
        
        
        if run_args['verbosity']>=2:
            data.set_z_display( [None, None, 'gamma', 1.0] )
            outfile = self.get_outfile('filtered', output_dir, ext='.png', ir=True)
            data.plot_image(save=outfile, zmin=0, zmax=255, cmap=mpl.cm.bone)    


#L0 = 27 # nm 'C48'
#L0 = 32 # nm 'C48'
#L0 = 34 # nm 'C67'
#L0 = 38.5 # nm 'L75'
#L0 = 43.6 # nm 'C99'
#L0 = 76 # nm 'C177'
#L0 = 79 # nm 'O184'
#L0 = 65 # nm 'L176'
#L0 = 128 # nm 'L570'

#L0 = 32 # nm 'SEO30'
#L0 = 30 # nm 'S2VP45'


# layering distance
#L0 = 34 # nm 'C67'
L0 = 36 # nm slightly distorted...
#L0 = 40 # nm distorted

q0 = 2*np.pi/(L0)




#load_args['scale'] = 3000.0/544 # nm/pixel (18k MAG)
#load_args['scale'] = 2000.0/454 # nm/pixel (25k MAG)
#load_args['scale'] = 1000.0/504 # nm/pixel (50k MAG)
#load_args['scale'] = 500.0/302 # nm/pixel (60k MAG)
#load_args['scale'] = 500.0/454 # nm/pixel (90k MAG)
#load_args['scale'] = 500.0/504 # nm/pixel (100k MAG)
#load_args['scale'] = 300.0/454 # nm/pixel (150k MAG)


process = Protocols.ProcessorIm()

run_args = { 'verbosity' : 3,
                'local_partition_image_size' : 30, # pixels
                'local_partition_step' : 5.0, # relative to image_size
                #'blur' : 1.0, # pixels
                'q0' : q0, # nm^-1
                'dq' : q0*0.6, # nm^-1
                'symmetry' : 2,
                #'rcParams': {'axes.labelsize': 45,},
                'NN_cutoff_distance_nm' : L0*1.36,
                'area_min' : L0*0.35,
                'correlation_edge_exclusion' : 20,
                'correlation_step_size_points' : 30,
                
                }

load_args = { 'format' : 'custom',
                #'scale' : 1000.0/1008, # nm/pixel
                'crop_edges' : [0, 0, 248, 0],
                #'crop_edges' : [0, 0, 124, 0],
                'load_param_file' : True,
                }
            

size_nm = 1000
thumb = thumb_view(name='thumb_view{:d}'.format(size_nm), size_nm=size_nm)
protocols = [ 
                Protocols.fft(blur=0.6, **run_args),
                #fft_custom(blur=0.8, size_nm=size_nm, **run_args),
                Protocols.grain_size(**run_args),
                #Protocols.thumbnails(resize=0.5, crop=0.5),
                #thumb ,
                ]






source_dir = './'
output_dir = './analysis/'
pattern = 'L75 on Dow 55_q02.tif'
tools.make_dir(output_dir)

infiles = glob.glob(source_dir + '/'+pattern)
infiles.sort()

print('{} infiles'.format(len(infiles)))
process.run(infiles, protocols, output_dir=output_dir, force=True, load_args=load_args, run_args=run_args)
