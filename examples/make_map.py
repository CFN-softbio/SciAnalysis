#!/usr/bin/python3
# -*- coding: utf-8 -*-


# This generates a map (2D false-color plot or 3D height plot) for a set of
# experiments (that are presumptively defined in some (x,y) space). The code
# assumes you've already used SciAnalysis to process your data, such that you
# have XML files in your "results" sub-folder with the analysis results of
# interest. This code then compiles that data and generates the plot.

# The code can also be used to generate an animation of the sequence of
# measurements during the experiment.


# Imports
########################################

import sys, os
#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH='/nsls2/xf11bm/software/SciAnalysis/'
SciAnalysis_PATH = '/home/etsai/BNL/Users/software/SciAnalysis'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
import numpy as np
import re

from SciAnalysis import tools
from SciAnalysis.Data import *
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols


# Helpers
########################################

def status_online(copy_from, copy_to='/home/kyager/software/statpage/current.png'):
    import shutil
    shutil.copyfile(copy_from, copy_to)



#AE_sample1_anneal1200_run1_x46.244_yy-39.304_SamX39.30_SamY46.24_10.00s_58036_saxs.tiff
#filename_re = re.compile('(.+)_x(-?\d+\.\d+)_yy(-?\d+\.\d+)_SamX(-?\d+\.\d+)_SamY(-?\d+\.\d+)_.+_(\d+)_saxs\.xml')
filename_re = re.compile('(.+)_x(-?\d+\.\d+)_y(-?\d+\.\d+)_.+_(\d+)_waxs\.xml')
def parse_filename(filename, verbosity=3):
    
    parsed = {'filename' : filename}
    m = filename_re.match(filename)
    if m:
        parsed['basename'] = m.groups()[0]
        parsed['x'] = float(m.groups()[1])
        parsed['y'] = float(m.groups()[2])

        #parsed['samx'] = float(m.groups()[3])
        #parsed['samy'] = float(m.groups()[4])

        parsed['scan_id'] = int(m.groups()[-1])
        
    else:
        if verbosity>=2:
            print("WARNING: RE doesn't match for {}".format(filename))
    
    return parsed
    

def val_stats(values, name='z'):
    span = np.max(values)-np.min(values)
    print("  {} = {:.2g} ± {:.2g} (span {:.2g}, from {:.2g} to {:.2g})".format(name, np.average(values), np.std(values), span, np.min(values), np.max(values)))
    
    
def power_N_list(N_max, N_min=3, num=40, exponent=3.0):
    '''Generates a list of integers that are spread out more logarithmically.
    That is, the list starts with small increments between integers, but finishes
    with large increments.'''
    
    #N_list = ( (np.exp( np.linspace(0, 1, num=40) ) - 1)/(np.exp(1)-1) )*len(z_vals)
    x = np.linspace(0, 1, num=num)
    N_list = np.power(x, exponent)*(N_max-N_min) + N_min
    N_list = np.unique(N_list.astype(int))
    #N_list = N_list[ (N_list>=N_min) & (N_list<=N_max) ] # Unnecessary
    
    return N_list
    


# Extract results from xml files
########################################
from SciAnalysis.Result import * # Results() object
def extract_results(infiles, extractions, outfile, verbosity=3):
    if verbosity>=3:
        print("Extracting results for {} infiles...".format(len(infiles)))
    
    results = Results().extract_multi_save_txt(outfile, infiles, extractions, verbosity=verbosity)
    
    return results




# Results
########################################

def load_file(infile, verbosity=3):
    if verbosity>=3:
        print(" Loading data from file: {}".format(infile))
    
    with open(infile, 'r') as fin:
        names = fin.readline().split()
        lines = fin.readlines()
        
    if verbosity>=4:
        print('Saved data has {} columns:'.format(len(names)))
        print(names)
        
    return names, lines

def load_results(names, lines, x_coord, y_coord, z_signal, sequence, verbosity=3):
    
    # Load results
    x_vals = []
    y_vals = []
    z_vals = []
    seq_vals = []
    
    signal_idx = names.index(z_signal)
    x_idx = names.index(x_coord)
    y_idx = names.index(y_coord)
    sequence_idx = names.index(sequence)
    
    if verbosity>=3:
        print(" Plot signal: {} (column index {})".format(z_signal, signal_idx))
        print(" Plot x: {} (column index {})".format(x_coord, x_idx))
        print(" Plot y: {} (column index {})".format(y_coord, y_idx))
        print(" Sorting: {} (column index {})".format(sequence, sequence_idx))
        
    
    for line in lines:
        els = line.split()
        if len(els)==len(names) and els[0][0]!='#' and els[signal_idx]!='-':
        #if len(els)==len(names) and els[0][0]!='#' and els[signal_idx]!='-' and float(els[sequence_idx])<=34809: # For sample4 (432s)

            parsed = parse_filename(els[0])


            #x_vals.append(float(els[x_idx]))
            #y_vals.append(float(els[y_idx]))

            x_vals.append(float( parsed['x'] ))
            y_vals.append(float( parsed['y'] ))


            z_vals.append(float(els[signal_idx]))
            seq_vals.append(int(float(els[sequence_idx])))
        else:
            if verbosity>=3:
                print('  Skipping line: {}'.format(line.strip()))
                
        
    if verbosity>=4:
        print(" Nx, Ny, Nz = {}, {}, {}".format(len(x_vals), len(y_vals), len(z_vals)))
        val_stats(z_vals, name='z')

        
    # Sort data
    indices = np.argsort(seq_vals)
    x_vals = np.asarray(x_vals)[indices]
    y_vals = np.asarray(y_vals)[indices]
    z_vals = np.asarray(z_vals)[indices]
    seq_vals = np.asarray(seq_vals)[indices]
    
    return x_vals, y_vals, z_vals, seq_vals


def trim_vals(vals_list, N_max, verbosity=3):
    
    trimmed_vals_list = []
    for i, vals in enumerate(vals_list):
        if N_max is not None and len(vals)>N_max:
            if verbosity>=4:
                print(' Reducing list {} to size N_max = {}'.format(i, N_max))
            vals = np.asarray(vals)[:N_max]
        trimmed_vals_list.append(vals)
        
    return trimmed_vals_list
    
    
# Plot
########################################
    
import SciAnalysis.colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.set_cmap(cmaps.viridis)    

class Data2D_current(Data2D):
    def _plot_extra(self, **plot_args):
        
        xi, xf, yi, yf = self.ax.axis()
        
        # Faded overlay
        rect = mpl.patches.Rectangle((xi,yi), xf-xi, yf-yi, linewidth=1, edgecolor='none', facecolor='white', alpha=self.alpha, zorder=10)
        #self.ax.add_patch(rect)            
        
        
        # Scatterplot
        cmap = plot_args['cmap'] if 'cmap' in plot_args else 'viridis'
        if plot_args['zmin'] is None:
            zmin = np.quantile(self.z_vals, 0.1)
        else:
            zmin = plot_args['zmin']
            
        if plot_args['zmin'] is None:
            zmax = np.quantile(self.z_vals, 0.9)
        else:
            zmax = plot_args['zmax']
        print('zmin, zmax = [{}, {}]'.format(zmin, zmax))
        #self.ax.scatter(self.x_vals, self.y_vals, s=100, c=self.z_vals, cmap=cmap, vmin=zmin, vmax=zmax, edgecolor='k', zorder=100)
                
        # Colorbar
        n = 5
        colorbar_labels = [ zmin + i*(zmax-zmin)/(n-1) for i in range(n) ]
        
        tick_positions = self._plot_z_transform(data=colorbar_labels, set_Z=True)
        cbar = self.fig.colorbar(self.im, ax=self.ax, ticks=tick_positions, fraction=0.03, pad=0.02)
        colorbar_labels = ["{:.3g}".format(c) for c in colorbar_labels]
        cbar.ax.set_yticklabels(colorbar_labels, size=17)
        
        # Titles
        size = plot_args['rcParams']['axes.labelsize']
        #size = plot_args['rcParams']['xtick.labelsize']
        plt.figtext(0, 1, '\n\nN = {:,d}'.format(len(self.z_vals)), size=size-20, weight='bold', verticalalignment='top', horizontalalignment='left')
        z_label = self.z_rlabel if self.z_rlabel is not None else self.z_label
        if z_label is not None:
            plt.figtext(1, 1, z_label, size=size, weight='bold', verticalalignment='top', horizontalalignment='right')
            
        # self.ax.set_aspect('equal', 'datalim')
        self.ax.set_aspect('equal')
        #self.ax.set_adjustable("datalim") 
        
        
    def _plot_extra3D(self, **plot_args):
        # Colorbar
        cbar = self.fig.colorbar(self.surf, ax=self.ax, aspect=40, fraction=0.02, pad=0.0)
        cbar.ax.yaxis.set_tick_params(labelsize=15)

        # Titles
        size = plot_args['rcParams']['axes.labelsize']
        #size = plot_args['rcParams']['xtick.labelsize']
        plt.figtext(0, 1, '$N = {:,d}$'.format(len(self.z_vals)), size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
        z_label = self.z_rlabel if self.z_rlabel is not None else self.z_label
        if z_label is not None:
            plt.figtext(1, 1, z_label, size=size, weight='bold', verticalalignment='top', horizontalalignment='right')
    
    
    
def plot_results(x_vals, y_vals, z_vals, outfile, plot2d=True, plot3d=False, grid=None, dgrid=None, title=None, cmap='viridis', alpha=0.5, x_label='$x$', y_label='$y$', z_label=None, interpolate_cyclic=None, verbosity=3):
    
    title = pattern 
    
    if verbosity>=3:
        print("Plotting for {}".format(outfile))
        val_stats(z_vals, name='z_vals')
    
    # Define grid for interpolation
    if grid is None:
        grid = [np.min(x_vals), np.max(x_vals), np.min(y_vals), np.max(y_vals)]
    if dgrid is None:
        dgrid = [ (grid[1]-grid[0])/200 , (grid[3]-grid[2])/200 ]
    if isinstance(dgrid, float):
        dgrid = [dgrid, dgrid]
        
    xi = np.arange(grid[0], grid[1]+dgrid[0], dgrid[0])
    yi = np.arange(grid[2], grid[3]+dgrid[1], dgrid[1])
    XI, YI = np.meshgrid(xi, yi)


    # Interpolation
    import scipy.interpolate
    POINTS = np.column_stack((x_vals,y_vals))
    VALUES = z_vals
    if verbosity>=4:
        print("Interpolating {:,} points to {:,}×{:,} = {:,} points".format(len(VALUES), len(xi), len(yi), len(xi)*len(yi)))    
        
        
    if interpolate_cyclic is not None:
        xlike = np.cos(VALUES*2*np.pi/interpolate_cyclic)
        ylike = np.sin(VALUES*2*np.pi/interpolate_cyclic)
        XLIKE = scipy.interpolate.griddata(POINTS, xlike, (XI, YI), method='linear')
        YLIKE = scipy.interpolate.griddata(POINTS, ylike, (XI, YI), method='linear')
        
        ZI = ( np.arctan2(YLIKE, XLIKE)/(2*np.pi) )*interpolate_cyclic
        ZI_mask = np.ma.masked_where( np.isnan(ZI), ZI)
        
    else:        
        ZI = scipy.interpolate.griddata(POINTS, VALUES, (XI, YI), rescale=True, method='linear') # method='nearest' 'linear' 'cubic'
        ZI_mask = np.ma.masked_where( np.isnan(ZI), ZI)

    if verbosity>=4:
        val_stats(ZI, name='ZI')
        val_stats(ZI_mask, name='ZI_mask')
    
    d = Data2D_current()
    d.data = ZI_mask
    d.x_axis = xi
    d.y_axis = yi
    
    d.x_vals = x_vals
    d.y_vals = y_vals
    d.z_vals = z_vals


    d.set_z_display([None, None, 'linear', 1.0])
    d.x_rlabel = x_label
    d.y_rlabel = y_label
    d.z_rlabel = z_label
    
    if plot2d:
        d.plot_args['rcParams'] = { 
                        'axes.labelsize': 45,
                        'xtick.labelsize': 40,
                        'ytick.labelsize': 40,  
                        }
        d.alpha = alpha

        
        d.plot(save=outfile, show=False, cmap=cmap, zmin=zmin, zmax=zmax, title=title, plot_buffers=[0.21, 0.12, 0.18, 0.10], plot_range=grid, plot_2D_type='pcolormesh', dpi=150, transparent=False)
        d.save_data(outfile[0:-4])
    
    
    if plot3d:
        d.plot_args['rcParams'] = { 
                        'axes.labelsize': 45,
                        'xtick.labelsize': 20,
                        'ytick.labelsize': 20,    
                        }    
        d.X = XI
        d.Y = YI
        
        elev = 30
        azim = 30
        azim = 30-90
        
        #outfile = outfile[:-4]+'-3D.png'
        outfile = tools.Filename(outfile).path_append('3D')
        d.plot3D(save=outfile, show=False, cmap=cmap, zmin=zmin, zmax=zmax, title=title, plot_buffers=[0.05, 0.10, 0.05, 0.05], plot_range=grid, elev=elev, azim=azim, dpi=150, transparent=False)
 

def plot_grids(results, N_list, grid=[None, None, None, None], grid_size=None, n_grid=200, plot2d=True, plot3d=False, cmap='viridis', x_label='$x \, (\mathrm{mm})$', y_label='$y \, (\mathrm{mm})$', z_label='$z$', interpolate_cyclic=None, verbosity=3):
    
    x_vals, y_vals, z_vals, seq_vals = results
    
    if grid[0] is None: grid[0] = np.min(x_vals)
    if grid[1] is None: grid[1] = np.max(x_vals)
    if grid[2] is None: grid[2] = np.min(y_vals)
    if grid[3] is None: grid[3] = np.max(y_vals)
    
    if grid_size is None:
        n_grid = [n_grid, n_grid]
    else:
        n1 = int((grid[1]-grid[0])/grid_size)
        n2 = int((grid[3]-grid[2])/grid_size)
        n_grid = [n1, n2]
    
    dgrid = [ (grid[1]-grid[0])/n_grid[0] , (grid[3]-grid[2])/n_grid[1] ]
    print("n_grid = {}, dgrid = {}".format(n_grid, dgrid))
    
    
    #N_list = [len(z_vals)] # Plot final value only
    
    for N in N_list:
        x_vals, y_vals, z_vals, seq_vals = trim_vals(results, N_max=N, verbosity=verbosity)
        if verbosity>=4:
            val_stats(z_vals, name='z_reduced')    
            
        outfile = os.path.join(output_dir, z_signal, '{}-{}-final.png'.format(pattern, z_signal, N))
        plot_results(x_vals, y_vals, z_vals, outfile=outfile, plot2d=plot2d, plot3d=plot3d, grid=grid, dgrid=dgrid, cmap=cmap, alpha=0.2, x_label=x_label, y_label=y_label, z_label=z_label, interpolate_cyclic=interpolate_cyclic, verbosity=verbosity)

        #if protocol=='circular_average_q2I_fit' and signal=='fit_peaks_d0': # TOCHANGE
            #status_online(outfile)
    
def save_grids(results):

    x_vals, y_vals, z_vals, seq_vals = results  

    grid[0] = np.min(x_vals)
    grid[1] = np.max(x_vals)
    grid[2] = np.min(y_vals)
    grid[3] = np.max(y_vals)
    
    grid_size0 = np.median(np.abs(np.diff(x_vals)))  
    grid_size1 = np.max(np.abs(np.diff(y_vals)))
    dgrid = [grid_size0, grid_size1]
    
    xi = np.arange(grid[0], grid[1]+dgrid[0], dgrid[0])
    yi = np.arange(grid[2], grid[3]+dgrid[1], dgrid[1])
    XI, YI = np.meshgrid(xi, yi)

    # Interpolation
    import scipy.interpolate
    POINTS = np.column_stack((x_vals,y_vals))
    VALUES = z_vals

    print("Interpolating {:,} points to {:,}×{:,} = {:,} points".format(len(VALUES), len(xi), len(yi), len(xi)*len(yi)))        
    
    ZI = scipy.interpolate.griddata(POINTS, VALUES, (XI, YI), rescale=True, method='linear') # method='nearest' 'linear' 'cubic'
    ZI_mask = np.ma.masked_where( np.isnan(ZI), ZI)   

    ### Save file  
    fn_tiff = '{}{}_{}.tiff'.format(output_dir, pattern, z_signal.replace('_fit_peaks_',''))
    img = PIL.Image.fromarray(np.float64(ZI))
    img.save(fn_tiff) 
    print('File saved to {}'.format(fn_tiff))
    
    
    
    
    

def animated_gif(source_dir='./', pattern='*.png', outfile=None, skip=None, verbosity=3):
    
    # If you get a 'cache resources exhausted' error, you can increase the cache sizes:
    # sudo nano /etc/ImageMagick-6/policy.xml    
    
    if verbosity>=3:
        print('Generating animated gif for {}/{}'.format(source_dir, pattern))
        
    # Select the files to animate
    infiles = glob.glob(os.path.join(source_dir, pattern))
    if verbosity>=3:
        print('    {} infiles'.format(len(infiles)))
        
        
    infiles.sort()
        
    if skip is not None:
        infiles = infiles[0::skip]
        
        
    if outfile is None:
        outfile = os.path.join(source_dir, 'anim.gif')
    elif outfile[-4:]!='.gif':
        outfile = outfile+'.gif'

    # Prepare command
    # (Animation is generated using imagemagick 'convert' bash command.)

    #cmd = "convert -delay 20 -resize 50% -fill white  -undercolor '#00000080'  -gravity NorthWest -annotate +0+5 ' Text ' "
    #cmd = "convert -delay 15 -loop 1 -resize 50% "
    #cmd = "convert -crop 450x450+60+220  +repage -delay 15 -loop 1 -resize 50% "
    cmd = "convert -dispose previous +repage -delay 30 -loop 1 -resize 30% "
            
    for infile in infiles:
        cmd += '{} '.format(infile)
    
    cmd += ' {}'.format(outfile)

    # Execute command
    print('  Saving {}'.format(outfile))
    os.system(cmd)
    
    # Add a white background
    #cmd = 'convert {} -coalesce -background white -alpha remove {}'.format(outfile, outfile[:-4]+'w.gif')
    #os.system(cmd)
    
    
def coord_transform(results, x_label, y_label):
    #x_vals, y_vals, z_vals, seq_vals = results
    #x_vals = (x_vals-2)*1000
    #y_vals = y_vals*1000
    #results = x_vals, y_vals, z_vals, seq_vals
    
    return results, x_label, y_label


def _coord_transform(results, x_label, y_label):
    # For sample 3 (1200s)

    x_vals, y_vals, z_vals, seq_vals = results
    
    x_vals = 50-x_vals
    y_vals = 50-y_vals
    
    results = x_vals, y_vals, z_vals, seq_vals

    return results, x_label, y_label

def _coord_transform(results, x_label, y_label):
    x_vals, y_vals, z_vals, seq_vals = results
    
    x_vals = (x_vals/50)*(200-140) + 140 # thickness
    y_vals = ((y_vals+50)/50)*(200-30) + 30 # Temperature
    
    results = x_vals, y_vals, z_vals, seq_vals
    return results, '$T \, (\mathrm{^{\circ}C})$', '$h \, (\mathrm{nm})$'    
    
    
# Run
########################################
 
# Settings
########################################
verbosity = 2
pattern = 'FN_B2_run3b'
#pattern = 'FN_B2_run6*_18*'
#pattern = 'FN_B2_run*_19*' #run6 (some) and run7

source_dir = '../B2_run3b' # The location of the SciAnalysis outputs
output_dir = './{}_offline/'.format(pattern)
tools.make_dir(output_dir)
       

# Extract results from XML files
results_dir = source_dir + '/results/' # Location of xml files
infiles = glob.glob(os.path.join(results_dir, '{}*.xml'.format(pattern)))
outfile = os.path.join(output_dir, '{}-extracted.txt'.format(pattern))

protocals = ['circular_fit_121', 'circular_fit_100', 'circular_fit_222', 'circular_fit_002','circular_fit_310']
protocals = ['circular_fit_310']


for protocal1 in protocals:
    print('============ {} ============'.format(protocal1))
    if '121' in protocal1: 
        num_curves = 3
    elif '310' in protocal1: 
        num_curves = 2
    else:
        num_curves = 1    
        
    names = []
    for nn in np.arange(1, num_curves+1):
        names.append('fit_peaks_prefactor{}'.format(nn))
        names.append('fit_peaks_x_center{}'.format(nn))
        names.append('fit_peaks_d0{}'.format(nn)) 
        names.append('fit_peaks_grain_size{}'.format(nn))
        names.append('fit_peaks_grain_size_G{}'.format(nn))
        names.append('fit_peaks_strain{}'.format(nn))         
        names.append('fit_peaks_gamma{}'.format(nn))
        names.append('fit_peaks_sigma{}'.format(nn))        
    names.append('fit_peaks_chi_squared')         
           
    extractions = [ 
            [ 'metadata_extract', ['x_position', 'y_position', 'sequence_ID'] ],
            [protocal1, names ]
            ]           
    
    
    results = extract_results(infiles, extractions, outfile=outfile, verbosity=verbosity)
    
    
    # Plot results
    x_coord, x_label = 'metadata_extract__x_position', '$x \, (\mathrm{mm})$'
    y_coord, y_label = 'metadata_extract__y_position', '$y \, (\mathrm{mm})$'
    sequence = 'metadata_extract__sequence_ID'    

  
    signals = []
    for nn in np.arange(1, num_curves+1):           
        signals.append(['{}__fit_peaks_prefactor{}'.format(protocal1,nn), 'p1', None, None, 'viridis'])
        signals.append([ '{}__fit_peaks_x_center{}'.format(protocal1,nn), '$d_0 \, (\mathrm{nm})$', None, None, 'inferno' ])
        signals.append([ '{}__fit_peaks_d0{}'.format(protocal1,nn), '$d_0 \, (\mathrm{nm})$', None, None, 'inferno' ])
        signals.append(['{}__fit_peaks_grain_size{}'.format(protocal1,nn), r'$\xi \, (\mathrm{nm})$', None, None, 'plasma' ])
        signals.append(['{}__fit_peaks_grain_size_G{}'.format(protocal1,nn), r'$\xi \, (\mathrm{nm})$', None, None, 'plasma' ])
        signals.append(['{}__fit_peaks_strain{}'.format(protocal1,nn), r'$\xi \, (\mathrm{nm})$', None, None, 'viridis' ])
        signals.append(['{}__fit_peaks_gamma{}'.format(protocal1,nn), '', None, None, 'viridis'])
        signals.append(['{}__fit_peaks_sigma{}'.format(protocal1,nn), '', None, None, 'viridis'])
    
    signals.append(['{}__fit_peaks_chi_squared'.format(protocal1), r'$\xi \, (\mathrm{nm})$', None, None, 'viridis' ])       
        


# Transform coordinates


    names, lines = load_file(outfile, verbosity=verbosity)
    print(names)
    
    
    z_vals_all = []
    flip_y_axis = True
    flag_mask = True
    
    for ii, signal in enumerate(signals):
        if verbosity>=3:
            print("========================================")
            print("Plotting signal {}...".format(signal))
            print("========================================")
    
        if ii>0:
            mask = z_vals_all[0]>0.05 # Make mask with p1 or chi?
        else: 
            mask = 1
        
        if ii<len(signals): #%-1:
            z_signal, z_label, zmin, zmax, cmap = signal
            protocol, signal = z_signal.split('__')
    
            grid = [-1.5, 2.2, -2, -0.2]
            #grid = [50, 0, None, None]
            #grid = [None, None, None, None]
            print(signal)
            tools.make_dir(os.path.join(output_dir,z_signal))
            #tools.make_dir(os.path.join(output_dir,signal,'3D'))
            
            results = load_results(names, lines, x_coord=x_coord, y_coord=y_coord, z_signal=z_signal, sequence=sequence, verbosity=verbosity)
            results, x_label, y_label = coord_transform(results, x_label, y_label)
            x_vals, y_vals, z_vals, seq_vals = results
            
            
            if flag_mask:
                z_vals = z_vals*mask
     
            if flip_y_axis: 
                y_vals = -y_vals
                #print('!!!--- Adjusting stage shift')
                #x_vals[y_vals>-1.0] = x_vals[y_vals>-1.0] + 0.356
                # x_vals[53]-x_vals[54*15+53]
                
                results = x_vals, y_vals, z_vals, seq_vals         
             
            #zmax = np.median(z_vals)*15
            z_vals_all.append(z_vals)
        
        else:
            z_signal, z_label, zmin, zmax, cmap = signal
            protocol, signal = z_signal.split('__')        
            tools.make_dir(os.path.join(output_dir,signal))
            
            mask = z_vals_all[1].copy()
            mask[mask<5] = np.nan
            mask[mask>=5] = 1.0
            z_vals_calc = z_vals_all[0] / z_vals_all[1] #* mask
            #zmin = np.nanmin(z_vals_calc)
            #zmax = np.nanmax(z_vals_calc)
            
            results = x_vals, y_vals, z_vals_calc, seq_vals            
            
            
        # Single result
        if True:
            plot_grids(results, [len(z_vals)], grid=grid, grid_size=60e-3, n_grid=300, plot2d=True, plot3d=False, cmap=cmap, x_label=x_label, y_label=y_label, z_label=z_label, verbosity=verbosity)
            
            save_grids(results)
            
            
            #plot_grids(results, [len(z_vals)], grid=grid, n_grid=40, plot2d=False, plot3d=True, cmap=cmap, x_label=x_label, y_label=y_label, z_label=z_label, verbosity=verbosity)
    
        '''
        # Sequence of results
        if False:
            # 2D plots
            N_spacing = 50
            N_list = np.arange(N_spacing, len(z_vals), N_spacing)
            N_list = power_N_list(len(z_vals), num=140, exponent=5.0)
            plot_grids(results, N_list, grid=grid, n_grid=200, plot2d=True, plot3d=False, cmap=cmap, x_label=x_label, y_label=y_label, z_label=z_label, verbosity=verbosity)
            
            # 3D plots
            plot_grids(results, N_list, grid=grid, n_grid=40, plot2d=False, plot3d=True, cmap=cmap, x_label=x_label, y_label=y_label, z_label=z_label, verbosity=verbosity)
            
            
        # Animated gif
        if False:
            outfile = os.path.join(output_dir, '{}-{}.gif'.format(pattern, signal))
            animated_gif(source_dir=os.path.join(output_dir, signal), outfile=outfile, verbosity=verbosity)
    
            outfile = os.path.join(output_dir, '{}-{}-3D.gif'.format(pattern, signal))
            animated_gif(source_dir=os.path.join(output_dir, signal, '3D'), outfile=outfile, verbosity=verbosity)
        '''










        
