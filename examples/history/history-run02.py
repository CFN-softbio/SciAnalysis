#!/usr/bin/python

#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.Data import *
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols


import re




source_dir = '../linecut_qr/'
#source_dir = '../circular_average/'

source_dir_T = '../temperature_curve/'


run_code = 'run04'

sample_code = 'run04_kinetics-120C_150nm-chip1'
sample_code = 'run04_kinetics-120C_65nm-chip1'
sample_code = 'run04_kinetics-120C_65nm-chip2'
sample_code = 'run04_kinetics-120C_15nm-chip1'
sample_code = 'run04_kinetics-120C_15nm-chip2'


run_code = 'run01_linear-ramp'

sample_code = 'run01_linear-ramp_150nm-chip1'
sample_code = 'run01_linear-ramp_150nm-chip2'
sample_code = 'run01_linear-ramp_65nm-chip1'
sample_code = 'run01_linear-ramp_65nm-chip2'
sample_code = 'run01_linear-ramp_15nm-chip1'
sample_code = 'run01_linear-ramp_15nm-chip2'


#run_code = 'run02_kinetics-100C'
#max_time = 10.1
run_code = 'run03_kinetics-150C'
max_time = 2
run_code = 'run05_kinetics-60C'
max_time = 10

sample_code = '{}_150nm-chip1'.format(run_code)
sample_code = '{}_150nm-chip2'.format(run_code)
sample_code = '{}_65nm-chip1'.format(run_code)
sample_code = '{}_65nm-chip2'.format(run_code)
sample_code = '{}_15nm-chip1'.format(run_code)
sample_code = '{}_15nm-chip2'.format(run_code)



file_re = re.compile('^(.+)_th.+_(\d+\.\d)s_T(\d+.\d\d\d)C_.+s_(\d+)_saxs.dat')

# Load time series
infile = ''.join([source_dir_T, run_code, '.dat'])
filenames = []
sIDs = []
times = []
temperatures = []
with open(infile) as fin:
    for line in fin.readlines():
        if line[0] != '#':
            filename, sID, time, temperature = line.split()
            filenames.append(filename)
            sIDs.append(int(sID))
            times.append(float(time))
            temperatures.append(float(temperature))






tT_sIDs = []
tT_times = []
tT_temperatures = []
with open(''.join((source_dir_T, run_code, '.dat'))) as fin:
    for line in fin.readlines():
        if line[0]!='#':
            els = line.strip().split()
            if len(els)>1:
                tT_sIDs.append(int(els[1]))
                tT_times.append(float(els[2]))
                tT_temperatures.append(float(els[3]))
    
def lookup_clock(sID):
    idx = np.where( tT_sIDs==sID )[0][0]
    #print(sID, idx, tT_sIDs[idx], tT_times[idx])
    return tT_times[idx]


def normalize_linecut(linecut, q_normalize=-0.18):
    
    idx = (np.abs(linecut[:,0]-q_normalize)).argmin()    
    val = linecut[idx,1]
    
    linecut[:,1] /= val
    
    return linecut



infiles = glob.glob( ''.join([source_dir, sample_code, '*_saxs.dat']) )

infiles_use = []
infile_sIDs = []
for infile in infiles:
    m = file_re.match(infile)
    if m:
        sample, clock, temperature, sID = m.groups()
        infiles_use.append(infile)
        infile_sIDs.append(int(sID))
    
    else:
        print("Warning: RE didn't match {}".format(infile))



infile_sIDs = np.asarray(infile_sIDs)
infiles = np.asarray(infiles_use)
idx = np.argsort(infile_sIDs)
infile_sIDs = infile_sIDs[idx]
infiles = infiles[idx]


# Prepare matrix for data

points = None
values = np.asarray([])



# Load individual linecuts
for i, (infile, sID) in enumerate(zip(infiles, infile_sIDs)):

    m = file_re.match(infile)
    if m:
        sample_current, clock, temperature, sID_current = m.groups()
        clock = float(clock)
        temperature = float(temperature)
        sID_current = int(sID_current)
        sample_current = sample_current[len(source_dir):]
        
        if sID!=sID_current:
            print("ERROR: sIDs don't match: {} != {}".format(sID, sID_current))
        
        if sample_current==sample_code:

            linecut = np.loadtxt(infile)
            
            #linecut = normalize_linecut(linecut)

            qs = linecut[:,0]
            times = np.ones(len(qs))*lookup_clock(sID)
            new_points = np.column_stack((qs,times))
            
            new_values = linecut[:,1]
            
            if points is None:
                points = new_points
            else:
                points = np.concatenate( (points, new_points) )
            values = np.concatenate( (values, new_values) )
            
        else:
            print("ERROR: sample codes don't match: {} != {}".format(sample_current, sample_code))
        
    else:
        print("Warning: RE didn't match {}".format(infile))
    



# Massage data
points[:,0] *= -1 # Convert negative qs to positive
values *= np.power( points[:,0], 2.0 ) # Kratky-style



# Generate grid of data
from scipy.interpolate import griddata

qs = np.linspace(0, 0.25, 300)
times = np.linspace(0, max_time, 2000)
Qs, TIMEs = np.meshgrid(qs, times)
Z = griddata( points, values, (Qs, TIMEs), method='cubic' )


Z = np.nan_to_num(Z)





# Plot

class Data2D_custom(Data2D):
    
    def xy_axes(self):
        return self.x_axis, self.y_axis
    
    
    def _plot(self, save=None, show=False, ztrim=[0.01, 0.01], size=10.0, plot_buffers=[0.1,0.1,0.1,0.1], top_plot=0.25, **kwargs):
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        self.fig = plt.figure( figsize=(size,size*0.75), facecolor='white' )
        
        
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height*(1-top_plot)] )
        
        
        
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
        
        # TODO: Handle 'origin' correctly. (E.g. allow it to be set externally.)
        self.im = plt.imshow(self.Z, vmin=0, vmax=1, cmap=cmap, interpolation='nearest', extent=extent, origin='lower', aspect='auto')
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
        
        
        
        
        # Top plot
        self.axt = self.fig.add_axes( [left_buf, bottom_buf+fig_height*(1-top_plot), fig_width, fig_height*top_plot] )
                
        self.axt.plot(self.top_x, self.top_y, 'k', linewidth=4.0)
        xi, xf, yi, yf = self.ax.axis()
        xit, xft, yit, yft = self.axt.axis()
        #self.axt.axis( [xi, xf, yit, yft] )
        self.axt.axis( [xi, xf, 20, 190] )
        self.axt.xaxis.set_ticklabels([])
        self.axt.set_ylabel('$T \, (\mathrm{^{\circ}C})$')
        #self.axt.yaxis.set_ticks( range(120, 200, 20) )
        self.axt.yaxis.set_ticks( range(20, 200, 40) )
        self.axt.tick_params(axis='both', which='major', labelsize=20)

        
        
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
        
        pass




# Output
data = Data2D_custom()
data.x_label = 'time'
data.x_rlabel = '$\mathrm{time \, (h)}$'
data.y_label = 'q'
data.y_rlabel = '$q_r \, (\mathrm{\AA^{-1}})$'
data.plot_args = { 'rcParams': {'axes.labelsize': 45,
                            'xtick.labelsize': 30,
                            'ytick.labelsize': 30,
                            },
                }


data.x_axis = times
data.y_axis = qs
data.data = Z.transpose()

data.top_x = tT_times
data.top_y = tT_temperatures

outfile = ''.join([sample_code, '.png'])
print(outfile)
data.plot(show=False, save=outfile, ztrim=[0.1,0.002], yticks=np.arange(0, 0.24, 0.05), cmap=cmap_vge_hdr, plot_buffers=[0.20,0.05,0.15,0.05],  )

