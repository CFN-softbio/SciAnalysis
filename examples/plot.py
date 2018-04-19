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



root_dir = './'
source_dir = os.path.join(root_dir, './')
output_dir = os.path.join(root_dir, './')


integration_time = 1*60 # seconds
rescale = 1e6/integration_time


class DataLines_current(DataLines):
    
    def _plot_extra(self, **plot_args):
        
        if True:
            # loc
            # http://matplotlib.org/users/legend_guide.html#legend-location
            # 2  9  1
            # 6 10  5
            # 3  8  4
            #leg = plt.legend( self.get_labels(), loc=2)
            
            if 'legend_location' in plot_args:
                loc_str = plot_args['legend_location']
                if loc_str=='NW':
                    loc = 2
                elif loc_str=='N':
                    loc = 9
                elif loc_str=='NE':
                    loc = 1
                elif loc_str=='E':
                    loc = 5
                elif loc_str=='SE':
                    loc = 4
                elif loc_str=='S':
                    loc = 8
                elif loc_str=='SW':
                    loc = 3
                elif loc_str=='W':
                    loc = 6
                elif loc_str=='C':
                    loc = 10
            else:
                loc = 1
                
                
            handles, labels = self.ax.get_legend_handles_labels()
            leg = self.ax.legend(handles, labels, loc=loc)
            #self.ax.legend()
            leg.draw_frame(False)
            
        
        
        line = self.lines[0]
        x, y = line.target_y(np.max(line.y))
        self.ax.text(x, y*1.02, '$q_0$', size=40, horizontalalignment='center', verticalalignment='bottom')
            
            
            
            
y_rlabel = '$q^2 I(q) \, (10^6 \AA^{-2} \mathrm{cts/pix}/s)$'
y_rlabel = '$q^2 I(q) \, (\mathrm{a.u.})$'
lines = DataLines_current(x_label='q', y_label='q2I', x_rlabel='$q \, (\AA^{-1})$', y_rlabel=y_rlabel )



x = np.linspace(0, 100, 100)
y = np.square(x) + np.random.normal(0,100, len(x))
line = DataLine(x=x, y=y)
line.y_rlabel = r'$\langle I \rangle$'
line.plot(save='output.png')


exit()



if True:
    
    infiles = glob.glob(source_dir+'*8p25*.dat')
    infile = infiles[0]
                          
    line = DataLine(infile, format='custom', skiplines=1, comment_char='#', xindex=0, yindex=2, name='8.25%')
    line.y *= rescale
    
    line.plot_args['color'] = 'purple'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None
    
    lines.add_line(line)
    
    
if True:
    
    infiles = glob.glob(source_dir+'*4p13*.dat')
    infile = infiles[0]
                          
    line = DataLine(infile, format='custom', skiplines=1, comment_char='#', xindex=0, yindex=2, name='4.13%')
    line.y *= rescale
    
    line.plot_args['color'] = 'b'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None
    
    lines.add_line(line)
    
    
if True:
    
    infiles = glob.glob(source_dir+'*1p65*.dat')
    infile = infiles[0]
                          
    line = DataLine(infile, format='custom', skiplines=1, comment_char='#', xindex=0, yindex=2, name='1.65%')
    line.y *= rescale
    
    line.plot_args['color'] = 'k'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None
    
    lines.add_line(line)
    

lines.plot_args = { 'color' : 'k',
                'marker' : 'o',
                'linewidth' : 3.0,
                'legend_location': 'NE'
                }  

lines.plot_args['rcParams'] = { 
                'axes.labelsize': 45,
                'xtick.labelsize': 35,
                'ytick.labelsize': 35,    
                #'legend.borderpad' : 0 ,
                'legend.fontsize' : 30 ,
                #'legend.numpoints' : 1 ,
                #'legend.handlelength' : 1.0 ,
                'legend.labelspacing' : 0.25 ,
                'legend.handletextpad' : 0.5 ,
                #'legend.columnspacing' : 0.0 ,
                }
    
outfile = os.path.join(output_dir, 'Iq2_concentration.png')
lines.plot(save=outfile, plot_range=[0, 0.13, 0, 110], plot_buffers=[0.22, 0.04, 0.22, 0.04], xticks=[0, 0.04, 0.08, 0.12], dpi=200)
    
