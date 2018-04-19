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



class DataLines_current(DataLines):
    
    def _plot_extra(self, **plot_args):
        pass
        

            


lines = DataLines_current()


x = np.linspace(-1, 1, num=500)
y = np.arccos(x)
line = DataLine(x=x, y=y)
line.plot_args['color'] = 'k'
line.plot_args['linestyle'] = '-'
line.plot_args['linewidth'] = 4.0
line.plot_args['marker'] = None
lines.add_line(line)
            

y = np.sqrt(2*(-x+1))
line = DataLine(x=x, y=y)
line.plot_args['color'] = 'purple'
line.plot_args['linestyle'] = '-'
line.plot_args['linewidth'] = 4.0
line.plot_args['marker'] = None
lines.add_line(line)


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
    
outfile = os.path.join(output_dir, 'compare.png')
lines.plot(save=outfile, _plot_range=[0, 0.13, 0, 110], plot_buffers=[0.22, 0.04, 0.22, 0.04], _xticks=[0, 0.04, 0.08, 0.12], dpi=200)
    
