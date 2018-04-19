#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys # To get commandline arguments
import re # Regular expressions

ANALYSIS_PATH='/home/kyager/current/code/SciAnalysis/main'
ANALYSIS_PATH in sys.path or sys.path.append(ANALYSIS_PATH)


from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *


if len(sys.argv)>1:
    # Use commandline argument
    source_dir = ''
    output_dir = './'
    infile = sys.argv[1]
    
else:
    # Use hardcoded file
    root_dir = '/media/extend/CHX/'
    data_dir = 'LDRD_Meso_correlation/2015_10Oct-ref_patterns/21/'
    
    protocol = 'view'
    source_dir = root_dir+'data/'+data_dir
    output_dir = root_dir+'analysis/'+data_dir+'/'+protocol+'/'
    
    tools.make_dir(output_dir)
    
    infile = 'AgBH_22_master.h5'
    infile = 'Coral_porous_glass_0.3s_30_master.h5'
    infile = 'chip411_found_star_x-0.356_y13.468_922_master.h5'




f = tools.Filename(source_dir+infile)

if f.get_ext()=='.h5':

    # Check if this is a HDF5 data (as opposed to master) file
    m = re.compile('(.+)_data_\d{5,20}\.h5').match(f.get_filename())
    if m:
        # Switch to corresponding master file
        f = tools.Filename(f.get_path() + m.groups()[0] + '_master.h5')
    
    data = Data2DScattering(f.get_filepath())
    data.threshold_pixels(4294967295-1) # Eiger inter-module gaps
    data.set_z_display( [None, None, 'gamma', 0.3] )
    data.plot(show=True)
    

elif f.get_ext()=='.txt' or f.get_ext()=='.dat' or f.get_ext()=='.data':
    
    data = np.loadtxt(infile)
    x = data[:,0]
    y = data[:,-1]
    line = DataLine(x=x, y=y)
    line.plot(show=True)







