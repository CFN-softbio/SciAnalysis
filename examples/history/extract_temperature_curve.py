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


run_code = 'run04'
run_code = 'run01_linear-ramp'
run_code = 'run02_kinetics-100C'
run_code = 'run03_kinetics-150C'

run_code = 'run05_kinetics-60C'

run_code = 'SAXS-'


source_dir = '../../'

infiles = glob.glob( ''.join([source_dir, run_code, '*_saxs.tiff']) )

file_re = re.compile('^(.+)_th0.150_(\d+\.\d)s_T(\d+.\d\d\d)C_.+s_(\d+)_saxs.tiff$')
file_re = re.compile('^(.+)_th0\.\d+_(\d+\.\d)s_T(\d+.\d\d\d)C_.+s_(\d+)_saxs.tiff$')

print('Considering {} infiles'.format(len(infiles)))

filenames = []
times = []
temperatures = []
sIDs = []
for infile in infiles:
    
    m = file_re.match(infile)
    if m:
        sample_code, clock, temperature, sID = m.groups()

        filenames.append(infile)
        times.append(float(clock))
        temperatures.append(float(temperature))
        sIDs.append(int(sID))
    else:
        print('WARNING: No RE match for {}'.format(infile))
        
        
print('Extracted {} infiles ({:.1f}%)'.format(len(times), 100.0*len(times)/len(infiles) ))


# Sort
filenames = np.asarray(filenames)
times = np.asarray(times)
temperatures = np.asarray(temperatures)
sIDs = np.asarray(sIDs)

indices = np.argsort(sIDs)
filenames = filenames[indices]
times = times[indices]
temperatures = temperatures[indices]
sIDs = sIDs[indices]





# Fix 'clock resets'
fixed_times = []
running_delta = 0
for i, (filename, time, sID) in enumerate(zip(filenames, times, sIDs)):
    if i>0:
        if time<times[i-1]:
            # There is a break in the sorted order of times
            print('Break in the sorting of data:')
            print('    i={}: sID={}; clock={:.1f}'.format(i-1, sIDs[i-1], times[i-1]))
            print('    i={}: sID={}; clock={:.1f}'.format(i, sID, time))
            if i<len(filenames):
                print('    i={}: sID={}; clock={:.1f}'.format(i+1, sIDs[i+1], times[i+1]))
            
            datetime_current_file = os.path.getmtime(filename)
            datetime_previous_file = os.path.getmtime(filenames[i-1])
            delta = datetime_current_file - datetime_previous_file
            if delta<0:
                print('WARNING: The time delta is negative.')
            else:
                running_delta += times[i-1] + delta
                
    fixed_times.append(time+running_delta)
        
times = np.asarray(fixed_times)

times /= 60.0*60.0 # convert: seconds to hours

line = DataLine(x=times, y=temperatures)
line.x_label = 'time'
line.x_rlabel = r'$\mathrm{time} \, (\mathrm{h})$'
line.y_label = 'temperature'
line.y_rlabel = r'$T \, (^{\circ}C)$'

line.plot(show=False, save=''.join([run_code, '.png']))


with open(''.join([run_code, '.dat']), 'w') as fout:
    fout.write('#filename\tsID\ttime\ttemperature\n')
    for filename, sID, time, temperature in zip(filenames, sIDs, times, temperatures):
        fout.write('{}\t{}\t{}\t{}\n'.format(filename, sID, time, temperature))
        fout