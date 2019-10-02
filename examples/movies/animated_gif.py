#!/usr/bin/python

# This simple script generates an animated GIF based on a sequence of input image files.

import os
import glob
import re
import numpy as np






for analysis in ['thumbnails2', 'circular_average_q2I']:
#for analysis in ['circular_average_q2I']:
    
    source_dir = '../{}/'.format(analysis)
    output_dir = './{}/'.format(analysis)
    
    for run_id in ['run4_']:
        #for film_id in ['A1', 'A2', 'A3']:
        for film_id in ['CB-aPP_1kDa', 'CB-aPP_500Da', 'GAL-aPP_1kDa', 'Lact-aPP_1kDa']:
            
            
            if analysis=='thumbnails2':
                pattern = '{}*{}*_saxs.jpg'.format(run_id, film_id)
                file_re = re.compile('^(.+)_x-?\d+\.\d+_.+_(\d+)_saxs.jpg$')
                
            elif analysis=='circular_average_q2I':
                pattern = '{}*{}*_saxs_q2I.png'.format(run_id, film_id)
                file_re = re.compile('^(.+)_x-?\d+\.\d+_.+_(\d+)_saxs_q2I.png$')
                
            else:
                pattern = '{}*{}*_saxs.png'.format(run_id, film_id)
                file_re = re.compile('^(.+)_x-?\d+\.\d+_.+_(\d+)_saxs.png$')
                
            print('Doing pattern {}'.format(pattern))
                        

            #run04_kinetics-120C_150nm-chip1_th0.150_7301.4s_T155.007C_10.00s_75771_saxs.jpg
            #pattern = 'run04_kinetics-120C_150nm-chip1*_saxs.jpg'
            #pattern = 'run04_kinetics-120C_65nm-chip1*_saxs.jpg'
            #pattern = 'run04_kinetics-120C_15nm-chip1*_saxs.jpg'

            outname = pattern[:-len('*_saxs.jpg')]
            outname = outname.replace('*', '_')
            outname = '{}{}'.format(output_dir, outname)


            # Select the files to animate
            infiles = glob.glob(source_dir+pattern)
            print('{} infiles'.format(len(infiles)))


            # Sort into correct order

            sort_array = []
            filenames_array = []
            sort_idx = 1
            for infile in infiles:
                
                filename = infile[len(source_dir):]
                m = file_re.match(filename)
                
                if m:
                    sort_array.append( int(m.groups()[sort_idx]) )
                    filenames_array.append(filename)
                
                
            sort_array = np.asarray(sort_array)
            filenames_array = np.asarray(filenames_array)

            idx = np.argsort(sort_array)
            sort_array = sort_array[idx]
            filenames_array = filenames_array[idx]
            
            print('    {} in filenames_array'.format(len(filenames_array)))



            # Don't include every frame
            filenames_fewer = []
            for i, filename in enumerate(filenames_array):
                
                if i%2==0:
                    filenames_fewer.append(filename)


            # Prepare command
            # (Animation is generated using imagemagick 'convert' bash command.)

            #cmd = "convert -delay 20 -resize 50% -fill white  -undercolor '#00000080'  -gravity NorthWest -annotate +0+5 ' Text ' "
            #cmd = "convert -delay 15 -loop 1 -resize 50% "
            #cmd = "convert -crop 450x450+60+220  +repage -delay 15 -loop 1 -resize 50% "
            cmd = "convert +repage -delay 15 -loop 1 -resize 80% "
                    
            print('    {} in filenames_fewer'.format(len(filenames_fewer)))
            for filename in filenames_fewer:
                
                cmd += '{}/{} '.format(source_dir, filename)

            cmd += ' {}.gif'.format(outname)


            # Execute command
            print('  Saving {}'.format(outname))
            os.system(cmd)


