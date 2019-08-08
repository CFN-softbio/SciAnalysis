#!/usr/bin/python3

import os
import glob
from shutil import copyfile


source_dir = './'

def make_dir(directory):
    if not os.path.isdir(directory):
        os.makedirs( directory )
 


ignore = ['collate']

for folder in glob.glob(source_dir+'*'):
    
    if folder not in ignore and '.py' not in folder and '.tif' not in folder:
        
        for protocol in [ ['histogram', 'aggregate.png'], ['tile_img', 'tile.jpg'] ]:
            
            protocol, filename = protocol
            
            source = os.path.join('./', folder, 'analysis', protocol, filename)
            output_dir = os.path.join('./', 'collate', protocol)
            destination = os.path.join(output_dir, '{}.png'.format(folder))
            
            if os.path.exists(source):
                print('Copying {}'.format(folder))
                print('   source: {}'.format(source))
                print('   destination: {}'.format(destination))
                
                make_dir(output_dir)
                copyfile(source, destination)
                
            else:
                print('Skipping {}'.format(folder))
