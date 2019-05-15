#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

from ..Protocols import *


class thumbnails_contrast(thumbnails):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : 0.5,
                        'blur' : 1.0,
                        'resize' : 0.5,
                        'cmap' : mpl.cm.bone,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'])
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize']) # Shrink
        
        data.set_z_display([None, None, 'gamma', 1.0])
        outfile = self.get_outfile(data.name, output_dir)
        data.plot_image(outfile, ztrim=[0.01, 0.005], **run_args)
        
        return results
        
