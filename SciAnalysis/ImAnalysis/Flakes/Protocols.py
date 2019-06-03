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
        
        outfile = self.get_outfile(data.name, output_dir)
        data.plot_image(save=outfile, size=10*run_args['resize'], **run_args)
        
        return results
        
