#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.ImAnalysis.Data` - Base objects for ImAnalysis
================================================
.. module:: SciAnalysis.ImAnalysis
   :synopsis: Provides base classes for doing analysis of images
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  This code defines some baseline objects for image analysis.
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################



#import sys
import re # Regular expressions

import numpy as np
import pylab as plt
import matplotlib as mpl
#from scipy.optimize import leastsq
#import scipy.special

import PIL # Python Image Library (for opening PNG, etc.)    

from .. import tools
from ..Data import *




   
    
# Data2DImage        
################################################################################    
class Data2DImage(Data2D):
    
    def __init__(self, infile=None, name=None, **kwargs):
        '''Creates a new Data2D object, which stores a scattering area detector
        image.'''
        
        super(Data2DImage, self).__init__(infile=infile, **kwargs)
        
        if name is not None:
            self.name = name
        elif infile is not None:
            self.name = tools.Filename(infile).get_filebase()            
        
        
        
    # Coordinate methods
    ########################################
                    
    def get_origin(self):
        
        x0 = 0
        y0 = 0
        
        return x0, y0


    def get_center(self):
        
        height, width = self.data.shape
        
        x0 = int( width/2 )
        y0 = int( height/2 )
        
        return x0, y0    
        
        
    # Data modification
    ########################################
        
    def invert(self, max_val=255.0):
        
        self.data = max_val-self.data
        
        
    def threshold(self, threshold, invert=False):

        if invert:
            self.data = np.where( self.data>threshold, 0, 1 )
        else:
            self.data = np.where( self.data>threshold, 1, 0 )
        
        
    def threshold_pixels(self, threshold, new_value=0.0):
        
        self.data[self.data>threshold] = new_value
        
        
    def maximize_intensity_spread(self, max_val=255.0):
        
        self.data = self.data - np.min(self.data)
        self.data = self.data*(max_val/np.max(self.data))
        
        
    def crop(self, size):
        '''Crop the image, relative to the image center.'''
        
        height, width = self.data.shape
        #self.data = self.data[ 0:height, 0:width ] # All the data
        
        x0, y0 = self.get_center()
        xi = max( int(x0 - size*width/2), 0 )
        xf = min( int(x0 + size*width/2), width )
        yi = max( int(y0 - size*height/2), 0 )
        yf = min( int(y0 + size*height/2), height )
        
        self.data = self.data[ yi:yf, xi:xf ]
        
        
    def equalize(self, num_bins=256, max_val=255):
        '''Change the histogram of intensities to be more 'even'.'''
        # Adapted from example at:
        # http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html

        try:
            data = self.data
            
            #get image histogram
            imhist, bins = np.histogram(data.flatten(), num_bins, normed=True)
            cdf = imhist.cumsum() #cumulative distribution function
            cdf = max_val * cdf / cdf[-1] #normalize

            #use linear interpolation of cdf to find new pixel values
            data_equalized = np.interp(data.flatten(), bins[:-1], cdf)
            data_equalized = data_equalized.reshape(data.shape)
            
            self.data = data_equalized
            
            
        except TypeError:
            # Avoid unexplained error:
            # File "/usr/lib/python2.7/dist-packages/numpy/lib/function_base.py", line 169, in histogram
            # mn, mx = [mi+0.0 for mi in range]
            # TypeError: unsupported operand type(s) for +: 'instance' and 'float'

            pass
        
        
    def enhance(self, contrast=1.5, contrast_passes=0, resharpen_passes=0):
        
        import PIL.ImageEnhance
        import PIL.ImageFilter
        
        img = PIL.Image.fromarray(np.uint8(self.data))
        
        # Contrast
        enhancer = PIL.ImageEnhance.Contrast(img)
        img = enhancer.enhance(contrast)
        
        for i in range(contrast_passes):
            img = img.filter(PIL.ImageFilter.BLUR)
            enhancer = PIL.ImageEnhance.Contrast(img)
            img = enhancer.enhance(contrast)
        
        for i in range(resharpen_passes):
            img = img.filter(PIL.ImageFilter.SMOOTH) # Smooth
            img = img.filter(PIL.ImageFilter.SHARPEN) # Sharpen
        
        self.data = np.asarray(img)
        
 
    
       
        

    # End class Data2DImage(Data2D)
    ########################################
    