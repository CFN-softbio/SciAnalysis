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

import matplotlib as mpl
from ..settings import *
if MATPLOTLIB_BACKEND is not None:
    mpl.use(MATPLOTLIB_BACKEND)
mpl.rcParams['mathtext.fontset'] = 'cm'
import pylab as plt

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
        
        
    def get_scale(self):
        return np.average([self.x_scale, self.y_scale])
        
        
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
        
    def crop_edges(self, left=0, right=0, bottom=0, top=0, relative=True):
        '''Crop the image, relative to the image center.'''
        
        height, width = self.data.shape
        
        if relative:
            # Convert from ratio (0 to 1) into pixels
            left *= width
            right *= width
            top *= height
            bottom *= bottom
        
        xi = max( int(left), 0 )
        xf = min( int(width-right), width )
        yi = max( int(top), 0 )
        yf = min( int(height-bottom), height )
        
        self.data = self.data[ yi:yf, xi:xf ]
        
                
        
    def equalize(self, num_bins=256, max_val=255):
        '''Change the histogram of intensities to be more 'even'.'''
        # Adapted from example at:
        # http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html

        try:
            data = self.data
            
            #get image histogram
            imhist, bins = np.histogram(data.flatten(), num_bins, density=True)
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
    


# Data2DImageRGB 
################################################################################    
class Data2DImageRGB(Data2DImage):
    
    def load_image(self, infile):
        
        img1 = PIL.Image.open(infile)
        img2 = img1.convert('I') # 'I' : 32-bit integer pixels
        
        self.data_rgb = np.array(img1)
        self.data = np.asarray(img2)
        del img1
        del img2
        
        
    def plot_image(self, save=None, show=False, size=10, ztrim=[0.01, 0.01], **plot_args):
        '''Generates a false-color image of the 2D data.'''
        
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
            
        h, w, c = self.data_rgb.shape
        
        data = self.data_rgb
        if 'image_contrast' in plot_args:
            in_range = ( plot_args['image_contrast'][0]*255, plot_args['image_contrast'][1]*255 )
        else:
            in_range = (zmin, zmax)
            
        import skimage
        data = skimage.exposure.rescale_intensity(data, in_range=in_range, out_range='dtype')
        
        self.fig = plt.figure( figsize=(size, size*h/w), facecolor='white' )
        self.ax = self.fig.add_axes([0,0,1,1])
        self.ax.imshow(data)
        
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
        


    # End class Data2DImageRGB(Data2DImage)
    ########################################
