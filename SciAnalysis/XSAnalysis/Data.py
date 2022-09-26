#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.Data` - Base objects for XSAnalysis
================================================
.. module:: SciAnalysis.XSAnalysis
   :synopsis: Provides base classes for doing analysis of x-ray scattering data
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  This code defines some baseline objects for x-ray analysis.
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
from SciAnalysis.settings import *
if MATPLOTLIB_BACKEND is not None:
    mpl.use(MATPLOTLIB_BACKEND)
mpl.rcParams['mathtext.fontset'] = 'cm'
import pylab as plt

#from scipy.optimize import leastsq
#import scipy.special

import PIL # Python Image Library (for opening PNG, etc.)    

from SciAnalysis import tools
from SciAnalysis.Data import *
try:
    from .Eiger import *
except ImportError:
    # Eiger support is optional
    pass



   
    
# Data2DScattering        
################################################################################    
class Data2DScattering(Data2D):
    '''Represents the data from a 2D (area) detector in a scattering measurement.'''
    
    def __init__(self, infile=None, format='auto', calibration=None, mask=None, name=None, **kwargs):
        '''Creates a new Data2D object, which stores a scattering area detector
        image.'''
        
        super(Data2DScattering, self).__init__(infile=None, **kwargs)
        
        self.set_z_display([None, None, 'gamma', 0.3])
        
        self.calibration = calibration
        self.mask = mask
        
        self.detector_data = None # Detector-specific object
        self.data = None # 2D data
        self.measure_time = 0.0
        
        if name is not None:
            self.name = name
        elif infile is not None:
            if 'full_name' in kwargs and kwargs['full_name']:
                self.name = tools.Filename(infile).get_filename()
            else:
                self.name = tools.Filename(infile).get_filebase()            
                
        if infile is not None:
            self.load(infile, format=format)
        

    # Data loading
    ########################################

    def load(self, infile, format='auto', **kwargs):
        '''Loads data from the specified file.'''
        
        if format=='eiger' or infile[-10:]=='_master.h5':
            self.load_eiger(infile, **kwargs)
            
        elif format=='hdf5' or infile[-3:]=='.h5' or infile[-4:]=='.hd5':
            self.load_hdf5(infile)
            
        elif format=='tiff' or infile[-5:]=='.tiff' or infile[-4:]=='.tif':
            self.load_tiff(infile)
            
        elif format=='BrukerASCII' or infile[-6:]=='.ascii' or infile[-4:]=='.dat':
            self.load_BrukerASCII(infile)
            
        else:
            super(Data2DScattering, self).load(infile=infile, format=format, **kwargs)


        # Masking is now applied in the Processor's load method
        #if self.mask is not None:
            #self.data *= self.mask.data
            
            
    def load_eiger(self, infile, frame='all'):
        
        self.detector_data = EigerImages(infile)
        
        self.measure_time = self.detector_data.exposuretime
        
        if frame=='all':
            # Sum all frames together
            self.data = np.zeros(self.detector_data.dims)
            num_frames = self.detector_data.__len__()
            for i in range(num_frames):
                if num_frames>1 and i%50==0:
                    print('    Adding frame %d of %d (%.1f%%).'%(i, num_frames, i*100./num_frames))
                self.data += self.detector_data.get_frame(i)
        else:
            self.data = self.detector_data.get_frame(frame)
            
        # TODO: Support slicing to select specific frame ranges.

        
    def load_hdf5(self, infile):
        
        # WARNING: This code doesn't work. Loading data form HDF5 depends on
        # the hiearchy. This code snippet can be copied and revised based on 
        # the actual keys in a given file.
        f = h5py.File(infile, 'r')
        self.data = np.asarray( f['entry']['data']['data_000001'] )
        # TODO: Make this code work in general case.
        
        
    def load_tiff(self, infile):
        
        img = PIL.Image.open(infile).convert('I') # 'I' : 32-bit integer pixels
        self.data = ( np.copy( np.asarray(img) ) ).astype(np.float)
        del img
        
        
    def load_BrukerASCII(self, infile):
        '''Brute-force loading of Bruker ASCII data.'''
        
        ascii_pair_re = re.compile( '(^.+):(.+)$')
        ascii_data_re = re.compile( '^[\d \r]{1,90}$')

        print( '    Opening BrukerASCII: %s' % (infile) )
        with open(infile, 'r') as fin:
            
            height = None
            width = None
            data = None
            
            for i, line in enumerate(fin.readlines()):
                if height is None or width is None:
                    # Still searching for the data size
                    m = ascii_pair_re.match(line)
                    if m:
                        category = m.groups()[0].strip()
                        if category=='NROWS':
                            height = int(m.groups()[1])
                        if category=='NCOLS':
                            width = int(m.groups()[1])
                            
                elif data is None:
                    data = np.zeros( (height*width) )
                    idx = 0
                    
                else:
                    m = ascii_data_re.match(line)
                    if m:
                        els = line.split()
                        nums = np.asarray(els).astype(np.int)
                        data[idx:idx+len(nums)] = nums
                        idx += len(nums)
                    

        
        #self.data = np.loadtxt(infile, skiprows=i).reshape(width, height) # Fails if lines of different length
        self.data = data.reshape(height, width)

                
                
    # Coordinate methods
    ########################################
                    
    def get_origin(self):
        
        x0 = self.calibration.x0
        y0 = self.calibration.y0
        
        return x0, y0

        
    def _xy_axes(self):
        # TODO: test, and integrate if it makes sense
        
        dim_y,dim_x = self.data.shape
        
        x_axis = (np.arange(dim_x) - self.calibration.x0)*self.x_scale
        y_axis = (np.arange(dim_y) - self.calibration.y0)*self.y_scale
        
        return x_axis, y_axis
        
        
    # Data modification
    ########################################
        
    def threshold_pixels(self, threshold, new_value=0.0):
        
        self.data[self.data>threshold] = new_value
        
        
    def crop(self, size, shift_crop_up=0.0, make_square=False):
        '''Crop the data, centered about the q-origin. I.e. this throws away 
        some of the high-q information. The size specifies the size of the new
        image (as a fraction of the original full image width).
        
        shift_crop_up forces the crop to be off-center in the vertical (e.g. a
        value of 1.0 will shift it up so the q-origin is at the bottom of the
        image, which is nice for GISAXS).'''
        
        
        height, width = self.data.shape
        #self.data = self.data[ 0:height, 0:width ] # All the data
        
        x0, y0 = self.get_origin()
        
        if make_square:
            yi = max( int(y0 - size*height*(0.5+0.5*shift_crop_up) ), 0 )
            yf = min( int(y0 + size*height*(0.5-0.5*shift_crop_up) ), height )

            yspan = abs(yf-yi)
            xi = max( int(x0 - yspan/2), 0 )
            xf = min( int(x0 + yspan/2), width )
            
        else:
            xi = max( int(x0 - size*width/2), 0 )
            xf = min( int(x0 + size*width/2), width )
            yi = max( int(y0 - size*height*(0.5+0.5*shift_crop_up) ), 0 )
            yf = min( int(y0 + size*height*(0.5-0.5*shift_crop_up) ), height )
        
        
        self.data = self.data[ yi:yf, xi:xf ]
        
        
    def dezinger(self, sigma=3, tol=100, mode='median', mask=True, fill=False):
        # NOTE: This could probably be improved.
        
        if mode=='median':
            avg = ndimage.filters.median_filter(self.data, size=(sigma,sigma))
            variation = ndimage.filters.maximum_filter(avg, size=(sigma,sigma)) - ndimage.filters.minimum_filter(avg, size=(sigma,sigma))
            variation = np.where(variation > 1, variation, 1)
            idx = np.where( (self.data-avg)/variation > tol )
            
        elif mode=='gauss':
            # sigma=3, tol=1e5
            avg = ndimage.filters.gaussian_filter( self.data, sigma )
            local = avg - self.data/np.square(sigma)
            
            #dy, dx = np.gradient(self.data)
            #var = np.sqrt( np.square(dx) + np.square(dy) )
            
            idx = np.where( (self.data > avg) & (self.data > local) & (self.data-avg > tol)  )
        
        
        #self.data[idx] = 0
        if fill:
            self.data[idx] = avg[idx]
            
        if mask:
            self.mask.data[idx] = 0

                       
        
        
    # Data reduction
    ########################################
    
    def circular_average(self, **kwargs):
        '''Returns a 1D curve that is a circular average of the 2D data. The
        data is average over 'chi', so that the resulting curve is as a function
        of q.'''
        
        #return self.circular_average_q(**kwargs)
        return self.circular_average_q_bin(**kwargs)
    
    
    def circular_average_pixel(self, **kwargs):
        '''Returns a 1D curve that is a circular average of the 2D data. The
        data is average over the angular direction around the origin. Data is
        returned in terms of pixel distance from origin.'''
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
            
        dim_y, dim_x = self.data.shape
        x0 = self.calibration.x0
        y0 = self.calibration.y0
        
        # .ravel() is used to convert the 2D grids into 1D arrays.
        # This is not strictly necessary, but improves speed somewhat.
        
        data = self.data.ravel()
        pixel_list = np.where(mask.ravel()==1) # Non-masked pixels
        
        # Generate map of distances-from-origin
        x = np.arange(dim_x) - x0
        y = np.arange(dim_y) - y0
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2).ravel()
        #R = self.calibration.r_map().ravel()
        Rd = (R + 0.5).astype(int) # Simplify the R pixel-distances to closest integers
        
        num_per_R = np.bincount(Rd[pixel_list])
        idx = np.where(num_per_R!=0) # R-distances that actually have data
        
        r_vals = np.bincount( Rd[pixel_list], weights=R[pixel_list] )[idx]/num_per_R[idx]
        #q_vals = r_vals*self.calibration.get_q_per_pixel()
        I_vals = np.bincount( Rd[pixel_list], weights=data[pixel_list] )[idx]/num_per_R[idx]
        
        line = DataLine( x=r_vals, y=I_vals, x_label='r', y_label='I', x_rlabel='$r \, (\mathrm{pixels})$', y_rlabel=r'$\langle I \rangle \, (\mathrm{counts/pixel})$' )
        
        return line
        
            
            
    def circular_average_q(self, error=True, **kwargs):
        '''Returns a 1D curve that is a circular average of the 2D data. The
        data is average over 'chi', so that the resulting curve is as a function
        of q.'''
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
            
        # .ravel() is used to convert the 2D grids into 1D arrays.
        # This is not strictly necessary, but improves speed somewhat.
        
        data = self.data.ravel()
        pixel_list = np.where(mask.ravel()==1) # Non-masked pixels
        
        Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()
        Qd = (Q/dq + 0.5).astype(int)
        
        num_per_bin = np.bincount(Qd[pixel_list])
        idx = np.where(num_per_bin!=0) # q-distances that actually have data
        
        if error:
            x_vals = np.bincount( Qd[pixel_list], weights=Q[pixel_list] )[idx]/num_per_bin[idx]
            x_err = np.ones(len(x_vals))*dq/2
            
            I_vals = np.bincount( Qd[pixel_list], weights=data[pixel_list] )[idx]
            I_err = np.sqrt(I_vals) # shot-noise is sqrt of TOTAL counts
            I_vals /= num_per_bin[idx]
            I_err /= num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_err=x_err, y_err=I_err, x_label='r', y_label='I', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
            
        else:
            x_vals = np.bincount( Qd[pixel_list], weights=Q[pixel_list] )[idx]/num_per_bin[idx]
            I_vals = np.bincount( Qd[pixel_list], weights=data[pixel_list] )[idx]/num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_label='q', y_label='I(q)', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
        
        
        return line
    
    
    def circular_average_q_bin(self, bins_relative=1.0, error=False, **kwargs):
        '''Returns a 1D curve that is a circular average of the 2D data. The
        data is average over 'chi', so that the resulting curve is as a function
        of q.
        
        'bins_relative' controls the binning (q-spacing in data).
            1.0 means the q-spacing is (approximately) a single pixel
            2.0 means there are twice as many bins (spacing is half a pixel)
            0.1 means there are one-tenth the number of bins (i.e. each data point is 10 pixels)
            
        'error' sets whether or not error-bars are calculated.
        '''
        
        # This version uses numpy.histogram instead of converting the binning
        # into an equivalent integer list (and then using numpy.bincount).
        # This code is slightly slower (30-40% longer to run). However, this
        # version is slightly more general inasmuch as one can be more
        # arbitrary about the number of bins to be used (doesn't have to match
        # the q-spacing).
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
            
        # .ravel() is used to convert the 2D grids into 1D arrays.
        # This is not strictly necessary, but improves speed somewhat.
        
        data = self.data.ravel()
        pixel_list = np.where(mask.ravel()==1) # Non-masked pixels
        
        Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()

        x_range = [np.min(Q[pixel_list]), np.max(Q[pixel_list])]
        bins = int( bins_relative * abs(x_range[1]-x_range[0])/dq )
        num_per_bin, rbins = np.histogram(Q[pixel_list], bins=bins, range=x_range)
        idx = np.where(num_per_bin!=0) # Bins that actually have data


        if error:
            # TODO: Include error calculations
            
            x_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=Q[pixel_list] )
            
            # Create array of the average values (mu), in the original array layout
            locations = np.digitize( Q[pixel_list], bins=rbins, right=True) # Mark the bin IDs in the original array layout
            mu = (x_vals/num_per_bin)[locations-1]
            
            weights = np.square(Q[pixel_list] - mu)
            
            x_err, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=weights )
            x_err = np.sqrt( x_err[idx]/num_per_bin[idx] )
            x_err[0] = dq/2 # np.digitize includes all the values less than the minimum bin into the first element
            
            x_vals = x_vals[idx]/num_per_bin[idx]
            
            I_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=data[pixel_list] )
            I_err_shot = np.sqrt(I_vals)[idx]/num_per_bin[idx]
            
            mu = (I_vals/num_per_bin)[locations-1]
            weights = np.square(data[pixel_list] - mu)
            I_err_std, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=weights )
            I_err_std = np.sqrt( I_err_std[idx]/num_per_bin[idx] )
            
            
            y_err = np.sqrt( np.square(I_err_shot) + np.square(I_err_std) )
            if True:
                # Student t-test formulation of error
                import scipy.stats
                confidence_2 = 0.95 # Two-tailed
                confidence = 0.5*(1+confidence_2)
                DF = num_per_bin[idx] - 1
                z = stats.t.ppf(confidence, DF)
                y_err = z*I_err_std/np.sqrt(num_per_bin[idx])
                
            
            I_vals = I_vals[idx]/num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_err=x_err, y_err=y_err, x_label='q', y_label='I(q)', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
            
            
        else:
            x_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=Q[pixel_list] )
            x_vals = x_vals[idx]/num_per_bin[idx]
            I_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=data[pixel_list] )
            I_vals = I_vals[idx]/num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_label='q', y_label='I(q)', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
        
        
        
        return line
    
    
    def circular_average_q_rich(self, bins_relative=1.0, error=False, **kwargs):
        # TODO:
        # Give q, q_err, I, I_err_shot, I_err_std, I_err_total, I_min, I_max
        # other measures of anisotropy? (e.g. eta, S)
        # measures of count statistics (coherence estimate)
        # 'average' from fitting count stats
        # 'average' from doing a 'smart smoothing' (width based on local gradient)
        pass
    
    
    def circular_average_q_range(self, q, dq, error=True, **kwargs):
        '''Returns a 1D curve that is a circular average of the 2D data. The
        data is average over 'chi', so that the resulting curve is as a function
        of q.'''
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
            
        # .ravel() is used to convert the 2D grids into 1D arrays.
        # This is not strictly necessary, but improves speed somewhat.
        
        data = self.data.ravel()
        pixel_list = np.where( (abs(self.calibration.q_map().ravel()-q)<dq) & (mask.ravel()==1) )
        
        Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()
        Qd = (Q/dq + 0.5).astype(int)
        
        num_per_bin = np.bincount(Qd[pixel_list])
        idx = np.where(num_per_bin!=0) # q-distances that actually have data
        
        if error:
            x_vals = np.bincount( Qd[pixel_list], weights=Q[pixel_list] )[idx]/num_per_bin[idx]
            x_err = np.ones(len(x_vals))*dq/2
            
            I_vals = np.bincount( Qd[pixel_list], weights=data[pixel_list] )[idx]
            I_err = np.sqrt(I_vals) # shot-noise is sqrt of TOTAL counts
            I_vals /= num_per_bin[idx]
            I_err /= num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_err=x_err, y_err=I_err, x_label='r', y_label='I', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
            
        else:
            x_vals = np.bincount( Qd[pixel_list], weights=Q[pixel_list] )[idx]/num_per_bin[idx]
            I_vals = np.bincount( Qd[pixel_list], weights=data[pixel_list] )[idx]/num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_label='q', y_label='I(q)', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
        
        
        return line    
    
    
    def sector_average_q_bin(self, angle=0, dangle=30, bins_relative=1.0, error=False, **kwargs):
        '''Returns a 1D curve that is a sector average of the 2D data. The
        data is average over 'chi' across the range specified by angle and 
        dangle, so that the resulting curve is as a function of q.
        
        'bins_relative' controls the binning (q-spacing in data).
            1.0 means the q-spacing is (approximately) a single pixel
            2.0 means there are twice as many bins (spacing is half a pixel)
            0.1 means there are one-tenth the number of bins (i.e. each data point is 10 pixels)
            
        'error' sets whether or not error-bars are calculated.
        '''
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
            
        # .ravel() is used to convert the 2D grids into 1D arrays.
        # This is not strictly necessary, but improves speed somewhat.
        
        data = self.data.ravel()
        
        
        Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()
        A = self.calibration.angle_map().ravel()
        
        
        pixel_list = np.where( (mask.ravel()==1) & (abs(A-angle)<dangle/2) )
        
        if 'show_region' in kwargs and kwargs['show_region']:
            region = np.ma.masked_where(abs(self.calibration.angle_map()-angle)>dangle/2, self.calibration.q_map())
            self.regions = [region]
        

        x_range = [np.min(Q[pixel_list]), np.max(Q[pixel_list])]
        bins = int( bins_relative * abs(x_range[1]-x_range[0])/dq )
        num_per_bin, rbins = np.histogram(Q[pixel_list], bins=bins, range=x_range)
        idx = np.where(num_per_bin!=0) # Bins that actually have data


        if error:
            # TODO: Include error calculations
            
            x_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=Q[pixel_list] )
            
            # Create array of the average values (mu), in the original array layout
            locations = np.digitize( Q[pixel_list], bins=rbins, right=True) # Mark the bin IDs in the original array layout
            mu = (x_vals/num_per_bin)[locations-1]
            
            weights = np.square(Q[pixel_list] - mu)
            
            x_err, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=weights )
            x_err = np.sqrt( x_err[idx]/num_per_bin[idx] )
            x_err[0] = dq/2 # np.digitize includes all the values less than the minimum bin into the first element
            
            x_vals = x_vals[idx]/num_per_bin[idx]
            
            I_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=data[pixel_list] )
            I_err_shot = np.sqrt(I_vals)[idx]/num_per_bin[idx]
            
            mu = (I_vals/num_per_bin)[locations-1]
            weights = np.square(data[pixel_list] - mu)
            I_err_std, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=weights )
            I_err_std = np.sqrt( I_err_std[idx]/num_per_bin[idx] )
                
            y_err = np.sqrt( np.square(I_err_shot) + np.square(I_err_std) )
            I_vals = I_vals[idx]/num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_err=x_err, y_err=y_err, x_label='q', y_label='I(q)', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
            
            
        else:
            x_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=Q[pixel_list] )
            x_vals = x_vals[idx]/num_per_bin[idx]
            I_vals, rbins = np.histogram( Q[pixel_list], bins=bins, range=x_range, weights=data[pixel_list] )
            I_vals = I_vals[idx]/num_per_bin[idx]
            
            line = DataLine( x=x_vals, y=I_vals, x_label='q', y_label='I(q)', x_rlabel='$q \, (\mathrm{\AA^{-1}})$', y_rlabel=r'$I(q) \, (\mathrm{counts/pixel})$' )
        
        
        
        return line    
    
    def overlay_ring(self, q0, dq, clear=False):
        '''Add an overlay region that is a ring of constant q.'''
        region = self.calibration.q_map()
        region = np.ma.masked_where(abs(region-q0)>dq, region)
        
        if clear or self.regions is None:
            self.regions = [region]
        else:
            self.regions.append(region)
        
    
    # Data extraction
    ########################################        
    
    def _linecut_test(self, q0, dq):
        '''Internal function for testing/prototyping.'''
        
        #region = np.ones(self.data.shape)
        region = self.calibration.q_map()
        region = np.ma.masked_where(abs(region-q0)>dq, region)
        
        self.region = region
        
        pass
        # TODO
        # numpy.histogram
        # scipy.ndimage.measurements.histogram
        
        
    def linecut_angle(self, q0, dq, x_label='angle', x_rlabel='$\chi \, (^{\circ})$', y_label='I', y_rlabel=r'$I (\chi) \, (\mathrm{counts/pixel})$', mask_fraction_cutoff=0, **kwargs):
        '''Returns the intensity integrated along a ring of constant q.'''
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
        
        
        data = self.data.ravel()
        pixel_list = np.where( (abs(self.calibration.q_map().ravel()-q0)<dq) & (mask.ravel()==1) )
        pixel_list_maskless = np.where( (abs(self.calibration.q_map().ravel()-q0)<dq)  )
        

        if 'show_region' in kwargs and kwargs['show_region']:
            region = np.ma.masked_where(abs(self.calibration.q_map()-q0)>dq, self.calibration.angle_map())
            self.regions = [region]

        #Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()
        
        # Generate map
        M = self.calibration.angle_map().ravel()
        scale = np.degrees( np.abs(np.arctan(1.0/(q0/dq))) ) # approximately 1-pixel
        
        Md = (M/scale + 0.5).astype(int) # Simplify the distances to closest integers
        Md -= np.min(Md)
        
        num_per_m = np.bincount(Md[pixel_list])


        idx = np.where(num_per_m!=0) # Old method: consider bins that have >0 pixels
        
        # TODO: Fix the new method
        # New method: Only include bins that don't have too many masked pixels
        # e.g. mask_fraction_cutoff=0.8 excludes bins where >20% of pixels were masked
        #num_per_m_maskless = np.bincount(Md[pixel_list_maskless])
        #mask_fractions = num_per_m/num_per_m_maskless
        #idx = np.where(mask_fractions>mask_fraction_cutoff)
        
        
        x_vals = np.bincount( Md[pixel_list], weights=M[pixel_list])[idx]/num_per_m[idx]
        I_vals = np.bincount( Md[pixel_list], weights=data[pixel_list])[idx]/num_per_m[idx]
        
        line = DataLineAngle( x=x_vals, y=I_vals, x_label=x_label, y_label=y_label, x_rlabel=x_rlabel, y_rlabel=y_rlabel )

        #line.mask_fractions = mask_fractions[idx]
        
        #x_vals_full = np.bincount(Md, weights=M)/num_per_m
        #line.f_chi = len(x_vals)/(len(x_vals_full)+1) # Fraction of full circle that we have actually sampled
        line.f_chi = len(x_vals)*scale/360
        line.dchi = scale

        
        return line         
    
    
    def linecut_qr(self, qz, dq, x_label='qr', x_rlabel='$q_r \, (\mathrm{\AA}^{-1})$', y_label='I', y_rlabel=r'$I (q_r) \, (\mathrm{counts/pixel})$', **kwargs):
        '''Returns the intensity integrated along a line of constant qz.'''

        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data

        data = self.data.ravel()
        pixel_list = np.where( (abs(self.calibration.qz_map().ravel()-qz)<dq) & (mask.ravel()==1) )


        if 'show_region' in kwargs and kwargs['show_region']:
            region = np.ma.masked_where(abs(self.calibration.qz_map()-qz)>dq, self.calibration.qr_map())
            self.regions = [region]

        #Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()

        # Generate map
        M = self.calibration.qr_map().ravel()
        scale = dq # approximately 1-pixel

        Md = (M/scale + 0.5).astype(int) # Simplify the distances to closest integers
        Md -= np.min(Md)

        num_per_m = np.bincount(Md[pixel_list])
        idx = np.where(num_per_m!=0) # distances that actually have data

        x_vals = np.bincount( Md[pixel_list], weights=M[pixel_list] )[idx]/num_per_m[idx]
        I_vals = np.bincount( Md[pixel_list], weights=data[pixel_list] )[idx]/num_per_m[idx]

        line = DataLineAngle( x=x_vals, y=I_vals, x_label=x_label, y_label=y_label, x_rlabel=x_rlabel, y_rlabel=y_rlabel )

        return line


    def linecut_qz(self, qr, dq, x_label='qz', x_rlabel='$q_z \, (\mathrm{\AA}^{-1})$', y_label='I', y_rlabel=r'$I (q_z) \, (\mathrm{counts/pixel})$', q_mode='qr', **kwargs):
        '''Returns the intensity integrated along a line of constant qr.'''

        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data

        data = self.data.ravel()
        
        if q_mode=='qr':
            pixel_list = np.where( (abs(self.calibration.qr_map().ravel()-qr)<dq) & (mask.ravel()==1) )
        elif q_mode=='qx':
            pixel_list = np.where( (abs(self.calibration.qx_map().ravel()-qr)<dq) & (mask.ravel()==1) )
        else:
            print('ERROR: q_mode {} not recognized in linecut_qz.'.format(q_mode))


        if 'show_region' in kwargs and kwargs['show_region']:
            
            if q_mode=='qr':
                map_use = self.calibration.qr_map()
            elif q_mode=='qx':
                map_use = self.calibration.qx_map()
            else:
                print('ERROR: q_mode {} not recognized in linecut_qz.'.format(q_mode))
            
            
            region = np.ma.masked_where(abs(map_use-qr)>dq, self.calibration.qz_map())
            self.regions = [region]

        #Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()

        # Generate map
        M = self.calibration.qz_map().ravel()
        scale = dq # approximately 1-pixel

        Md = (M/scale + 0.5).astype(int) # Simplify the distances to closest integers
        Md -= np.min(Md)

        num_per_m = np.bincount(Md[pixel_list])
        idx = np.where(num_per_m!=0) # distances that actually have data

        x_vals = np.bincount( Md[pixel_list], weights=M[pixel_list] )[idx]/num_per_m[idx]
        I_vals = np.bincount( Md[pixel_list], weights=data[pixel_list] )[idx]/num_per_m[idx]

        line = DataLineAngle( x=x_vals, y=I_vals, x_label=x_label, y_label=y_label, x_rlabel=x_rlabel, y_rlabel=y_rlabel )

        return line

        
    def linecut_q(self, chi0, dq, x_label='q', x_rlabel='$q \, (\mathrm{\AA}^{-1})$', y_label='I', y_rlabel=r'$I (q) \, (\mathrm{counts/pixel})$', **kwargs):
        '''Returns the intensity integrated along a radial line with linewidth = 2 * dq.'''
        
        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data
        
        
        data = self.data.ravel()
        
        #QR = self.calibration.qr_map()
        QX = self.calibration.qx_map()
        QZ = self.calibration.qz_map()

        if np.isclose(chi0,0): # edge case for a vertical line
            
            pixel_list = np.where( (abs(QX.ravel()) < dq) & (mask.ravel()==1) )
            
            if 'show_region' in kwargs and kwargs['show_region']:
                region = np.ma.masked_where(abs(QX) > dq, self.calibration.q_map())
                self.regions = [region]
        
        elif np.isclose(chi0, np.pi/2) or np.isclose(chi0, -np.pi/2): # edge case for a horizontal line
            
            pixel_list = np.where( (abs(QZ.ravel()) < dq) & (mask.ravel()==1) )
            
            if 'show_region' in kwargs and kwargs['show_region']:
                region = np.ma.masked_where(abs(QZ) > dq, self.calibration.q_map())
                self.regions = [region]
                
        else:
            SLOPE = - np.tan(np.pi/2 + np.radians(chi0))
            INTCPT = dq / abs(np.sin(np.radians(chi0)))
            pixel_list = np.where( (abs(QZ.ravel() - SLOPE * QX.ravel()) < INTCPT) & (mask.ravel()==1) )
                    
            if 'show_region' in kwargs and kwargs['show_region']:
                region = np.ma.masked_where(abs(QZ - SLOPE * QX) > INTCPT, self.calibration.q_map())
                self.regions = [region]

        #Q = self.calibration.q_map().ravel()
        dq = self.calibration.get_q_per_pixel()
        
        # Generate map
        M = self.calibration.q_map().ravel()
        scale = dq # approximate 1-pixel
                
        Md = (M/scale + 0.5).astype(int) # Simplify the distances to closest integers
        Md -= np.min(Md)
        
        num_per_m = np.bincount(Md[pixel_list])
        idx = np.where(num_per_m!=0) # distances that actually have data
        
        x_vals = np.bincount( Md[pixel_list], weights=M[pixel_list] )[idx]/num_per_m[idx]
        I_vals = np.bincount( Md[pixel_list], weights=data[pixel_list] )[idx]/num_per_m[idx]
        
        line = DataLineAngle( x=x_vals, y=I_vals, x_label=x_label, y_label=y_label, x_rlabel=x_rlabel, y_rlabel=y_rlabel )
        
        return line
    
    
    def roi_q(self, qx, dqx, qz, dqz, prepend='stats_', **kwargs):
        '''Returns the intensity integrated in a box around (qx, qz).'''

        if self.mask is None:
            mask = np.ones(self.data.shape)
        else:
            mask = self.mask.data

        data = self.data.ravel()
        
        pixel_list = np.where( (abs(self.calibration.qx_map().ravel()-qx)<dqx) & (abs(self.calibration.qz_map().ravel()-qz)<dqz) & (mask.ravel()==1) )


        if 'show_region' in kwargs and kwargs['show_region']:
            
            map_use = self.calibration.q_map()
            region = np.ma.masked_where( (abs(self.calibration.qx_map()-qx)>dqx) | (abs(self.calibration.qz_map()-qz)>dqz), self.calibration.q_map())
            self.regions = [region]


        values = data[pixel_list]
        
        results = {}
        results['qx'] = qx
        results['dqx'] = dqx
        results['qz'] = qz
        results['dqz'] = dqz

        results[prepend+'max'] = np.max(values)
        results[prepend+'min'] = np.min(values)
        results[prepend+'average'] = np.average(values)
        results[prepend+'std'] = np.std(values)
        results[prepend+'N'] = len(values)
        results[prepend+'total'] = np.sum(values)
        
        results[prepend+'skew'] = stats.skew(values)
        
        results[prepend+'spread'] = results[prepend+'max'] - results[prepend+'min']
        results[prepend+'std_rel'] = results[prepend+'std'] / results[prepend+'average']        
        

        return results
    
    
        
    # Data remeshing
    ########################################
    
    def remesh_q_interpolate(self, bins_relative=1.0, method='linear', **kwargs):
        '''Converts the data from detector-space into reciprocal-space. The returned
        object has a regular grid in reciprocal-space.
        The data is converted into a (qx,qz) plane (qy contribution ignored).'''
        
        # Determine limits
        dq = self.calibration.get_q_per_pixel()/bins_relative
        qx_min = np.min(self.calibration.qx_map())
        qx_max = np.max(self.calibration.qx_map())
        qz_min = np.min(self.calibration.qz_map())
        qz_max = np.max(self.calibration.qz_map())
        
        qx = np.arange(qx_min, qx_max+dq, dq)
        qz = np.arange(qz_min, qz_max+dq, dq)
        QX, QZ = np.meshgrid(qx, qz)
        
        from scipy.interpolate import griddata

        points = np.column_stack((self.calibration.qx_map().ravel(), self.calibration.qz_map().ravel()))
        values = self.data.ravel()
        
        remesh_data = griddata(points, values, (QX, QZ), method=method)
        
        q_data = Data2DReciprocal()
        q_data.data = remesh_data
        
        return q_data


    def remesh_q_interpolate_explicit(self, qx_min=0, qx_max=1, qz_min=0, qz_max=1, method='linear', **kwargs):
        '''Converts the data from detector-space into reciprocal-space. The returned
        object has a regular grid in reciprocal-space.
        The data is converted into a (qx,qz) plane (qy contribution ignored).'''
        
        # Determine limits
        dq = kwargs['dq']
        qx = np.arange(qx_min, qx_max+dq, dq)
        qz = np.arange(qz_min, qz_max+dq, dq)
        QX, QZ = np.meshgrid(qx, qz)
        
        from scipy.interpolate import griddata

        points = np.column_stack((self.calibration.qx_map().ravel(), self.calibration.qz_map().ravel()))
        values = self.data.ravel()
        
        remesh_data = griddata(points, values, (QX, QZ), method=method)
        num_per_pixel = griddata(points, self.mask.data.ravel(), (QX, QZ), method=method)
        
        return remesh_data, num_per_pixel
    

    def remesh_q_bin(self, bins_relative=1.0, **kwargs):
        '''Converts the data from detector-space into reciprocal-space. The returned
        object has a regular grid in reciprocal-space.
        The data is converted into a (qx,qz) plane (qy contribution ignored).'''
        
        # TODO: Account for masking
        
        # Determine limits
        dq = self.calibration.get_q_per_pixel()/bins_relative
        QZ = self.calibration.qz_map().ravel()
        QX = self.calibration.qx_map().ravel()
        D = self.data.ravel()
        

        qz_min = np.min(QZ)
        qz_max = np.max(QZ)
        qx_min = np.min(QX)
        qx_max = np.max(QX)
        
        bins = [ int( abs(qz_max-qz_min)/dq ) , int( abs(qx_max-qx_min)/dq ) ]
        
        remesh_data, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=[[qz_min,qz_max], [qx_min,qx_max]], normed=False, weights=D)

        # Normalize by the binning
        num_per_bin, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=[[qz_min,qz_max], [qx_min,qx_max]], normed=False, weights=None)
        remesh_data = np.nan_to_num( remesh_data/num_per_bin )
        
        q_data = Data2DReciprocal()
        q_data.data = remesh_data
        
        q_data.x_scale = (xbins[1]-xbins[0])
        q_data.y_scale = (zbins[1]-zbins[0])
        q_data.x_axis = xbins[:-1] + (xbins[1]-xbins[0]) # convert from bin edges to bin centers
        q_data.y_axis = zbins[:-1] + (zbins[1]-zbins[0]) # convert from bin edges to bin centers
        
        
        return q_data
    
    
    def remesh_qr_bin(self, bins_relative=1.0, **kwargs):
        '''Converts the data from detector-space into reciprocal-space. The returned
        object has a regular grid in reciprocal-space.
        The data is converted into a (qx,qz) plane (qy contribution ignored).'''
        
        # TODO: Account for masking
        
        # Determine limits
        dq = self.calibration.get_q_per_pixel()/bins_relative
        
        QZ = self.calibration.qz_map().ravel()
        QX = self.calibration.qr_map().ravel()
        D = self.data.ravel()
        
        qz_min = np.min(QZ)
        qz_max = np.max(QZ)
        qx_min = np.min(QX)
        qx_max = np.max(QX)
        
        bins = [ int( abs(qz_max-qz_min)/dq ) , int( abs(qx_max-qx_min)/dq ) ]
        
        remesh_data, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=[[qz_min,qz_max], [qx_min,qx_max]], normed=False, weights=D)

        # Normalize by the binning
        num_per_bin, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=[[qz_min,qz_max], [qx_min,qx_max]], normed=False, weights=None)
        remesh_data = np.nan_to_num( remesh_data/num_per_bin )
        
        q_data = Data2DReciprocal()
        q_data.data = remesh_data
        
        q_data.x_scale = (xbins[1]-xbins[0])
        q_data.y_scale = (zbins[1]-zbins[0])
        q_data.x_axis = xbins[:-1] + (xbins[1]-xbins[0]) # convert from bin edges to bin centers
        q_data.y_axis = zbins[:-1] + (zbins[1]-zbins[0]) # convert from bin edges to bin centers
        
        
        return q_data
        
    
    def remesh_q_phi(self, bins_relative=1.0, bins_phi=None, **kwargs):
        '''Converts the data from detector-space into a (q,phi) map.'''

        # TODO: Account for masking
        
        # Determine limits
        dq = self.calibration.get_q_per_pixel()/bins_relative
        
        
        Q = self.calibration.q_map().ravel()
        q_min = np.min(np.abs(Q))
        q_max = np.max(np.abs(Q))
        q_mid = q_max-q_min
        
        if bins_phi is None:
            dphi = np.degrees(np.arctan(dq/q_mid))
            bins_phi = 360.0/dphi
            
        else:
            dphi = 360.0/bins_phi
        
        #QZ = self.calibration.qz_map().ravel()
        #QX = self.calibration.qx_map().ravel()
        PHI = self.calibration.angle_map().ravel() # degrees
        #phi_min = np.min(PHI)
        #phi_max = np.max(PHI)
        phi_min = -180.0
        phi_max = +180.0
        
        D = self.data.ravel()
        

        
        bins = [ int( abs(phi_max-phi_min)/dphi ) , int( abs(q_max-q_min)/dq ) ]
        
        remesh_data, zbins, xbins = np.histogram2d(PHI, Q, bins=bins, range=[[phi_min,phi_max], [q_min,q_max]], normed=False, weights=D)
        #num_per_bin, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=[[qz_min,qz_max], [qx_min,qx_max]], normed=False, weights=None)
        #remesh_data = np.nan_to_num( remesh_data/num_per_bin )
        
        q_phi_data = Data2DQPhi()
        q_phi_data.data = remesh_data
        
        q_phi_data.x_scale = (xbins[1]-xbins[0])
        q_phi_data.y_scale = (zbins[1]-zbins[0])
        q_phi_data.x_axis = xbins[:-1] + (xbins[1]-xbins[0]) # convert from bin edges to bin centers
        q_phi_data.y_axis = zbins[:-1] + (zbins[1]-zbins[0]) # convert from bin edges to bin centers
                
        
        return q_phi_data



    def remesh_q_bin_explicit(self, qx_min, qx_max, num_qx, qz_min, qz_max, num_qz, **kwargs):
        '''Converts the data from detector-space into reciprocal-space.
        
        This version allows explicit control of the remesh matrix, and returns
        the corresponding mask. (This can be useful, e.g. for stitching/tiling 
        images together into a combined/total reciprocal-space.)'''
        
        
        
        pixel_list = np.where( self.mask.data.ravel()==1 ) # Only consider non-masked pixels
        
        QZ = self.calibration.qz_map().ravel()[pixel_list]
        QX = self.calibration.qx_map().ravel()[pixel_list]
        D = self.data.ravel()[pixel_list]

        bins = [num_qz, num_qx]
        range = [ [qz_min,qz_max], [qx_min,qx_max] ]
        
        remesh_data, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=range, normed=False, weights=D)

        # Normalize by the binning
        num_per_bin, zbins, xbins = np.histogram2d(QZ, QX, bins=bins, range=range, normed=False, weights=None)
        #remesh_data = np.nan_to_num( remesh_data/num_per_bin )
        
        return remesh_data, num_per_bin
        
        
    # Plotting
    ########################################
        
    def plot(self, save=None, show=False, ztrim=[0.01, 0.001], **kwargs):
        
        super(Data2DScattering, self).plot(save=save, show=show, ztrim=ztrim, **kwargs)
        

    # Plot interaction
    ########################################

    def _format_coord(self, x, y):
        
        h, w = self.data.shape
        
        #xp = (x-self.calibration.x0)/self.x_scale
        #yp = (y-self.calibration.y0)/self.y_scale
        xp = (x)/self.x_scale
        yp = (h-y)/self.y_scale
        
        col = int(xp+0.5)
        row = int(yp+0.5 - 1)
        
        numrows, numcols = self.data.shape
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = self.data[row,col]
            #z = self.Z[row,col]
            #return 'x=%1.1f, y=%1.1f, z=%1.1f'%(x, y, z)
            return 'x=%g, y=%g, z=%g'%(x, y, z)
        else:
            return 'x=%g, y=%g'%(x, y)        

    
    # End class Data2DScattering(Data2D)
    ########################################





# Mask
################################################################################    
class Mask(object):
    '''Stores the matrix of pixels to be excluded from further analysis.'''
    
    def __init__(self, infile=None, format='auto'):
        '''Creates a new mask object, storing a matrix of the pixels to be 
        excluded from further analysis.'''
        
        self.data = None
        
        if infile is not None:
            self.load(infile, format=format)
        
        
    def load(self, infile, format='auto', invert=False):
        '''Loads a mask from a a file. If this object already has some masking
        defined, then the new mask is 'added' to it. Thus, one can load multiple
        masks to exlude various pixels.'''
        
        if format=='png' or infile[-4:]=='.png':
            self.load_png(infile, invert=invert)
            
        elif format=='hdf5' or infile[-3:]=='.h5' or infile[-4:]=='.hd5':
            self.load_hdf5(infile, invert=invert)
            
        else:
            print("Couldn't identify mask format for %s."%(infile))
            
            
    def load_blank(self, width, height):
        '''Creates a null mask; i.e. one that doesn't exlude any pixels.'''
        
        # TODO: Confirm that this is the correct order for x and y.
        self.data = np.ones((height, width))
        
            
    def load_png(self, infile, threshold=127, invert=False):
        '''Load a mask from a PNG image file. High values (white) are included, 
        low values (black) are exluded.'''
        
        # Image should be black (0) for excluded pixels, white (255) for included pixels
        img = PIL.Image.open(infile).convert("L") # black-and-white
        img2 = img.point(lambda p: p > threshold and 255)
        data = np.asarray(img2)/255
        data = data.astype(int)
        
        if invert:
            data = -1*(data-1)
        
        if self.data is None:
            self.data = data
        else:
            self.data *= data
        
        
    def load_hdf5(self, infile, invert=False):
        
        with h5py.File(infile, 'r') as f:
            data = np.asarray( f['mask'] )

        if invert:
            data = -1*(data-1)
        
        if self.data is None:
            self.data = data
        else:
            self.data *= data

        
    def invert(self):
        '''Inverts the mask. Can be used if the mask file was written using the
        opposite convention.'''
        self.data = -1*(self.data-1)


    # End class Mask(object)
    ########################################
    
    
    
    
    
# Calibration
################################################################################    
class Calibration(object):
    '''Stores aspects of the experimental setup; especially the calibration
    parameters for a particular detector. That is, the wavelength, detector
    distance, and pixel size that are needed to convert pixel (x,y) into
    reciprocal-space (q) value.
    
    This class may also store other information about the experimental setup
    (such as beam size and beam divergence).
    '''
    
    def __init__(self, wavelength_A=None, distance_m=None, pixel_size_um=None, incident_angle=0):
        
        self.wavelength_A = wavelength_A
        self.distance_m = distance_m
        self.pixel_size_um = pixel_size_um
        
        self.incident_angle = incident_angle
        self.sample_normal = None
        self._beam_positions = {}
        
        
        # Data structures will be generated as needed
        # (and preserved to speedup repeated calculations)
        self.clear_maps()
    
    
    # Experimental parameters
    ########################################
    
    def set_wavelength(self, wavelength_A):
        '''Set the experimental x-ray wavelength (in Angstroms).'''
        
        self.wavelength_A = wavelength_A
    
    
    def get_wavelength(self):
        '''Get the x-ray beam wavelength (in Angstroms) for this setup.'''
        
        return self.wavelength_A
    
        
    def set_energy(self, energy_keV):
        '''Set the experimental x-ray beam energy (in keV).'''
        
        energy_eV = energy_keV*1000.0
        energy_J = energy_eV/6.24150974e18
        
        h = 6.626068e-34 # m^2 kg / s
        c = 299792458 # m/s
        
        wavelength_m = (h*c)/energy_J
        self.wavelength_A = wavelength_m*1e+10
    
    
    def get_energy(self):
        '''Get the x-ray beam energy (in keV) for this setup.'''
        
        h = 6.626068e-34 # m^2 kg / s
        c = 299792458 # m/s
        
        wavelength_m = self.wavelength_A*1e-10 # m
        E = h*c/wavelength_m # Joules
        
        E *= 6.24150974e18 # electron volts
        
        E /= 1000.0 # keV
        
        return E
    
    
    def get_k(self):
        '''Get k = 2*pi/lambda for this setup, in units of inverse Angstroms.'''
        
        return 2.0*np.pi/self.wavelength_A
    
    
    def set_distance(self, distance_m):
        '''Sets the experimental detector distance (in meters).'''
        
        self.distance_m = distance_m
        
    
    def set_pixel_size(self, pixel_size_um=None, width_mm=None, num_pixels=None):
        '''Sets the pixel size (in microns) for the detector. Pixels are assumed
        to be square.'''
        
        if pixel_size_um is not None:
            self.pixel_size_um = pixel_size_um
            
        else:
            if num_pixels is None:
                num_pixels = self.width
            pixel_size_mm = width_mm*1./num_pixels
            self.pixel_size_um = pixel_size_mm*1000.0
        
        
    def set_beam_position(self, x0, y0, named=None):
        '''Sets the direct beam position in the detector images (in pixel 
        coordinates).'''
        
        if named is not None:
            self._beam_positions[named] = [x0, y0]
        else:
            self._beam_positions['default'] = [x0, y0]
            self.x0 = x0
            self.y0 = y0
        
        
    def use_beam_position(self, name):
        self.x0, self.y0 = self._beam_positions[name]
        
        
    def set_image_size(self, width, height=None):
        '''Sets the size of the detector image, in pixels.'''
        
        self.width = width
        if height is None:
            # Assume a square detector
            self.height = width
        else:
            self.height = height
    
    
    def get_q_per_pixel(self):
        '''Gets the delta-q associated with a single pixel. This is computed in
        the small-angle limit, so it should only be considered approximate.
        For instance, wide-angle detectors will have different delta-q across
        the detector face.'''
        
        if self.q_per_pixel is not None:
            return self.q_per_pixel
        
        c = (self.pixel_size_um/1e6)/self.distance_m
        twotheta = np.arctan(c) # radians
        
        self.q_per_pixel = 2.0*self.get_k()*np.sin(twotheta/2.0)
        
        return self.q_per_pixel
    
    
    def set_incident_angle(self, incident_angle=0, sample_normal=None):
        
        self.clear_maps() # Presumptively invalidate cached maps
        self.incident_angle = incident_angle
        if sample_normal is not None:
            self.sample_normal = sample_normal
    
    
    def set_angles(self, sample_normal=0, incident_angle=None):
        
        self.clear_maps() # Presumptively invalidate cached maps
        self.sample_normal = sample_normal
        if incident_angle is not None:
            self.incident_angle = incident_angle
    
    
    # Convenience methods
    ########################################
    def q_to_angle(self, q):
        '''Convert from q to angle (full scattering angle, 2theta, in degrees).'''
        kpre = 2.0*self.get_k()
        return np.degrees( 2.0*np.arcsin(q/kpre) )
    
    def angle_to_q(self, angle):
        '''Convert from scattering angle (full scattering angle, in degrees)
        to q-value (in inverse angstroms).'''
        kpre = 2.0*self.get_k()
        return kpre*np.sin(np.radians(angle/2))    
    
    def compute_qy(self, QX, QZ):
        '''Compute the (sometimes ignored!) qy component for the given (QX,QZ) matrices.'''
        
        k = self.get_k()
        
        
        alpha_f = np.arcsin(QZ/k)
        theta_f = np.arcsin( QX/(k*np.cos(alpha_f)) )
        
        QY = k*( np.cos(theta_f)*np.cos(alpha_f) - 1 )
        
        # Alternate computation:
        #QZk2 = (1 - np.square(QZ/k))
        #QY = k*( np.sqrt( 1 - np.square(QX/k)*(1/QZk2) )*np.sqrt(QZk2) - 1 )
        
        return QY
    
    # Maps
    ########################################
    
    def clear_maps(self):
        self.r_map_data = None
        self.q_per_pixel = None
        self.q_map_data = None
        self.angle_map_data = None
        
        self.qx_map_data = None
        self.qy_map_data = None
        self.qz_map_data = None
        self.qr_map_data = None

    
    def r_map(self):
        '''Returns a 2D map of the distance from the origin (in pixel units) for
        each pixel position in the detector image.'''
        
        if self.r_map_data is not None:
            return self.r_map_data

        x = np.arange(self.width) - self.x0
        y = np.arange(self.height) - self.y0
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        self.r_map_data = R
        
        return self.r_map_data
        
    
    def q_map(self):
        '''Returns a 2D map of the q-value associated with each pixel position
        in the detector image.'''

        if self.q_map_data is not None:
            return self.q_map_data
        
        c = (self.pixel_size_um/1e6)/self.distance_m
        twotheta = np.arctan(self.r_map()*c) # radians
        
        self.q_map_data = 2.0*self.get_k()*np.sin(twotheta/2.0)
        
        return self.q_map_data
        
    
    def angle_map(self):
        '''Returns a map of the angle for each pixel (w.r.t. origin).
        0 degrees is vertical, +90 degrees is right, -90 degrees is left.'''

        if self.angle_map_data is not None:
            return self.angle_map_data
        
        x = (np.arange(self.width) - self.x0)
        y = (np.arange(self.height) - self.y0)
        X,Y = np.meshgrid(x,y)
        #M = np.degrees(np.arctan2(Y, X))
        # Note intentional inversion of the usual (x,y) convention.
        # This is so that 0 degrees is vertical.
        #M = np.degrees(np.arctan2(X, Y))

        # TODO: Lookup some internal parameter to determine direction
        # of normal. (This is what should befine the angle convention.)
        M = np.degrees(np.arctan2(X, -Y))

        
        self.angle_map_data = M

        if self.sample_normal is not None:
            self.angle_map_data += self.sample_normal
        
        
        return self.angle_map_data
    
    
    def qx_map(self):
        if self.qx_map_data is not None:
            return self.qx_map_data
        
        self._generate_qxyz_maps()
        
        return self.qx_map_data    

    def qy_map(self):
        if self.qy_map_data is not None:
            return self.qy_map_data
        
        self._generate_qxyz_maps()
        
        return self.qy_map_data    

    def qz_map(self):
        if self.qz_map_data is not None:
            return self.qz_map_data
        
        self._generate_qxyz_maps()
        
        return self.qz_map_data    
    
    def qr_map(self):
        if self.qr_map_data is not None:
            return self.qr_map_data

        self._generate_qxyz_maps()

        return self.qr_map_data



    def _generate_qxyz_maps(self):
        
        # Conversion factor for pixel coordinates
        # (where sample-detector distance is set to d = 1)
        c = (self.pixel_size_um/1e6)/self.distance_m
        
        x = np.arange(self.width) - self.x0
        y = np.arange(self.height) - self.y0
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        #twotheta = np.arctan(self.r_map()*c) # radians
        theta_f = np.arctan2( X*c, 1 ) # radians
        #alpha_f_prime = np.arctan2( Y*c, 1 ) # radians
        alpha_f = np.arctan2( Y*c*np.cos(theta_f), 1 ) # radians
        
        cos_inc = np.cos(np.radians(self.incident_angle))
        sin_inc = np.sin(np.radians(self.incident_angle))
        self.qx_map_data = self.get_k()*np.sin(theta_f)*np.cos(alpha_f)
        self.qy_map_data = self.get_k()*( np.cos(theta_f)*np.cos(alpha_f) - cos_inc ) # TODO: Check sign
        self.qz_map_data = -1.0*self.get_k()*( np.sin(alpha_f) + sin_inc ) 
        

        if self.sample_normal is not None:
            s = np.sin(np.radians(self.sample_normal))
            c = np.cos(np.radians(self.sample_normal))
            self.qx_map_data, self.qz_map_data = c*self.qx_map_data - s*self.qz_map_data, s*self.qx_map_data + c*self.qz_map_data
        
        self.qr_map_data = np.sign(self.qx_map_data)*np.sqrt(np.square(self.qx_map_data) + np.square(self.qy_map_data))




        
    
    
    # End class Calibration(object)
    ########################################
    
    

    
    
# Data2DReciprocal
################################################################################
class Data2DReciprocal(Data2D):
    '''Represents a 2D slice in reciprocal space.
    Note that Data2DScattering requires a corresponding Calibration object (to
    convert from pixel to q), whereas Data2DReciprocal is implicitly 'flat', and
    contains the required values for etablishing the q-axes.'''
    
    def __init__(self, infile=None, format='auto', name=None, **kwargs):
        
        iargs = {
                'x_label' : 'qx',
                'x_rlabel' : '$q_x \, (\mathrm{\AA^{-1}})$',
                'y_label' : 'qz',
                'y_rlabel' : '$q_z \, (\mathrm{\AA^{-1}})$',
                 }
        iargs.update(kwargs)
        
        super(Data2DReciprocal, self).__init__(infile=None, format=format, name=name, **iargs)
    
        self.plot_args = { 'rcParams': {'axes.labelsize': 48,
                                        'xtick.labelsize': 33,
                                        'ytick.labelsize': 33,
                                        },
                            }    
        
        
    
    def xy_axes(self):
        
        return self.x_axis, self.y_axis


    def plot(self, save=None, show=False, ztrim=[0.01, 0.01], size=10.0, plot_buffers=[0.25,0.05,0.25,0.05], **kwargs):
        '''Plots the data.
        
        Parameters
        ----------
        save : str
            Set to 'None' to avoid saving to disk. Provide filename to save.
        show : bool
            Set to true to open an interactive window.
        ztrim : [float, float]
            Specify how to auto-set the z-scale. The floats indicate how much of
            the z-scale to 'trim' (relative units; i.e. 0.05 indicates 5%).
        '''  
        
        self._plot(save=save, show=show, ztrim=ztrim, size=size, plot_buffers=plot_buffers, **kwargs)
        
        
        
    def _plot_extra(self, **plot_args):
        self.ax.get_yaxis().set_tick_params(which='both', direction='out')
        self.ax.get_xaxis().set_tick_params(which='both', direction='out')       
        
    def _format_coord(self, x, y):
        
        
        
        # TODO: Avoid IndexError
        xp = np.where(self.x_axis>x)[0][0]
        yp = np.where(self.y_axis>y)[0][0]
        
        # TODO: Not pixel-perfect
        col = int(xp)
        row = int(yp)
        
        numrows, numcols = self.data.shape
        row = numrows-row-1
        if col>=0 and col<numcols and row>=0 and row<numrows:
            # TODO: For some reason, this accesses data outside the array.
            
            z = self.data[row,col]
            #z = self.Z[row,col]
            #return 'x=%1.1f, y=%1.1f, z=%1.1f'%(x, y, z)
            #return 'x=%g, y=%g, z=%g'%(x, y, z)
            return 'x=%g, y=%g, z=%g'%(col, row, z)
        else:
            return 'x=%g, y=%g'%(x, y)        
        
    # End class Data2DReciprocal(Data2D)
    ########################################
        
        
        
# Data2DQPhi
################################################################################
class Data2DQPhi(Data2D):
    '''Represents a 2D (q,phi) map.
    '''
    
    def __init__(self, infile=None, format='auto', name=None, **kwargs):
        
        iargs = {
                'x_label' : 'q',
                'x_rlabel' : '$q \, (\mathrm{\AA^{-1}})$',
                'y_label' : 'phi',
                'y_rlabel' : '$\phi \, (^{\circ})$',
                 }
        iargs.update(kwargs)
        
        super(Data2DQPhi, self).__init__(infile=None, format=format, name=name, **iargs)
    
        self.plot_args = { 'rcParams': {'axes.labelsize': 48,
                                        'xtick.labelsize': 33,
                                        'ytick.labelsize': 33,
                                        },
                            }    
    
    def xy_axes(self):
        
        return self.x_axis, self.y_axis


    def plot(self, save=None, show=False, ztrim=[0.01, 0.01], size=10.0, plot_buffers=[0.25,0.05,0.25,0.05], **kwargs):
        '''Plots the data.
        
        Parameters
        ----------
        save : str
            Set to 'None' to avoid saving to disk. Provide filename to save.
        show : bool
            Set to true to open an interactive window.
        ztrim : [float, float]
            Specify how to auto-set the z-scale. The floats indicate how much of
            the z-scale to 'trim' (relative units; i.e. 0.05 indicates 5%).
        '''  
        
        self._plot(save=save, show=show, ztrim=ztrim, size=size, plot_buffers=plot_buffers, **kwargs)
        
        
    def _plot_extra(self, **plot_args):
        '''This internal function can be over-ridden in order to force additional
        plotting behavior.'''
        
        self.ax.set_aspect('auto')
        
    def _format_coord(self, x, y):
        
        
        
        # TODO: Avoid IndexError
        xp = np.where(self.x_axis>x)[0][0]
        yp = np.where(self.y_axis>y)[0][0]
        
        # TODO: Not pixel-perfect
        col = int(xp)
        row = int(yp)
        
        numrows, numcols = self.data.shape
        row = numrows-row-1
        if col>=0 and col<numcols and row>=0 and row<numrows:
            # TODO: For some reason, this accesses data outside the array.
            
            z = self.data[row,col]
            #z = self.Z[row,col]
            #return 'x=%1.1f, y=%1.1f, z=%1.1f'%(x, y, z)
            #return 'x=%g, y=%g, z=%g'%(x, y, z)
            return 'x=%g, y=%g, z=%g'%(col, row, z)
        else:
            return 'x=%g, y=%g'%(x, y)        
        
    # End class Data2DQPhi(Data2D)
    ########################################        
