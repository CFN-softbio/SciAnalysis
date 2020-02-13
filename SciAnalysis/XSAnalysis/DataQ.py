#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.DataQ` - Handling data that has already been saved in q-space
================================================
.. module:: SciAnalysis.XSAnalysis
   :synopsis: Provides base classes for doing analysis of x-ray data with data
   that has already been converted to reciprocal-space (q-space)
'''

import pickle
from .Data import *


# CalibrationQ
################################################################################    
class CalibrationQ(Calibration):
    """
    This object (and variants thereof) is intended for opening up data that has already
    been mapped into q-space. This is typically obtained by some other conversion operation.
    For instance, if different detector images are merged into q-space (accounting for 
    detector position), then the resultant saved matrix should be "calibrated" using
    and calibration object like this.
    
    """    
    
    
    def __init__(self, infile=None):
        
        self.clear_maps()
        
        if infile is not None and 'axes' in infile:
            self.load_axes(infile)
            self._load_mode = 'axes'
        else:
            self.load_maps(infile)
            
            
    def load_axes(self, infile):
        
        #self.qxs, self.qzs = pickle.load(open(infile,'rb'))
        #self.qxs, self.qzs = np.load(infile)
        
        npzfile = np.load(infile)
        self.qxs = npzfile['qxs']
        self.qzs = npzfile['qzs']
        
        dqx = np.average(abs(self.qxs[:-1] - self.qxs[1:]))
        dqz = np.average(abs(self.qzs[:-1] - self.qzs[1:]))
        
        self.q_per_pixel = (dqx + dqz)*0.5
        
    def load_maps(self, infile):
        npzfile = np.load(infile)
        self.qx_map_data = npzfile['QX']
        self.qy_map_data = npzfile['QY']
        self.qz_map_data = npzfile['QZ']
        self._generate_qxyz_maps()
        
        self.q_per_pixel = abs(self.qx_map_data[0,0] - self.qx_map_data[1,1])

    
    def get_q_per_pixel(self):
        
        if self.q_per_pixel is not None:
            return self.q_per_pixel
        
        c = (self.pixel_size_um/1e6)/self.distance_m
        twotheta = np.arctan(c) # radians
        
        self.q_per_pixel = 2.0*self.get_k()*np.sin(twotheta/2.0)
        
        return self.q_per_pixel
    
    
    # Maps
    ########################################
    
    def q_map(self):
        if self.q_map_data is None:
            self._generate_qxyz_maps()
        
        return self.q_map_data

    
    def angle_map(self):
        if self.angle_map_data is not None:
            self._generate_qxyz_maps()
        
        return self.angle_map_data

    
    def _generate_qxyz_maps(self):
        
        if self._load_mode=='axes':
        
            QX, QZ = np.meshgrid(self.qxs, self.qzs)
            
            self.qx_map_data = QX
            self.qz_map_data = QZ
            
            self.q_map_data = np.sqrt( np.square(QX) + np.square(QZ) )
            
            #self.qy_map_data = np.zeros_like(QX)
            #self.qr_map_data = np.zeros_like(QX)
            
        else:
            
            self.q_map_data = np.sqrt( np.square(QX) + np.square(QY) + np.square(QZ) )
            self.qr_map_data = np.sqrt( np.square(QX) + np.square(QY) )
            
            
            
        
    
