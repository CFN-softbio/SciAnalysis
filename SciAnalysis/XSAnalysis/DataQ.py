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
    The geomatric claculations used here are described in Yang, J Synch Rad (2013) 20, 211–218
    http://dx.doi.org/10.1107/S0909049512048984
    
    """    
    
    
    def __init__(self, infile=None):
        
        self.clear_maps()
        
        if infile is not None:
            self.load_axes(infile)
            
    def load_axes(self, infile):
        
        self.qxs, self.qzs = pickle.load(open(infile,'rb'))
        
        dqx = np.average(abs(self.qxs[:-1] - self.qxs[1:]))
        dqz = np.average(abs(self.qzs[:-1] - self.qzs[1:]))
        
        self.q_per_pixel = (dqx + dqz)*0.5
        
    
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
        
        QX, QZ = np.meshgrid(self.qxs, self.qzs)
        
        self.qx_map_data = QX
        self.qz_map_data = QZ
        
        self.q_map_data = np.sqrt( np.square(QX) + np.square(QZ) )
        
        #self.qy_map_data = np.zeros_like(QX)
        #self.qr_map_data = np.zeros_like(QX)
        
    