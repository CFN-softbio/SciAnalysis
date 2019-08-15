#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.DataGonio` - Handling wide-angle area detectors
================================================
.. module:: SciAnalysis.XSAnalysis
   :synopsis: Provides base classes for doing analysis of x-ray data with tilted detectors
'''

################################################################################
#  This code defines some baseline objects for analysis of area detectors with
# that are swept about an arc centered on the sample.
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################

from .Data import *


# CalibrationGonio
################################################################################    
class CalibrationGonio(Calibration):
    """
    The geometric claculations used here are described:
    http://gisaxs.com/index.php/Geometry:WAXS_3D
    
    """    
    
    # Experimental parameters
    ########################################
    
    def set_angles(self, det_phi_g=0., det_theta_g=0.):
        
        self.clear_maps() # Any change to the detector position will presumptively invalidate cached maps
        
        self.det_phi_g = det_phi_g
        self.det_theta_g = det_theta_g
    
    def get_ratioDw(self):
        
        width_mm = self.width*self.pixel_size_um/1000.
        return self.distance_m/(width_mm/1000.)


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
        """
        The geometric claculations used here are described:
        http://gisaxs.com/index.php/Geometry:WAXS_3D
        
        """    
        
        d = self.distance_m
        pix_size = self.pixel_size_um/1e6
        phi_g = np.radians(self.det_phi_g)
        theta_g = np.radians(self.det_theta_g)
        
        xs = (np.arange(self.width) - self.x0)*pix_size
        ys = (np.arange(self.height) - self.y0)*pix_size
        #ys = ys[::-1]

        X_c, Y_c = np.meshgrid(xs, ys)
        Dprime = np.sqrt( np.square(d) + np.square(X_c) + np.square(Y_c) )
        k_over_Dprime = self.get_k()/Dprime
        
        
        qx_c = k_over_Dprime*( X_c*np.cos(phi_g) - np.sin(phi_g)*(d*np.cos(theta_g) - Y_c*np.sin(theta_g)) )
        qy_c = k_over_Dprime*( X_c*np.sin(phi_g) + np.cos(phi_g)*(d*np.cos(theta_g) - Y_c*np.sin(theta_g)) - Dprime )
        qz_c = -1*k_over_Dprime*( d*np.sin(theta_g) + Y_c*np.cos(theta_g) )

        qr_c = np.sqrt(np.square(qx_c) + np.square(qy_c))        
        q_c = np.sqrt(np.square(qx_c) + np.square(qy_c) + np.square(qz_c))
        
        
        
        
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
        
        
        self.qx_map_data = self.get_k()*np.sin(theta_f)*np.cos(alpha_f)
        self.qy_map_data = self.get_k()*( np.cos(theta_f)*np.cos(alpha_f) - 1 ) # TODO: Check sign
        self.qz_map_data = -1.0*self.get_k()*np.sin(alpha_f)
        
        self.qr_map_data = np.sign(self.qx_map_data)*np.sqrt(np.square(self.qx_map_data) + np.square(self.qy_map_data))

        
        self.qx_map_data = qx_c
        self.qy_map_data = qy_c
        self.qz_map_data = qz_c
        self.q_map_data = q_c
        
        
        
        
