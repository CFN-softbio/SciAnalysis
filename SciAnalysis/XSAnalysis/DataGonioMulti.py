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
from .DataGonio import *


# CalibrationGonio
################################################################################    
class CalibrationGonioMulti(CalibrationGonio):
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
        if self.angle_map_data is None:
            self._generate_qxyz_maps()
        
        return self.angle_map_data
    
    
    def combine_waxs(self, img, img0, img1, img2):
        # print('##### IN combine_waxs........')       
        img = img*0.0  
        print(img.shape)
        # print(img0.shape)
        if 0:
            img[0:195, :] = img0
            img[212:407, :] = img1
            img[-195::, :] = img2
        else:
            img[:, 0:195] = img0
            img[:, 212:407] = img1
            img[:, -195::] = img2

        return img

    def _generate_qxyz_maps_0(self):
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

        self.angle_map_data = np.degrees(np.arctan2(X, Y))
        
        


    def _generate_qxyz_maps(self):
        """
        The geometric claculations used here are described:
        http://gisaxs.com/index.php/Geometry:WAXS_3D
        
        """    

        print('##### IN Multi _generate_qxyz_maps........')     

        self._generate_qxyz_maps_0() 


        phi_center = 0.6
        phi_tilt = 7.6

        cali0 = CalibrationGonio(wavelength_A=self.wavelength_A)
        # cali0.set_image_size(1475, height=195) 
        cali0.set_image_size(195, height=1475)
        cali0.set_beam_position(self.x0, self.y0)
        cali0.set_distance(self.distance_m)
        cali0.set_pixel_size(self.pixel_size_um)
        cali0.set_angles(det_phi_g=self.det_phi_g + phi_center+ phi_tilt, det_theta_g=self.det_theta_g  ) # 6.1, 7.37, 7.8, 7.59
        cali0._generate_qxyz_maps() 
  

        cali1 = CalibrationGonio(wavelength_A=self.wavelength_A)
        cali1.set_image_size(195, height=1475)
        cali1.set_beam_position(self.x0, self.y0)
        cali1.set_distance(self.distance_m)
        cali1.set_pixel_size(self.pixel_size_um)
        cali1.set_angles(det_phi_g=self.det_phi_g + phi_center, det_theta_g=self.det_theta_g  )
        cali1._generate_qxyz_maps() 


        cali2 = CalibrationGonio(wavelength_A=self.wavelength_A)
        cali2.set_image_size(195, height=1475)
        cali2.set_beam_position(self.x0, self.y0)
        cali2.set_distance(self.distance_m)
        cali2.set_pixel_size(self.pixel_size_um)
        cali2.set_angles(det_phi_g=self.det_phi_g + phi_center - phi_tilt , det_theta_g=self.det_theta_g ) #7.8
        cali2._generate_qxyz_maps() 
        

        # self.q_map_data = self.q_map_data*0.0
        # self.q_map_data[0:195, :] = cali0.q_map_data
        # self.q_map_data[212:407, :] = cali1.q_map_data
        # self.q_map_data[-195::, :] = cali2.q_map_data

        # self.q_map_data = self.combine_waxs(img=self.q_map_data, img0=cali0.q_map_data, img1=cali1.q_map_data, img2=cali2.q_map_data)
        # self.qr_map_data = self.combine_waxs(img=self.qr_map_data, img0=cali0.qr_map_data,  img1=cali1.qr_map_data, img2=cali2.qr_map_data)


        qx_c = self.combine_waxs(img=self.qx_map_data, img0=cali0.qx_map_data,  img1=cali1.qx_map_data, img2=cali2.qx_map_data)
        qy_c = self.combine_waxs(img=self.qy_map_data, img0=cali0.qy_map_data,  img1=cali1.qy_map_data, img2=cali2.qy_map_data)
        qz_c = self.combine_waxs(img=self.qz_map_data, img0=cali0.qz_map_data,  img1=cali1.qz_map_data, img2=cali2.qz_map_data)
 
        self.qx_map_data = qx_c
        self.qy_map_data = qy_c
        self.qz_map_data = qz_c

        self.qr_map_data = np.sqrt(np.square(qx_c) + np.square(qy_c))        
        self.q_map_data = np.sqrt(np.square(qx_c) + np.square(qy_c) + np.square(qz_c))

        if 0:
            x = np.arange(1475) - self.x0
            y = np.arange(619) - self.y0
        else:
            x = np.arange(619) - self.x0
            y = np.arange(1475) - self.y0
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)

        self.angle_map_data = np.degrees(np.arctan2(X, Y))

        
        # self.qr_map_data         
        # self.qx_map_data = qx_c
        # self.qy_map_data = qy_c
        # self.qz_map_data = qz_c
        # self.q_map_data = q_c       
