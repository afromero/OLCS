#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 17:09:54 2018

@author: romerowo
"""

import numpy as np
from geometry import Geometry
from scipy.interpolate import Akima1DInterpolator

class Xmax_calc:
    def __init__(self, d_lam=0.1, detector_altitude_km=1.222):
        self.Earth_radius   = 6371.      # km
        self.detector_altitude_km = detector_altitude_km
        self.d_lam = d_lam
        # Initialize US standard atmosphere 
        self.h_km       =             [-1, 0, 1, 2,  3, 4,  5,  6,  7,  8,  9,  10,  15,  20,  25,  30,  40,  50,  60,  70,  80, 81., 1.e20] 
        self.dens_kg_m3 =0.1*np.array([13.47, 12.25, 11.12, 10.07, 9.093, 8.194,7.364, 6.601, 5.900, 5.258, 4.671, 4.135, 1.948, 0.8891, 0.4008, 0.1841, 0.03996,  0.01027, 0.003097, 0.0008283, 0.0001846,0.,0.])
        # self.get_density = interp1d(self.h_km, self.dens_kg_m3)
        self.get_density = Akima1DInterpolator(self.h_km, self.dens_kg_m3)
        
        self.lam_array = np.arange(0.,80., self.d_lam)[::-1]
        
        self.proton_log10_energy = np.array([14.,   15.,   16.,   17.,   18.,   19.,   20.,   21.])
        self.proton_Xmax         = np.array([5000., 5800., 6400., 6900., 7600., 8000., 8450., 8800.]) # in kg/m^2
        self.get_proton_Xmax     = Akima1DInterpolator(self.proton_log10_energy, self.proton_Xmax)
        
        #XM_log10_energy = np.array()


    def get_path_altitude(self, lam, x_init, y_init, z_init, k_x, k_y, k_z):
        x = x_init + lam * k_x
        y = y_init + lam * k_y
        z = z_init + lam * k_z
        return np.sqrt(x**2 + y**2 + z**2) - self.Earth_radius
        
    def get_X(self, x_init, y_init, z_init, k_x, k_y, k_z, plots=False):
        self.alt_array = self.get_path_altitude(self.lam_array, x_init, y_init, z_init, k_x, k_y, k_z)
        self.density   = self.get_density(self.alt_array)
        self.depth     = np.cumsum( self.density ) * (self.d_lam) * 1.e3 #factor of 1000. for atmosphere in kg/m^3 but distances are in km
        # in kg/m^2 = 0.1 g/cm^2
        if(plots):
            plt.figure()
            plt.subplot(221)
            plt.plot(self.lam_array, self.density, '.')
            plt.subplot(222)
            plt.plot(self.lam_array, self.alt_array, '.')
            plt.subplot(223)
            plt.plot(self.lam_array, self.depth, '.')
            plt.subplot(224)
            plt.plot(self.alt_array, self.depth, '.')
        
    def get_lam_Xmax(self, Xmax_kg_m2, x_init, y_init, z_init, k_x, k_y, k_z, plots=False):
        depth = self.get_X(x_init, y_init, z_init, k_x, k_y, k_z, plots=plots)
        idx = np.argmin(np.abs(self.depth - Xmax_kg_m2))
        #print self.lam_array[idx], self.depth[idx]
        return self.lam_array[idx]
    
    def get_position(self, lam, x_pos, y_pos, z_pos, shower_k_x, shower_k_y, shower_k_z):
        x = x_pos + lam * shower_k_x
        y = y_pos + lam * shower_k_y
        z = z_pos + lam * shower_k_z
        return x,y,z
    
    def get_Xmax_position(self, log10_proton_energy, geom):
        proton_Xmax = self.get_proton_Xmax(log10_proton_energy)
        self.Xmax_dist = np.zeros((len(geom.x_pos)))
        self.x_max     = np.zeros(len(geom.x_pos))
        self.y_max     = np.zeros(len(geom.x_pos))
        self.z_max     = np.zeros(len(geom.x_pos))
        for k in range(0,len(geom.x_pos)):
            self.Xmax_dist[k] = self.get_lam_Xmax(proton_Xmax, 
                                    geom.x_pos[k], geom.y_pos[k], self.Earth_radius+self.detector_altitude_km, # array altitude is 1222 m a.s.l. 
                                    geom.k_x[k], geom.k_y[k], geom.k_z[k], 
                                    plots=False)
        self.x_max, self.y_max, self.z_max = self.get_position(self.Xmax_dist,
                                                          geom.x_pos, geom.y_pos, self.Earth_radius+self.detector_altitude_km,
                                                          geom.k_x, geom.k_y, geom.k_z)
