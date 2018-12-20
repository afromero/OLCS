#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 17:06:16 2018

@author: romerowo
"""

import numpy as np

class Geometry:
    def __init__(self, radius_km, num_events, zen_angle_deg=None):
        self.num_events = num_events
        self.radius_km = radius_km
        self.Area = np.pi*radius_km**2
        self.AOmega = np.pi * self.Area
        self.get_random_interaction_point()
        self.zen_angle_deg = zen_angle_deg
        self.get_random_direction(zen_angle_deg=self.zen_angle_deg)
        
    ####################################################################################

    def get_random_interaction_point(self):
        self.R_pos   = np.sqrt( np.random.uniform(0., self.radius_km**2, self.num_events) )
        self.phi_pos = np.random.uniform(0., 2.*np.pi, self.num_events)
        self.x_pos   = self.R_pos*np.cos(self.phi_pos)
        self.y_pos   = self.R_pos*np.sin(self.phi_pos)
        self.z_pos   = 1.222*np.ones(len(self.R_pos))
        
    ####################################################################################

    def get_random_direction(self, zen_angle_deg=None):
        if zen_angle_deg == None:
            self.cos_th_CR = np.sqrt(np.random.uniform(0.,1., self.num_events) )
            self.th_CR = np.arccos(self.cos_th_CR)
        else:
            self.th_CR = zen_angle_deg * np.pi/180. * np.ones(self.num_events) 
            self.cos_th_CR = np.cos(self.th_CR)
        self.phi_CR    = np.random.uniform(0.,2.*np.pi, self.num_events)
        self.k_x       = np.sqrt(1.-self.cos_th_CR**2)*np.cos(self.phi_CR)
        self.k_y       = np.sqrt(1.-self.cos_th_CR**2)*np.sin(self.phi_CR)
        self.k_z       = self.cos_th_CR
        
        