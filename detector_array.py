#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 12:50:02 2018

@author: romerowo
"""

import pandas as pd
import numpy as np

class Detector_Array:
    def __init__(self, mode='file', fnm='OV-LWA_256_32_64_ant_config_B_2018JAN08.xls'):
        self.df = None
        if mode=='file':
            self.df = pd.read_excel(fnm, sheet_name='Sheet1', header = None)
            self.x = self.df[0][0:256]*1.e-3
            self.y = self.df[1][0:256]*1.e-3
            #self.x = self.df[0][0:251]*1.e-3
            #self.y = self.df[1][0:251]*1.e-3
            self.z = 1.222*np.ones(len(self.x)) # elevation of OVRO-LWA
        if mode=='Ryan':
            self.x = []
            self.y = []
            self.z = []
            self.z0 = 1.222
            for line in file('OVRO-LWA_256_ant_pos_Ryan.txt'):
                #print line
                self.x.append(float(line.split(',')[0]))
                self.y.append(float(line.split(',')[1]))
                self.z.append(float(line.split(',')[2]))
            self.x = np.array(self.x)*1.e-3
            self.y = np.array(self.y)*1.e-3
            self.z = self.z0 + np.array(self.z)*1.e-3
        
    ####################################################################################

    def get_distances_and_view_angles(self, XmaxC, Geom, event=0):
        self.th_view   = np.zeros(len(self.x))
        self.dist_Xmax = np.zeros(len(self.x))
        dx1 = XmaxC.x_max[event] - Geom.x_pos[event]
        dy1 = XmaxC.y_max[event] - Geom.y_pos[event]
        dz1 = (XmaxC.z_max[event]- XmaxC.Earth_radius) - Geom.z_pos[event]
        r1 = np.sqrt( dx1**2 + dy1**2 + dz1**2 )
        dx2 = XmaxC.x_max[event] - self.x
        dy2 = XmaxC.y_max[event] - self.y
        dz2 = (XmaxC.z_max[event]- XmaxC.Earth_radius) - self.z
        r2  = np.sqrt( dx2**2 + dy2**2 + dz2**2 )
        dx3 = self.x - Geom.x_pos[event]
        dy3 = self.y - Geom.y_pos[event]
        dz3 = self.z - Geom.z_pos[event]
        r3 = np.sqrt( dx3**2 + dy3**2 + dz3**2 )
        #cos_view_angle    = (r1**2 + r2**2 - r3**2) / (2.*r1*r2)
        cos_view_angle = (dx1*dx2 + dy1*dy2 + dz1*dz2) / (r1*r2)
        #print cos_view_angle
        self.th_view   = np.arccos(cos_view_angle)
        self.dist_Xmax = np.array(r2) 
