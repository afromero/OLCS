#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 14:09:53 2018

@author: romerowo
"""

import numpy as np
from xmax_calc import Xmax_calc
from pylab import *


class Radio_Emission:
    def __init__(self, plots=False):
        self.XMC = Xmax_calc()
        '''
        Digitized data points from arXiv:1402.3504v2
        '''
        d_core_km_25 = [-0.4235, -0.2948, -0.2093, -0.1669, -0.146, -0.1182, -0.0685, -0.0114, 0.0681, 0.1119, 0.1412, 0.1778, 0.2359, 0.3443, 0.5023, 0.5023]
        E_field_25 = [2.572e-06, 2.0588e-05, 4.9099e-05, 7.7979e-05, 0.000104607, 0.00013911, 0.000164991, 0.000176249, 0.000151884,0.000120764, 9.7143e-05, 6.7147e-05, 3.7904e-05, 9.043e-06, 1.938e-06, 1.938e-06]
        
        d_core_km_45  = [-0.9689, -0.725, -0.5959, -0.5244, -0.4459, -0.4105, -0.3534, -0.2894, -0.2325, -0.2002, -0.1539, -0.1236, -0.0601, -0.0174, 0.0529, 0.107, 0.1523, 0.1841, 0.212, 0.2528, 0.2897, 0.3607, 0.4322, 0.5253, 0.6273, 0.8543, 1.0624]
        E_field_45 = [1.754e-06, 3.659e-06, 9.676e-06, 2.0934e-05, 4.0444e-05, 6.0699e-05, 7.2331e-05, 9.7089e-05, 0.000118846, 0.000127679, 0.000133106, 0.000136193, 0.00013124, 0.000127498, 0.000131727, 0.000136458, 0.000133865, 0.000126546, 0.000120061, 0.000104558, 8.1819e-05, 5.7187e-05, 3.6001e-05, 1.7998e-05, 7.983e-06, 1.468e-06, 1.564e-06 ]
        
        d_core_km_60  = [-1.9665, -1.4212, -1.0627, -0.8623, -0.7192, -0.5762, -0.3897, -0.2244, -0.0593, 0.199, 0.3854, 0.5292, 0.6155, 0.7452, 0.9251, 1.0906, 1.3779, 1.8445]
        E_field_60 = [1.629e-06, 5.072e-06, 1.9742e-05, 4.3017e-05, 6.216e-05, 8.3928e-05, 8.8076e-05, 8.1722e-05, 7.7618e-05, 8.29e-05, 8.9298e-05, 8.1441e-05, 7.1702e-05, 5.1093e-05, 3.0866e-05, 1.4761e-05, 4.297e-06, 1.356e-06 ]
        
        d_core_km_70 = [-2.8165, -2.6043, -2.3777, -2.1799, -1.9748, -1.759, -1.5216, -1.259, -0.8957, -0.5432, -0.2338, 0.018, 0.3165, 0.5935, 0.8669, 1.1223, 1.3705, 1.5971, 1.8525, 1.9928, 2.1619, 2.4065, 2.7986]
        E_field_70=[7.204e-06, 1.1773e-05, 1.8419e-05, 2.5815e-05, 3.3778e-05, 4.1366e-05, 4.877e-05, 5.3159e-05, 5.3416e-05, 5.0841e-05, 4.9201e-05, 4.906e-05, 4.9682e-05, 5.2376e-05, 5.3937e-05, 5.4929e-05, 5.1391e-05, 4.4264e-05, 3.45e-05, 2.7546e-05, 2.2295e-05, 1.4039e-05, 6.566e-06]
        
        # x is (+N|-S), y (+E|-W), z (+D|-U)
        # units are Teslas
        B_field_x = +22506.5e-9
        B_field_y =  +5012.5e-9
        B_field_z = +42700.9e-9
        
        # I am assuming the antenna positions are x (+E|-W), y (+N|-S), z(+D|-U)
        # NOTE: This needs to be checked
        self.B_x = +B_field_y
        self.B_y = +B_field_x
        self.B_z = -B_field_z
        
        '''
        I want to parametrize this as a double-Gaussian beam model with strength 
        depending on distance to Xmax
        '''
        
        
        if(plots):
            x = np.arange(-3.5, 3.51, 0.01)
            
            figure(1)
            plot(d_core_km_25, E_field_25, 'b-', lw=3, alpha=0.5)
            plot(x, self.double_gaussian(x, 0.0, 0.15, 9.e-5), 'b--', label=r'$\theta_z=25^\circ$')
            
            plot(d_core_km_45, E_field_45, 'r-', lw=3, alpha=0.5)
            plot(x, self.double_gaussian(x, 0.175, 0.155, 1.22e-4), 'r--', label=r'$\theta_z=45^\circ$')
            
            plot(d_core_km_60, E_field_60, 'g-', lw=3, alpha=0.5)
            plot(x, self.double_gaussian(x, 0.45, 0.37, 8.3e-5), 'g--', label=r'$\theta_z=60^\circ$')
            
            plot(d_core_km_70, E_field_70, 'k-', lw=3, alpha=0.5)
            plot(x, self.double_gaussian(x, 1.15, 0.95, 5.1e-5), 'k--', label=r'$\theta_z=70^\circ$')
            legend(loc=0)
            xlim(-2.5, 2.5)
            xlabel('Distance to Core, km')
            ylabel('$E_\mathrm{pk}$ (30-80 MHz), V/m')
            

            figure(2)
            th_smooth = np.arange(0., 10.*pi/180., 0.001)
            th_view_25 = self.distances_angles(d_core_km_25, 25.*pi/180., 1.4)
            sep_25, d25 = self.angle_sep(th_smooth, 25.*pi/180., 1.4)
            plot(th_view_25*180./pi, E_field_25, 'b-', lw=3, alpha=0.4)
            plot(th_smooth*180./pi, self.double_gaussian(sep_25, 0.0, 0.15, 9.e-5), 'b--', lw=2, label=r'$\theta_z=25^\circ$')
            
            th_view_45 = self.distances_angles(d_core_km_45, 45.*pi/180., 1.4)
            sep_45, d45 = self.angle_sep(th_smooth, 45.*pi/180., 1.4)
            plot(th_view_45*180./pi, E_field_45, 'r-', lw=3, alpha=0.4)
            plot(th_smooth*180./pi, self.double_gaussian(sep_45, 0.175, 0.155, 1.22e-4), 'r--', lw=2, label=r'$\theta_z=45^\circ$')
            
            th_view_60 = self.distances_angles(d_core_km_60, 60.*pi/180., 1.4)
            sep_60, d60 = self.angle_sep(th_smooth, 60.*pi/180., 1.4)
            plot(th_view_60*180./pi, E_field_60, 'g-', lw=3, alpha=0.4)
            plot(th_smooth*180./pi, self.double_gaussian(sep_60, 0.45, 0.37, 8.3e-5), 'g--', lw=2, label=r'$\theta_z=60^\circ$')
            
            th_view_70 = self.distances_angles(d_core_km_70, 70.*pi/180., 1.4)
            sep_70, d70 = self.angle_sep(th_smooth, 70.*pi/180., 1.4)
            plot(th_view_70*180./pi, E_field_70, 'k-', lw=3, alpha=0.4)
            plot(th_smooth*180./pi, self.double_gaussian(sep_70, 1.15, 0.95, 5.1e-5), 'k--', lw=2, label=r'$\theta_z=70^\circ$')
            xlabel('View Angle, deg')
            ylabel('$E_\mathrm{pk}$ (30-80 MHz), V/m')
        
            figure()
            t = np.arange(0.,90.)
            th_vals = [25., 45., 60., 70.]
            sep_vals = [0.01, 0.175, 0.45, 1.15]
            sig_vals = [0.15, 0.155, 0.37, 0.95]
            norm_vals = [9.e-5, 1.22e-4, 8.3e-5, 5.1e-5]
            plot(th_vals, sep_vals, 'o-')
            plot(t, 1.1*(t/70.)**4, lw=2)
            xlabel('Zenith Angle, deg')
            ylabel('Gaussian Peak, km')
            xticks(np.arange(0.,91., 10.))
            
            figure()
            semilogy(th_vals, sig_vals, 'o-')
            plot(t, 0.15+0.85*(t/70.)**8, lw=2)
            xlabel('Zenith Angle, deg')
            ylabel('Gaussian $\sigma$, km')
            xticks(np.arange(0.,91., 10.))
            ylim(0.1, 10.)
            
            figure()
            plot(th_vals, norm_vals, 'o-')
            plot(t[t<=45.], 5.5e-5+6.6e-5*(t[t<=45.]/45.), lw=2)
            plot(t[t>=45.], 1.2e-4*(90.-t[t>=45.])/45., color='orange', lw=2)
            xlim(0.,90.)
            ylim(0.,1.5e-4)
            xlabel('Zenith Angle, deg')
            ylabel('Gaussian Amplitude, V/m')
            xticks(np.arange(0.,91., 10.))
            #plot(x, 1.8e-4*np.exp(-x**2/2./0.15**2))
            #xlim(-1.,1.)
        
    def distances_angles(self, sep_vals, zenith_angle_rad, altitude_km):
        # get distance from ground intersection point to Xmax.
        log10_proton_energy = 17. 
        proton_Xmax = self.XMC.get_proton_Xmax(log10_proton_energy)
    
        k_x = np.sin(zenith_angle_rad)
        k_y = 0.
        k_z = np.cos(zenith_angle_rad)
        Xmax_dist = self.XMC.get_lam_Xmax(proton_Xmax, 
                                    0., 0., self.XMC.Earth_radius+altitude_km, # array altitude is 1222 m a.s.l. 
                                    k_x, k_y, k_z, 
                                    plots=False)
        #x = np.arange(-1.0, 1.0, 0.05) # km
        x = np.array(sep_vals)
        print Xmax_dist # km
        y = np.sqrt(x**2 + Xmax_dist**2 - 2*x*Xmax_dist*np.cos(pi/2-zenith_angle_rad))
        thv = np.arccos((y**2 + Xmax_dist**2 - x**2)/(2.*y*Xmax_dist))
        '''
        figure()
        subplot(211)
        plot(x,y)
        subplot(212)
        plot(x,thv*180./pi)
        '''
        return thv

    def angle_sep(self, view_angles_rad, zenith_angle_rad, altitude_km):
        # get distance from ground intersection point to Xmax.
        log10_proton_energy = 17. 
        proton_Xmax = self.XMC.get_proton_Xmax(log10_proton_energy)
    
        k_x = np.sin(zenith_angle_rad)
        k_y = 0.
        k_z = np.cos(zenith_angle_rad)
        Xmax_dist = self.XMC.get_lam_Xmax(proton_Xmax, 
                                    0., 0., self.XMC.Earth_radius+altitude_km, # array altitude is 1222 m a.s.l. 
                                    k_x, k_y, k_z, 
                                    plots=False)
        #x = np.arange(-1.0, 1.0, 0.05) # km
        #y = np.sqrt()
        x = np.sin(view_angles_rad) * Xmax_dist / (np.sin(pi/2. - zenith_angle_rad - view_angles_rad))
        y = np.sin(pi/2. + zenith_angle_rad) * Xmax_dist / (np.sin(pi/2. - zenith_angle_rad - view_angles_rad))
        '''
        figure()
        subplot(211)
        plot(x,y)
        subplot(212)
        plot(x,thv*180./pi)
        '''
        return x, y
    
    def radio_beam_model(self, theta_CR_deg, dist_core_km):
        sep = 1.1*(theta_CR_deg/70.)**4
        sig = 0.15+0.85*(theta_CR_deg/70.)**8
        norm = 5.5e-5+6.6e-5*(theta_CR_deg/45.)
        if(theta_CR_deg>45.):
            norm = 1.2e-4*(90.-theta_CR_deg)/45.
        print sep, sig, norm
        return self.double_gaussian(dist_core_km, sep, sig, norm)
    
    def radio_beam_model2(self, theta_CR_deg, view_angle, altitude_km):
        sep = 1.1*(theta_CR_deg/70.)**4
        sig = 0.15+0.85*(theta_CR_deg/70.)**8
        norm = 5.5e-5+6.6e-5*(theta_CR_deg/45.)
        if(theta_CR_deg>45.):
            norm = 1.2e-4*(90.-theta_CR_deg)/45.
        #print sep, sig, norm
        # get sep from view_angle
        dist_core_km, d_xmax = self.angle_sep(view_angle, theta_CR_deg*pi/180., altitude_km)
        return self.double_gaussian(dist_core_km, sep, sig, norm), d_xmax
        
        #figure()
        #plot(x,radio_beam_model(25., x))
        #plot(x,radio_beam_model(45., x))
        #plot(d_core_km_60, E_field_60)
        #plot(x,radio_beam_model(60., x))
        #plot(d_core_km_70, E_field_70)
        #plot(x,radio_beam_model(70., x))

        #NEED TO CONVERT THIS TO A FUNCTION OF VIEW ANGLE (NOT DISTANCE TO THE CORE)
    def double_gaussian(self, x, sep, sigma, norm):
        return norm * (np.exp(-(x-sep)**2/2./sigma**2) + np.exp(-(x+sep)**2/2./sigma**2))
    
    def get_pol(self, x_Xmax, y_Xmax, z_Xmax, x_pos, y_pos, z_pos):
        #shower axis vector
        #print '* %+1.2f %+1.2f %+1.2f %+1.2f %+1.2f %+1.2f '%(x_pos, y_pos, z_pos, x_Xmax, y_Xmax, z_Xmax)
        x_ax = (x_pos - x_Xmax)
        y_ax = (y_pos - y_Xmax) 
        z_ax = (z_pos - z_Xmax)
        
        x_pol = y_ax * self.B_z - z_ax * self.B_y
        y_pol = z_ax * self.B_x - x_ax * self.B_z
        z_pol = x_ax * self.B_y - y_ax * self.B_x

        '''
        ph = np.arctan2(y_ax, x_ax)
        th = np.arccos(z_ax/np.sqrt(x_ax**2 + y_ax**2 + z_ax**2))
        k_x = np.sin(th)*np.cos(ph)
        k_y = np.sin(th)*np.sin(ph)
        k_z = np.cos(th)
        print '* in get_pol: k_x %+1.2f k_y %+1.2f k_z %+1.2f'%(k_x, k_y, k_z)
        print '* in get_pol: B_x %+1.2e B_y %+1.2e B_z %+1.2e'%(self.B_x, self.B_y, self.B_z)        
        print '* in get_pol: x_pol %+1.2e y_pol %+1.2e z_pol %+1.2e'%(x_pol, y_pol, z_pol)
        '''
        mag = np.sqrt(x_pol**2 + y_pol**2 + z_pol**2)
        return x_pol/mag, y_pol/mag, z_pol/mag
