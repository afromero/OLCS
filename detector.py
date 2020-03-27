#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 13:14:16 2018

@author: romerowo
"""

import numpy as np
from scipy.interpolate import Akima1DInterpolator

class Detector:
    def __init__(self):
        #print 'IMPORTING DETECTOR'
        self.kB = 1.3806488e-17 #Watt / MHz / Kelvin 
        self.c = 299.792458 # in m MHz
        self.Z0 = 120.*np.pi# Ohms
        f = np.load('filter_function_20190115.npz')
        self.filter = f['filter_array'] # made for self.fr_fine = np.arange(30., 80.1, 0.1)
        #self.filter = 1.
        #self.filter[250] = 1.
        #self.filter *= 1.e-8
        #self.filter[250] = 1.
        #self.filter[250]
        '''Antenna zenithal gain pattern'''
        zen =      np.array([0.,  15.,  30.,  45.,  60.,  75., 80.,  85.,   90. ])
        zen_gain = np.array([8.6, 8.03, 7.28, 5.84, 4.00,  0., -2.9, -8.3, -40.])
        self.zen_gain_interp = Akima1DInterpolator(zen, zen_gain-3.) # -3 is the conversion from dBic to dBi in linear polarization

        '''Antenna impeance'''
        fr =   [24., 34.,    44.,  54.,  64.,  74.,   84.,   94.]
        Z_re = [1.8,  10.,   40.,  128., 317., 271.,  158.,  129.]
        Z_im = [-68.5, 12.6, 91.5, 169., 88.2, -80.3, -89.4, -58]
        
        self.Z_re_interp = Akima1DInterpolator(fr, Z_re)
        self.Z_im_interp = Akima1DInterpolator(fr, Z_im)
        self.Z_in = 100. # Ohms, the impedance seen by the terminals of the antenna.
        
        
        '''Noise'''
        df = 0.1
        self.fr_fine = np.arange(30., 80.1, df)
        
        '''P_div is the power from the voltage divider'''
        P_div = np.abs(self.Z_in)**2/np.abs(self.Z_in + self.Z_re_interp(self.fr_fine)+1j*self.Z_im_interp(self.fr_fine))**2
        self.Noise = (4.*(self.kB*1.e-6)*self.galactic_noise_temp(self.fr_fine)) * self.Z_re_interp(self.fr_fine) 
        self.Noise *= P_div
        self.Noise *= self.filter
        self.Noise += self.kB*1.e-6*250.*np.real(self.Z_in) # 250. Kelvin internal noise
        

        self.f_c = 55. # MHz
        lam_c = self.c/self.f_c
        
        self.h_eff = 4. * self.Z_re_interp(self.fr_fine) / self.Z0 * lam_c**2 / 4. / np.pi 
        self.h_eff *= np.abs(self.Z_in)**2 / np.abs(self.Z_re_interp(self.fr_fine)+1j*self.Z_im_interp(self.fr_fine)+self.Z_in)**2
        self.h_eff *= self.filter
        self.h_eff = np.sqrt(self.h_eff)
        
        self.h_0 = np.mean(self.h_eff) # assume flat spectrum for CR pulse
        #self.A_0 = 4. * self.Z_re_interp(self.f_c) / self.Z0 * lam_c**2 / 4. / np.pi 
        #self.A_0 *= np.abs(self.Z_in)**2 / np.abs(self.Z_re_interp(self.f_c)+1j*self.Z_im_interp(self.f_c)+self.Z_in)**2

        
        self.V_rms = np.sqrt(np.sum(self.Noise * df*1.e6))
        
    ####################################################################################

    def Efield_2_Voltage(self, E_field, theta_zenith_deg): # input in V/m
        # Get directivity
        D = self.zen_gain_interp(theta_zenith_deg)
        return E_field * self.h_0 * np.sqrt( 10**(D/10.))
        
    ''' Frequencies are in MHz'''
    def galactic_noise(self, freq):
        Ig  = 2.48e-20
        Ieg = 1.06e-20
        tau = 5.0*freq**-2.1    
        return Ig*freq**-0.52*(1.-np.exp(-tau))/tau + Ieg*freq**-0.8*np.exp(-tau) # in Watts / M^2 / Hz / sr
    
    def galactic_noise_temp(self, freq):
        return 1./2./self.kB * self.c**2 / freq**2 * self.galactic_noise(freq)*1.e6
    
'''


print kB*1.e-6*250.*Z_in
figure()
plot(fr_fine, galactic_noise_temp(fr_fine))
grid(True)

figure()
plot(fr_fine, 10*np.log10(Noise))
plot(fr_fine, 10*np.log10(kB*1.e-6*250.*Z_in*np.ones(len(fr_fine))), 'r--')
#plot(fr_fine, Noise)
grid(True, which='both')
ylabel('dB(V$^2$/Hz)')

df = np.diff(fr_fine)[0]
print np.sum(Noise)*df, np.sqrt(np.sum(Noise)*df), np.sqrt(np.sum(Noise)*df)*1.e6
'''