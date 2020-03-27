#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 13:25:21 2018

@author: romerowo
"""

#import matplotlib
#matplotlib.use('Agg')
import sys
from pylab import *

from geometry import Geometry
from xmax_calc import Xmax_calc
from detector_array import Detector_Array
from radio_emission import Radio_Emission
from detector import Detector

def acceptance_sim(radius_km = 10., 
                   num_particles = 10000, 
                   zenith_angle_deg = None, 
                   log10_energy=17.,
                   SNR_thresh = 5.,
                   N_trig = 16,
                   verbose = False):
    '''
    SETUP
    '''
    geom = Geometry(radius_km, num_particles, zenith_angle_deg)
    
    XMC = Xmax_calc()
    # I want to provide the energy, cosmic ray ground position and direction, and get the Xmax position
    XMC.get_Xmax_position(log10_energy, geom)
    Xmax_altitude = np.sqrt(XMC.x_max**2 + XMC.y_max**2 + XMC.z_max**2) - XMC.Earth_radius - XMC.detector_altitude_km
    
    det_arr = Detector_Array(mode='Ryan')
    
    rad_em = Radio_Emission(plots=False)
    
    det = Detector()
    
    # loop through simulated events:
    trigger = 0
    trig_evn = []
    thz_array = []
    flat_residual = [] # a proxy for Ryan's 1 - % power in Gaussian.
    for evn in range(0, num_particles): 
        if evn%500 == 0: 
            if verbose: 
                print '%d of %d'%(evn, num_particles)
    
        det_arr.get_distances_and_view_angles(XMC, geom, event=evn)
        # first cut by minimum view angle here.
        if np.min(det_arr.th_view*180./pi)> 10.: continue
    
        th_z = np.arccos(geom.k_z[evn])*180./pi
        E_field,dist = rad_em.radio_beam_model2(th_z, det_arr.th_view, 1.2) # 1.2 is the altitude
        E_field *= 10**(log10_energy-17.)
        x_pol, y_pol, z_pol = rad_em.get_pol(XMC.x_max[evn],  XMC.y_max[evn],  XMC.z_max[evn]-XMC.Earth_radius, geom.x_pos[evn], geom.y_pos[evn], geom.z_pos[evn])
        max_val = np.max(np.array(E_field)*1.e6)
        # then cut by maximum voltage here
        
        V_x = det.Efield_2_Voltage(np.array(E_field*x_pol), th_z)
        V_y = det.Efield_2_Voltage(np.array(E_field*y_pol), th_z)
        V_z = det.Efield_2_Voltage(np.array(E_field*z_pol), th_z)
    
        V_mag = np.sqrt(V_x**2 + V_y**2 + V_z**2)
    
        max_val_V = np.max(V_mag)*1.e6
        
        #SNR_x = np.abs(V_x + np.random.normal(0., det.V_rms, len(V_x)))/det.V_rms
        #SNR_y = np.abs(V_y + np.random.normal(0., det.V_rms, len(V_x)))/det.V_rms
        #SNR_z = np.abs(V_z + np.random.normal(0., det.V_rms, len(V_x)))/det.V_rms
        SNR_x = np.abs(V_x)/det.V_rms
        SNR_y = np.abs(V_y)/det.V_rms
        SNR_z = np.abs(V_z)/det.V_rms
    
        cut_x = SNR_x>SNR_thresh
        cut_y = SNR_y>SNR_thresh
    
        #if evn%100 == 0:
        #    print evn, '%1.1f, %1.2e, %d, %d'%(np.min(det_arr.th_view*180./pi), max_val, np.sum(cut_x), np.sum(cut_y))
        #    print '\t %1.1f %1.1f'%(np.max(SNR_x), np.max(SNR_y))
        #if( np.sum(cut_x)>=N_trig or np.sum(cut_y) >= N_trig):
        event_trigger = 0
        for FPGA in range(16):
            if( np.sum(cut_x[FPGA*16:(FPGA+1)*16]) + np.sum(cut_y[FPGA*16:(FPGA+1)*16]) >= N_trig):
            #if( np.sum(cut_x[FPGA*16:(FPGA+1)*16])>=N_trig or np.sum(cut_y[FPGA*16:(FPGA+1)*16]) >= N_trig ):
        
                event_trigger = 1
                if verbose:
                    print '\t d_core  %1.2f'%(np.sqrt(geom.x_pos[evn]**2 + geom.y_pos[evn]**2))
                    print '\t th_view %1.2f %1.2f'%(np.min(det_arr.th_view)*180./pi, np.max(det_arr.th_view)*180./pi) 
                    print '\t th_zen %1.2f'%(np.arccos(geom.k_z[evn])*180./pi)
                    print '\t=================\n'
        # Estimate Guassian power residual if the event triggers
        if event_trigger == 1:
            thz_array.append(th_z)
            trig_evn.append(evn)
            SNR_cut = 1.
            min_x = np.max([np.min(1.+SNR_x**2),SNR_cut])
            min_y = np.max([np.min((1.+SNR_y**2)),SNR_cut])
            #print 'min_x, min_y', min_x, min_y, float(len(SNR_x[SNR_x**2>SNR_cut])), float(len(SNR_y[SNR_y**2>SNR_cut]))
            res_x = 0.
            res_y = 0.
            if len(SNR_x[SNR_x**2>SNR_cut]) > N_trig:
                P_tot_x = np.sum((1.+SNR_x[SNR_x**2>SNR_cut])**2)
                P_min_x = min_x*float(len(SNR_x[SNR_x**2>SNR_cut]))
                #res_x = np.sum(((1.+SNR_x) - np.mean((1.+SNR_x)))**2)/np.sum((1.+SNR_x)**2)
                res_x = np.sum(((1.+SNR_x**2) - np.mean((1.+SNR_x**2)))**2)/np.sum((1.+SNR_x**2)**2)
                #res_x = P_min_x / P_tot_x
                #res_x = np.sum(SNR_x[SNR_x**2>SNR_cut]**2 - min_x)/min_x/float(len(SNR_x[SNR_x**2>SNR_cut]))
            if len(SNR_y[SNR_y**2>SNR_cut]) > N_trig:
                P_tot_y = np.sum((1.+SNR_y[SNR_y**2>SNR_cut])**2)
                P_min_y = min_y*float(len(SNR_y[SNR_y**2>SNR_cut]))
                #res_y = np.sum(((1.+SNR_y) - np.mean((1.+SNR_y)))**2)/np.sum((1.+SNR_y)**2)
                res_y = np.sum(((1.+SNR_y**2) - np.mean((1.+SNR_y**2)))**2)/np.sum((1.+SNR_y**2)**2)
                #res_y = P_min_y / P_tot_y
                #res_y = np.sum(SNR_y[SNR_y**2>SNR_cut]**2 - min_y)/min_y/float(len(SNR_y[SNR_y**2>SNR_cut]))
            flat_residual.append([res_x, res_y])
        trigger += event_trigger
            # will want to save events
        #print '\t %1.2e %1.2e %1.2e'%(x_pol, y_pol, z_pol)
        # then cute by maximum V_mag
        #plot([np.min(det_arr.th_view)*180./pi], [max_val], 'k.')
        #semilogy([np.min(det_arr.th_view)*180./pi], [max([max(SNR_x), max(SNR_y), max(SNR_z)])], 'k.')
        
    
    #show()
    #print 'trigger', trigger
    return float(trigger)/float(num_particles), np.array(thz_array), np.array(flat_residual)
    #print 10**(log10_energy-17.)
    #print 'det.V_rms', det.V_rms
    '''
    figure(figsize=(4,4))
    plot(geom.x_pos, geom.y_pos, '.', alpha=0.1)
    xlabel('x, km')
    ylabel('y, km')
    figure()
    hist(geom.x_pos)
    figure()
    hist(geom.th_CR*180./pi, bins=np.arange(0.,91., 5.))
    xticks(np.arange(0.,91., 10.))
    xlabel('CR Zenith Angle, deg')
    show()
    '''
        