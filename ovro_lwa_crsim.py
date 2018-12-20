#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 13:25:21 2018

@author: romerowo
"""

import matplotlib
matplotlib.use('Agg')
import sys
from pylab import *

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='OVRO-LWA Acceptance simulation')

    parser.add_argument("-E",              "--log10_energy",      default=17.,  help="CR energy in log10 eV basis", type=float)
    parser.add_argument("-z",              "--zenith_angle_deg",  default=None,  help="zenith angle", type=float)
    parser.add_argument("-r",              "--radius_km",         default=10.,   help="radius of area in km", type=float)
    parser.add_argument("-nparticles",     "--num_particles",     default=10000,help="number of particles to simulate", type=int)
    parser.add_argument("-o",              "--output_tag",        default=None, help="Output tag, in quotes",type=str)
    parser.add_argument("-od",             "--outputdir",         default='./', help="Output directory, in quotes",type=str)
    parser.add_argument("-olcs_dir",      "--olcs_dir",           default='./',help="directory where simulations resided (for imports)",type=str)
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    #for i, arg in enumerate(sys.argv):
    #      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    #PARSE ARGUMENTS
    args=parser.parse_args()
    # ensure the output director ends with '/' to keep the format consisent.
    if(args.outputdir[-1]!='/'):
        args.outputdir = args.outputdir+'/'


    if(args.output_tag==None):
        if args.zenith_angle_deg == None:
            args.output_tag = 'run_E_%1.1f_r_%1.1f'%(args.log10_energy, args.radius_km)
        else:
            args.output_tag = 'run_E_%1.1f_thz_%1.1f_r_%1.1f'%(args.log10_energy, args.zenith_angle_deg, args.radius_km)
    print args.output_tag

    print ''
    print 'SIMULATION PARAMETERS' 
    print args
    print ''


    # point to the OVRO-LWA_CRSIM installation dir
    sys.path.insert(0, args.olcs_dir)

    from geometry import Geometry
    from xmax_calc import Xmax_calc
    from detector_array import Detector_Array
    from radio_emission import Radio_Emission
    from detector import Detector

    '''
    SETUP
    '''
    geom = Geometry(args.radius_km, args.num_particles, args.zenith_angle_deg)

    XMC = Xmax_calc()
    # I want to provide the energy, cosmic ray ground position and direction, and get the Xmax position
    XMC.get_Xmax_position(args.log10_energy, geom)
    Xmax_altitude = np.sqrt(XMC.x_max**2 + XMC.y_max**2 + XMC.z_max**2) - XMC.Earth_radius - XMC.detector_altitude_km

    det_arr = Detector_Array()

    rad_em = Radio_Emission(plots=False)

    det = Detector()

    # loop through simulated events:
    trigger = 0
    figure(1)
    for evn in range(0, args.num_particles): 
        if evn%500 == 0: print evn

        det_arr.get_distances_and_view_angles(XMC, geom, event=evn)
        # first cut by minimum view angle here.
        if np.min(det_arr.th_view*180./pi)> 10.: continue
    
        th_z = np.arccos(geom.k_z[evn])*180./pi
        E_field,dist = rad_em.radio_beam_model2(th_z, det_arr.th_view, 1.2) # 1.2 is the altitude
        E_field *= 10**(args.log10_energy-17.)
        x_pol, y_pol, z_pol = rad_em.get_pol(XMC.x_max[evn],  XMC.y_max[evn],  XMC.z_max[evn], geom.x_pos[evn], geom.y_pos[evn], geom.z_pos[evn])
        max_val = vmax=np.max(np.array(E_field)*1.e6)
        # then cut by maximum voltage here
        
        V_x = det.Efield_2_Voltage(np.array(E_field*x_pol), th_z)
        V_y = det.Efield_2_Voltage(np.array(E_field*y_pol), th_z)
        V_z = det.Efield_2_Voltage(np.array(E_field*z_pol), th_z)

        V_mag = np.sqrt(V_x**2 + V_y**2 + V_z**2)

        max_val_V = np.max(V_mag)*1.e6
        
        SNR_x = np.abs(V_x)/det.V_rms
        SNR_y = np.abs(V_y)/det.V_rms
        SNR_z = np.abs(V_z)/det.V_rms

        cut_x = SNR_x>5.
        cut_y = SNR_y>5.

        #if evn%100 == 0:
        #    print evn, '%1.1f, %1.2e, %d, %d'%(np.min(det_arr.th_view*180./pi), max_val, np.sum(cut_x), np.sum(cut_y))
        #    print '\t %1.1f %1.1f'%(np.max(SNR_x), np.max(SNR_y))
        if( np.sum(cut_x)>10 or np.sum(cut_y)>10):
            trigger += 1
            # will want to save events
        #print '\t %1.2e %1.2e %1.2e'%(x_pol, y_pol, z_pol)
        # then cute by maximum V_mag
        #plot([np.min(det_arr.th_view)*180./pi], [max_val], 'k.')
        #semilogy([np.min(det_arr.th_view)*180./pi], [max([max(SNR_x), max(SNR_y), max(SNR_z)])], 'k.')
        

    #show()
    print 'trigger', trigger
    print 10**(args.log10_energy-17.)
    print 'det.V_rms', det.V_rms
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
    