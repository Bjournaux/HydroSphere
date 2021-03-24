#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 18:59:06 2020

@author: daiwei
"""


"""
This is the seperated Hydrosphere script
Runs once every call
"""

import numpy as np
import seafreeze as sf
from pynverse import inversefunc as inv
from scipy.optimize import root


# Ice Ih conductivity (Andersson and Akira, 2005) 
def K_Ih(P,T):
    return -2.04104873*P+(632/T+0.38-0.00197*T)


def Calculate_Cp(P_s,T_s):
    P, T = P_s*1e-03, T_s # P in GPa; T in K
    
    '''Tchijov2003'''
    a = 72.49
    b = 3.024e-3
    c = 1.442e06
    T0 = 300
    
    '''Bezacier2014'''
    V0 = 12.49 # cm3mol-1
    K0 = 20.15 # GPa
    alpha0 = 11.58e-5 # K-1
    
    molar_mass = 18.01528 # molar mass for water/ice in grams
    
    P0 = 0.1*1e-03 # GPa
    scale = 100
    pressure = np.linspace(P0,P,scale)
    step = np.abs(P-P0)/scale
    
    '''First half of formula(2) from Tchijov2003'''
    C_p7_0 = a + b*T - c*T**-2
    '''
    Nested functions 
    func(x) is formula(1) from Bezacier2014
    
    interg(invers, step) calculates the integral part for formula(2) from Tchijov2003 
    uses Trapezoidal Rules to approximate the area
    '''
    def func(x):
        return 1.5*K0*((V0/x)**(7/3)-(V0/x)**(5/3))
    
    def interg(sol, step):
        size = len(sol)-1
        area = 0
        # Uses Trapezoidal Rules
        for i in range(size):
            area = area + (sol[i]+sol[i+1])*step/2
        return area
    
    
    '''This parts calculates the values for corresponding volumes for given pressure
       Produces a set of solutions on the x-axis(pressure)'''
    v_sol = np.zeros(scale)
    for i in range(scale):
        v_sol[i] = root(lambda x : func(x) - pressure[i], 0.1).x * np.exp(alpha0*(T-T0))
        
    '''This part calculates the density at the given pressure and tempurature'''
    V7 =  root(lambda x : func(x) - P, 0.1).x * np.exp(alpha0*(T-T0)) # in cm3mol-1
    rho = (molar_mass*1e+03/V7)[0] # in kgm-3
    
    '''integral part in formula(2) from Tchijov2003 '''
    integral = interg(v_sol, step)
    
    '''Formula(2) from Tchijov2003'''
    result = C_p7_0 + T*alpha0**2*np.exp(alpha0*(T-T0))*integral
    
    return rho, result



def general(Ps, Ts, MassWi, rb, rhocore):
    P_s = Ps   # Surface Pressure (MPa);
    T_s = Ts   # Surface Temperature (K);
    Mass_W_i = MassWi    # Mass of water in kg;
    r_b = rb              #1*6370*1e3  # Radius rocky core (m);
    rho_core = rhocore    # density of the rocky core (kg/m3);
    
    
    ##### Resolution
    res = 100;  # z grid
    Mass_it = 5 # Mass convergence iteration: 3-5 is enough (check % Mass difference)
    g_it = 2;  # Gravity convergence iteration (3 is enough)
    
    """
    ##### Plots and table
    plot_flg = 0 # plot flag: 0=no plot; 1=plot
    tab_flg = 0 # Table flag: 0=no table; 1=table
    gravity_flg = 0 # gravity profile flag: 0=no plot; 1=plot
    writefile_flg = 0 # Output file flag: 0=no output file; 1=output file
    IO_flg = 0 # Output input/output file flag: 0=no output file; 1=output file
    rawdata_flg = 1 # Output the raw data in columns; 1=output
    piechart_flg = 0 # plot a pie chart: 0=no plot; 1=plot
    scatterplt_flg = 0 # plot a scatter point between P and T: 0=no plot; 1=plot
    """
    
    #### Threshold 
    # Mass covergence loop threshhold, in percentage
    mass_thrshd = 1
    
    ###############################################################################
    
    # Ice VII approximation
    # rho_VII = 1.6169e+03
    # Cp_VII = 3400
    # alpha_VII = 2.45e-4
    '''new approimation'''
    alpha_VII = 11.58e-5
    rho_VII, Cp_VII = Calculate_Cp(P_s,T_s)
    
    
    
    # Rocky/metal core calculation
    Mass_core = 4/3*np.pi*r_b**3*rho_core
    g_s = 6.67430e-11*Mass_core/r_b**2 # Gravity at the Hydrosphere Mantle Boundary
    depth = (Mass_W_i/800/(4/3*np.pi)+r_b**3)**(1/3)-r_b # Depth approximation in m
    
    # Testing factors
    massDiff = 100
    EnableTest = 0 #0=disable, 1=enable
    
    
    # initializing the grids
    z = np.linspace(0, depth, num=res)  # depth grid
    rho = np.zeros(z.size)  # Density grid
    alpha = np.zeros(z.size)  # thermal epansivity grid
    Cp = np.zeros(z.size)  # Heat capacity grid
    dT_dz = np.zeros(z.size)  # thermal gradient grid
    phase = np.zeros(z.size)  # phase grid
    T = np.zeros(z.size)  # Temperature grid
    P = np.zeros(z.size)  # Pressure grid
    grav = np.zeros(z.size) # gravity grid
    M_L = np.zeros(z.size) # Mass grid
    
    Mass_WL = 0
    
    while (massDiff > mass_thrshd):
        """
        if EnableTest == 1:
            print("depth before " + str(depth))
        """
            
        # For mass loop the factor being iterated is /depth/
        
        # initializing the grids
        z = np.linspace(0, depth, num=res)  # depth grid
        
            
        grav[:]=g_s # Constant gravity to start with ## set all elements to g_s
        
        massDiff = np.abs(100*(Mass_W_i-Mass_WL)/Mass_W_i)
        #print(massDiff)
        
        # Gravity conversion loop
        for k in range(g_it) if (massDiff==100 or massDiff<mass_thrshd) else range(1): 
    
            # For gravity loop the factor being iterated is /grav/
    
            g = np.flip(grav,0)
            PT = np.empty((1,), np.object)
            PT[0] = (P_s, T_s)
            if P_s > 2200:
                out.rho = rho_VII
                out.alpha = alpha_VII
                out.Cp = Calculate_Cp(P_s,T_s)
                phase_s[0] = 7
            else:
                phase_s = sf.whichphase(PT)
                out = sf.seafreeze(PT,sf.phasenum2phase[phase_s[0]]) 
    
            rho_s = out.rho  # Density at the surface
            alpha_s = out.alpha  # Thermal expansivity at the surface
            Cp_s = out.Cp  # Heat capacity at the surface
            dT_dz_s = alpha_s*g[0]*T_s/Cp_s  # Thermal gradient at the surface
            T[0] = T_s
            P[0] = P_s
            rho[0] = rho_s
            alpha[0] = alpha_s
            Cp[0] = Cp_s
            dT_dz[0] = dT_dz_s
            phase[0] = phase_s[0]
            
    
            for i in range(z.size-1):  # Integration with depth
                T[i+1] = T[i] + dT_dz[i] * (z[i+1]-z[i]);
                P[i+1] = P[i] + rho[i] * g[i] * (z[i+1]-z[i])*1e-6;
                #print(rho[i])
                PT[0] = (P[i+1],T[i+1])
                '''
                print('T'+str(T))
                print('P'+str(P))
                '''
                if P[i+1] > 2200:
                    Tm=((P[i+1]*1e-3-2.17)/1.253+1)**(1/3)*354.8
                    if T[i+1] < Tm:
                        '''implemetation'''
                        out.rho = rho_VII
                        out.alpha = alpha_VII
                        out.Cp = Cp_VII
                        phase[i+1] = 7
                    else:
                        phase[i+1] = 0
                        out = sf.seafreeze(PT, 'water_IAPWS95')
                        print("water_IAPWS95 is used")
                        
                else:
                    phase[i+1] = sf.whichphase(PT)
                    out = sf.seafreeze(PT,sf.phasenum2phase[phase[i]])
                    
                rho[i+1] = out.rho;
                alpha[i+1] = out.alpha;
                Cp[i+1] = out.Cp;
                dT_dz[i+1] = alpha[i+1]*g[i+1]*T[i+1]/Cp[i+1];
    
            # Gravity in the hydrosphere
            for i in range(1,len(rho)):
                M_L[i]=rho[i]*4/3*np.pi*((r_b+z[i-1]+(depth/res))**3-(r_b+z[i-1])**3)
            Mass_Shells = np.cumsum(np.flip(M_L,0))
    
            for i in range(len(rho)):    
                grav[i] = 6.67430e-11*(Mass_core+Mass_Shells[i])/(r_b+z[i])**2  
    
        # Compute Mass
        Mass_WL = np.sum(M_L)
        Mass_diff = Mass_W_i-Mass_WL
        
        # depth difference for Mass convergence
        depth_diff = (np.abs(Mass_diff)/(np.mean(rho)*1.8)/(4/3*np.pi)+r_b**3)**(1/3)-r_b
        if   Mass_diff > 0:  
            depth = depth + depth_diff
        else:
            depth = depth - depth_diff
        
        """
        if EnableTest == 1:
            print("depth after " + str(depth))
            print()
        """
        
    # Compute Mass
    E_M_WL = Mass_WL/5.9722e24 # Mass Water Layer in Earth mass
    O_M_WL = Mass_WL/1.4e21 # Mass Water Layer in Ocean mass (Earth)
    
    # Boundary of each layer and their thickness
    bd = []
    phase_diff = phase[0]
    count = 1
    phasenum = 1
    phasediffstat = []
    phasediffstat.append(phase_diff)
    for i in range(phase.size-1):
        if phase_diff == phase[i+1]:
            count += 1
        else:
            bd.append(count)
            phase_diff = phase[i+1]
            count = 1
            phasenum += 1
            phasediffstat.append(phase_diff)
    bd.append(count)
    boundary = [[0, bd[0]-1]]
    for i in range(len(bd)-1):
        boundary.append([boundary[i][1],bd[i+1]-1+boundary[i][1]+1])
    
    sumdepth = [0]
    for i in range(len(boundary)):
        a = boundary[i][0]
        b = boundary[i][1]
        
        '''
        print('Layer '+str(i+1)+' is in phase '+str(int(phase[a])),\
                ' of depth from '+str(z[a])+'km to '+str(z[b])+'km')
        '''
        sumdepth.append(z[b])
        
    sumdepth.append(z[b]+r_b)
        
    Ra = np.zeros(len(boundary))
    Conductivity = np.zeros(len(boundary))
    # Compute Rayleigh number
    # Temperature dependent k
    # Ice Ih, II, III, V, VI, VII from top to bottom
    # k = a/T + b + c*T
    # [[a, b, c]]
    Temp_dataset = [[ 6.39162953e+02,  7.51104318e-02, -1.96063855e-04],
                    [ 4.86357760e+02, -5.05830073e-01,  7.59695073e-04],
                    [ 1.81045287e+02,  3.69540967e-01, -3.90757776e-04],
                    [ 1.60798154e+02,  8.00941009e-01, -5.96848739e-04],
                    [ 1.89382439e+02,  1.30683834e+00, -1.32469061e-03],
                    [ 6.58847367e+02,  1.00020949e+00, -7.84685730e-04]]
    # Pressure dependent k
    # Ice Ih, II, III, V, VI, VII from left to right
    # ln(k) = E + F*p - k = np.exp(E + F*p)
    # e^{E + F*p} = a/T + b + c*T
    E = [1.60, 1.25, -0.02, 0.16, 0.37, 0.65]
    F = [-0.44, 0.2, 0.2, 0.2, 0.16, 0.2]
    
    direc1 = {1:0, 2:1, 3:2, 5:3, 6:4, 7:5} # {phase:position in array}
    
    Rg = 8.314
    c1 = 1.43
    c2 = -0.03
    # The following arrays follow 'water','Ih','II','III','V','VI'
    nu0 = [0, 5e13, 1e18, 5e12, 5e14, 5e14]
    E_Jmol = [0, 60e3, np.mean([98, 55])*1e3, np.mean([103, 151])*1e3, 136e3, 110e3]
    direc2 = {0:0, 1:1, 2:2, 3:3, 5:4, 6:5, 7:4} # {phase:position in array} using ice V values for ice VII
    
    # Kappa fix
    # I, II, III, V, VI, VII
    T_intcpt = [130, 120, 240, 246, 246, 286]   # In K
    P_intcpt = [0.1, 0.24, 0.24, 0.53, 1, 2.4]  # In GPa
    k_Tfix = np.zeros(6)
    k_Pfix = np.zeros(6)
    for i in range(6):
        k_Tfix[i] = Temp_dataset[i][0]/T_intcpt[i] + Temp_dataset[i][1] + Temp_dataset[i][2]*T_intcpt[i]
        k_Pfix[i] = np.exp(E[i] + F[i]*P_intcpt[i])
        
    # print(k_Tfix)
    # print(k_Pfix)
    
    for i in range(len(boundary)):
        upper = boundary[i][0]
        lower = boundary[i][1]
        if (phase[upper] == 0):
            Ra[i] = -1
            Conductivity[i] = -1
        else:
            dir1 = direc1[phase[upper]]
            dataset = Temp_dataset[dir1]
            k = dataset[0]/T[lower]+dataset[1]+dataset[2]*T[lower]+np.exp(E[dir1]+F[dir1]*P[lower]*1e-3)
            Kappa = k/rho[lower]/Cp[lower]
            
    
            dir2 = direc2[phase[upper]]
            A = E_Jmol[dir2]/Rg/T[lower]
            B = E_Jmol[dir2]/2/Rg/c1   # B zero case?
            C = c2*(T[lower]-T[i])
            Tc = B*(np.sqrt(1+2/B*(T[lower]-C))-1)
            nu = nu0[dir2]*np.exp(A*(T[lower]/Tc-1))
    
            Ra[i] = alpha[lower]*rho[lower]*grav[lower]*(T[lower]-T[upper])*z[lower]**3/Kappa/nu
            Conductivity[i] = -k*np.abs(T[lower]-T[upper])
            
    # End of Rayleigh number computation
    
    return z,rho,alpha,Cp,dT_dz,phase,T,P,grav,M_L





