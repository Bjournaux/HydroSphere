#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:55:18 2020

@author: daiwei
"""

"""
This file writes the raw data in forms of arrays to different txt files 
"""


#%%
#######################################

# This part generates the output file #

#######################################


import numpy as np
import seafreeze as sf
import Hydrosph as Hyd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy.optimize import root
from matplotlib import cm
import random
import os

##### Inputs

current_version = "Version 0p8.4\n"

# Data set: Ganymede; Europa; Titan; Earth; 50x ocean Earth; 1.8xR 100x ocean Earth
P_s_set = [0.7e-12, 0.1e-12, 0.1467, 0.1, 0.1, 0.1]
#T_s_set = [130, 130, 90, 300, 300, 300]
#Mass_W_i_set = [1.5e+22, 3e+18, 7.2e+22, 1.4e+21, 7e+22, 1.4e+23]
#r_b_set = [1.63e+06, 1.46e+06, 2.58e+06, 6.36e+06, 6.36e+06, 1.14e+07]
rho_core_set = [5500, 3040, 1880, 5514, 5514, 5514]

# 0:Ganymede; 1:Europa; 2:Titan; 3:Earth; 4:50x ocean Earth; 5:1.8xR 100x ocean Earth
target_planet = 2
'''
#Ps = 0.1
#Ts = 263 - 373      5steps
#MassWi = 1.4e+21 - 100*1.4e+21   5steps
#rb = 0.6*6.36e+06 - 1.8*6.36e+06 5steps
#rhocore = 5514
'''
T_s_set = np.linspace(263, 373, num=15)
#Mass_W_i_set = np.linspace(1.4e+21, 100*1.4e+21, num=15)
r_b_set = np.linspace(0.6*6.36e+06, 1.8*6.36e+06, num=5)

rho_core = 5514    # density of the rocky core (kg/m3);

rawdata_flg = 1 # flag to write the txt file 1=write, 0=perform nothing

'''all the cases'''

'''
path = "seventh/"
os.makedirs(os.path.dirname(path))

for count1 in range(15):    #Temperature
    for count3 in range(5): #rocky core radius
        r_b = r_b_set[count3] #1*6370*1e3  # Radius rocky core (m);
        Mass_core = 4/3*np.pi*r_b**3*rho_core
        Mass_W_i_set = np.linspace(Mass_core/999, Mass_core/9, num=15)
        
        for count2 in range(15):    #Mass of Water
            P_s = 0.1   # Surface Pressure (MPa);
            T_s = T_s_set[count1]   # Surface Temperature (K);                        
            Mass_W_i = Mass_W_i_set[count2] # Mass of water in kg;
            
            z,rho,alpha,Cp,dT_dz,phase,T,P,grav,M_L = Hyd.general(P_s,T_s,Mass_W_i,r_b,rho_core)
            
            print(str(count1)+str(count2)+str(count3)+"\n")
            # Output the raw data to seperated files
            if rawdata_flg == 1:
                filename = "T"+str(count1+1)+"M"+str(count2+1)+"Rb"+str(count3+1)+".txt"
                with open(path+filename, mode = "w") as f:
                    f.write(current_version+"INPUT:\n")
                    f.write("Temperature(K) Surface_Pressure(MPa) Water_Mass(kg) Rocky_Core_Radius(m) Rocky_Core_Density(kg/m3)\n")
                    f.write(str(T_s)+" "+str(P_s)+" "+str(Mass_W_i)+" "+str(r_b)+" "+str(rho_core)+" \n")
                    f.write("OUTPUT:\n")
                    f.write("depth(m) density(kg/m3) thermal_expansivity(1/K) heat_capacity(J/K) thermal_gradient phase temperature(K) pressure(MPa) gravity(m/s2) mass(kg)\n")
                    for i in range(rho.size):
                        f.write(str(z[i])+" "+str(rho[i])+" "+str(alpha[i])+" "+str(Cp[i])+" "+str(dT_dz[i])+" ")
                        f.write(str(phase[i])+" "+str(T[i])+" "+str(P[i])+" "+str(grav[i])+" "+str(M_L[i])+"\n")

'''
    
'''
# This part is discarded

for count1 in range(15):
    for count3 in range(10):
        P_s = 0.1   # Surface Pressure (MPa);
        T_s = T_s_set[count1]   # Surface Temperature (K);
        #Mass_W_i = Mass_W_i_set[count2]    # Mass of water in kg;
        r_b = r_b_set[count3]              #1*6370*1e3  # Radius rocky core (m);
        rho_core = 5514    # density of the rocky core (kg/m3);
        
        # new way of calculating the water mass (kg)
        Mass_core = 4/3*np.pi*r_b**3*rho_core
        Mass_W_i = Mass_core/9
        
        z,rho,alpha,Cp,dT_dz,phase,T,P,grav,M_L = Hyd.general(P_s,T_s,Mass_W_i,r_b,rho_core)
        
        print(count3)
        # Output the raw data to seperated files 
        if rawdata_flg == 1:
            filename = "T"+str(count1+1)+"Rb"+str(count3+1)+".txt"
            with open(path+filename, mode = "w") as f:
                f.write("Version 0p8.1\n")
                f.write("Temperature: "+str(T_s)+"K; Water Mass:"+str(Mass_W_i)+"kg; Rocky Core Radius: "+str(r_b)+"km\n")
                f.write("depth density thermal_expansivity heat_capacity thermal_gradient ")
                f.write("phase temperature pressure gravity mass\n")
                for i in range(rho.size):
                    f.write(str(z[i])+" "+str(rho[i])+" "+str(alpha[i])+" "+str(Cp[i])+" "+str(dT_dz[i])+" ")
                    f.write(str(phase[i])+" "+str(T[i])+" "+str(P[i])+" "+str(grav[i])+" "+str(M_L[i])+"\n")
'''


#%%


'''running only once ver.'''
rawdata_flg = 0
plot_flag = 1
test_plot = 0

c1 = 11 #temp: 0-14
c2 = 7 #water mass: 0-14  8
c3 = 3 # rb radius: 0-4

T_s_set = np.linspace(263, 373, num=15)
r_b_set = np.linspace(0.6*6.36e+06, 1.8*6.36e+06, num=5)
rho_core = 5514    # density of the rocky core (kg/m3);

r_b = r_b_set[c3] # 1*6370*1e3 # Radius rocky core (m);
Mass_core = 4/3*np.pi*r_b**3*rho_core
# Mass_W_i_set = np.linspace(Mass_core/999, Mass_core/4, num=15)
Mass_W_i_set = np.linspace(Mass_core/999, Mass_core/9, num=15)


P_s = 0.1   # Surface Pressure (MPa);
T_s = T_s_set[c1]   # Surface Temperature (K);
#Mass_W_i = Mass_W_i_set[count2]    # Mass of water in kg;


Mass_W_i = Mass_W_i_set[c2]
# print(str(T_s)+" "+str(Mass_W_i))

z,rho,alpha,Cp,dT_dz,phase,T,P,grav,M_L = Hyd.general(P_s,T_s,Mass_W_i,r_b,rho_core)

# def reform(Cp):
#     for i in range(len(Cp)-2):
#         if Cp[i+1] < 80:
#             Cp[i+1] = (Cp[i]+Cp[i+2])/2
#     return Cp

# Cp = reform(Cp)

# Output the raw data to seperated files
if rawdata_flg == 1:
    filename = "T"+str(c1+1)+"M"+str(c2+1)+"Rb"+str(c3+1)+".txt"
    with open(filename, mode = "w") as f:
        f.write(current_version+"INPUT:\n")
        f.write("Temperature(K) Surface_Pressure(MPa) Water_Mass(kg) Rocky_Core_Radius(m) Rocky_Core_Density(kg/m3)\n")
        f.write(str(T_s)+" "+str(P_s)+" "+str(Mass_W_i)+" "+str(r_b)+" "+str(rho_core)+" \n")
        f.write("OUTPUT:\n")
        f.write("depth(m) density(kg/m3) thermal_expansivity(1/K) heat_capacity(J/K) thermal_gradient phase temperature(K) pressure(MPa) gravity(m/s2) mass(kg)\n")
        for i in range(rho.size):
            f.write(str(z[i])+" "+str(rho[i])+" "+str(alpha[i])+" "+str(Cp[i])+" "+str(dT_dz[i])+" ")
            f.write(str(phase[i])+" "+str(T[i])+" "+str(P[i])+" "+str(grav[i])+" "+str(M_L[i])+"\n")

#%%

if plot_flag == 1:

    font = {'family' : 'normal',
                'size'   : 15}
    
    plt.rc('font', **font)
    plt.rc('xtick', labelsize=13) 
    plt.rc('ytick', labelsize=13)
    
    
    fig = plt.figure(figsize=(12,12))
    plt.subplots_adjust(wspace=0.2 , hspace=0.2)
    plt.subplot(221) 
    plt.xlabel('Depth (km)')
    plt.ylabel('Pressure (GPa)')
    plt.plot(z*1e-3, P*1e-3, '-')
    
    plt.subplot(222)
    plt.xlabel('Depth (km)')
    plt.ylabel('density')
    plt.plot(z*1e-3, rho, '-')
    
    plt.subplot(223)  
    plt.xlabel('Depth (km)')
    plt.ylabel('Cp')
    plt.plot(z*1e-3, Cp, '-')
    
    plt.subplot(224)  
    plt.xlabel('Depth (km)')
    plt.ylabel('Phase')
    plt.plot(z*1e-3, phase, '*')
    
    
    plt.savefig('T'+str(c1)+'WM'+str(c2)+'Rb'+str(c3)+'.png', format='png', dpi=600)

#%%

if test_plot == 1:
    
    step = 20
    alpha0 = 11.58e-5
    
    T = np.linspace(250, 600, num = step) # K
    P = np.linspace(0, 50000, num = step) # MPa
    
    Cp_grid = np.zeros((step, step)) # x-axis:T y-axis:P
    
    for i in range(step):
        for j in range(step):
            rho, Cp = Hyd.Calculate_rho_Cp(P[i], T[j], alpha0)
            Cp_grid[step-i-1][j] = Cp
    
    P = P*1e-03
    TT, PP = np.meshgrid(T, P)
    
    fig = plt.figure()
    ax = Axes3D(fig)
    fig.add_axes(ax)
    ax.plot_surface(TT, PP, Cp_grid, rstride=1, cstride=1, cmap=cm.viridis)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Pressure (GPa)")
    
    ax.view_init(30,120)
    
    plt.show()
        

        
    



#%%

# def Calculate_rho_Cp(P_s,T_s):
#     P, T = P_s*1e-03, T_s # P in GPa; T in K
    
#     '''Tchijov2003'''
#     a = 72.49
#     b = 3.024e-3
#     c = 1.442e06
#     T0 = 300
    
#     '''Bezacier2014'''
#     V0 = 12.49 # cm3mol-1
#     K0 = 20.15 # GPa
#     alpha0 = 11.58e-5 # K-1
    
#     molar_mass = 18.01528 # molar mass for water/ice in grams
    
#     P0 = 0.1*1e-03 # GPa
#     scale = 100
#     pressure = np.linspace(P0,P,scale)
#     step = np.abs(P-P0)/scale
    
#     '''First half of formula(2) from Tchijov2003'''
#     C_p7_0 = a + b*T - c*T**-2
#     '''
#     Nested functions 
#     func(x) is formula(1) from Bezacier2014
    
#     interg(invers, step) calculates the integral part for formula(2) from Tchijov2003 
#     uses Trapezoidal Rules to approximate the area
#     '''
#     def func(x):
#         return 1.5*K0*((V0/x)**(7/3)-(V0/x)**(5/3))
    
#     def interg(sol, step):
#         size = len(sol)-1
#         area = 0
#         # Uses Trapezoidal Rules
#         for i in range(size):
#             area = area + (sol[i]+sol[i+1])*step/2
#         return area
    
    
#     '''This parts calculates the values for corresponding volumes for given pressure
#         Produces a set of solutions on the x-axis(pressure)'''
#     v_sol = np.zeros(scale)
#     for i in range(scale):
#         v_sol[i] = root(lambda x : func(x) - pressure[i], 0.1).x * np.exp(alpha0*(T-T0))
        
#     '''This part calculates the density at the given pressure and tempurature'''
#     V7 =  root(lambda x : func(x) - P, 0.1).x * np.exp(alpha0*(T-T0)) # in cm3mol-1
#     rho = (molar_mass*1e+03/V7)[0] # in kgm-3
    
#     '''integral part in formula(2) from Tchijov2003 '''
#     integral = interg(v_sol, step)
    
#     '''Formula(2) from Tchijov2003'''
#     result = C_p7_0 + T*alpha0**2*np.exp(alpha0*(T-T0))*integral
    
#     #print("cp70 is " + str(C_p7_0))
    
#     return rho, result


# # for i in range(len(Cp)):
# #     if Cp[i] < 100:
# #         print(str(i)+" "+str(Cp[i]))

# for i in [13,16,35]:
#     print(str(P[i])+" "+str(T[i]))
#     rho, cp = Calculate_rho_Cp(P[i+1],T[i+1])
#     print(cp)
# cp = np.zeros(len(P))
# for i in range(len(P)):
#     rho, c = Calculate_rho_Cp(P[i],T[i])
#     cp[i] = c

# print(cp)

#%%

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
