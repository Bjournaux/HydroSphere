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
"""

'''running only once ver.'''
c1 = 7 #temp
c2 = 8 #water mass
c3 = 4 # rb radius 

T_s_set = np.linspace(263, 373, num=15)
r_b_set = np.linspace(0.6*6.36e+06, 1.8*6.36e+06, num=5)
rho_core = 5514    # density of the rocky core (kg/m3);


r_b = r_b_set[c3] #1*6370*1e3  # Radius rocky core (m);
Mass_core = 4/3*np.pi*r_b**3*rho_core
Mass_W_i_set = np.linspace(Mass_core/999, Mass_core/4, num=15)
        
        
P_s = 0.1   # Surface Pressure (MPa);
T_s = T_s_set[c1]   # Surface Temperature (K);
#Mass_W_i = Mass_W_i_set[count2]    # Mass of water in kg;
            
Mass_W_i = Mass_W_i_set[c2]
#print(str(T_s)+" "+str(Mass_W_i))

z,rho,alpha,Cp,dT_dz,phase,T,P,grav,M_L = Hyd.general(P_s,T_s,Mass_W_i,r_b,rho_core)

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


    
"""
    
    
    
    
