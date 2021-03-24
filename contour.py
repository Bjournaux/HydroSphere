#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:04:44 2020

@author: daiwei
"""


########################################

# This part generates the contour plot #

########################################

    
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import os
import matplotlib
matplotlib.use('TkAgg')

#%%
'''Function part'''
def grid(step):
    return np.meshgrid(np.zeros(steps), np.zeros(steps))


def find_Thickness(data, phase):
    data = data[6:]
    start = 0
    end = 0
    for i in range(len(data)):
        currentphase = int(float(data[i].split()[5]))
        if currentphase == phase:
            if i == 0 or int(float(data[i-1].split()[5])) != currentphase:
                start = float(data[i].split()[0])
            end = float(data[i].split()[0])
    return np.abs(end - start)


def find_largest(tar):
    largest = 0
    for i in range(len(tar)):
        for j in range(len(tar[i])):
            if tar[i][j] > largest:
                largest = tar[i][j]
    return largest
            
#%%
'''constant'''
rho_core = 5514 # Constant density in kg/m3
nametag = ['Hydrosphere', 'High Pressure Ice', 'Liquid Water', 'Ice I', 'Ice II', 
           'Ice III', 'Ice IV', 'Ice V', 'IceVI', 'IceVII']
nametag2 = ['Bottom Temperature', 'Bottom Pressure']
os_path = "~/Desktop/document/Python/Test folder/dataset/sixth(1125)/"

#%%
'''control panel'''
steps = 15 # resolution
r_b_number = 5 # cases of radius to be chosen: from 1 to 5
limitset = np.array((0.11, 10)) # situable range to show the contour plots

target_phase = 3 # starts from 1 NOT 0; refer to nametag above
name = nametag[target_phase-1] # used for generating titles on the contour plots

special_plot = False # True=plot bottom temperature or bottom pressure
target_phase_s = 1 # 0=temperature 1=pressure

#%%
'''initializing grids'''
ttg, ttg = grid(steps) # hydrosphere
wg, wg = grid(steps) # liquid water
ICE, ICE = grid(steps) # HP ice
Ig, Ig = grid(steps) # phase I
IIg, IIg = grid(steps) # phase II
IIIg, IIIg = grid(steps) # phase III
IVg, IVg = grid(steps) # phase IV
Vg, Vg = grid(steps) # phase V
VIg, VIg = grid(steps) # phase VI
VIIg, VIIg = grid(steps) # phase VII

Tg, Tg = grid(steps)
Pg, Pg = grid(steps)

grid_set = [ttg, ICE, wg, Ig, IIg, IIIg, IVg, Vg, VIg, VIIg]
grid_set2 = [Tg, Pg]

'''initializing some contants through calculation'''
r_b = np.linspace(0.6*6.36e+06, 1.8*6.36e+06, num=5)
r_b_earth = r_b/6.36e+06 # represent in earth radius

'''x-axis: scaling temperature'''
x = np.array(np.linspace(263, 373, steps))
# y needs to change accoding to radius

'''calculation of mass of the core and water based on chosen rocky core radius'''
r_b_current = r_b[r_b_number-1] #1*6370*1e3  # Radius rocky core (m);
Mass_core = 4/3*np.pi*r_b_current**3*rho_core
water_mass = np.linspace(Mass_core/999, Mass_core/9, num=steps)
y = np.array(100*water_mass/(water_mass+Mass_core))
xg, yg = np.meshgrid(x,y)

#%%
'''filling the grids from data files'''
for i in range(steps): # temp
    for j in range(steps): # mass of water
        filename = "T"+str(i+1)+"M"+str(j+1)+"Rb"+str(r_b_number)+".txt"
        path = os_path + filename
        with open(os.path.expanduser(path)) as f:
            lines = f.readlines()
        # thickness of hydrosphere
        ttg[j][i] = float(lines[-1].split()[0])
        # bottom temperature
        Tg[j][i] = float(lines[-1].split()[-4])
        # bottom pressure
        Pg[j][i] = float(lines[-1].split()[-3])
        # remaining grids
        for k in range(len(grid_set)-2):
            grid_set[k+2][j][i] = find_Thickness(lines, k)
        
        '''
        # liquid water
        wg[j][i] = find_Thickness(lines, 0)
        # phase I ice
        Ig[j][i] = find_Thickness(lines, 1)
        # phase II ice
        IIg[j][i] = find_Thickness(lines, 2)
        # phase III ice
        IIIg[j][i] = find_Thickness(lines, 3)
        # phase IV ice
        IVg[j][i] = find_Thickness(lines, 4)
        # phase V ice
        Vg[j][i] = find_Thickness(lines, 5)
        # phase VI ice
        VIg[j][i] = find_Thickness(lines, 6)
        # phase VII ice
        VIIg[j][i] = find_Thickness(lines, 7)
        '''
        
'''Avoid zero cases'''
for grids in grid_set:
    if grids.any() != 0:
        grids[grids==0] = 1

'''
if wg.any() != 0:
    wg[wg==0] = 1
    #print('water')
    #print(wg)
if Ig.any() != 0:
    #print('I')
    Ig[Ig==0] = 1
    #print(Ig)
if IIg.any() != 0:
    #print('II')
    IIg[IIg==0] = 1
    #print(IIg)
if IIIg.any() != 0:
    #print('III')
    IIIg[IIIg==0] = 1
    #print(IIIg)
if IVg.any() != 0:
    #print('IV')
    IVg[IVg==0] = 1
    #print(IVg)
if Vg.any() != 0:
    #print('V')
    Vg[Vg==0] = 1
    #print(Vg)
if VIg.any() != 0:
    #print('VI')
    VIg[VIg==0] = 1
    #print(VIg)
if VIIg.any() != 0:
    #print('VII')
    VIIg[VIIg==0] = 1
    #print(VIIg)
'''
    
#print(str(y)+"\n")

# Total thickness of high pressure ices
grid_set[1] = Vg + VIg + VIIg
# suppress to 1
grid_set[1][grid_set[1]==2] = 1


if special_plot:
    target = grid_set2[target_phase_s]
    name = nametag2[target_phase_s]
else:
    target = grid_set[target_phase-1]
    

'''
with open("depth_data.txt", mode = "w") as f:
    f.write("x-axis\n"+str(x)+" \n")
    f.write("y-axis\n"+str(y)+" \n")
    f.write("meshgrid x-axis\n"+str(xg)+" \n")
    f.write("meshgrid y-axis\n"+str(yg)+" \n")
    f.write("depth grid\n"+str(target)+" \n")
'''


#%%
'''plotting part'''
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
# gnuplot
# viridis
# ocean
# Blues_r
colormap = ax.contourf(xg, yg, target/1000,levels=100, #target/1000
                       cmap=plt.cm.viridis, alpha=0.7
                       )

cbar=fig.colorbar(colormap,
             orientation='vertical',  # horizontal colour bar
             shrink=0.85,
             )
scale = np.linspace(0,find_largest(target),10)/1000#'''/1000'''
scale_string = str(scale.astype(int))[1:-1].split()
cbar.ax.set_yticklabels(scale_string)

'''plot color line if it is HP ice plot'''
if target_phase==2:
    colorline = ax.contour(xg,yg,target/1000,levels=[3.9999,4],linewidths=1,linestyles='solid',
                       cmap=plt.cm.viridis)
    

'''
cl = fig.colorbar(colorline,orientation='horizontal',shrink=0.85)
scale = np.linspace(0,find_largest(target),10)/1000
scale_string = str(scale.astype(int))[1:-1].split()
cl.ax.set_yticklabels(scale_string)
'''

'''giving level line labels by hand'''
#ax.clabel(colorline, inline=1,fontsize=12,fmt='%1.0f',inline_spacing=1,manual=True)


ax.set_xlim((263, 373))
ax.set_ylim(limitset)
txt = "rocky core radius = "+str(r_b[r_b_number-1])+"km"

if special_plot:
    title_name = name + ' '
else:
    title_name = name + ' Thickness '


unit = ['(K)', '(MPa)']
label = title_name+unit[target_phase_s] if special_plot else title_name+'(km)'
plt.title(label=label,fontsize=20)

plt.xlabel('Surface Temperature (K)',fontsize=17)
plt.ylabel('Mass of Water in wt% H2O',fontsize=17)

txt = "rocky core radius = "+str(r_b_earth[r_b_number-1])+" Earth rocky core radius\nsurface pressure = 1bar"
fig.text(.5, .01, txt, ha='center')

#plt.savefig('Hydrosphere Thickness Case '+str(r_b_number)+'.png')
plt.savefig(title_name + 'Case '+str(r_b_number)+'.png', format='png', dpi=600)






























