# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 22:50:26 2025

@author: HP
"""

import numpy as np

chl = np.array([0.35,0.40,0.50,0.55,0.48,0.30]) 
z = np.array([5,10,20,40,70,100])

def chlorophyll_values(chl,z):
    
    
    chl_midpoint =[]
    depth_z =[]
    chlorophyll_sum =[]
    midpoint_nz =[]
    delta_z_new=[]
    
    for i in range(len(z)-1):
        midpoint_z = ((z[i+1]+z[i])/2)
        midpoint_nz.append(midpoint_z)
    
    grid=np.array(0)
    grid=np.append(grid,midpoint_nz)
    
    for i in range(len(grid)-1):
       delta_z_ne = ((grid[i+1]-grid[i]))
       delta_z_new.append(delta_z_ne)
    
    delta_z_new.append(abs(midpoint_nz[-1]-z[-1]))
    
    for i in range(len(delta_z_new)):
        chlorophyll = (chl[i]*delta_z_new[i])
        chlorophyll_sum.append(chlorophyll) 
    
    chlorophyll_value_new = np.sum(chlorophyll_sum)
    return chlorophyll_value_new 

chl_v = chlorophyll_values(chl,z)
print('value of chlorophyll Concentration=',chl_v)   
       
