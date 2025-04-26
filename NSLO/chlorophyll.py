# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 12:14:12 2025

@author: user
"""

import numpy as np

chl = np.array([0.35,0.40,0.50,0.55,0.48,0.30]) 
z = np.array([5,10,20,40,70,100])
delta_x =0.5
z1 = np.arange(5,105,5)

z1_mn =[]
z1_midpoint=[]
rate_chl= []
    
for i in range(len(z1)-1):
          z1_m = (z1[i+1]+z1[i])/2
          z1_mn.append(z1_m)
          z1_midpoint.append(z1_mn[i]-5)
for i in range(len(chl)-1):
        rate = ((chl[i+1]-chl[i])/5)
        rate_chl.append(rate)          
        
for i in range(len(z1_midpoint)-1):
    (rate_chl[i])*(z1_midpoint[i])
    
        
    
    
        

    
       
        
    