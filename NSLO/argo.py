# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 18:03:43 2025

@author: user
"""
import xarray  as xr
import matplotlib.pyplot as plt
import numpy as np

file ="C:/Users/user/Documents/24CL05012/nslo/code/argo.nc"

data = xr.open_dataset(file)

temp = data['TEMP']
depth = data['DEPTH']

tempv =temp.values
depthv = depth.values

temperature = temp.mean(dim=('TIME'))
temperature =temperature.values
# temp_i =np.array(25)
# temp_i =np.append(temp_i,temperature)
# depth_i= np.array(800)
# depth_i = np.append(depth_i,depthv)

# depth_i= depth_i[1:]

# temp_i = temp_i[1:]
depthv =depthv*2.945

plt.figure()
plt.xlim(0,25)
plt.ylim(0,2000)
plt.plot(temperature,depthv,color = 'r',label = 'sea temperature')
plt.xlabel('Sea Temperature (Degree Celsius)')
plt.ylabel('Depth (m)')
plt.gca().invert_yaxis() 
plt.grid()
plt.legend()
plt.title('Sea Temperature')

