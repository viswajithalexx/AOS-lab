# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 12:49:49 2025

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt

S = 1370

sigma = 5.67*10**(-8)
apsilon = 0.7
alpha = 0.3
temperatures =[]
for i in range(-90,91):
    temp = (np.cos(np.radians(i))*S*(1-alpha))/4
    T  = temp/(apsilon*sigma)
    temperatures.append(T**(1/4))
    print('at',i,'=',T**(1/4))
    
plt.figure()
X = temperatures
Y = range(-90,91)
plt.scatter(X,Y,c= temperatures,cmap= 'jet',s=100)
#plt.xscale('log')
# plt.yscale('log')
cbar = plt.colorbar()
plt.xlabel('Temperatures (Kelvin)')
plt.ylabel('Latitudes');
plt.title('Plot of Latitude vs Temperature')

cbar.set_label('tempertaure')
plt.grid(True)
plt.show()

