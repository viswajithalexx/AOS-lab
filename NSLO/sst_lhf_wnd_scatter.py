# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:28:40 2025

@author: user
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file1 = "C:/Users/user/Documents/24CL05012/nslo/data/OAFLUX_1deg_LHF_1990_2008.nc"
file2 = 'C:/Users/user/Documents/24CL05012/nslo/data/OAFLUX_1deg_SST_1990_2008.nc'
file3 =  "C:/Users/user/Documents/24CL05012/nslo/data/OAFLUX_1deg_WSP10m_1990_2008.nc" 

data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)
data3 = xr.open_dataset(file3)

lon3 = data1['LONN179_180']
lat3 = data1['LAT']
time3 = data1['TIME']
sst3  = data2['TMPSF']
lh = data1['LHTFL']
ws = data3['WND10']

#%%


sst_bob = sst3.sel(LAT = slice(7,30),LONN179_180 = slice(75,100))
lhf_bob = lh.sel(LAT = slice(7,30),LONN179_180 = slice(75,100))
wnd_bob = ws.sel(LAT = slice(7,30),LONN179_180 = slice(75,100))



sst_bob  = sst_bob.values
lhf_bob  = lhf_bob.values
wnd_bob  = wnd_bob.values


sst = sst_bob.flatten()
lhf = lhf_bob.flatten()
wnd = wnd_bob.flatten()




#%%

plt.figure()
plt.scatter(sst,wnd,s=4,alpha=0.5)
plt.xlabel('SST') 
plt.ylabel('WINDSPEED')  


plt.figure()
plt.scatter(sst,lhf, color = 'red',s=4,alpha=0.5)
plt.xlabel('SST') 
plt.ylabel('LHF')  

       
    
plt.figure()
plt.scatter(wnd,lhf,color = 'green',s=4,alpha=0.5)
plt.xlabel('WINDSPEED') 
plt.ylabel('LHF')  





#%%



