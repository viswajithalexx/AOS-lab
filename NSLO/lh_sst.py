# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 17:01:24 2025

@author: HP
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



#%% Selecting region of interest


sst_bob = sst3.sel(LAT = slice(-10,25),LONN179_180 = slice(50,110),TIME= slice('1990-01-15','1990-12-15')).mean(dim=('LAT','LONN179_180'))
lhf_bob = lh.sel(LAT = slice(-10,25),LONN179_180 = slice(50,110),TIME= slice('1990-01-15','1990-12-15')).mean(dim=('LAT','LONN179_180'))
wnd_bob = ws.sel(LAT = slice(-10,25),LONN179_180 = slice(50,110),TIME= slice('1990-01-15','1990-12-15')).mean(dim=('LAT','LONN179_180'))


sst_values = np.unique(sst_bob)
mask1 = ~np.isnan(sst_values)
sst_values = sst_values[mask1]


lhf = []
wnd = []

for i in range(len(sst_values)):
    
    sst = sst_values[i]
    
    lhfv = lhf_bob.where(sst_values ==sst).mean().values
    lhf.append(lhfv)
    wndv = wnd_bob.where(sst_values == sst).mean().values
    wnd.append(wndv)
    
sst_kelvin = sst_values+273.15    

#%%
plt.figure(figsize = (10,5),dpi =100)

plt.scatter(sst_kelvin, lhf, c='blue', label='Mean LHF')
plt.xlabel('Sea Surface Temperature (SST) [K]')
plt.ylabel('Latent Heat Flux (LHF) [W/m²]')
plt.title('SST vs Mean Latent Heat Flux')
plt.ylim(0, 200)
plt.grid(True)
plt.legend()
plt.show()    


plt.figure(figsize = (10,5),dpi =100)
plt.scatter(sst_kelvin, wnd, c='red', label='Mean WS')
plt.xlabel('Sea Surface Temperature (SST) [K]')
plt.ylabel('Wind Speed [m/s]')
plt.title('SST vs Mean Wind Speed')
plt.ylim(0,10)
plt.grid(True)
plt.legend()
plt.show()

#%% lhf variability over sst

lhf = lh.sel(LAT = slice(-30,25),LONN179_180 = slice(30,110))

# lhf_g = lhf.groupby('TIME.year')

lhf_t = np.std(lhf,axis =0)


LON = np.linspace(-30,25,55)
LAT =np.linspace(30,110,80)

plt.figure(figsize = (15,10),dpi = 100)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels =np.linspace(0,50)
im =  plt.contourf(LAT,LON,lhf_t,levels,cmap ='seismic')
ax.add_feature(cf.LAND,color = 'grey')
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)

#%%


seasons = ['DJF','MAM','JJA','SON']

for i in range(len(seasons)):
    lhf = lh.sel(TIME = time3.dt.season == seasons[i],LAT = slice(-30,25),LONN179_180 = slice(30,110))

    # lhf_g = lhf.groupby('TIME.year')

    lhf_t = np.std(lhf,axis =0)


    LON = np.linspace(-30,25,55)
    LAT =np.linspace(30,110,80)

    plt.figure(figsize = (15,10),dpi = 100)
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    levels =np.linspace(0,50)
    im =  plt.contourf(LAT,LON,lhf_t,levels,cmap ='seismic')
    ax.add_feature(cf.LAND,color = 'grey')
    ax.gridlines(visible=True,draw_labels=True)
    cbar = plt.colorbar(im)

    
