#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 16:48:29 2025

@author: bobco-08
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import cmaps


file = '/home/bobco-08/24cl05012/CO2/data/oc_v2025.flux.nc'

data = xr.open_dataset(file)


#CO2 data variables

lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2_ocean = data['co2flux_ocean']
area = data['dxyp']


co2_flux = co2_ocean/area

esio_flux = co2_ocean.sel(lat = slice(-5,10),lon = slice(92,110))*(10**15)


#%%

esio_on = esio_flux.sel(mtime = esio_flux.mtime.dt.month.isin([10,11]))
esio_yr = esio_on.resample(mtime = 'Y').mean('mtime')
esio_flux  = esio_yr.mean(dim ='mtime')
#%%
ax=plt.figure(1, figsize=(10,8),dpi = 500)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im=ax.contourf(esio_flux.lon, esio_flux.lat, esio_flux,levels =20,cmap= cmaps.cmp_b2r ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)

plt.title('DJF climatology of  pCO$_2$ Flux (1957–2023)', pad = 20,fontsize = 16)
plt.grid(True)

#%%
nwio_flux = co2_ocean.sel(lat = slice(-0,24),lon = slice(30,70))*(10**15)
nwio_jjas = nwio_flux.sel(mtime = nwio_flux.mtime.dt.month.isin([6,7,8,9]))
nwio_yr = nwio_jjas.resample(mtime = 'Y').mean('mtime')
nwio_flux  = nwio_yr.mean(dim ='mtime')

#%%
ax=plt.figure(1, figsize=(15,12),dpi = 300)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im=ax.contourf(nwio_flux.lon, nwio_flux.lat, nwio_flux,levels =20,cmap= cmaps.cmp_b2r ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)

plt.title('pCO$_2$ Flux in the Northwestern Indian Ocean (0°-24°N, 30°-70°E) (1957–2023)', pad = 20,fontsize = 16)
plt.grid(True)
#%%
eio_flux = co2_ocean.sel(lat = slice(-5,5),lon = slice(30,96))*(10**15)
eio_flux  = eio_flux.mean(dim ='mtime')

#%%
ax=plt.figure(1, figsize=(15,12),dpi = 300)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im=ax.contourf(eio_flux.lon, eio_flux.lat, eio_flux,levels =20,cmap= cmaps.cmp_b2r ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)

plt.title('pCO$_2$ Flux in the Equatorial Indian Ocean (5°S-5°N, 30°E-96°E)(1957–2023)', pad = 20,fontsize = 16)
plt.grid(True)
#%%
io_flux = co2_ocean.sel(lat = slice(-30,30),lon = slice(40,110))*(10**15)
io_flux  = io_flux.mean(dim ='mtime')
ax=plt.figure(1, figsize=(15,12),dpi = 300)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(zorder = 15)
im=ax.contourf(io_flux.lon, io_flux.lat, io_flux,levels =20,cmap= cmaps.cmp_b2r ,transform=ccrs.PlateCarree(), extend ='both')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im,label= r'g C m$^{-2}$ yr$^{-1}$',pad = 0.13)

plt.title('Indian Ocean (1957–2023)', pad = 20,fontsize = 16)
plt.grid(True)