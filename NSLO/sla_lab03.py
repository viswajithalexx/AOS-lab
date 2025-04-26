# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 22:47:50 2025

@author: HP
"""

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import cmocean


file = "C:/Users/user/Documents/24CL05012/nslo/data/seasurfaceheight_1993_2020.nc"
data = xr.open_dataset(file)

lon = data['longitude']
lat = data['latitude']
time = data['time']
sla = data['sla']

lonv =lon.values
latv = lat.values
timev =time.values


sla_clim = sla.mean(dim = 'time')
sla_climv = sla_clim.values
# Group by the time.season and compute the mean for each season
sla_mean = sla.groupby('time.season').mean(dim='time')

# Now select the seasonal means
sla_mean_djf = sla_mean.sel(season='DJF')
sla_mean_mam = sla_mean.sel(season='MAM')
sla_mean_jja = sla_mean.sel(season='JJA')
sla_mean_son = sla_mean.sel(season='SON')
season = [sla_mean_djf,sla_mean_mam,sla_mean_jja,sla_mean_son]
season_n = ['DJF','MAM','JJA','SON']

#%% Gradient
import math
import numpy as np

g = 9.8
omega = 7.2921150 * (10**-5)
dx = (lonv[1] - lonv[0]) *1000*110
dy = (latv[1] - latv[0]) *1000*110


def gradient(a, latv):
    # Calculate gradients
    sla_dy = np.gradient(a, axis=0) / dy
    sla_dx= np.gradient(a, axis=1) / dx
    
    coriolis = []
    for i in range(len(latv)):
        f = 2 * omega * math.sin(math.radians(latv[i]))
        coriolis.append(f)          
    
    cor_f = np.array(coriolis)
    cor_f[0] = cor_f[0] + 0.000001
    cor_f = cor_f.reshape(80,1)
    
    # Calculate surface_u and surface_v
    surface_v = (g / cor_f) * sla_dx
    surface_u = -(g / cor_f) * sla_dy
    
    return surface_v, surface_u
#%% climatology
sla_climv =sla_clim.values

ax=plt.figure(figsize=(20,18),dpi =100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(visible=True,draw_labels=True)
im = ax.contourf(lonv,latv,sla_clim,cmap = cmocean.cm.balance,levels =22,vmin =-0.004,vmax =0.2)
ax.add_feature(cfeature.LAND,color="black",zorder=5)
surface_v, surface_u = gradient(sla_clim, latv)
var = 2
ax.streamplot(lon[::var], lat[::var],surface_u[::var,::var],surface_v[::var,::var],color='black')
cbar = plt.colorbar(im)
cbar.set_label('SLA (m)')
plt.title('Sea level anomaly of climatological mean', fontsize=15, weight='bold')

plt.show()
#%% seasons
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20, 18), dpi=100, subplot_kw={'projection': ccrs.PlateCarree()})
fig.subplots_adjust(hspace=0.4, wspace=0.4)

for i, ax in enumerate(axes.flat):
    ax.coastlines()
    im = ax.contourf(lonv, latv, season[i], 50,cmap =cmocean.cm.balance, vmin=-0.2, vmax=0.24)
    ax.gridlines(visible=True, draw_labels=True)
    ax.add_feature(cfeature.LAND, color="black", zorder=11)
    surface_v, surface_u = gradient(season[i], latv)
    var = 2
    ax.streamplot(lon[::var], lat[::var],surface_u[::var,::var],surface_v[::var,::var],color='black')
    ax.set_title(season_n[i], fontsize=15, weight='bold')

cbar = fig.colorbar(im, ax=axes, orientation='vertical', fraction=0.05, pad=0.1)
cbar.set_label('SLA (m)')

plt.show()    
#%%
ax=plt.figure(figsize=(20,18),dpi =100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(visible=True,draw_labels=True)
im = ax.contourf(lonv,latv,sla_mean_djf,cmap = cmocean.cm.balance,levels =22,vmin = -0.2,vmax =0.24)
ax.add_feature(cfeature.LAND,color="black",zorder=5)
surface_v, surface_u = gradient(sla_mean_djf, latv)
var =3
ax.streamplot(lon[::var], lat[::var],surface_u[::var,::var],surface_v[::var,::var],color='black')
cbar = plt.colorbar(im)
cbar.set_label('SLA (m)')
plt.title('Sea level anomaly of DJF', fontsize=15, weight='bold')

plt.show()


ax=plt.figure(figsize=(20,18),dpi =100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(visible=True,draw_labels=True)
im = ax.contourf(lonv,latv,sla_mean_mam,cmap = cmocean.cm.balance,levels =22,vmin = -0.2,vmax =0.24)
ax.add_feature(cfeature.LAND,color="black",zorder=5)
surface_v, surface_u = gradient(sla_mean_mam, latv)
var = 3
ax.streamplot(lon[::var], lat[::var],surface_u[::var,::var],surface_v[::var,::var],color='black')
cbar = plt.colorbar(im)
cbar.set_label('SLA (m)')
plt.title('Sea level anomaly of MAM', fontsize=15, weight='bold') 

plt.show()


ax=plt.figure(figsize=(20,18),dpi =100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(visible=True,draw_labels=True)
im = ax.contourf(lonv,latv,sla_mean_jja,cmap = cmocean.cm.balance,levels =22,vmin = -0.2,vmax =0.24)
ax.add_feature(cfeature.LAND,color="black",zorder=5)
surface_v, surface_u = gradient(sla_mean_jja, latv)
var =3
ax.streamplot(lon[::var], lat[::var],surface_u[::var,::var],surface_v[::var,::var],color='black')
cbar = plt.colorbar(im)
cbar.set_label('SLA (m)')
plt.title('Sea level anomaly of JJA', fontsize=15, weight='bold')

plt.show()


ax=plt.figure(figsize=(20,18),dpi =100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines(visible=True,draw_labels=True)
im = ax.contourf(lonv,latv,sla_mean_son,cmap = cmocean.cm.balance,levels =22,vmin = -0.2,vmax =0.24)
ax.add_feature(cfeature.LAND,color="black",zorder=5)
surface_v, surface_u = gradient(sla_mean_son, latv)
var =3
ax.streamplot(lon[::var], lat[::var],surface_u[::var,::var],surface_v[::var,::var],color='black')
cbar = plt.colorbar(im)
cbar.set_label('SLA (m)')
plt.title('Sea level anomaly of SON', fontsize=15, weight='bold')

plt.show()
