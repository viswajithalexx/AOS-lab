# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 09:50:20 2025

@author: user
"""

import xarray as xr
import numpy as np

file1 ="C:/Users/user/Documents/24CL05012/nsla/WIND_DATA/U_V_Wind_June_1979_2021_ERA5_Hr_BOB.nc"
data1 = xr.open_dataset(file1)


lon = data1['longitude']
lat = data1['latitude']
time = data1['time']
u = data1[ 'u10']
v = data1[ 'v10']

latv = lat.values
lonv =lon.values

u_10 = u.sel(time = slice('2021-6-1','2021-6-10'))
v_10 = v.sel(time = slice('2021-6-1','2021-6-10')) 

u_10v =u_10.values
v_10v =v_10.values

dx = (lonv[1]-lonv[0])*110*1000

dy = (latv[0]-latv[1])*110*1000

for n in range(0,len(u_10v)-1):
    du_divergence = (u_10v[n+1,:,:]-u_10v[n,:,:])/dx
    # print(du)

for n in range(0,len(u_10v)-1):
    dv_divergence = (v_10v[n+1,:,:]-v_10v[n,:,:])/dy
    #print(dv)
    
diverg = (du_divergence+dv_divergence)*(10**4)


#vortcity
for s in range(0,len(u_10v)-1):
    du_vort = (u_10v[s+1,:,:]-u_10v[n,:,:])/dy
    
for s in range(0,len(v_10v)-1):
    dv_vort = (v_10v[n+1,:,:]-v_10v[n,:,:])/dx
    
vorticity = (du_vort-dv_vort)*(10**4)

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#divergence
xtick = np.linspace(75,100,len(lon))
ytick = np.linspace(20,5,len(lat))

plt.figure(figsize=(15,10),dpi=100)
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.LAND)
im=ax1.contourf(xtick,ytick,diverg,cmap ='jet')
sp= 3 
ax1.quiver(xtick[::sp],ytick[::sp],u_10v[0,::sp,::sp],v_10v[0,::sp,::sp],scale = 200,width = 0.002)      
ax1.gridlines(visible=True,draw_labels=True)
plt.colorbar(im,label ='Divergence (e^-4)')
ax1.set_title('Divergence of wind at June (2021-6-1 to 2021-6-10)')
plt.show()


#vorticty
xtick = np.linspace(75,100,len(lon))
ytick = np.linspace(20,5,len(lat))

plt.figure(figsize=(15,10),dpi=100)
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE)
ax2.add_feature(cfeature.LAND)
ax2.add_feature(cfeature.BORDERS)
iu = ax2.contourf(xtick,ytick,vorticity,cmap= 'jet')
sp=3
ax2.quiver(xtick[::sp],ytick[::sp],u_10v[0,::sp,::sp],v_10v[0,::sp,::sp],scale = 200,width = 0.002)
ax2.gridlines(visible=True,draw_labels=True)
plt.colorbar(iu,label ='Vorticity (e^-4)')
ax2.set_title('Vorticity of wind at June (2021-6-1 to 2021-6-10)')