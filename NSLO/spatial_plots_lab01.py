# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:45:54 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 23:55:32 2025

@author: HP
"""

import numpy as np
import xarray as xr

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/u10_daily.nc"
file2 ="C:/Users/user/Documents/24CL05012/nslo/data/v10_daily.nc"


data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)

lat = data1['latitude']
lon = data1['longitude']
v = data2['v10']
u = data1['u10']
time = data1['valid_time']

latv = lat.values
lonv = lon.values

ws = np.sqrt((u**2)+(v**2))


import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature


# xticks = (30,120,lon)
# yticks = (30,0,lat)
# ws_sio = ws.sel(latitude = slice(30,0),longitude = slice(30,120))
ws_siom = ws.mean(dim=('valid_time')) 
# plt.figure(figsize=(15,10),dpi=100)
# ws_siom.plot(cmap ='jet')

ax=plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(lon,lat,ws_siom)          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im)
plt.show()


plt.title('Mean Wind Speed over Selected Region')
plt.xlabel('Longitude')
plt.ylabel('Wind Speed (m/s)')
plt.grid(True)
plt.show()

#%% sst


import numpy as np
import xarray as xr

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/sst-daily.nc"


data3 = xr.open_dataset(file1)


lat = data3['latitude']
lon = data3['longitude']
sst = data3['sst']-273
time = data3['valid_time']



import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature


# sst_sio = sst.sel(latitude = slice(30,0),longitude = slice(30,120))
sst_siom = sst.mean(dim=('valid_time')) 
# plt.figure()
# sst_siom.plot(cmap ='jet')


ax=plt.figure(figsize=(15,10),dpi=100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(lon,lat,sst_siom,cmap ='jet')          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im)
plt.show()

plt.title('Sea surface Temperature over NIO')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True)
plt.show()
#%% sp

import numpy as np
import xarray as xr

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/sp_daily.nc"


data4 = xr.open_dataset(file1)


lat = data4['latitude']
lon = data4['longitude']
sp = data4['sp']/100
time = data4['valid_time']



import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature



sp_siom = sp.mean(dim=('valid_time')) 

ax=plt.figure(figsize=(15,10),dpi=100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(lon,lat,sp_siom,cmap ='jet',vmin = 900,vmax =1040,levels = 1000)
# im.clim(900,1040)         
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="white",zorder=11)
plt.colorbar(im)
plt.show()


plt.title('Surface pressure over NIO')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True)
plt.show()
#%% t2m

import numpy as np
import xarray as xr

file1 ="C:/Users/user/Documents/24CL05012/nslo/data/t2m_daily.nc"


data5 = xr.open_dataset(file1)


lat = data5['latitude']
lon = data5['longitude']
temp = data5['t2m']-273
time = data5['valid_time']



import matplotlib.pyplot as plt


t2m_siom = temp.mean(dim=('valid_time')) 


ax=plt.figure(figsize=(15,10),dpi=100)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(lon,lat,t2m_siom,cmap ='jet')
# im.clim(900,1040)         
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="white",zorder=11)
plt.colorbar(im)
plt.show()

plt.title('2-meter temperature over NIO')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True)
plt.show()
