# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 16:21:18 2025

@author: HP
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature


file = "/home/bobco-08/Desktop/24cl05012/CO2/data/oc_v2024E.flux.nc"
file2 = "/home/bobco-08/Desktop/24cl05012/CO2/data/monthly_sst_io.nc"

data = xr.open_dataset(file,decode_timedelta=True)
data2 = xr.open_dataset(file2,decode_timedelta = True)

#SST data variables

lon_sst = data2['longitude']
lat_sst = data2['latitude']
time_sst = data2['valid_time']
sst  = data2['sst']


#CO2 data variables

lat =  data['lat']
lon = data['lon']
time = data['mtime']
co2_ocean = data['co2flux_ocean']

co2_io = co2_ocean.sel(lat = slice(-40,30),lon = slice(30,120),mtime = slice("2000-1-1","2023-12-31")).mean(dim=('mtime'))

lon_io = lon.sel(lon =slice(30,120))
lat_io = lat.sel(lat = slice(-40,30))

#%% spatial

# Create figure and axis with PlateCarree projection
x = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot coastlines and land
ax.coastlines()
ax.add_feature(cfeature.LAND, color='lightgray', zorder=11)

# Plot the data with 'coolwarm' colormap
im = ax.contourf(lon_io, lat_io, co2_io, transform=ccrs.PlateCarree())

# Add gridlines with labels
gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10}
gl.ylabel_style = {'size': 10}

# Add colorbar with label
cbar = plt.colorbar(im, orientation='vertical', pad=0.05, aspect=30)
cbar.set_label("CO₂ Flux (PgC/yr)", fontsize=12)

# Show plot
plt.title("CO₂ Flux over Indian Ocean Climatology", fontsize=14)
plt.tight_layout()
plt.show()

#%% IOD and CO2 flux

co2_flux = co2_ocean.sel(lat = slice(2,30),lon = slice(30,100))

# Remove all Feb 29s from the time axis
co2_flux = co2_flux.sel(mtime=~((co2_flux['mtime.month'] == 2) & (co2_flux['mtime.day'] == 29)))


co2_monthly = co2_flux.resample(mtime = '1MS').mean()

co2_clim = co2_monthly.groupby('mtime.month').mean('mtime')

co2_anomaly = co2_monthly.groupby('mtime.month') - co2_clim #no mean so they will be mtime = 804

co2a_t = co2_anomaly.mean(dim=('lat','lon'))

co2a_3mo = co2a_t.rolling(mtime=3, center=True).mean()


#IOD index

sst_modified = sst.sel(valid_time = slice("1957-1-1","2023-12-31"))

sst_clim = sst_modified.groupby('valid_time.month').mean()

sst_anomaly = sst_modified.groupby('valid_time.month')- sst_clim

wbox = sst_anomaly.sel(latitude=slice(10, -10), longitude=slice(50, 70)).mean(dim=('latitude', 'longitude'))
ebox = sst_anomaly.sel(latitude=slice(0, -10), longitude=slice(90, 110)).mean(dim=('latitude', 'longitude'))

iod = wbox-ebox

iod_3mo = iod.rolling(valid_time=3, center=True).mean()


#%% plotting iod and co2


plt.subplot(2,1,2)
plt.plot(iod.valid_time,co2a_t, label='CO2 flux anomaly over IO ', color='black')

plt.grid()
plt.legend()
plt.title("Indian Ocean CO2 anomaly from 1957 to 2023")
plt.ylabel("PgC/yr")
plt.xlabel("Time")
plt.tight_layout()
plt.show()

plt.subplot(2,1,1)
plt.plot(iod.valid_time, iod, label='IOD Index', color='black')

# Add horizontal lines
plt.axhline(0.4, color='red', linestyle='--', linewidth=1)
plt.axhline(0, color='gray', linestyle='-', linewidth=1)
plt.axhline(-0.4, color='blue', linestyle='--', linewidth=1)

# Shade regions
plt.fill_between(iod.valid_time,iod, 0.4, where=(iod > 0.4), color='red', alpha=0.5,label ='PIOD')
plt.fill_between(iod.valid_time,iod, -0.4, where=(iod <-0.4), color='b', alpha=0.5,label ='NIOD')

plt.grid()
plt.legend()
plt.title("Indian Ocean Dipole (IOD) Phases")
plt.ylabel("IOD Index")
plt.xlabel("Time")
plt.tight_layout()
plt.show()


#%% 3 months runnig mean

plt.figure(figsize=(10, 8))
plt.subplot(2,1,2)
plt.plot(iod.valid_time,co2a_3mo, label='CO2 flux anomaly over IO ', color='black')
plt.title("Indian Ocean CO2 anomaly from 1957 to 2023")
plt.ylabel("PgC/yr")
plt.xlabel("Time")

plt.subplot(2,1,1)
plt.plot(iod.valid_time,iod_3mo, label='IOD Index', color='black')

# Add horizontal lines
plt.axhline(0.4, color='red', linestyle='--', linewidth=1)
plt.axhline(0, color='gray', linestyle='-', linewidth=1)
plt.axhline(-0.4, color='blue', linestyle='--', linewidth=1)

# Shade regions
plt.fill_between(iod_3mo.valid_time,iod_3mo, 0.4, where=(iod_3mo > 0.4), color='red', alpha=0.5,label ='PIOD')
plt.fill_between(iod_3mo.valid_time,iod_3mo, -0.4, where=(iod_3mo <-0.4), color='b', alpha=0.5,label ='NIOD')

plt.title("Indian Ocean Dipole (IOD) Phases")
plt.ylabel("IOD Index")
plt.xlabel("Time")


plt.legend()
plt.tight_layout()
plt.show()

#%% co2 monthly value
 

co2mon_mean = co2_flux.groupby("mtime.month").mean(dim=('mtime'))

co2mon_mean =  co2mon_mean.where(co2mon_mean != 0)

fig, axs = plt.subplots(3, 4, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

for i, ax in enumerate(axs.flat):
    
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2mon_mean[i], 
                     levels = 20 , cmap='coolwarm', transform=ccrs.PlateCarree())
    
    ax.set_title(months[i])
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND,color="gray",zorder=11)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5)
    gl.top_labels = gl.right_labels = False
    gl.left_labels = (i % 4 == 0)
    gl.bottom_labels = (i >= 8)  # 8–11 are bottom row in 3x4 grid

    
# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(label ='CO₂ Flux (PgC/yr)',labelpad = 5)


plt.suptitle('Monthly Mean CO₂ Flux (1957–2023)', fontsize=16)
plt.tight_layout(rect=[0, 0, 0.9, 0.95])
plt.show()

#%% 

co2sen_mean = co2_flux.groupby('mtime.season').mean(dim=("mtime"))

fig, axs = plt.subplots(2, 2, figsize=(24,24), subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw={'wspace': -0.3, 'hspace': 0.3})
season = ['DJF','MAM','JJAS',"SON"]

for i, ax in enumerate(axs.flat):
    im = ax.contourf(co2_flux.lon, co2_flux.lat, co2sen_mean[i], 
                     levels= np.linspace(-0.00225,0.00300,9) , cmap='coolwarm', transform=ccrs.PlateCarree(), extend ='both')
    ax.set_title(season[i])
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND,color="gray",zorder=11)
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
    gl.top_labels =  False
    gl.right_labels =  False

# Add colorbar
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.25, 0.50, 0.5, 0.015]) # left,bottom,width,height
cbar = fig.colorbar(im, cax=cbar_ax, orientation = 'horizontal')



plt.suptitle('Monthly Mean CO₂ Flux (1957–2023)', fontsize=16, y= 0.95)
plt.show()
#%%

co2mon_mean = co2_flux.groupby("mtime.month").mean(dim=('mtime'))
co2mon_mean =  co2mon_mean.where(co2mon_mean != 0)

ax=plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im=ax.contourf(co2_flux.lon, co2_flux.lat, co2mon_mean[1], 
                 levels=20 , cmap='viridis', transform=ccrs.PlateCarree())          
ax.gridlines(visible=True,draw_labels=True)
ax.add_feature(cfeature.LAND,color="gray",zorder=11)
plt.colorbar(im)
plt.show()


plt.title('Mean Wind Speed over Selected Region')
plt.xlabel('Longitude')
plt.ylabel('Wind Speed (m/s)')
plt.grid(True)
plt.show()

