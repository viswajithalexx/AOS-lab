# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 14:39:29 2025

@author: HP
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = '/home/bobco-08/Desktop/24cl05012/ADWD/daily_sst.nc'

data = xr.open_dataset(file)

lat = data['latitude']
lon = data['longitude']
time = data['valid_time']

sst = data['sst']-273.15

sst_io = sst.sel(latitude = slice(30,-40),longitude = slice(40,100))

sst_a = sst.sel(latitude = 16.57,longitude = 88.10, method ='nearest')

corr = xr.corr(sst_io,sst_a,  dim='valid_time')

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im = ax.contourf(sst_io.longitude, sst_io.latitude, corr,
                   cmap='rainbow' ,vmin=-1, vmax=1)
# highlight the reference point
ax.scatter(88.10, 16.57, color='yellow', s=60, marker='*',
           transform=ccrs.PlateCarree(), label="Reference Point")

ax.gridlines(draw_labels=True)
ax.add_feature(cf.LAND, color="gray", zorder=11)
plt.colorbar(im, ax=ax, label="Correlation", pad =0.1)
# legend for reference point
plt.legend(loc='lower left')

plt.show()

#%%


fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ia = ax.contourf(sst_io.longitude, sst_io.latitude,sst_io.mean(dim=('valid_time')),
                   cmap='rainbow' ,levels =20)
# highlight the reference point
ax.scatter(88.10, 16.57, color='black', s=60, marker='o',
           transform=ccrs.PlateCarree(), label="Reference Point")

ax.gridlines(draw_labels=True)
ax.add_feature(cf.LAND, color="gray", zorder=11)
plt.colorbar(ia, ax=ax, label="SST($^\circ$C)", pad = 0.1)
# legend for reference point
plt.legend(loc='lower left')

plt.show()

#%%

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = '/home/bobco-08/Desktop/24cl05012/ADWD/daily_sst.nc'

data = xr.open_dataset(file)

lat = data['latitude']
lon = data['longitude']
time = data['valid_time']

# Convert from Kelvin to Celsius
sst = data['sst'] - 273.15  

# Select Indian Ocean region
sst_io = sst.sel(latitude=slice(22,-40), longitude=slice(20, 100))

# Reference point (Bay of Bengal)
sst_a = sst.sel(latitude=15, longitude=90, method='nearest')

# Correlation of reference point with rest of IO
corr = xr.corr(sst_io, sst_a, dim='valid_time')

# Plot
fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im = ax.contourf(sst_io.longitude, sst_io.latitude, corr,
                 cmap='seismic', vmin=-1, vmax=1)

# Add contour lines on top
contours = ax.contour(sst_io.longitude, sst_io.latitude, corr,
                      colors='black', linewidths=0.6, levels= 20)
ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours


# Add coastlines, land, grid
ax.coastlines()
ax.add_feature(cf.LAND, color="lightgray", zorder=11)
gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.5)
gl.top_labels = False
gl.right_labels = False

# Highlight reference point
ax.scatter(90,15, color='yellow', s=80, marker='*',
           transform=ccrs.PlateCarree(), label="Reference Point")

# Colorbar
cbar = plt.colorbar(im, ax=ax, orientation="vertical", pad=0.05, shrink=0.8)
cbar.set_label("Correlation Coefficient")

# Title and labels
plt.title("Correlation of Daily SST (2022)\nReference Point (16.57°N, 88.10°E) vs Indian Ocean",
          fontsize=14, weight='bold')
plt.legend(loc='lower left')

plt.show()
#%%

plt.figure()
plt.plot(sst_io.longitude,corr.sel(latitude = 5 ))

plt.grid()
plt.show()

#%%
plt.figure()
plt.plot(sst_io.latitude,corr.sel(longitude = 90))
plt.grid()
plt.show()
#%%
# Your base scatter plot
plt.figure()
plt.scatter(sst_io.longitude, corr.mean('latitude'), color="blue", label="All points")

# Choose the longitude you want to highlight
target_lon = 90  # for example
target_idx = (abs(sst_io.longitude - target_lon)).argmin().item()  # nearest index

# Plot the highlighted point
plt.scatter(sst_io.longitude[target_idx],
            corr.mean('latitude')[target_idx],
            color="red", s=100, marker="o", label=f"Lon {target_lon}")

plt.grid()
plt.legend()
plt.show()
#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

# -----------------------------
# Load data
# -----------------------------
path = "/home/bobco-08/Desktop/24cl05012/ADWD/daily_sst.nc"
data = xr.open_dataset(path)

sst = data["sst"]
lat = data["latitude"]
lon = data["longitude"]
time = data["valid_time"]

# -----------------------------
# Select target location
# -----------------------------
lat0, lon0 = 15, 90
target = sst.sel(latitude=lat0, longitude=lon0, method="nearest") - 273



# -----------------------------
# Annual correlation map
# -----------------------------
corr_annual = xr.corr(sst, target, dim="valid_time")

plt.figure(figsize=(10,4))
ax = plt.axes(projection=ccrs.PlateCarree())
cs = ax.contourf(lon, lat, corr_annual,
                 levels=np.linspace(-1,1,10),
                 cmap="bwr", transform=ccrs.PlateCarree())
ax.add_feature(cf.LAND, facecolor="lightgrey") # land in grey
ax.coastlines()
ax.gridlines(draw_labels=True, linestyle="--", alpha=0.5)
ax.plot(lon0, lat0, "k*", markersize=8)
plt.colorbar(cs, ax=ax, orientation="vertical", label="Correlation")
ax.set_title("Annual Correlation with SST at (15N, 65E)")
plt.show()

# -----------------------------
# Seasonal correlation maps
# -----------------------------
seasons = {"DJF":[12,1,2], "MAM":[3,4,5], "JJA":[6,7,8], "SON":[9,10,11]}
fig, axs = plt.subplots(2,2, figsize=(12,8),
                        subplot_kw={"projection": ccrs.PlateCarree()})

for ax, (season, months) in zip(axs.flat, seasons.items()):
    target_season = target.sel(valid_time=target['valid_time.month'].isin(months))
    sst_season = sst.sel(valid_time=sst['valid_time.month'].isin(months))
    corr_season = xr.corr(sst_season, target_season, dim="valid_time")

    cs = ax.contourf(lon, lat, corr_season,
                     levels=np.linspace(-1,1,10),
                     cmap="bwr", transform=ccrs.PlateCarree())
    ax.add_feature(cf.LAND, facecolor="lightgrey")
    ax.coastlines()
    ax.gridlines(draw_labels=True, linestyle="--", alpha=0.5)
    ax.plot(lon0, lat0, "k*", markersize=6)
    ax.set_title(f"{season} Correlation")

# shared colorbar
cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
fig.colorbar(cs, cax=cbar_ax, label="Correlation")
plt.suptitle("Seasonal Correlation with SST at (15N, 65E)", fontsize=14)
plt.tight_layout(rect=[0,0,0.9,0.95])
plt.show()


#%%

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Load data
# -----------------------------
path = "/home/bobco-08/Desktop/24cl05012/ADWD/daily_sst.nc"
data = xr.open_dataset(path)

sst = data["sst"] - 273  # K -> °C
lat = data["latitude"]
lon = data["longitude"]
time = data["valid_time"]

# -----------------------------
# Select target location
# -----------------------------
lat0, lon0 = 15, 90
target = sst.sel(latitude=lat0, longitude=lon0, method="nearest")

# -----------------------------
# Find grid index of target
# -----------------------------
ilat = np.argmin(np.abs(lat.values - lat0))
ilon = np.argmin(np.abs(lon.values - lon0))

# -----------------------------
# Parameters
# -----------------------------
step_deg = 0.25   # grid step in degrees
max_steps = 100   # how far south to go

# -----------------------------
# Compute correlation map (full climatology)
# -----------------------------
corr_map = xr.corr(sst, target, dim="valid_time")

# -----------------------------
# Extract southward profile
# -----------------------------
distances = []
south_corr = []
for step in range(1, max_steps+1):
    d = step * step_deg
    if ilat+step < len(lat):
        cval = float(corr_map.isel(latitude=ilat+step, longitude=ilon))
        distances.append(d)
        south_corr.append(cval)

distances = np.array(distances)
south_corr = np.array(south_corr)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(distances, south_corr, "o-", color="blue", label="Climatology Corr")

plt.axhline(0, color="k", lw=0.8)
plt.ylim(-1,1)
plt.xlabel("Distance southwards (°)")
plt.ylabel("Correlation with target SST")
plt.title(f"Southward Climatological Correlation\nfrom ({lat0}°N, {lon0}°E)")
plt.legend()
plt.grid(alpha=0.3)
plt.show()
#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf

file = '/home/bobco-08/Desktop/24cl05012/ADWD/daily_sst.nc'

data = xr.open_dataset(file)

lat = data['latitude']
lon = data['longitude']
time = data['valid_time']

# Convert from Kelvin to Celsius
sst = data['sst'] - 273.15  

# Select Indian Ocean region
sst_io = sst.sel(latitude=slice(22,-40), longitude=slice(20, 100))

# Reference point (Bay of Bengal)
sst_a = sst.sel(latitude= -22.77, longitude=90, method='nearest')

# Correlation of reference point with rest of IO
corr = xr.corr(sst_io, sst_a, dim='valid_time')

# Plot
fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
im = ax.contourf(sst_io.longitude, sst_io.latitude, corr,
                 cmap='seismic', vmin=-1, vmax=1)

# Add contour lines on top
contours = ax.contour(sst_io.longitude, sst_io.latitude, corr,
                      colors='black', linewidths=0.6, levels= 20)
ax.clabel(contours, inline=True, fontsize=8, fmt="%.1f")  # label the contours


# Add coastlines, land, grid
ax.coastlines()
ax.add_feature(cf.LAND, color="lightgray", zorder=11)
gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.5)
gl.top_labels = False
gl.right_labels = False

# Highlight reference point
ax.scatter(90,-22.77, color='yellow', s=80, marker='*',
           transform=ccrs.PlateCarree(), label="Reference Point")

# Colorbar
cbar = plt.colorbar(im, ax=ax, orientation="vertical", pad=0.05, shrink=0.8)
cbar.set_label("Correlation Coefficient")

# Title and labels
plt.title("Correlation of Daily SST (2022)\nReference Point (16.57°N, 88.10°E) vs Indian Ocean",
          fontsize=14, weight='bold')
plt.legend(loc='lower left')

plt.show()