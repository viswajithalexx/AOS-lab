# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 09:12:10 2026

@author: user
"""

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cf
import xarray as xr
import cmaps

file = "/media/bobco-08/valex/oisst_global_monthly_mean_1982_2020.nc"
data = xr.open_dataset(file)

lon = data['lon']
lat = data['lat']
t = data['time']
temp = data['sst']

temp_m = temp.groupby('time.month').mean('time')



sstio_j = temp_m.sel(month = 1)
sstio_a = temp_m.sel(month = 4)
sstio_s = temp_m.sel(month = 9)
sstio_o = temp_m.sel(month = 10)

months = [sstio_j,sstio_a,sstio_s,sstio_o]
mn =['JAN','APR','SEP','OCT']
#%%
fig, axs = plt.subplots(
    2, 2,
    figsize=(10,8),dpi= 300,
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'hspace':-0.42, 'wspace': 0.2}
)

for i, ax in enumerate(axs.flat):
    
    im = ax.contourf(sstio_j.lon,sstio_j.lat,months[i],levels = 20 ,cmap = cmaps.MPL_RdYlBu_r,transform=ccrs.PlateCarree(),extend ='max')
    ax.coastlines(zorder =15)
    ax.add_feature(cf.LAND, color = 'grey',zorder = 10)
    
    gl = ax.gridlines(draw_labels = True,visible = False,linewidth = 0.5 , color = 'grey', alpha =0.2)
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'fontsize':10}
    gl.ylabel_style = {'fontsize':10}
    
    ax.set_title(f'{mn[i]}')
    
cbar_ax = fig.add_axes([0.92, 0.28, 0.01, 0.45])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im, cax=cbar_ax ,orientation='vertical')
cbar.set_label(label='°C',fontsize= 10)
plt.savefig("/home/bobco-08/24cl05012/npdl/sst_global.png", dpi=500, bbox_inches="tight")
plt.show()


#%% tropical io

temp_io = temp_m.sel(lat = slice(-30,30),lon= slice(30,120))

mn =['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

fig1, axs = plt.subplots(
    3,4,
    figsize=(15,10),
    subplot_kw={'projection': ccrs.PlateCarree()},
    gridspec_kw={'hspace':-0.40, 'wspace': 0.28}
)

plot_level = np.arange(25,30,0.1)

for i, ax in enumerate(axs.flat):
    
    im_io = ax.contourf(temp_io.lon,temp_io.lat,temp_io[i],levels = plot_level,cmap = cmaps.MPL_RdYlBu_r,transform=ccrs.PlateCarree(),extend ='both')
    ax.coastlines(zorder =15)
    ax.add_feature(cf.LAND, color = 'grey',zorder = 10)
    
    gl = ax.gridlines(draw_labels = True,visible = False,linewidth = 0.5 , color = 'grey', alpha =0.2)
    gl.right_labels = False
    gl.top_labels = False
    gl.xlabel_style = {'fontsize': 10}
    gl.ylabel_style = {'fontsize': 10}
    
    ax.set_title(f'{mn[i]}')
    
cbar_ax = fig1.add_axes([0.92, 0.28, 0.01, 0.45])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(im_io, cax=cbar_ax ,orientation='vertical')
cbar.set_label(label='°C',fontsize= 10)
plt.savefig("/home/bobco-08/24cl05012/npdl/itcz_I.png", dpi=500, bbox_inches="tight")
plt.show()



#%% nino
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr

# --- Your preprocessing (unchanged) ---
sst_pacific = temp.sel(lat=slice(-5,5), lon=slice(190,240))

cof = sst_pacific.polyfit(dim='time', deg=1)
trends = xr.polyval(sst_pacific['time'], cof.polyfit_coefficients)
dtrend = sst_pacific - trends

sst_clim = dtrend.groupby('time.month').mean()
sst_monthly = dtrend.groupby('time.month')
sst_anomalies = sst_monthly - sst_clim

nino_index = sst_anomalies.mean(dim=('lat','lon'))

# 3-month running mean
index = nino_index.to_series().rolling(window=5, center=True).mean()

# --- Plot using actual datetime ---
time = index.index   # pandas DatetimeIndex
y = index.values

plt.figure(figsize=(15,5))

plt.plot(time, y, color='black')
plt.axhline(0.4, linestyle='--', color='r')
plt.axhline(0, linestyle='-', color='black')
plt.axhline(-0.4, linestyle='--', color='b')

plt.fill_between(time, y, 0.4, where=(y > 0.5), color='red', alpha=0.5)
plt.fill_between(time, y, -0.4, where=(y < -0.5), color='b', alpha=0.5)

plt.xlabel('Year')
plt.ylabel('Detrended SST (°C)')
plt.title('SST Anomaly in Niño 3.4 Region (5°N–5°S, 170°–120°W)')

# ---- Year ticks formatting ----
ax = plt.gca()
ax.xaxis.set_major_locator(mdates.YearLocator(1))     # every year
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xlim(time.min(), time.max())
plt.xticks(rotation=45)

# --- Save ---
plt.savefig("/home/bobco-08/24cl05012/npdl/nino_3.4.png",
            dpi=500, bbox_inches="tight")

plt.show()

#%%
# Boolean mask for values > 0.5
index_peak = index > 0.5

# Get actual values where True
peak_values = index[index_peak]

# Print
print(peak_values)
#%%
import numpy as np
import pandas as pd

# Your Series of values > 0.5
peak_values = index[index < -0.5]

# Group by YEAR instead of index//12
groups = peak_values.groupby(peak_values.index.year)

# Loop through each year to find the peak
for year, group in groups:
    min_val = group.min()
    min_time = group.idxmin()   # datetime of peak
    print(f"Year {year}: Peak value = {min_val:.3f} at {min_time.strftime('%Y-%m')}")


#%%


r =  [11, 13, 59, 68, 72, 119, 121, 136, 155, 156, 191, 192, 251, 252, 273, 276, 299, 300, 335, 336, 395, 407, 408, 443, 448
]

# r = array of peak indices
# t = your time array (e.g., nino_index['time'].values)

peak_months = []  # store the corresponding times

for i in range(len(r)):
    peak_months.append(t[r[i]])

# If you want as a numpy array
#%%

# 1. Detrend full SST (grid-point wise)
cof = temp.polyfit(dim='time', deg=1)
trend = xr.polyval(temp.time, cof.polyfit_coefficients)
dtrend = temp - trend

sst_dec = dtrend.sel(time=dtrend.time.dt.month == 12)

# 3. One December per year
dec_year = sst_dec.groupby('time.year').mean('time')

# 4. El Niño years
yr = [1982,1986,1991,1994,1997,2006,2009,2014,2015]

# El Niño December fields
el_nino_dec = dec_year.sel(year=yr)

# December climatology (mean of El Niño Decembers)
dec_clim = dec_year.mean(dim='year')

# December anomalies
dec_anom = el_nino_dec - dec_clim

# El Niño December composite anomaly
el_ano = dec_anom.mean(dim='year')

#%%

fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=150)) 
plot_levels = np.linspace(-2.5, 2.5, 21)  
im2 = ax.pcolormesh(
    el_ano.lon,
    el_ano.lat,
    el_ano.values,        
    cmap=cmaps.cmp_b2r,
    transform=ccrs.PlateCarree(),
    vmin=-2.5,  # set min color value
    vmax=2.5    # set max color value
)

ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='grey', zorder=10)

# gridlines (kept but corrected)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.2)
gl.top_labels = False
gl.right_labels = False

cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8,ticks =plot_levels)
cbar.set_label('SST Anomaly (°C)', fontsize=12)

ax.set_title('Composite December SST Anomaly during El Niño Years', fontsize=16)

plt.savefig("/home/bobco-08/24cl05012/npdl/december_el_nino_sst_anomaly_cbar.png",
            dpi=500, bbox_inches="tight")

plt.show()

#%%
sst_jan = dtrend.sel(time=dtrend.time.dt.month == 1)

# 3. One December per year
jan_year = sst_jan.groupby('time.year').mean('time')

# 4. El Niño years
yr = [1985, 1989, 1996, 2000, 2001, 2006, 2008, 2009, 2011, 2012, 2013, 2018]


# El Niño December fields
la_nina_jan = jan_year.sel(year=yr)

# December climatology (mean of El Niño Decembers)
jan_clim = jan_year.mean(dim='year')

# December anomalies
jan_anom = la_nina_jan - jan_clim

# El Niño December composite anomaly
la_ano = jan_anom.mean(dim='year')
#%%
fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=150))  
plot_levels = np.linspace(-2.5, 2.5, 21)
im2 = ax.pcolormesh(
    la_ano.lon,
    la_ano.lat,
    la_ano.values,          
    cmap=cmaps.MPL_RdYlBu_r,
    transform=ccrs.PlateCarree(),
    vmin=-2.5,  # set min color value
    vmax=2.5 
)

ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='grey', zorder=10)

# gridlines (kept but corrected)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.2)
gl.top_labels = False
gl.right_labels = False

cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8,ticks =plot_levels)
cbar.set_label('SST Anomaly (°C)', fontsize=12)

ax.set_title('Composite January SST Anomaly during La Nina Years', fontsize=16)

plt.savefig("/home/bobco-08/24cl05012/npdl/january_el_nino_sst_anomaly_cbar.png",
            dpi=500, bbox_inches="tight")

plt.show()