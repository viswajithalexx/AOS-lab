# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 14:40:13 2026

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

#%%
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
index = nino_index.to_series().rolling(window=3, center=True).mean()
# --- Plot using actual datetime ---
time = index.index # pandas DatetimeIndex

y = index.values
plt.figure(figsize=(15,5))
plt.plot(time, y, color='black')
plt.axhline(0.5, linestyle='--', color='r')
plt.axhline(0, linestyle='-', color='black')
plt.axhline(-0.5, linestyle='--', color='b')
plt.fill_between(time, y, 0.5, where=(y > 0.5), color='red', alpha=0.5)
plt.fill_between(time, y, -0.5, where=(y < -0.5), color='b', alpha=0.5)
plt.xlabel('Year')
plt.ylabel('Detrended SST (°C)')
plt.title('SST Anomaly in Niño 3.4 Region (5°N–5°S, 120°–170°W)')
# ---- Year ticks formatting ----
ax = plt.gca()
ax.xaxis.set_major_locator(mdates.YearLocator(1)) # every year
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xlim(time.min(), time.max())
plt.xticks(rotation=45)
# --- Save --
plt.show()

#%% correlation

# Chunk inputs BEFORE computing correlation
nino_index = nino_index.chunk({'time': 468})

temp = temp.chunk({
    'time': 468,
    'lat': 90,
    'lon': 180
})

# Compute correlation
nino_sst_corr = xr.corr(nino_index, temp, dim='time')

#%%

fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude= 150))
plot_levels = np.arange(-1,1.2,0.2)
im2 = ax.contourf(nino_sst_corr.lon,nino_sst_corr.lat,nino_sst_corr,cmap=cmaps.MPL_RdYlBu_r,levels =plot_levels,transform=ccrs.PlateCarree())
ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='grey', zorder=10)
# gridlines (kept but corrected)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.2)
gl.top_labels = False
gl.right_labels = False
cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8,ticks = plot_levels)
cbar.set_label('correlation', fontsize=12)
ax.set_title('SST Nino Index correlation', fontsize=16)
plt.show()
#%% 

nino_index_lead = nino_index.shift(time = -3)

temp_lag = temp

nino_index_lead =nino_index_lead.chunk({'time': 468})

temp_lag = temp_lag.chunk({
    'time': 468,
    'lat': 90,
    'lon': 180
})

corr_3months_lead = xr.corr(nino_index_lead,temp_lag , dim='time')


fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude= 150))
plot_levels = np.arange(-1,1.2,0.2)
im2 = ax.contourf(corr_3months_lead.lon,corr_3months_lead.lat,corr_3months_lead,cmap=cmaps.MPL_RdYlBu_r,levels =plot_levels,transform=ccrs.PlateCarree())
ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='grey', zorder=10)
# gridlines (kept but corrected)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.2)
gl.top_labels = False
gl.right_labels = False
cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8,ticks = plot_levels)
cbar.set_label('correlation', fontsize=12)
ax.set_title('Correlation between Niño Index (3-month lead) and SST', fontsize=16)
plt.show()

#%%

nino_index_lag = nino_index.shift(time = 3)

temp_3_lead = temp

nino_index_lag  = nino_index_lag .chunk({'time': 468})

temp_3_lead = temp_3_lead.chunk({
    'time': 465,
    'lat': 90,
    'lon': 180
})

corr_3months_lag = xr.corr(nino_index_lag,temp_3_lead, dim='time')

fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude= 150))
plot_levels = np.arange(-1,1.2,0.2)
im2 = ax.contourf(corr_3months_lag.lon,corr_3months_lag.lat,corr_3months_lag,cmap=cmaps.MPL_RdYlBu_r,levels =plot_levels,transform=ccrs.PlateCarree())
ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='grey', zorder=10)
# gridlines (kept but corrected)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.2)
gl.top_labels = False
gl.right_labels = False
cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8,ticks = plot_levels)
cbar.set_label('correlation', fontsize=12)
ax.set_title('Correlation between Niño Index (3-month lag) and SST', fontsize=16)
plt.show()
#%%

wio = temp.sel(lat=slice(-10, 10), lon=slice(50, 70)).mean(dim=('lat','lon'))
eio = temp.sel(lat=slice(-10, 0),  lon=slice(90,110)).mean(dim=('lat','lon'))

wio_clim = wio.groupby('time.month').mean('time')
eio_clim = eio.groupby('time.month').mean('time')

wio_ano = wio.groupby('time.month') - wio_clim
eio_ano = eio.groupby('time.month') - eio_clim

dmi = wio_ano - eio_ano

dmi_norm = (dmi - dmi.mean('time')) / dmi.std('time')



#%%


import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# convert xarray time to matplotlib-friendly values
t = dmi['time'].values

# --- Figure setup ---
fig, ax = plt.subplots(figsize=(16, 5))

# --- Zero line ---
ax.axhline(0, color='k', linewidth=1)

# --- Main line plot ---
ax.plot(t, dmi, color='black', linewidth=1.5, label='IOD (DMI)')

# --- Positive & negative shading (based on DMI itself) ---
ax.fill_between(t, dmi, 0, where=(dmi > 0), interpolate=True,
                alpha=0.5, label='Positive IOD')

ax.fill_between(t, dmi, 0, where=(dmi < 0), interpolate=True,
                alpha=0.5, label='Negative IOD')

# --- Labels & title ---
ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('DMI (°C)', fontsize=12)
ax.set_title('Dipole Mode Index (IOD)', fontsize=14, fontweight='bold')

# --- Year ticks formatting ---
ax.xaxis.set_major_locator(mdates.YearLocator(1))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.set_xlim(t.min(), t.max())
plt.xticks(rotation=45)

# --- Grid & legend ---
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(loc='upper right', frameon=False)

# --- Layout ---
plt.tight_layout()
plt.show()

#%%
dmi_neg = dmi.where(dmi > 0)

dmi_pd = dmi_neg.to_series()


for year, group in dmi_pd.groupby(dmi_pd.index.year):
    if len(group.dropna()) == 0:
        continue
    max_val = group.max()
    max_time = group.idxmax()
    print(f"Year {year}: max DMI = {max_val:.2f} at {max_time.strftime('%Y-%m')}")


#%%

io = temp.sel(lat = slice(-30,30), lon = slice(35,120))

io_aug = io.sel(time = time.month == 8)

io_aug = io_aug.groupby('time.year').mean('time')

yr = [1989, 1991, 1992, 1998, 2009, 2014, 2020]

io_aug_yr = io_aug.sel(year=yr)

io_aug_yr = io_aug_yr.mean('year')

io_aug_clim = io_aug.mean(dim =('year'))

iod_negative = io_aug_yr - io_aug_clim

#%%

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np

# --- Figure ---
fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=150))

plot_levels = np.linspace(-1, 1, 21)

# --- SST anomaly map ---
im2 = ax.contourf(
    iod_negative.lon,
    iod_negative.lat,
    iod_negative,
    cmap=cmaps.MPL_RdYlBu_r,
    levels =plot_levels,
    extend='both',
    transform=ccrs.PlateCarree()
)

# --- Coast & land ---
ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='lightgrey', zorder=10)

# --- Gridlines ---
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.3)
gl.top_labels = False
gl.right_labels = False

# --- Colorbar ---
cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8, ticks =plot_levels)
cbar.set_label('SST Anomaly (°C)', fontsize=12)

# --- Title ---
ax.set_title('Composite August SST Anomaly during Negative IOD Years', fontsize=16)

# =========================
#        IOD BOXES
# =========================

# West Indian Ocean (WIO)
wio_box = mpatches.Rectangle(
    (50, -10), 20, 20,
    fill=False, edgecolor='black', linewidth=2,
    transform=ccrs.PlateCarree(), zorder=20
)

# East Indian Ocean (EIO)
eio_box = mpatches.Rectangle(
    (90, -10), 20, 10,
    fill=False, edgecolor='black', linewidth=2,
    transform=ccrs.PlateCarree(), zorder=20
)

ax.add_patch(wio_box)
ax.add_patch(eio_box)

# --- Box labels (optional but nice) ---
ax.text(60, 10, 'WIO', transform=ccrs.PlateCarree(),
        ha='center', va='bottom', fontsize=15, fontweight='bold')

ax.text(95, 0, 'EIO', transform=ccrs.PlateCarree(),
        ha='center', va='bottom', fontsize=15, fontweight='bold',zorder = 15)

# --- Show ---
plt.tight_layout()
plt.show()

#%%

io_mar = io.sel(time=io['time.month'] == 3)


io_mar_yr = io_mar.groupby('time.year').mean('time')


yr_mar = [1986, 1999, 2000, 2001, 2004, 2009]
io_mar_pos = io_mar_yr.sel(year=yr_mar).mean(dim=('year'))


io_mar_clim = io_mar_yr.mean('year')


iod_positive_mar = io_mar_pos - io_mar_clim

#%%

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np

# --- Figure ---
fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=150))

plot_levels = np.linspace(-1, 1, 21)

# --- SST anomaly map ---
im2 = ax.contourf(
    iod_positive_mar.lon,
    iod_positive_mar.lat,
    iod_positive_mar,
    cmap=cmaps.MPL_RdYlBu_r,
    levels =plot_levels,
    extend='both',
    transform=ccrs.PlateCarree()
)

# --- Coast & land ---
ax.coastlines(zorder=15)
ax.add_feature(cf.LAND, color='lightgrey', zorder=10)

# --- Gridlines ---
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='grey', alpha=0.3)
gl.top_labels = False
gl.right_labels = False

# --- Colorbar ---
cbar = fig.colorbar(im2, ax=ax, orientation='vertical', shrink=0.8,ticks =plot_levels)
cbar.set_label('SST Anomaly (°C)', fontsize=12)

# --- Title ---
ax.set_title('Composite March SST Anomaly during Positive IOD Years', fontsize=16)

# =========================
#        IOD BOXES
# =========================

# West Indian Ocean (WIO)
wio_box = mpatches.Rectangle(
    (50, -10), 20, 20,
    fill=False, edgecolor='black', linewidth=2,
    transform=ccrs.PlateCarree(), zorder=20
)

# East Indian Ocean (EIO)
eio_box = mpatches.Rectangle(
    (90, -10), 20, 10,
    fill=False, edgecolor='black', linewidth=2,
    transform=ccrs.PlateCarree(), zorder=20
)

ax.add_patch(wio_box)
ax.add_patch(eio_box)

# --- Box labels (optional but nice) ---
ax.text(60, 10, 'WIO', transform=ccrs.PlateCarree(),
        ha='center', va='bottom', fontsize=20, fontweight='bold')

ax.text(95, 0, 'EIO', transform=ccrs.PlateCarree(),
        ha='center', va='bottom', fontsize=20, fontweight='bold',zorder = 15)

# --- Show ---
plt.tight_layout()
plt.show()


#%%

import numpy as np
import matplotlib.pyplot as plt

fs = 22
maxlag = 12

# --- Align and remove NaNs ---
mask = (~np.isnan(nino_index.values)) & (~np.isnan(dmi.values))
time_v = time[mask]
nino = nino_index.values[mask]
iod  = dmi.values[mask]

# --- Normalize (z-score) ---
nino_n = (nino - np.mean(nino)) / np.std(nino)
iod_n  = (iod  - np.mean(iod))  / np.std(iod)

# =========================
# Time series plot
# =========================
thr = 1.0
fig, ax = plt.subplots(figsize=(16,7))

ax.plot(time_v, iod_n,  lw=3, label="DMI (normalized)")
ax.plot(time_v, nino_n, lw=3, label="Nino (normalized)")

ax.axhline(0, color="k", lw=1.5)
ax.axhline(+thr, color="k", lw=1.5, ls="--", alpha=0.7)
ax.axhline(-thr, color="k", lw=1.5, ls="--", alpha=0.7)

ax.set_title("ENSO (Nino) and IOD (DMI) Indices", fontsize=fs)
ax.set_xlabel("Time", fontsize=fs)
ax.set_ylabel("Normalized Index (σ)", fontsize=fs)
ax.tick_params(labelsize=fs-2)
ax.grid(True, ls="--", alpha=0.4)
ax.legend(fontsize=fs-4, frameon=False)
plt.show()

# =========================
# Lead–lag correlation
# =========================
lags = np.arange(-maxlag, maxlag+1)
cc = np.zeros(len(lags))

for i, L in enumerate(lags):
    if L > 0:     # ENSO leads
        cc[i] = np.corrcoef(nino_n[:-L], iod_n[L:])[0,1]
    elif L < 0:   # IOD leads
        k = abs(L)
        cc[i] = np.corrcoef(nino_n[k:], iod_n[:-k])[0,1]
    else:
        cc[i] = np.corrcoef(nino_n, iod_n)[0,1]

# --- Max correlation ---
imax = np.nanargmax(cc)
best_lag = lags[imax]
best_corr = cc[imax]

if best_lag > 0:
    txt = f"ENSO leads IOD by {best_lag} months"
elif best_lag < 0:
    txt = f"IOD leads ENSO by {abs(best_lag)} months"
else:
    txt = "Simultaneous"

print(f"Max correlation = {best_corr:.2f} at lag = {best_lag} months → {txt}")

# =========================
# Plot correlation
# =========================
fig, ax = plt.subplots(figsize=(12,5))
ax.plot(lags, cc, lw=2, marker="o")
ax.axhline(0, color="k", lw=1)
ax.axvline(0, color="k", lw=1)
ax.axvline(best_lag, color="r", ls="--", lw=2, label=txt)

ax.set_title("Lead–Lag Correlation: ENSO vs IOD", fontsize=fs)
ax.set_xlabel("Lag (months)   [ +ve = ENSO leads IOD ]", fontsize=fs)
ax.set_ylabel("Correlation", fontsize=fs)
ax.tick_params(labelsize=fs-2)
ax.grid(True, ls="--", alpha=0.4)
ax.legend(fontsize=fs-4)
plt.show()
