#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 15:35:22 2026

@author: bobco-08
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmaps
# ===================== LOAD DATA =====================
path = "/home/bobco-08/24cl05012/npdl/oisst_global_monthly_mean_1982_2020.nc"
ds = xr.open_dataset(path)
sst = ds['sst'] # [time, lat, lon]
t = ds['time']

weights = np.cos(np.deg2rad(ds['lat']))

weight = weights/weights.mean()


sst_global = sst.weighted(weight).mean(dim=('lat', 'lon'))
#%%
# ===================== TIME AXIS (FRACTIONAL YEARS)
#=====================
years = (sst_global['time'].dt.year+ (sst_global['time'].dt.month - 0.5) / 12)
t = years.values
sst_vals = sst_global.values
# ===================== LINEAR TREND =====================
mask = np.isfinite(sst_vals)
slope, intercept, r, p, stderr = linregress(t[mask], sst_vals[mask])
trend_line = intercept + slope * t
trend_decade = slope * 10

#%%
# ===================== PLOT: TIME SERIES =====================
fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
ax.plot(
t, sst_global,
color='k',
linewidth=1.5,
label='Global Mean SST'
)
ax.plot(
t, trend_line,
color='red',
linestyle='--',
linewidth=2,
label=f'Trend = {trend_decade:.2f} °C / decade'
)
ax.set_xlabel('Year', fontsize=13)
ax.set_ylabel('Sea Surface Temperature (°C)', fontsize=13)
ax.set_title(
'Global Mean Sea Surface Temperature (1982–2020)',
fontsize=15,
weight='bold'
)
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(fontsize=11, frameon=False)
plt.tight_layout()
plt.show()

#%%  Regime in Global SST
import matplotlib.pyplot as plt
import ruptures as rpt



signal = sst_vals.astype(float).reshape(-1,+1)
time = sst['time'].values



algo = rpt.Pelt(model ='l2').fit(signal)
result = algo.predict(pen=1)

plt.figure(figsize=(12,5))
plt.plot(time, signal[:,0], 'k', lw=1.8, label="Global Mean SST")

for bp in result[:-1]:
    plt.axvline(time[bp], color='r', ls='--', lw=1)

plt.title("Change Point Detection in Global Mean SST", fontsize=16)
plt.xlabel("Time")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

#%%
regimes = []
start = 0

for bp in result:
    regimes.append(signal[start:bp, 0])
    start = bp
    
plt.figure(figsize=(12,5))

start = 0
for i, bp in enumerate(result):
    plt.plot(
        time[start:bp],
        signal[start:bp, 0],
        lw=2,
        label=f"Regime {i+1}"
    )
    start = bp

plt.title("Global Mean SST Regimes", fontsize=16)
plt.xlabel("Time")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.show()
#%%

regime_1 = sst.sel(time = slice('1982-01-31','2001-03-31')).mean(dim=('time'))
regime_2 = sst.sel(time = slice('2001-03-31','2014-07-31')).mean(dim=('time'))
regime_3 = sst.sel(time = slice('2014-07-31','2020-12-31')).mean(dim=('time'))

#%% Global SST regimes

import matplotlib.ticker as mticker

# ---------- Compute scalar color limits ----------
vmin = min(regime_1.min().item(), regime_2.min().item(), regime_3.min().item())
vmax = max(regime_1.max().item(), regime_2.max().item(), regime_3.max().item())

# ---------- Create 1x3 subplots ----------
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(18,15),
    dpi=300,
    subplot_kw=dict(projection=ccrs.Robinson())
)

fig.subplots_adjust(wspace=0.1)

regimes = [regime_1, regime_2, regime_3]
titles = [
    '(a) Mean SST (1982–2001)',
    '(b) Mean SST (2001–2014)',
    '(c) Mean SST (2014–2020)'
]

# ---------- Loop ----------
for ax, reg, title in zip(axes, regimes, titles):

    pcm = reg.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,
        levels=20,
        vmin=vmin,
        vmax=vmax,
        add_colorbar=False
    )

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)

    # ---------- Gridlines WITH labels ----------
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.4,
        color='gray',
        alpha=0.6,
        linestyle='--'
    )

    # Label placement (important for Robinson)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True

    # Tick locations
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 60))
    gl.ylocator = mticker.FixedLocator(range(-90, 91, 30))

    # Label style
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    ax.set_title(title, fontsize=12, weight='bold')

# ---------- Shared colorbar ----------
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.5,
    pad=0.05
)
cbar.set_label('Mean SST (°C)', fontsize=11)

plt.show()

#%% Global Regime Difference

import matplotlib.ticker as mticker

diff_1 = regime_2 - regime_1
diff_2 = regime_3 - regime_1
diff_3 = regime_3 - regime_2

diff_list = [diff_1, diff_2, diff_3]
titles = [
    '(a) Regime 2 (2001–2014) − Regime 1 (1982–2001)',
    '(b) Regime 3 (2014–2020) − Regime 1 (1982–2001)',
    '(c) Regime 3 (2014–2020) − Regime 2 (2001–2014)'
]

# ---------- Create 1x3 subplot ----------
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(22, 6),
    dpi=300,
    subplot_kw=dict(projection=ccrs.Robinson())
)
fig.subplots_adjust(wspace=0.1)

# ---------- Loop ----------
for ax, diff, title in zip(axes, diff_list, titles):

    pcm = diff.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,
        levels=21,
        vmin=-1,
        vmax=1,
        add_colorbar=False
    )

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)

    # ---------- Gridlines WITH lat/lon labels ----------
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.4,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )

    # Label placement (important for Robinson)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True

    # Tick locations
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 60))
    gl.ylocator = mticker.FixedLocator(range(-90, 91, 30))

    # Label styling
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    ax.set_title(title, fontsize=12, weight='bold')

# ---------- Shared colorbar ----------
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    fraction=0.055,
    pad=0.05
)
cbar.set_label('ΔSST (°C)', fontsize=11)

plt.show()


#%%   Region of Interest - Indian Ocean

sst_io = sst.sel(lat=slice(-30, 30), lon=slice(35, 110)).weighted(weight).mean(dim=('lat','lon'))


fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
ax.plot(
sst_io.time, sst_io,
color='b',
linewidth=1.5,
label='Global Mean SST'
)

ax.set_xlabel('Year', fontsize=13)
ax.set_ylabel('Sea Surface Temperature (°C)', fontsize=13)
ax.set_title(
'Mean SST over the Indian Ocean (30°S–30°N, 35°E–110°E)',
fontsize=15,
weight='bold'
)
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(fontsize=11, frameon=False)
plt.tight_layout()
plt.show()
#%% Change Point Detection Indian Ocean

signal_io = sst_io.values.astype(float).reshape(-1,+1)
time_io = sst_io['time'].values



algo = rpt.Pelt(model ='l2').fit(signal_io)
result_io = algo.predict(pen=3)

plt.figure(figsize=(12,5))
plt.plot(time_io, signal_io[:,0], 'b', lw=1.8, label="Global Mean SST")

for bp in result_io[:-1]:
    plt.axvline(time_io[bp], color='r', ls='--', lw=1)

plt.title("Change Point Detection in  SST over the Indian Ocean", fontsize=16)
plt.xlabel("Time")
plt.ylabel("SST (°C)")
plt.xlim(time_io.min(),time_io.max())
plt.legend()
plt.grid(alpha=0.3)
plt.show()
#%%  Regimes Indian Ocean


regio_1 = sst.sel(time=slice('1982-01-31', '1997-11-30'),lat=slice(-30, 30),lon=slice(35, 110)).mean(dim='time')
regio_2 = sst.sel(time=slice('1997-11-30', '2009-02-28'),lat=slice(-30, 30),lon=slice(35, 110)).mean(dim='time')
regio_3 = sst.sel(time=slice('2009-02-28', '2020-12-31'),lat=slice(-30, 30),lon=slice(35, 110)).mean(dim='time')

regions = [regio_1, regio_2, regio_3]
titles = [
    'Indian Ocean (a) Mean SST (1982–1997)',
    'Indian Ocean (b) Mean SST (1997–2009)',
    'Indian Ocean (c) Mean SST (2009–2020)'
]

vmin = min(r.min().item() for r in regions)
vmax = max(r.max().item() for r in regions)

sst_io_levels = np.linspace(18, 30, 9)

# -----------------------------
# Create 1x3 subplot
# -----------------------------
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(20, 10),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# -----------------------------
# Loop over regions
# -----------------------------
for ax, region, title in zip(axes, regions, titles):
    pcm = region.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,
        levels=20,
        vmin=vmin,
        vmax=vmax,
        add_colorbar=False
    )
    
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    
    # Set map extent for Indian Ocean
    ax.set_extent([35, 110, -30, 30], crs=ccrs.PlateCarree())
    
    # Add gridlines
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
  
    ax.set_title(title, fontsize=12, weight='bold')

# -----------------------------
# Shared horizontal colorbar
# -----------------------------
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.5,
    pad=0.05,
    ticks=sst_io_levels
)
cbar.set_label('Mean SST (°C)', fontsize=11)

plt.show()

#%% Difference Between Regimes

diff_1 = regio_2 - regio_1
diff_2 = regio_3 - regio_1
io_list = [diff_1, diff_2]

titles = [
    'Indian Ocean (a) 1997–2009 minus 1982–1997',
    'Indian Ocean (b) 2009–2020 minus 1982–1997'
]

# ---------- Create 1×2 subplot ----------
fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(18, 6),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)
fig.subplots_adjust(wspace=-0.3)

# ---------- Plot each panel ----------
for ax, sst_diff, title in zip(axes, io_list, titles):

    pcm = sst_diff.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,  # diverging colormap for differences
        vmin=-1,
        vmax=1,
        levels=21,extend = 'both',
        add_colorbar=False
    )

    # Coastlines & land
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)

    # Gridlines with labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    ax.set_title(title, fontsize=12, weight='bold')

# ---------- Shared horizontal colorbar ----------
# Use shrink to control length; fraction only affects thickness
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.6,  
    pad=0.08,
)
cbar.set_label('SST difference (°C)', fontsize=11)

plt.show()
#%% Labrador Sea - Highest  Positive trend 
sst_ls = sst.sel( lat=slice(53,63), lon=slice(302.08,317.2)).weighted(weight).mean(dim=('lat','lon'))


fig, ax1 = plt.subplots(figsize=(12, 5), dpi=300)
ax1.plot(
sst_ls.time, sst_ls,
color='b',
linewidth=1.5,
label='Mean SST'
)

ax1.set_xlabel('Year', fontsize=13)
ax1.set_ylabel('Sea Surface Temperature (°C)', fontsize=13)
ax1.set_title(
'Mean SST over the Labradore Sea region',
fontsize=15,
weight='bold'
)
ax1.grid(True, linestyle='--', alpha=0.4)
ax1.legend(fontsize=11, frameon=False)
plt.tight_layout()
plt.show()
#%% 

signal_ls = sst_ls.values.astype(float).reshape(-1,1)
time_ls = sst_ls['time'].values

# Binary segmentation with l2 model
algo = rpt.Binseg(model="l2", min_size=24).fit(signal_ls)

# Force exactly 3 change points
result_gs = algo.predict(n_bkps=3)
gs_yrs = np.array(result_gs)



for i in range(len(gs_yrs)-1):
    regime_years  = time_ls[gs_yrs[i]]
    print(regime_years)
    
plt.figure(figsize=(12,5))
plt.plot(time_ls, signal_ls[:,0], 'b', lw=1.8, label="Global Mean SST")

for bp in result_gs[:-1]:
    plt.axvline(time_ls[bp], color='r', ls='--', lw=1.8, alpha=1)
plt.title(f"Change Point Detection in SST\n Labrador Sea (53°N–63°N, 57.92°W – 42.80°W)", fontsize=16)

plt.xlabel("Year")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
#%% Regimes Labrador Sea

ls_1 = sst.sel(time=slice('1982-01-31', '1984-07-31'),lat=slice(53,63), lon=slice(302.08,317.2)).mean(dim='time')

ls_2 = sst.sel(time=slice('1984-07-31', '1997-06-30'),lat=slice(53,63), lon=slice(302.08,317.2)).mean(dim='time')

ls_3 = sst.sel(time=slice('1997-06-30', '2003-04-30'),lat=slice(53,63), lon=slice(302.08,317.2)).mean(dim='time')

ls_4 = sst.sel(time=slice('2003-04-30', '2020-12-31'),lat=slice(53,63), lon=slice(302.08,317.2)).mean(dim='time')


ls_list = [ls_1, ls_2, ls_3, ls_4]
titles = [
    'Labrador Sea (a) Mean SST (1982–1984)',
    'Labrador Sea (b) Mean SST (1984–1997)',
    'Labrador Sea (c) Mean SST (1997–2003)',
    'Labrador Sea (d) Mean SST (2003–2020)'
]

# -----------------------------
# Common color limits
# -----------------------------
vmin = min(arr.min().item() for arr in ls_list)
vmax = max(arr.max().item() for arr in ls_list)

# -----------------------------
# Create 2×2 subplots
# -----------------------------
fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(17, 14),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

fig.subplots_adjust(hspace=0.2, wspace=0.01)

axes_flat = axes.flatten()

# -----------------------------
# Plot each panel
# -----------------------------
for ax, sst_regime, title in zip(axes_flat, ls_list, titles):

    pcm = sst_regime.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,
        vmin=vmin,
        vmax=vmax,
        levels=20,
        add_colorbar=False
    )

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)

    # Gridlines with labels
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )

    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    ax.set_title(title, fontsize=12, weight='bold')

cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.6,
    pad=0.06
)
cbar.set_label('Mean SST (°C)', fontsize=11)
plt.show()
#%% Regime Difference Labradore sea

ls_diff1 = ls_2 - ls_1 
ls_diff2 = ls_3 - ls_1 
ls_diff3 = ls_4 - ls_1 


ls_list = [ls_diff1, ls_diff2, ls_diff3]
titles = [
    'Labrador Sea (a) 1984–1997 minus 1982–1984',
    'Labrador Sea (b) 1997–2003 minus 1982–1984',
    'Labrador Sea(c) 2003–2020 minus 1982–1984'
]

# -----------------------------
# Common color limits (SCALARS)
# -----------------------------
vmin = min(arr.min().item() for arr in ls_list)
vmax = max(arr.max().item() for arr in ls_list)

# Force symmetric limits for difference plots (recommended)
abs_max = max(abs(vmin), abs(vmax))
vmin, vmax = -abs_max, abs_max

# -----------------------------
# Create 1×3 subplots
# -----------------------------
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(22, 6),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Adjust spacing (horizontal only)
fig.subplots_adjust(wspace=0.15)

# -----------------------------
# Plot each panel
# -----------------------------
for ax, sst_diff, title in zip(axes, ls_list, titles):

    pcm = sst_diff.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,          # better for differences than jet
        vmin=vmin,
        vmax=vmax,
        levels=21,
        add_colorbar=False
    )

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)

    # Gridlines with labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    ax.set_title(title, fontsize=12, weight='bold')

# -----------------------------
# Shared colorbar
# -----------------------------
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.8,fraction = 0.08,
    pad=0.08,ticks =np.round(np.linspace(-3,3,10),2)
)
cbar.set_label('SST difference (°C)', fontsize=11)

plt.show()
#%% Antarctic region Lowest negative Trend

sst_ant = sst.sel( lat=slice(-65,-59), lon=slice(261.73,288.28)).weighted(weight).mean(dim=('lat','lon'))

fig, ax1 = plt.subplots(figsize=(12, 5), dpi=300)
ax1.plot(
sst_ant.time, sst_ant,
color='b',
linewidth=1.5,
label='Mean SST'
)

ax1.set_xlabel('Year', fontsize=13)
ax1.set_ylabel('Sea Surface Temperature (°C)', fontsize=13)
ax1.set_title(
'Mean SST over the Antartic region',
fontsize=15,
weight='bold'
)
ax1.grid(True, linestyle='--', alpha=0.4)
ax1.legend(fontsize=11, frameon=False)
plt.tight_layout()
plt.show()
#%%

signal_ant = sst_ant.values.astype(float).reshape(-1,1)
time_ant = sst_ant['time'].values

# Binary segmentation with l2 model
algo = rpt.Binseg(model="l2", min_size=24).fit(signal_ant)

# Force exactly 3 change points
result_ant = algo.predict(n_bkps=3)
ant_yrs = np.array(result_ant)



for i in range(len(ant_yrs)-1):
    regime_years  = time_ls[ant_yrs[i]]
    print(regime_years)
    
plt.figure(figsize=(12,5))
plt.plot(time_ant, signal_ant[:,0], 'b', lw=1.8, label="Global Mean SST")

for bp in result_ant[:-1]:
    plt.axvline(time_ant[bp], color='r', ls='--', lw=1)
plt.title(f"Change Point Detection in SST\n Antarctic region 34°N-40°N, 74°W – 68°W"
, fontsize=16)

plt.xlabel("Year")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
#%% Regimes Antarctic region

ant_1 = sst.sel(
    time=slice('1982-01-31', '1995-05-31'),
    lat=slice(-65,-59), lon=slice(261.73,288.28)
).mean(dim='time')

ant_2 = sst.sel(
    time=slice('1995-10-31', '2002-11-30'),
    lat=slice(-65,-59), lon=slice(261.73,288.28)
).mean(dim='time')

ant_3 = sst.sel(
    time=slice('2002-11-30', '2008-04-30'),
    lat=slice(-65,-59), lon=slice(261.73,288.28)
).mean(dim='time')

ant_4 = sst.sel(
    time=slice('2008-04-30', '2020-12-31'),
    lat=slice(-65,-59), lon=slice(261.73,288.28)
).mean(dim='time')


ant_list = [ant_1, ant_2, ant_3, ant_4]
titles = [
    'Antarctic Region (a) Mean SST (1982–1995)',
    'Antarctic Region (b) Mean SST (1995–2002)',
    'Antarctic Region (c) Mean SST (2002–2008)',
    'Antarctic Region (d) Mean SST (2008–2020)'
]

# -----------------------------
# Lock color limits (SCALARS)
# -----------------------------
vmin = min(arr.min().item() for arr in ant_list)
vmax = max(arr.max().item() for arr in ant_list)

# -----------------------------
# Create 2×2 subplot
# -----------------------------
fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(14,13),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())  # Antarctic projection
)

# Adjust spacing
fig.subplots_adjust(hspace=-0.6, wspace=0.15)

axes_flat = axes.flatten()

# -----------------------------
# Plot each panel
# -----------------------------
for ax, sst_regime, title in zip(axes_flat, ant_list, titles):

    pcm = sst_regime.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,
        vmin=vmin,
        vmax=vmax,
        levels=20,
        add_colorbar=False
    )

    # Coastlines & land
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)


    # Gridlines with labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    # Title
    ax.set_title(title, fontsize=12, weight='bold')

# -----------------------------
# Single shared colorbar
# -----------------------------
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.6,
    pad=0.08
)
cbar.set_label('Sea Surface Temperature (°C)', fontsize=11)

plt.show()
#%% Regime Difference Antarctic region

ant_diff1 = ant_2 - ant_1 
ant_diff2 = ant_3 - ant_1 
ant_diff3 = ant_4 - ant_1 

ls_list = [ant_diff1,ant_diff2, ant_diff3]
titles = [
    'Antarctic Region (a) 1995–2002 minus 1982–1995',
    'Antarctic Region (b) 2002–2008 minus 1982–1995',
    'Antarctic Region (c) 2008–2020 minus 1982–1995'
]

# -----------------------------
# Common color limits (SCALARS)
# -----------------------------
vmin = min(arr.min().item() for arr in ls_list)
vmax = max(arr.max().item() for arr in ls_list)

# Force symmetric limits for difference plots (recommended)
abs_max = max(abs(vmin), abs(vmax))
vmin, vmax = -abs_max, abs_max

# -----------------------------
# Create 1×3 subplots
# -----------------------------
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(22,8),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Adjust spacing (horizontal only)
fig.subplots_adjust(wspace=0.15)

# -----------------------------
# Plot each panel
# -----------------------------
for ax, sst_diff, title in zip(axes, ls_list, titles):

    pcm = sst_diff.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmaps.cmp_b2r,          # better for differences than jet
        vmin=vmin,
        vmax=vmax,
        levels=21,
        add_colorbar=False
    )

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)

    # Gridlines with labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    ax.set_title(title, fontsize=12, weight='bold')

# -----------------------------
# Shared colorbar
# -----------------------------
cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation='horizontal',
    shrink=0.7,fraction = 0.04,
    pad=0.08,ticks =np.round(np.linspace(-1,1,11),2)
)
cbar.set_label('SST difference (°C)', fontsize=11)

plt.show()