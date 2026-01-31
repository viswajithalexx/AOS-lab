# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 14:40:55 2026

@author: user
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.stats import t as t_dist
# ===================== LOAD DATA =====================
file_path = "C:/Users/user/Documents/24CL05012/4th/NPDL2026/oisst_global_monthly_mean_1982_2020.nc"
ds_sst = xr.open_dataset(file_path, engine="netcdf4")
sst_field = ds_sst["sst"]
lat_arr = ds_sst["lat"]
lon_arr = ds_sst["lon"]
time_arr = ds_sst["time"]
# ===================== TIME TO FLOAT YEARS =====================
time_float = ((time_arr - time_arr.isel(time=0)) / np.timedelta64(1, "D")).values / 365.25
time_float = time_float.astype(np.float64)
sst_float_time = sst_field.assign_coords(time=("time", time_float))
# ===================== SIMPLE TREND =====================
trend_model = sst_float_time.polyfit(dim="time", deg=1)
coef_slope = trend_model["polyfit_coefficients"].sel(degree=1) # °C/year
trend_10yr = coef_slope * 10 # °C/decade
# ===================== PLOT 1: RAW TREND =====================
#%%
proj = ccrs.Robinson()
fig1, ax1 = plt.subplots(
figsize=(12, 5),
subplot_kw={"projection": proj}
)
levels = np.linspace(-0.5, 0.5, 11)
map1 = ax1.contourf(
lon_arr, lat_arr, trend_10yr,
levels=levels,
cmap="bwr",
extend="both",
transform=ccrs.PlateCarree()
)
ax1.add_feature(cfeature.LAND, facecolor="lightgrey", zorder=10)
ax1.coastlines(linewidth=0.8)
ax1.set_global()
ax1.set_title("SST Trend (°C per decade)", fontsize=18)
cb1 = plt.colorbar(map1, ax=ax1, orientation="horizontal", fraction=0.045,pad=0.05,
aspect=40)
cb1.set_label("°C / decade", fontsize=14)
cb1.ax.tick_params(labelsize=12)
plt.show()

#%%
# ===================== TREND + P-VALUE FUNCTION
#=====================
def linear_trend_pvalue(series, time_x):
    series = np.asarray(series, dtype=float)
    time_x = np.asarray(time_x, dtype=float)

    good = np.isfinite(series) & np.isfinite(time_x)
    series = series[good]
    time_x = time_x[good]

    n = series.size
    if n < 10:
        return np.nan, np.nan

    xm = time_x.mean()
    ym = series.mean()

    Sxx = np.sum((time_x - xm) ** 2)
    Sxy = np.sum((time_x - xm) * (series - ym))

    beta = Sxy / Sxx
    alpha = ym - beta * xm

    residuals = series - (beta * time_x + alpha)
    dof = n - 2

    std_err = np.sqrt(np.sum(residuals ** 2) / dof)
    beta_se = std_err / np.sqrt(Sxx)

    t_stat = beta / beta_se
    p_val = 2 * (1 - t_dist.cdf(np.abs(t_stat), dof))

    return beta, p_val

    
#%%    
# ===================== APPLY TREND + SIGNIFICANCE
#=====================
trend_coef, p_values = xr.apply_ufunc(
linear_trend_pvalue,
sst_field,
xr.DataArray(time_float, dims=["time"], coords={"time": time_arr}),
input_core_dims=[["time"], ["time"]],
output_core_dims=[[], []],
vectorize=True,
dask="parallelized",
output_dtypes=[float, float],
)
trend_dec = trend_coef * 10
trend_sig = trend_dec.where(p_values < 0.05)

#%%
# ===================== PLOT 2: SIGNIFICANT TREND
#=====================
fig2, ax2 = plt.subplots(
figsize=(12, 5),
subplot_kw={"projection": proj}
)
map2 = ax2.contourf(
lon_arr, lat_arr, trend_sig,

levels=levels,
cmap="bwr",
extend="both",
transform=ccrs.PlateCarree()
)
cbar = fig2.colorbar(map2, ax=ax2, orientation="horizontal", fraction=0.045,pad=0.05,
aspect=40)
ax2.add_feature(cfeature.LAND, facecolor="lightgrey", zorder=10)
ax2.coastlines(linewidth=0.8)
ax2.set_global()
ax2.set_title("SST Trend (°C/decade), Significant at p < 0.05", fontsize=18)

#%%
# ---- STAR MARKERS FOR SIGNIFICANCE ----
lon_grid, lat_grid = np.meshgrid(lon_arr, lat_arr)
skip = 10
sig_mask = (p_values < 0.05).values
sig_sub = sig_mask[::skip, ::skip]
cb2 = plt.colorbar(map2, ax=ax2, orientation="horizontal", fraction=0.045, pad=0.05,
aspect=40)
cb2.set_label("°C / decade", fontsize=14)
cb2.ax.tick_params(labelsize=12)
plt.show()
#%%
# ===================== PLOT 3: INSIGNIFICANT MARKED
#=====================
fig3, ax3 = plt.subplots(
figsize=(12, 5),
subplot_kw={"projection": proj}
)
lvl2 = np.linspace(-0.2, 0.5, 11)
map3 = ax3.contourf(
lon_arr, lat_arr, trend_dec,
levels=lvl2,
cmap="jet",
extend="both",
transform=ccrs.PlateCarree()
)
ax3.add_feature(cfeature.LAND, facecolor="lightgrey", zorder=10)
ax3.coastlines(linewidth=0.8)

ax3.set_global()
ax3.set_title("SST Trend (°C/decade): * = Insignificant", fontsize=18)
bad_mask = (p_values >= 0.05).values
bad_sub = bad_mask[::skip, ::skip]
ax3.scatter(
lon_grid[::skip, ::skip][bad_sub],
lat_grid[::skip, ::skip][bad_sub],
s=8,
marker="*",
color="k",
transform=ccrs.PlateCarree(),
zorder=20
)
cb3 = plt.colorbar(
map3,
ax=ax3,
orientation="horizontal",
fraction=0.045, # ↓ decrease this to shorten length
pad=0.05,
aspect=40
)
cb3.set_label("°C / decade", fontsize=14)
cb3.ax.tick_params(labelsize=12)
plt.show()

#%%

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import cartopy.crs as ccrs
import cartopy.feature as cfeature
file_loc = "C:/Users/user/Documents/24CL05012/4th/NPDL2026/oisst_global_monthly_mean_1982_2020.nc"
dataset = xr.open_dataset(file_loc)
sst_data = dataset['sst']
sst_data = sst_data.assign_coords(lon=((sst_data.lon + 360) % 360))
sst_data = sst_data.sortby('lon')
#%%
mean_period1 = sst_data.sel(time=slice('1982-01-01', '1996-12-31')).mean('time')
mean_period2 = sst_data.sel(time=slice('1997-01-01', '2020-12-31')).mean('time')
mean_diff = mean_period2 - mean_period1
#%%
fig = plt.figure(figsize=(14,6), dpi=300)

ax = plt.axes(projection=ccrs.Robinson())
mean_period1.plot(
ax=ax,
transform=ccrs.PlateCarree(),
cmap='turbo',
cbar_kwargs=dict(label='Mean SST (°C)', shrink=0.85)
)
ax.coastlines(linewidth=0.8)
ax.set_title('Mean SST (1982–1996)', fontsize=14, weight='bold')
plt.tight_layout()
plt.show()
#%%
fig = plt.figure(figsize=(14,6), dpi=300)
ax = plt.axes(projection=ccrs.Robinson())
mean_period2.plot(
ax=ax,
transform=ccrs.PlateCarree(),
cmap='turbo',
cbar_kwargs=dict(label='Mean SST (°C)', shrink=0.85)
)
ax.coastlines(linewidth=0.8)
ax.set_title('Mean SST (1997–2020)', fontsize=14, weight='bold')
plt.tight_layout()
plt.show()
#%%
fig = plt.figure(figsize=(14,6), dpi=300)
ax = plt.axes(projection=ccrs.Robinson())
diff_levels = np.arange(-1.0, 1.05, 0.1)
mean_diff.plot.contourf(
ax=ax,
transform=ccrs.PlateCarree(),
levels=diff_levels,
cmap='RdBu_r',
extend='both',

cbar_kwargs=dict(
label='ΔSST (°C)',
shrink=0.85,
pad=0.05
)
)
ax.coastlines(linewidth=0.8)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
ax.set_title(
'Change in Mean SST (1997–2020 minus 1982–1996)',
fontsize=14,
weight='bold'
)
plt.tight_layout()
plt.show()
#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# ===================== LOAD DATA =====================
path = "C:/Users/user/Documents/24CL05012/4th/NPDL2026/oisst_global_monthly_mean_1982_2020.nc"
ds = xr.open_dataset(path)
sst = ds['sst'] # [time, lat, lon]
t = ds['time']

weights = np.cos(np.deg2rad(ds['lat']))

weight = weights/weights.mean()


sst_global = sst.weighted(weight).mean(dim=('lat', 'lon'))
#%%
# ===================== TIME AXIS (FRACTIONAL YEARS)
#=====================
years = (

sst_global['time'].dt.year
+ (sst_global['time'].dt.month - 0.5) / 12
)
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

#%%
import matplotlib.pyplot as plt
import ruptures as rpt



signal = sst_vals.astype(float).reshape(-1,+1)
time = sst['time'].values



algo = rpt.Pelt(model ='l2').fit(signal)
result = algo.predict(pen=10)

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

#%%

fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(18, 6),
    dpi=300,
    subplot_kw=dict(projection=ccrs.Robinson())
)

# Optional: lock color limits for fair comparison
vmin = min(regime_1.min(), regime_2.min(), regime_3.min())
vmax = max(regime_1.max(), regime_2.max(), regime_3.max())

# ---------- Panel (a) ----------
regime_1.plot(
    ax=axes[0],
    transform=ccrs.PlateCarree(),
    cmap='turbo',
    vmin=vmin,
    vmax=vmax,
    cbar_kwargs=dict(label='Mean SST (°C)', shrink=0.75)
)
axes[0].coastlines(linewidth=0.8)
axes[0].set_title('(a) Mean SST (1982–2001)', fontsize=12, weight='bold')

# ---------- Panel (b) ----------
regime_2.plot(
    ax=axes[1],
    transform=ccrs.PlateCarree(),
    cmap='turbo',
    vmin=vmin,
    vmax=vmax,
    cbar_kwargs=dict(label='Mean SST (°C)', shrink=0.75)
)
axes[1].coastlines(linewidth=0.8)
axes[1].set_title('(b) Mean SST (2001–2014)', fontsize=12, weight='bold')

# ---------- Panel (c) ----------
regime_3.plot.contourf(
    ax=axes[2],
    transform=ccrs.PlateCarree(),
    cmap='turbo',
    levels = 20,
    vmin=vmin,
    vmax=vmax,
    cbar_kwargs=dict(label='Mean SST (°C)', shrink=0.75, pad=0.05)
)
axes[2].coastlines(linewidth=0.8)
axes[2].add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
axes[2].set_title('(c) Mean SST (2014–2020)', fontsize=12, weight='bold')

plt.tight_layout()
plt.show()
#%%
diff_1 = regime_2 - regime_1
diff_2 = regime_3 - regime_1
diff_3 = regime_3 - regime_2

fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(18,8),
    dpi=300,
    subplot_kw=dict(projection=ccrs.Robinson()),
    
)

# ---------- Panel (a) ----------
diff_1.plot(
    ax=axes[0],
    transform=ccrs.PlateCarree(),
    cmap='RdBu_r',
    vmin = -1,
    vmax = 1,
    cbar_kwargs=dict(label='ΔSST (°C)', shrink=0.35)
)
axes[0].coastlines(linewidth=0.8)
axes[0].add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
axes[0].set_title(
    '(a) Regime 2 − Regime 1',
    fontsize=12,
    weight='bold'
)

# ---------- Panel (b) ----------
diff_2.plot(
    ax=axes[1],
    transform=ccrs.PlateCarree(),
    cmap='RdBu_r',
    vmin = -1,
    vmax = 1,
    cbar_kwargs=dict(label='ΔSST (°C)', shrink=0.35)
)
axes[1].coastlines(linewidth=0.8)
axes[1].add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
axes[1].set_title(
    '(b) Regime 3 − Regime 1',
    fontsize=12,
    weight='bold'
)
diff_3.plot(
    ax=axes[2],
    transform=ccrs.PlateCarree(),
    cmap='RdBu_r',
    vmin = -1,
    vmax = 1,
    cbar_kwargs=dict(label='ΔSST (°C)', shrink=0.35)
)
axes[2].coastlines(linewidth=0.8)
axes[2].add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
axes[2].set_title(
    '(a) Regime 3 − Regime 2',
    fontsize=12,
    weight='bold'
)


plt.tight_layout()
plt.show()

#%%
# Week 4:
# USE the given OISST data. 
# 1. Get the time series plot of (i) the Global mean (spatial) temperature, (ii) the spatial mean temperature in any region of interest, 
# (iii) identify a region with the highest positive trend, and (iv) a region with the negative trend. 
# 2. Plot the SST in different regimes for the above four regions, and the differences from the initial regime. The number of regimes could
#  exceed two.

#%%

fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
ax.plot(t, sst_global,color='k',linewidth=1.5,label='Global Mean SST')
ax.set_xlabel('Year', fontsize=13)
ax.set_ylabel('Sea Surface Temperature (°C)', fontsize=13)
ax.set_title('Global Mean Sea Surface Temperature (1982–2020)',fontsize=15,weight='bold')
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(fontsize=11, frameon=False)
plt.tight_layout()
plt.show()


#%%

sst_io = sst.sel(lat=slice(-30, 30), lon=slice(35, 110)).weighted(weight).mean(dim=('lat','lon'))


fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
ax.plot(
sst_io.time, sst_io,
color='k',
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
#%%

signal_io = sst_io.values.astype(float).reshape(-1,+1)
time_io = sst_io['time'].values



algo = rpt.Pelt(model ='l2').fit(signal_io)
result_io = algo.predict(pen=3)

plt.figure(figsize=(12,5))
plt.plot(time_io, signal_io[:,0], 'k', lw=1.8, label="Global Mean SST")

for bp in result_io[:-1]:
    plt.axvline(time_io[bp], color='r', ls='--', lw=1)

plt.title("Change Point Detection in  SST over the Indian Ocean", fontsize=16)
plt.xlabel("Time")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.show()
#%%


regio_1 = sst.sel(time = slice('1982-01-31','1997-11-30'),lat=slice(-30, 30), lon=slice(35, 110)).mean(dim=('time'))
regio_2 = sst.sel(time = slice('1997-11-30','2009-02-28'),lat=slice(-30, 30), lon=slice(35, 110)).mean(dim=('time'))
regio_3 = sst.sel(time = slice('2009-02-28','2020-12-31'),lat=slice(-30, 30), lon=slice(35, 110)).mean(dim=('time'))

# List of datasets and corresponding titles
regions = [regio_1, regio_2, regio_3]
titles = [
    'Indian Ocean (a) Mean SST (1982–1997)',
    'Indian Ocean (b) Mean SST (1997–2009)',
    'Indian Ocean (c) Mean SST (2009–2020)'
]

# Lock color limits for fair comparison
vmin = min(regio_1.min(), regio_2.min(), regio_3.min())
vmax = max(regio_1.max(), regio_2.max(), regio_3.max())

# Create subplots
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(18, 6),
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Loop through datasets and axes
for ax, region, title in zip(axes, regions, titles):
    region.plot.contourf(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='turbo',
        levels=20,       # You can increase or decrease levels for smoother contours
        vmin=vmin,
        vmax=vmax,
        cbar_kwargs=dict(label='Mean SST (°C)', shrink=0.75, pad=0.05)
    )
    
    # Add coastlines and land
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.set_title(title, fontsize=12, weight='bold')

plt.tight_layout()
plt.show()


#%%
sst_ls = sst.sel( lat=slice(34,40), lon=slice(286,292)).weighted(weight).mean(dim=('lat','lon'))


fig, ax1 = plt.subplots(figsize=(12, 5), dpi=300)
ax1.plot(
sst_ls.time, sst_ls,
color='k',
linewidth=1.5,
label='Mean SST'
)

ax1.set_xlabel('Year', fontsize=13)
ax1.set_ylabel('Sea Surface Temperature (°C)', fontsize=13)
ax1.set_title(
'Mean SST over the Gulf Stream region',
fontsize=15,
weight='bold'
)
ax1.grid(True, linestyle='--', alpha=0.4)
ax1.legend(fontsize=11, frameon=False)
plt.tight_layout()
plt.show()
#%%
sst_ls = sst_ls

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
plt.plot(time_ls, signal_ls[:,0], 'k', lw=1.8, label="Global Mean SST")

for bp in result_gs[:-1]:
    plt.axvline(time_ls[bp], color='r', ls='--', lw=1, alpha=0.6)
plt.title(f"Change Point Detection in SST\n54N-58N, 305E-318E", fontsize=16)

plt.xlabel("Year")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
#%%
# New time slices based on PELT change points
gs_1 = sst.sel(
    time=slice('1982-01-31', '1995-10-31'),
    lat=slice(34,40), lon=slice(286,292)
).mean(dim='time')

gs_2 = sst.sel(
    time=slice('1995-10-31', '1999-07-31'),
    lat=slice(34,40), lon=slice(286,292)
).mean(dim='time')

gs_3 = sst.sel(
    time=slice('1999-07-31', '2010-05-31'),
    lat=slice(34,40), lon=slice(286,292)
).mean(dim='time')

gs_4 = sst.sel(
    time=slice('2010-05-31', '2020-12-31'),
    lat=slice(34,40),lon=slice(286,292)
).mean(dim='time')

#%%
# List of SST arrays and titles
ls_list = [gs_1, gs_2, gs_3, gs_4]
titles = [
    'Gulf Stream (a) Mean SST (1982–1995)',
    'Gulf Stream (b) Mean SST (1995–1999)',
    'Gulf Stream (c) Mean SST (1999–2010)',
    'Gulf Stream (d) Mean SST (2010–2020)'
]

# Lock color limits
vmin = min([arr.min() for arr in ls_list])
vmax = max([arr.max() for arr in ls_list])

# Create 2x2 subplot
fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(18, 12),  # taller figure for 2 rows
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Flatten axes for easy iteration
axes_flat = axes.flatten()

# Loop over SST arrays and axes
for ax, sst_regime, title in zip(axes_flat, ls_list, titles):
    sst_regime.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='jet',
        vmin=vmin,
        vmax=vmax,
        cbar_kwargs=dict(label='SST (°C)', shrink=0.70)
    )
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.set_title(title, fontsize=12, weight='bold')

plt.tight_layout()
plt.show()

#%%

gs_diff1 = gs_2 - gs_1 
gs_diff2 = gs_3 - gs_1 
gs_diff3 = gs_4 - gs_1 

#%%


# List of SST difference arrays and titles
ls_list = [gs_diff1, gs_diff2, gs_diff3]
titles = [
    'Gulf Stream diff: region 2 (1995–1999) minus region 1 (1982–1995)',
    'Gulf Stream diff: region 3 (1999–2010) minus region 1 (1982–1995)',
    'Gulf Stream diff: region 3 (2010–2020) minus region 1 (1982–1995)'
]

# Lock color limits across all panels
vmin = min([arr.min() for arr in ls_list])
vmax = max([arr.max() for arr in ls_list])

# Create 1x3 subplot
fig, axes = plt.subplots(
    nrows=1,
    ncols=3,
    figsize=(24, 6),  # wider for horizontal layout
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Flatten axes in case it's 1D
axes_flat = axes.flatten() if hasattr(axes, 'flatten') else [axes]

# Loop over SST arrays and axes
for ax, sst_diff, title in zip(axes_flat, ls_list, titles):
    sst_diff.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='jet',
        vmin=vmin,
        vmax=vmax,
        cbar_kwargs=dict(label='SST difference (°C)', shrink=0.7)
    )
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.set_title(title, fontsize=12, weight='bold')

plt.tight_layout()
plt.show()

#%%

sst_ant = sst.sel( lat=slice(-65,-59), lon=slice(261.73,288.28)).weighted(weight).mean(dim=('lat','lon'))

fig, ax1 = plt.subplots(figsize=(12, 5), dpi=300)
ax1.plot(
sst_ant.time, sst_ant,
color='k',
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
plt.plot(time_ant, signal_ant[:,0], 'k', lw=1.8, label="Global Mean SST")

for bp in result_ant[:-1]:
    plt.axvline(time_ant[bp], color='r', ls='--', lw=1, alpha=0.6)
plt.title(f"Change Point Detection in SST\n34N-40N, 286E-292E"
, fontsize=16)

plt.xlabel("Year")
plt.ylabel("SST (°C)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
#%%
# New time slices based on PELT change points
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

#%%
# List of SST arrays and titles
ant_list = [ant_1, ant_2, ant_3, ant_4]
titles = [
    'Antartic Region (a) Mean SST (1982–1995)',
    'Antartic Region (b) Mean SST (1995–2002)',
    'Antartic Region (c) Mean SST (2002–2008)',
    'Antartic Region (d) Mean SST (2008–2020)'
]

# Lock color limits
vmin = min([arr.min() for arr in ls_list])
vmax = max([arr.max() for arr in ls_list])

# Create 2x2 subplot
fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(18, 12),  # taller figure for 2 rows
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Flatten axes for easy iteration
axes_flat = axes.flatten()

# Loop over SST arrays and axes
for ax, sst_regime, title in zip(axes_flat, ant_list, titles):
    sst_regime.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='jet',
        vmin=vmin,
        vmax=vmax,
        cbar_kwargs=dict(label='SST (°C)', shrink=0.50)
    )
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.set_title(title, fontsize=12, weight='bold')

plt.tight_layout()
plt.show()

#%%

ant_diff1 = ant_2 - ant_1 
ant_diff2 = ant_3 - ant_1 
ant_diff3 = ant_4 - ant_1 


#%% SST differences list and titles
gs_list = [ant_diff1, ant_diff2, ant_diff3]
titles = [
    'Antarctic Region diff: region 2 (1995–2002) minus region 1 (1982–1995)',
    'Antarctic Region diff: region 3 (2002–2008) minus region 1 (1982–1995)',
    'Antarctic Region diff: region 4 (2008–2020) minus region 1 (1982–1995)'
]

# Lock color limits for fair comparison
vmin = min([arr.min() for arr in gs_list])
vmax = max([arr.max() for arr in gs_list])

fig, axes = plt.subplots(
    nrows=3,
    ncols=1,
    figsize=(10, 18),  # taller figure for 3 rows
    dpi=300,
    subplot_kw=dict(projection=ccrs.PlateCarree())
)

# Flatten axes if needed
axes_flat = axes.flatten() if hasattr(axes, 'flatten') else [axes]

# Loop over SST differences and axes
for ax, sst_diff, title in zip(axes_flat, gs_list, titles):
    sst_diff.plot.contourf(   # use contourf for smoother diff visualization
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='coolwarm',      # good for differences (positive/negative)
        levels=21,            # more levels for smoother gradient
        vmin=vmin,
        vmax=vmax,
        cbar_kwargs=dict(label='SST difference (°C)', shrink=0.7, pad=0.05)
    )
    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.set_title(title, fontsize=12, weight='bold')

plt.tight_layout()
plt.show()

#%%
