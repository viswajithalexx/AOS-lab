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
file_path = ""
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
