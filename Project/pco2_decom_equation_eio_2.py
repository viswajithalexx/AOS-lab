#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 10:16:06 2026

@author: bobco-08
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


file = '/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

data = xr.open_dataset(file)
data = data.rename({'latitude':'lat', 'longitude': 'lon'})
lat = data['lat']
lon =data['lon']
time = data['time']
dic = data['tco2'].sel(lat = slice(-6.5,5),lon = slice(49,92))
alk = data['talk'].sel(lat = slice(-6.5,5),lon = slice(49,92))
pco2 = data['spco2'].sel(lat = slice(-6.5,5),lon = slice(49,92))

data2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc')
data2 = data2.rename({'latitude':'lat', 'longitude': 'lon'})
sst = data2['thetao_oras'].isel(depth=0)
sss = data2['so_oras'].isel(depth=0)
sst = sst.sel(lat = slice(-6.5,5),lon = slice(49,92))
sss = sss.sel(lat = slice(-6.5,5),lon = slice(49,92))
#%% pco2-term


pco2_reg = pco2.mean(dim=('lat', 'lon'))

pco2_annual = pco2_reg.groupby('time.year').mean(dim = 'time')



# Select JJAS months
pco2_jja = pco2_reg.sel(time=pco2_reg['time'].dt.month.isin([6,7,8]))

# Group by year
pco2_jja_grouped = pco2_jja.groupby('time.year')
pco2_jja_y = pco2_jja_grouped.mean()


pco2_mam = pco2_reg.sel(time=pco2_reg['time'].dt.month.isin([2,3,4]))

# Group by year
pco2_mam_grouped = pco2_mam.groupby('time.year')
pco2_mam_y = pco2_mam_grouped.mean()

# DJF using quarterly resample starting in December
# pco2_mam = pco2_reg.resample(time='QS-JAN').mean(dim='time')  # Dec–Jan–Feb mean
# pco2_mam = pco2_mam.isel(time=slice(1,-1))                 # Skip first incomplete DJF if needed
# pco2_mam = pco2_mam.sel(time=pco2_jfm.time.dt.month == )

# pco2_mam_y = pco2_mam.groupby('time.year').mean()

# JJAS – DJF difference
delta_pco2 = (pco2_mam_y - pco2_jja_y)


#%% t-term

# Average SST over latitude and longitude
sst_reg = sst.mean(dim=['lat', 'lon'])

# Select JJAS months
sst_jja = sst_reg.sel(time=sst_reg['time'].dt.month.isin([6,7,8]))

# Group by year
sst_jja_grouped = sst_jja.groupby('time.year')

sst_jja_y = sst_jja_grouped.mean()


sst_mam = sst_reg.sel(time=sst_reg['time'].dt.month.isin([2,3,4]))

# Group by year
sst_mam_grouped = sst_mam.groupby('time.year')

sst_mam_y = sst_mam_grouped.mean()


# sst_jfm = sst_reg.resample(time='QS-JAN').mean(dim='time')  # Quarterly starting in Dec
# sst_jfm = sst_jfm.isel(time=slice(1, -1))                 # Skip first incomplete DJF if needed
# sst_jfm = sst_jfm.sel(time=sst_jfm.time.dt.month == 1)  

# sst_jfm_y = sst_mam.groupby('time.year').mean()

delta_T = (sst_mam_y-sst_jja_y)



T_term = 2 * pco2_annual * (np.exp(0.0423 * (delta_T / 2)) - 1)


#%% dic_term

# Average over lat/lon
dic_reg = dic.mean(dim=['lat', 'lon'])

# Annual mean
dic_annual = dic_reg.groupby('time.year').mean(dim='time')

# Select JJAS months (June, July, August, September if needed)
dic_jja = dic_reg.sel(time=dic_reg['time'].dt.month.isin([6,7,8]))

# Group by year and take mean over time
dic_jja_grouped = dic_jja.groupby('time.year')

dic_jja_y = dic_jja_grouped.mean()


dic_mam = dic_reg.sel(time=dic_reg['time'].dt.month.isin([2,3,4]))

# Group by year and take mean over time
dic_mam_grouped = dic_mam.groupby('time.year')

dic_mam_y = dic_mam_grouped.mean()


# dic_mam = dic_reg.resample(time='QS-JAN').mean(dim='time')  # Quarterly starting in Dec
# dic_mam  = dic_mam.isel(time=slice(1,-1))                 # Skip first incomplete DJF if needed
# dic_jfm  = dic_jfm.sel(time=dic_jfm.time.dt.month == 1)  

# dic_jfm_y = dic_jfm.groupby('time.year').mean()

delta_dic = (dic_mam_y - dic_jja_y)

dic_term = (9.5 * (pco2_annual/dic_annual))* delta_dic


#%%

#%% alk_term

# Average over lat/lon
alk_reg = alk.mean(dim=['lat', 'lon'])

# Annual mean
alk_annual = alk_reg.groupby('time.year').mean(dim='time')

# Select JJAS months (June, July, August, September)
alk_jja = alk_reg.sel(time=alk_reg['time'].dt.month.isin([6,7,8]))

# Group by year and take mean over time
alk_jja_grouped = alk_jja.groupby('time.year')
alk_jja_y = alk_jja_grouped.mean()

alk_mam = alk_reg.sel(time=alk_reg['time'].dt.month.isin([2,3,4]))

# Group by year and take mean over time
alk_mam_grouped = alk_mam.groupby('time.year')
alk_mam_y = alk_mam_grouped.mean()



# DJF using quarterly resample starting in December
# alk_jfm = alk_reg.resample(time='QS-JAN').mean(dim='time')   # Dec–Jan–Feb mean
# alk_jfm = alk_jfm.isel(time=slice(1,-1))                  # Skip first incomplete DJF
# alk_jfm = alk_jfm.sel(time=alk_jfm.time.dt.month == 1)      # Keep DJF labeled by Dec

# alk_jfm_y = alk_jfm.groupby('time.year').mean()

# JJAS – DJF difference
delta_alk = (alk_mam_y - alk_jja_y)


alk_term = (- 8.9 * (pco2_annual/alk_annual))* delta_alk

#%% sss_term

#%% sal_term

# Average over lat/lon
sal_reg = sss.mean(dim=['lat', 'lon'])

# Annual mean
sal_annual = sal_reg.groupby('time.year').mean(dim='time')

# Select JJAS months (June, July, August, September)
sal_jja = sal_reg.sel(time=sal_reg['time'].dt.month.isin([6,7,8]))

# Group by year and take mean over time
sal_jja_grouped = sal_jja.groupby('time.year')
sal_jja_y = sal_jja_grouped.mean()

sal_mam = sal_reg.sel(time=sal_reg['time'].dt.month.isin([2,3,4]))

# Group by year and take mean over time
sal_mam_grouped = sal_mam.groupby('time.year')
sal_mam_y = sal_mam_grouped.mean()




# DJF using quarterly resample starting in December
# sal_jfm = sal_reg.resample(time='QS-JAN').mean(dim='time')   # Dec–Jan–Feb mean
# sal_jfm = sal_jfm.isel(time=slice(1,-1))                  # Skip first incomplete DJF
# sal_jfm = sal_jfm.sel(time=sal_jfm.time.dt.month == 1)      # Keep DJF labeled by Dec

# sal_jfm_y = sal_jfm.groupby('time.year').mean()

# JJAS – DJF difference
delta_sal = (sal_mam_y - sal_jja_y)


sal_term = (0.026* pco2_annual)* delta_sal

#%%
 
delta_pco2_cal = T_term + dic_term + alk_term + sal_term

#%%

#%% Plot observed vs calculated seasonal ΔpCO2 (1994–2024)

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# Align observed and calculated ΔpCO2 in time
delta_pco2, delta_pco2_cal = xr.align(delta_pco2, delta_pco2_cal)

# Extract years and values
years = delta_pco2['year'].values
dpco2_obs = delta_pco2.values
dpco2_cal = delta_pco2_cal.values

# Bar plot settings
x = np.arange(len(years))
width = 0.35

plt.figure(figsize=(18, 6))

plt.bar(
    x - width/2,
    dpco2_obs,
    width,
    label='Observed ΔpCO$_2$',
    color='gray'
)

plt.bar(
    x + width/2,
    dpco2_cal,
    width,
    label='Calculated ΔpCO$_2$',
    color='steelblue'
)

plt.axhline(0, color='k', linewidth=0.8)

plt.xticks(x[::2], years[::2], rotation=45)
plt.xlabel('Year')
plt.ylabel('ΔpCO$_2$ (µatm)')
plt.title('Observed vs Calculated Year-wise pCO$_2$ Amplitude in EIO (1994–2024)')
plt.legend(frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig("/home/bobco-08/24cl05012/CO2/plot/equation_pco2/valid_plots_eqn/eio/eio_obsvscal",dpi=500,bbox_inches="tight"
)
plt.show()


#%% Line plot of observed, calculated seasonal ΔpCO2 and residual

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# Align observed and calculated ΔpCO2
delta_pco2, delta_pco2_cal = xr.align(delta_pco2, delta_pco2_cal)

# Extract years and values
years = delta_pco2['year'].values
dpco2_obs = delta_pco2.values
dpco2_cal = delta_pco2_cal.values

# Residual (misfit)
residual = dpco2_cal - dpco2_obs

plt.figure(figsize=(16, 6))

# Observed and calculated ΔpCO2
plt.plot(
    years,
    dpco2_obs,
    marker='o',
    linewidth=2,
    label='Observed ΔpCO$_2$'
)

plt.plot(
    years,
    dpco2_cal,
    marker='s',
    linewidth=2,
    label='Calculated ΔpCO$_2$'
)

# Residual
plt.plot(
    years,
    residual,
    marker='^',
    linestyle='--',
    linewidth=1.8,
    label='Residual (Calc − Obs)'
)

plt.axhline(0, color='k', linewidth=0.8)

plt.xlabel('Year')
plt.ylabel('ΔpCO$_2$ (µatm)')
plt.title('Observed vs Calculated pCO$_2$ and Residualin EIO (1994–2024)')
plt.legend(frameon=False, ncol=3)
plt.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig(
    "/home/bobco-08/24cl05012/CO2/plot/equation_pco2/valid_plots_eqn/eio/eio_obs_cal_res",
    dpi=500,
    bbox_inches="tight"
)
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# ---- Align all data in time ----
delta_pco2, T_term, dic_term, alk_term, sal_term, delta_pco2_cal = xr.align(
    delta_pco2, T_term, dic_term, alk_term, sal_term, delta_pco2_cal
)

years = delta_pco2['year'].values

obs = delta_pco2.values
T   = T_term.values
DIC = dic_term.values
ALK = alk_term.values
SAL = sal_term.values
TOT = delta_pco2_cal.values

x = np.arange(len(years))

# width of each bar
w = 0.12

plt.figure(figsize=(24, 7))

plt.bar(x - 2.5*w, T,   w, label='Temperature')
plt.bar(x - 1.5*w, DIC, w, label='DIC')
plt.bar(x - 0.5*w, ALK, w, label='Alkalinity')
plt.bar(x + 0.5*w, SAL, w, label='Salinity')
plt.bar(x + 1.5*w, obs, w, label='Observed ΔpCO$_2$', color='k')
plt.bar(x + 2.5*w, TOT, w, label='Reconstructed ΔpCO$_2$', color='gray')

plt.axhline(0, color='k', linewidth=0.8)

plt.xticks(x[::2], years[::2], rotation=45)
plt.ylabel('ΔpCO$_2$ (µatm)')
plt.title('pCO$_2$ Decomposition (FMA − JJA)\nEIO')
plt.legend(ncol=6, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig(
    "/home/bobco-08/24cl05012/CO2/plot/equation_pco2/valid_plots_eqn/eio/eio_drivers",
    dpi=500,
    bbox_inches="tight"
)

plt.show()
#%%
import matplotlib.pyplot as plt
import numpy as np

pco2_clim = pco2_reg.groupby('time.month').mean()

months = pco2_clim['month'].values
values = pco2_clim.values

# --- Find peak and minimum month indices ---
imax = int(values.argmax())
imin = int(values.argmin())

# --- Indices for ±1 month (with circular wrap) ---
def wrap_idx(i):
    return i % 12

peak_window = [wrap_idx(imax-1), imax, wrap_idx(imax+1)]
min_window  = [wrap_idx(imin-1), imin, wrap_idx(imin+1)]

plt.figure(figsize=(10, 5))

# Main line
plt.plot(months, values, marker='o', linewidth=2, color='k')

# Highlight peak ±1 months (red)
plt.scatter(months[peak_window], values[peak_window],
            color='red', s=90, zorder=5, label='Peak ±1 month')

# Highlight minimum ±1 months (blue)
plt.scatter(months[min_window], values[min_window],
            color='blue', s=90, zorder=5, label='Min ±1 month')

# Optional labels
for i in peak_window:
    plt.text(months[i], values[i]-1, f'{values[i]:.0f}', ha='center', color='red')

for i in min_window:
    plt.text(months[i], values[i]+1, f'{values[i]:.0f}', ha='center', color='blue')

plt.xticks(range(1, 13),
           ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])

plt.xlabel('Month')
plt.ylabel('pCO$_2$ (µatm)')
plt.title('Monthly Climatology of Surface pCO$_2$ in EIO\n(Peak and Trough ±1 Month Highlighted)')

plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

plt.savefig(
    "/home/bobco-08/24cl05012/CO2/plot/equation_pco2/valid_plots_eqn/eio/mon_clim_eio_peak_trough.png",
    dpi=500,
    bbox_inches="tight"
)
plt.show()

#%%
plt.figure(figsize=(18,7),dpi = 500)

plt.plot(years, dpco2_obs, marker='o', linewidth=2.5, label='Observed ΔpCO$_2$', color='k')
plt.plot(years, dpco2_cal, marker='s', linewidth=2.5, label='Reconstructed ΔpCO$_2$', color='gray')

plt.plot(years, T_term,  marker='^', label='Temp-driven')
plt.plot(years, dic_term, marker='d',  label='DIC-driven')
plt.plot(years, alk_term, marker='v',  label='ALK-driven')
plt.plot(years, sal_term, marker='x',  label='Salinity-driven')

plt.axhline(0, color='k', linewidth=0.8)

plt.xlabel('Year')
plt.ylabel('ΔpCO$_2$ (µatm)')
plt.title('Annual Time Series of pCO$_2$ and Drivers (FMA − JJA)\nEIO (1994–2024)')
plt.legend(ncol=3, frameon=False)
plt.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig(
    "/home/bobco-08/24cl05012/CO2/plot/equation_pco2/valid_plots_eqn/eio/eio_annual_timeseries_drivers.png",
    dpi=500, bbox_inches="tight"
)
plt.show()
