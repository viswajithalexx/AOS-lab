#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 19:11:21 2026

@author: bobco-08
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 16:58:48 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# --- Open datasets ---
dic  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/esio/pCO2_dic_esio.nc')
alk  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/esio/pCO2_alk_esio.nc')
sss  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/esio/pCO2_sss_esio.nc')
sst  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/esio/pCO2_temp_esio.nc')
pco2_ref = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/esio/pCO2_ref_esio.nc')

# --- Rename variables MANUALLY ---
dic = dic.rename({list(dic.data_vars)[0]: "pco2_dic_esio"})
alk = alk.rename({list(alk.data_vars)[0]: "pco2_alk_esio"})
sss = sss.rename({list(sss.data_vars)[0]: "pco2_sss_esio"})
sst = sst.rename({list(sst.data_vars)[0]: "pco2_temp_esio"})
pco2_ref = pco2_ref.rename({list(pco2_ref.data_vars)[0]: "pco2_ref_esio"})

# --- Extract DataArrays (IMPORTANT STEP) ---
pco2_dic_esio  = dic["pco2_dic_esio"]
pco2_alk_esio  = alk["pco2_alk_esio"]
pco2_sss_esio  = sss["pco2_sss_esio"]
pco2_temp_esio = sst["pco2_temp_esio"]
pco2_ref_esio  = pco2_ref["pco2_ref_esio"]


#%%

file = '/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

data = xr.open_dataset(file)
data = data.rename({'latitude':'lat', 'longitude': 'lon'})
dic = data['tco2']
alk = data['talk']

data2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc')
data2 = data2.rename({'latitude':'lat', 'longitude': 'lon'})
sst = data2['thetao_oras'].isel(depth=0)
sss = data2['so_oras'].isel(depth=0)


#%% esio

dic_esio = dic.sel(lat = slice(-6.6,8),lon = slice(92,109))

alk_esio = alk.sel(lat = slice(-6.6,8),lon = slice(92,109))

sst_esio = sst.sel(lat = slice(-6.6,8),lon = slice(92,109), time=slice('1994-01-01','2024-12-01'))

sss_esio = sss.sel(lat = slice(-6.6,8),lon = slice(92,109), time=slice('1994-01-01','2024-12-01'))

#%%
sst_esio = sst_esio.interp(
    lat=dic_esio.lat,
    lon=dic_esio.lon
)

sss_esio = sss_esio.interp(
    lat=dic_esio.lat,
    lon=dic_esio.lon
)

alk_esio = alk_esio.interp(
    lat=dic_esio.lat,
    lon=dic_esio.lon
)

#%% individual drivers 

# data = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/oc_v2025.pCO2.nc')

# pco2_jena = data['pCO2'].sel(mtime = slice('2021-11-12','2024-12-31'),lat = slice(-31,31),lon = slice(35,110) )


# DIC = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/DIC_2deg.nc", engine="netcdf4")
# ALK = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/ALK_2deg.nc", engine="netcdf4")
# SST = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SST_2deg.nc", engine="netcdf4")
# SSS = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SSS_2deg.nc", engine="netcdf4")

#%%  perturbation added to the drivers

dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1

#%% pco2 sensitivity of the drivers

sen_dic_esio = (pco2_dic_esio - pco2_ref_esio) / dDIC
sen_alk_esio = (pco2_alk_esio - pco2_ref_esio) / dALK
sen_T_esio   = (pco2_temp_esio - pco2_ref_esio) / dT
sen_s_esio   = (pco2_sss_esio - pco2_ref_esio) / dS

#%%

sen_alkv = sen_alk_esio.values

#%% time derivative of t,s,dic,alk

dT   = sst_esio.diff('time')
dDIC = dic_esio.diff('time')
dALK = alk_esio.diff('time')
dS   = sss_esio.diff('time')

#%% each terms

T_term   = sen_T_esio * dT
dic_term = sen_dic_esio * dDIC
alk_term = sen_alk_esio * dALK
S_term   = sen_s_esio * dS

dic_termv = dic_term.values
T_termv   = T_term.values
sen_T_esiov = sen_T_esio.values

#%% calculating constructed pco2

pco2_recon_esio = T_term + dic_term + alk_term + S_term  

pco2_reconv = pco2_recon_esio.values

#%%

pco2_recon = pco2_recon_esio.sel(time =slice('1994-01-01','2024-12-01'))

#%% anomalies of drivers

uatm_temp_esio = T_term.sel(time =slice('1994-01-01','2024-12-01'))
uatm_dic_esio  = dic_term.sel(time =slice('1994-01-01','2024-12-01'))
uatm_alk_esio  = alk_term.sel(time =slice('1994-01-01','2024-12-01'))
uatm_sal_esio  = S_term.sel(time =slice('1994-01-01','2024-12-01'))

#%%

time = uatm_temp_esio['time']

plt.figure(figsize=(12, 5))

plt.plot(time, pco2_recon.mean(dim=('lat','lon')), label='Reconstructed pCO$_2$', linewidth=2)
plt.plot(time, uatm_temp_esio.mean(dim=('lat','lon')), label='temperature-driven pCO$_2$', linewidth=2)
plt.plot(time, uatm_alk_esio.mean(dim=('lat','lon')), label='Alkalinity-driven pCO$_2$', linewidth=2)
plt.plot(time, uatm_dic_esio.mean(dim=('lat','lon')), label='DIC-driven pCO$_2$', linewidth=2)
plt.plot(time, uatm_sal_esio.mean(dim=('lat','lon')), label='Salinity-driven pCO$_2$', linewidth=2)

plt.xlabel('Year')
plt.ylabel('pCO$_2$ (µatm)')
plt.title('Contribution of pCO$_2$ by different variables in Eastern South Indian Ocean')
plt.xlim(time.min(), time.max())
plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

#%%
pco2_bar = pco2_recon.mean(dim=('lat','lon'))
temp_bar = uatm_temp_esio.mean(dim=('lat','lon'))
dic_bar  = uatm_dic_esio.mean(dim=('lat','lon'))
alk_bar  = uatm_alk_esio.mean(dim=('lat','lon'))
sal_bar  = uatm_sal_esio.mean(dim=('lat','lon'))

#%%
pco2_bar = pco2_bar.resample(time='1Y').sum()
temp_bar = temp_bar.resample(time='1Y').sum()
dic_bar  = dic_bar.resample(time='1Y').sum()
alk_bar  = alk_bar.resample(time='1Y').sum()
sal_bar  = sal_bar.resample(time='1Y').sum()

time = pco2_bar['time'].dt.year

#%%
plt.figure(figsize=(12,6))

plt.bar(time, temp_bar, label='Temperature')
plt.bar(time, dic_bar,  bottom=temp_bar, label='DIC')
plt.bar(time, alk_bar,  bottom=temp_bar+dic_bar, label='Alkalinity')
plt.bar(time, sal_bar,  bottom=temp_bar+dic_bar+alk_bar, label='Salinity')

plt.plot(time, pco2_bar, color='k', marker='o', linewidth=2, label='Reconstructed pCO$_2$')

plt.xlabel('Year')
plt.ylabel('Annual pCO₂ change (µatm yr$^{-1}$)')
plt.title('Annual accumulated pCO₂ contributions (ESIO)')
plt.legend(ncol=2, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()

#%%


years = pco2_bar['time'].dt.year.values.astype(float)

width = 0.15

plt.figure(figsize=(14,6))

plt.bar(years - 1.5*width, temp_bar.values, width, label='Temperature')
plt.bar(years - 0.5*width, dic_bar.values,  width, label='DIC')
plt.bar(years + 0.5*width, alk_bar.values,  width, label='Alkalinity')
plt.bar(years + 1.5*width, sal_bar.values,  width, label='Salinity')

plt.plot(years, pco2_bar.values, color='k', marker='o', linewidth=2, label='Reconstructed pCO$_2$')

plt.xlabel('Year')
plt.ylabel('Annual pCO₂ change (µatm yr$^{-1}$)')
plt.title('Annual pCO₂ contributions by individual drivers (ESIO)')
plt.xticks(years, rotation=45)
plt.axhline(0, color='k', linewidth=0.8)
plt.legend(ncol=3, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
