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
dic  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nwio/pCO2_dic_nwio.nc')
alk  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nwio/pCO2_alk_nwio.nc')
sss  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nwio/pCO2_sss_nwio.nc')
sst  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nwio/pCO2_temp_nwio.nc')
pco2_ref = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nwio/pCO2_ref_nwio.nc')

# --- Rename variables MANUALLY ---
dic = dic.rename({list(dic.data_vars)[0]: "pco2_dic_nwio"})
alk = alk.rename({list(alk.data_vars)[0]: "pco2_alk_nwio"})
sss = sss.rename({list(sss.data_vars)[0]: "pco2_sss_nwio"})
sst = sst.rename({list(sst.data_vars)[0]: "pco2_temp_nwio"})
pco2_ref = pco2_ref.rename({list(pco2_ref.data_vars)[0]: "pco2_ref_nwio"})

# --- Extract DataArrays (IMPORTANT STEP) ---
pco2_dic_nwio  = dic["pco2_dic_nwio"]
pco2_alk_nwio  = alk["pco2_alk_nwio"]
pco2_sss_nwio  = sss["pco2_sss_nwio"]
pco2_temp_nwio = sst["pco2_temp_nwio"]
pco2_ref_nwio  = pco2_ref["pco2_ref_nwio"]


#%%

file = '/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/cmems_obs-mob_glo_bgc-car_my_irr-i_multi-vars_35.12E-119.88E_29.88S-29.88N_1994-01-01-2024-12-01.nc'

data = xr.open_dataset(file)
data = data.rename({'latitude':'lat', 'longitude': 'lon'})
dic = data['tco2']
alk = data['talk']

data2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/cmems_mod_glo_phy-all_my_0.25deg_P1M-m_multi-vars_35.00E-120.00E_30.00S-30.00N_0.51m_1993-01-01-2024-12-01.nc')
data2 = data2.rename({'latitude':'lat', 'longitude': 'lon'})
sst = data2['thetao_oras'].isel(depth =0)
sss = data2['so_oras'].isel(depth =0)



#%% nwio

dic_nwio = dic.sel(lat = slice(5,22.5),lon = slice(45,65))

alk_nwio = alk.sel(lat = slice(5,22.5),lon = slice(45,65))

sst_nwio = sst.sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01'))

sss_nwio = sss.sel(lat = slice(5,22.5),lon = slice(45,65),time = slice('1994-01-01','2024-12-01'))

#%%
sst_nwio = sst_nwio.interp(
    lat=dic_nwio.lat,
    lon=dic_nwio.lon
)

sss_nwio = sss_nwio.interp(
    lat=dic_nwio.lat,
    lon=dic_nwio.lon
)

alk_nwio = alk_nwio.interp(
    lat=dic_nwio.lat,
    lon=dic_nwio.lon
)

#%% individual drivers 

# data = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/oc_v2025.pCO2.nc')

# pco2_jena = data['pCO2'].sel(mtime = slice('2021-11-12','2024-12-31'),lat = slice(-31,31),lon = slice(35,110) )


# DIC = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/DIC_2deg.nc", engine="netcdf4")
# ALK = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/ALK_2deg.nc", engine="netcdf4")
# SST = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SST_2deg.nc", engine="netcdf4")
# SSS = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SSS_2deg.nc", engine="netcdf4")

# #variables

# dic = DIC['DIC']
# alk = ALK['ALK']
# sss = SSS['SSS']
# sst = SST['SST']
#%%  perturbation added to the drivers

dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1
#%% pco2 sensitivity of the drivers

sen_dic_nwio = (pco2_dic_nwio-pco2_ref_nwio)/ dDIC
sen_alk_nwio = (pco2_alk_nwio-pco2_ref_nwio)/ dALK
sen_T_nwio = (pco2_temp_nwio-pco2_ref_nwio)/ dT
sen_s_nwio = (pco2_sss_nwio-pco2_ref_nwio)/ dS
#%%

sen_alkv = sen_alk_nwio.values

#%% time derivative of t,s,dic,alk

dT = sst_nwio.diff('time')
dDIC = dic_nwio.diff('time')
dALK = alk_nwio.diff('time')
dS = sss_nwio.diff('time')

#%%
# dT_delta = sst_nwio.groupby('time.month') - sst_nwio.groupby('time.month').mean('time')
# ddic_delta = dic_nwio.groupby('time.month') - dic_nwio.groupby('time.month').mean('time')
# dalk_delta = alk_nwio.groupby('time.month') - alk_nwio.groupby('time.month').mean('time')
# dS_delta = sss_nwio.groupby('time.month') - sss_nwio.groupby('time.month').mean('time')
#%% each terms
T_term = sen_T_nwio * dT
dic_term = sen_dic_nwio * dDIC
alk_term = sen_alk_nwio *dALK
S_term =  sen_s_nwio * dS

dic_termv = dic_term.values
T_termv = T_term.values
sen_T_nwiov = sen_T_nwio.values

#%% calculating constructed pco2


pco2_recon_nwio = T_term + dic_term + alk_term + S_term  

pco2_reconv = pco2_recon_nwio.values

# pco2_recon_nwio = pco2_recon_nwio.sel(time = slice('2020-01-01','2020-12-31'))

#%%

pco2_recon = pco2_recon_nwio
#%% anomalies of drivers

uatm_temp_nwio =T_term
uatm_dic_nwio = dic_term
uatm_alk_nwio = alk_term
uatm_sal_nwio = S_term

#%%

time = uatm_temp_nwio['time']

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon.mean(dim =('lat','lon')),
    label='Reconstructed pCO$_2$ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_nwio.mean(dim=('lat','lon')),
    label='temperature-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_nwio.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_nwio.mean(dim=('lat','lon')),
    label='DIC-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_nwio.mean(dim=('lat','lon')),
    label='Salinity-driven pCO₂',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in Northwestern Indian Ocean')
plt.xlim(time.min(), time.max())
plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
#%%
pco2_bar = pco2_recon.mean(dim=('lat','lon'))
temp_bar = uatm_temp_nwio.mean(dim=('lat','lon'))
dic_bar  = uatm_dic_nwio.mean(dim=('lat','lon'))
alk_bar  = uatm_alk_nwio.mean(dim=('lat','lon'))
sal_bar  = uatm_sal_nwio.mean(dim=('lat','lon'))

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

plt.plot(
    time,
    pco2_bar,
    color='k',
    marker='o',
    linewidth=2,
    label='Reconstructed pCO$_2$'
)

plt.xlabel('Year')
plt.ylabel('Annual pCO₂ change (µatm yr$^{-1}$)')
plt.title('Annual accumulated pCO₂ contributions (NWIO)')
plt.legend(ncol=2, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
#%%

# --- Annual SUM ---
pco2_bar = pco2_bar.resample(time='1Y').sum()
temp_bar = temp_bar.resample(time='1Y').sum()
dic_bar  = dic_bar.resample(time='1Y').sum()
alk_bar  = alk_bar.resample(time='1Y').sum()
sal_bar  = sal_bar.resample(time='1Y').sum()

years = pco2_bar['time'].dt.year.values.astype(float)

# --- Bar settings (fraction of a year) ---
width = 0.15

plt.figure(figsize=(14,6))

plt.bar(years - 1.5*width, temp_bar.values, width, label='Temperature')
plt.bar(years - 0.5*width, dic_bar.values,  width, label='DIC')
plt.bar(years + 0.5*width, alk_bar.values,  width, label='Alkalinity')
plt.bar(years + 1.5*width, sal_bar.values,  width, label='Salinity')

# --- Reconstructed pCO2 line ---
plt.plot(
    years,
    pco2_bar.values,
    color='k',
    marker='o',
    linewidth=2,
    label='Reconstructed pCO$_2$'
)

# --- Axis formatting ---
plt.xlabel('Year')
plt.ylabel('Annual pCO₂ change (µatm yr$^{-1}$)')
plt.title('Annual pCO₂ contributions by individual drivers (NWIO)')
plt.xticks(years, rotation=45)
plt.axhline(0, color='k', linewidth=0.8)
plt.legend(ncol=3, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
