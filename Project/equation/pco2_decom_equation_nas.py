#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 21:51:33 2026

@author: bobco-08
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# --- Open datasets ---
dic  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nas/pCO2_dic_nas.nc')
alk  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nas/pCO2_alk_nas.nc')
sss  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nas/pCO2_sss_nas.nc')
sst  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nas/pCO2_temp_nas.nc')
pco2_ref = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/nas/pCO2_ref_nas.nc')

# --- Rename variables MANUALLY ---
dic = dic.rename({list(dic.data_vars)[0]: "pco2_dic_nas"})
alk = alk.rename({list(alk.data_vars)[0]: "pco2_alk_nas"})
sss = sss.rename({list(sss.data_vars)[0]: "pco2_sss_nas"})
sst = sst.rename({list(sst.data_vars)[0]: "pco2_temp_nas"})
pco2_ref = pco2_ref.rename({list(pco2_ref.data_vars)[0]: "pco2_ref_nas"})

# --- Extract DataArrays (IMPORTANT STEP) ---
pco2_dic_nas  = dic["pco2_dic_nas"]
pco2_alk_nas  = alk["pco2_alk_nas"]
pco2_sss_nas  = sss["pco2_sss_nas"]
pco2_temp_nas = sst["pco2_temp_nas"]
pco2_ref_nas  = pco2_ref["pco2_ref_nas"]


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



#%% nas

dic_nas = dic.sel(lat = slice(22.5,28),lon = slice(56,70))

alk_nas = alk.sel(lat = slice(22.5,28),lon = slice(56,70))

sst_nas = sst.sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01'))

sss_nas = sss.sel(lat = slice(22.5,28),lon = slice(56,70),time = slice('1994-01-01','2024-12-01'))

#%%
sst_nas = sst_nas.interp(
    lat=dic_nas.lat,
    lon=dic_nas.lon
)

sss_nas = sss_nas.interp(
    lat=dic_nas.lat,
    lon=dic_nas.lon
)

alk_nas = alk_nas.interp(
    lat=dic_nas.lat,
    lon=dic_nas.lon
)

#%% individual drivers 

#%%  perturbation added to the drivers

dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1
#%% pco2 sensitivity of the drivers

sen_dic_nas = (pco2_dic_nas-pco2_ref_nas)/ dDIC
sen_alk_nas = (pco2_alk_nas-pco2_ref_nas)/ dALK
sen_T_nas = (pco2_temp_nas-pco2_ref_nas)/ dT
sen_s_nas = (pco2_sss_nas-pco2_ref_nas)/ dS
#%%

sen_alkv = sen_alk_nas.values

#%% time derivative of t,s,dic,alk

dT = sst_nas.diff('time')
dDIC = dic_nas.diff('time')
dALK = alk_nas.diff('time')
dS = sss_nas.diff('time')

#%%
#%% each terms
T_term = sen_T_nas * dT
dic_term = sen_dic_nas * dDIC
alk_term = sen_alk_nas *dALK
S_term =  sen_s_nas * dS

dic_termv = dic_term.values
T_termv = T_term.values
sen_T_nasv = sen_T_nas.values

#%% calculating constructed pco2


pco2_recon_nas = T_term + dic_term + alk_term + S_term  

pco2_reconv = pco2_recon_nas.values

#%%

pco2_recon = pco2_recon_nas.sel(time = slice('2010-01-01','2020-12-31'))
#%% anomalies of drivers

uatm_temp_nas =T_term.sel(time = slice('2010-01-01','2020-12-31'))
uatm_dic_nas = dic_term.sel(time = slice('2010-01-01','2020-12-31'))
uatm_alk_nas = alk_term.sel(time = slice('2010-01-01','2020-12-31'))
uatm_sal_nas = S_term.sel(time = slice('2010-01-01','2020-12-31'))

#%%

time = uatm_temp_nas['time']

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon.mean(dim =('lat','lon')),
    label='Reconstructed pCO$_2$ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_nas.mean(dim=('lat','lon')),
    label='temperature-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_nas.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_nas.mean(dim=('lat','lon')),
    label='DIC-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_nas.mean(dim=('lat','lon')),
    label='Salinity-driven pCO₂',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in (nas)')
plt.xlim(time.min(), time.max())
plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
#%%
pco2_bar = pco2_recon.mean(dim=('lat','lon'))
temp_bar = uatm_temp_nas.mean(dim=('lat','lon'))
dic_bar  = uatm_dic_nas.mean(dim=('lat','lon'))
alk_bar  = uatm_alk_nas.mean(dim=('lat','lon'))
sal_bar  = uatm_sal_nas.mean(dim=('lat','lon'))

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
plt.title('Annual accumulated pCO₂ contributions (NAS)')
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

plt.plot(
    years,
    pco2_bar.values,
    color='k',
    marker='o',
    linewidth=2,
    label='Reconstructed pCO$_2$'
)

plt.xlabel('Year')
plt.ylabel('Annual pCO₂ change (µatm yr$^{-1}$)')
plt.title('Annual pCO₂ contributions by individual drivers (NAS)')
plt.xticks(years, rotation=45)
plt.axhline(0, color='k', linewidth=0.8)
plt.legend(ncol=3, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
