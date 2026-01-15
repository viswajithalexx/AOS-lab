#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 22:09:17 2026

@author: bobco-08
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# --- Open datasets ---
dic  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/eio/pCO2_dic_eio.nc')
alk  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/eio/pCO2_alk_eio.nc')
sss  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/eio/pCO2_sss_eio.nc')
sst  = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/eio/pCO2_temp_eio.nc')
pco2_ref = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/eio/pCO2_ref_eio.nc')

# --- Rename variables MANUALLY ---
dic = dic.rename({list(dic.data_vars)[0]: "pco2_dic_eio"})
alk = alk.rename({list(alk.data_vars)[0]: "pco2_alk_eio"})
sss = sss.rename({list(sss.data_vars)[0]: "pco2_sss_eio"})
sst = sst.rename({list(sst.data_vars)[0]: "pco2_temp_eio"})
pco2_ref = pco2_ref.rename({list(pco2_ref.data_vars)[0]: "pco2_ref_eio"})

# --- Extract DataArrays (IMPORTANT STEP) ---
pco2_dic_eio  = dic["pco2_dic_eio"]
pco2_alk_eio  = alk["pco2_alk_eio"]
pco2_sss_eio  = sss["pco2_sss_eio"]
pco2_temp_eio = sst["pco2_temp_eio"]
pco2_ref_eio  = pco2_ref["pco2_ref_eio"]


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



#%% eio

dic_eio = dic.sel(lat = slice(-6.5,5),lon = slice(49,92))

alk_eio = alk.sel(lat = slice(-6.5,5),lon = slice(49,92))

sst_eio = sst.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01'))

sss_eio = sss.sel(lat = slice(-6.5,5),lon = slice(49,92),time = slice('1994-01-01','2024-12-01'))

#%%
sst_eio = sst_eio.interp(
    lat=dic_eio.lat,
    lon=dic_eio.lon
)

sss_eio = sss_eio.interp(
    lat=dic_eio.lat,
    lon=dic_eio.lon
)

alk_eio = alk_eio.interp(
    lat=dic_eio.lat,
    lon=dic_eio.lon
)

#%% individual drivers 

#%%  perturbation added to the drivers

dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1
#%% pco2 sensitivity of the drivers

sen_dic_eio = (pco2_dic_eio-pco2_ref_eio)/ dDIC
sen_alk_eio = (pco2_alk_eio-pco2_ref_eio)/ dALK
sen_T_eio = (pco2_temp_eio-pco2_ref_eio)/ dT
sen_s_eio = (pco2_sss_eio-pco2_ref_eio)/ dS
#%%

sen_alkv = sen_alk_eio.values

#%% time derivative of t,s,dic,alk

dT = sst_eio.diff('time')
dDIC = dic_eio.diff('time')
dALK = alk_eio.diff('time')
dS = sss_eio.diff('time')

#%%
#%% each terms
T_term = sen_T_eio * dT
dic_term = sen_dic_eio * dDIC
alk_term = sen_alk_eio *dALK
S_term =  sen_s_eio * dS

dic_termv = dic_term.values
T_termv = T_term.values
sen_T_eiov = sen_T_eio.values

#%% calculating constructed pco2


pco2_recon_eio = T_term + dic_term + alk_term + S_term  

pco2_reconv = pco2_recon_eio.values

#%%

pco2_recon = pco2_recon_eio.sel(time = slice('2010-01-01','2020-12-31'))
#%% anomalies of drivers

uatm_temp_eio =T_term.sel(time = slice('2010-01-01','2020-12-31'))
uatm_dic_eio = dic_term.sel(time = slice('2010-01-01','2020-12-31'))
uatm_alk_eio = alk_term.sel(time = slice('2010-01-01','2020-12-31'))
uatm_sal_eio = S_term.sel(time = slice('2010-01-01','2020-12-31'))

#%%

time = uatm_temp_eio['time']

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon.mean(dim =('lat','lon')),
    label='Reconstructed pCO$_2$ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_eio.mean(dim=('lat','lon')),
    label='temperature-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_eio.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_eio.mean(dim=('lat','lon')),
    label='DIC-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_eio.mean(dim=('lat','lon')),
    label='Salinity-driven pCO₂',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in (eio)')
plt.xlim(time.min(), time.max())
plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
#%%
pco2_bar = pco2_recon.mean(dim=('lat','lon'))
temp_bar = uatm_temp_eio.mean(dim=('lat','lon'))
dic_bar  = uatm_dic_eio.mean(dim=('lat','lon'))
alk_bar  = uatm_alk_eio.mean(dim=('lat','lon'))
sal_bar  = uatm_sal_eio.mean(dim=('lat','lon'))

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
plt.title('Annual accumulated pCO₂ contributions (EIO)')
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
plt.title('Annual pCO₂ contributions by individual drivers (EIO)')
plt.xticks(years, rotation=45)
plt.axhline(0, color='k', linewidth=0.8)
plt.legend(ncol=3, frameon=False)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
