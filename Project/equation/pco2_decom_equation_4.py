#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 19:42:37 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

dic = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/pCO2_dic.nc')
alk = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/pCO2_alk.nc')
sss = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/pCO2_S.nc')
sst = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/pCO2_T.nc')
pco2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/pCO2_ref.nc')


pco2_dic = dic['pCO2_dic'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
pco2_alk = alk['pCO2_alk'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
pco2_t = sst['pCO2_T'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
pco2_s = sss['pCO2_S'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
pco2_ref = pco2['pCO2'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))

pco2_dic_m = pco2_dic.resample(time='1MS').mean()
pco2_alk_m = pco2_alk.resample(time='1MS').mean()
pco2_t_m   = pco2_t.resample(time='1MS').mean()
pco2_s_m   = pco2_s.resample(time='1MS').mean()
pco2_ref_m = pco2_ref.resample(time='1MS').mean()


data = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/oc_v2025.pCO2.nc')

pco2_jena = data['pCO2'].sel(mtime = slice('2021-11-12','2024-12-31'),lat = slice(-31,31),lon = slice(35,110) )


DIC = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/DIC_2deg.nc", engine="netcdf4")
ALK = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/ALK_2deg.nc", engine="netcdf4")
SST = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SST_2deg.nc", engine="netcdf4")
SSS = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SSS_2deg.nc", engine="netcdf4")

#variables

dic = DIC['DIC'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
alk = ALK['ALK'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
sss = SSS['SSS'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
sst = SST['SST'].sel(time = slice('2022-01-01','2022-12-31'),lat = slice(5,22.5),lon = slice(45,65))
#%%

import PyCO2SYS as pyco2
import gsw as gsw

# Broadcast lon/lat
lon2d, lat2d = xr.broadcast(sss.lon, sss.lat)

# Absolute Salinity
SA = gsw.SA_from_SP(sss, 0, lon2d, lat2d)

CT = gsw.CT_from_pt(SA,sst)
# Density (kg/m3)
rho = gsw.rho(SA,CT, 0)

# Convert mol/m3 → µmol/kg
dic_umolkg = dic * 1e6 / rho
alk_umolkg = alk * 1e6 / rho

dic_ukg = dic_umolkg.where(dic_umolkg > 0)
alk_ukg = alk_umolkg.where(alk_umolkg > 0)

dic_m = dic_ukg.resample(time='1MS').mean()
alk_m = alk_ukg.resample(time='1MS').mean()
sss_m = sss.resample(time='1MS').mean()
sst_m = sst.resample(time='1MS').mean()

#%%  perturbation

dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1
#%% sensitivity of variables

sen_dic = (pco2_dic_m-pco2_ref_m)/ dDIC
sen_alk = (pco2_alk_m-pco2_ref_m)/ dALK
sen_T = (pco2_t_m-pco2_ref_m)/ dT
sen_s = (pco2_s_m-pco2_ref_m)/ dS

#%% derivative of terms

dT_dt = sst_m.diff('time')
ddic_dt = dic_m.diff('time')
dalk_dt = alk_m.diff('time')
dS_dt = sss_m.diff('time')

#%% each terms
T_term = sen_T * dT_dt
dic_term = sen_dic * ddic_dt
alk_term = sen_alk * dalk_dt
S_term =  sen_s * dS_dt

#%% calculating constructed pco2



dpco2_dt_cal = T_term + dic_term + alk_term + S_term  #tendency

pco2_recon = dpco2_dt_cal



#%%

# pco2_recon = pco2_cal.isel(time=0)+ dpdt.cumsum(dim='time')*dt

#%% anomalies of drivers

uatm_temp =T_term
uatm_dic = dic_term
uatm_alk = alk_term
uatm_sal = S_term

#%%

time = uatm_temp['time']

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    dpco2_dt_cal.mean(dim =('lat','lon')),
    label='Reconstructed pCO$_2$ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp.mean(dim=('lat','lon')),
    label='temperature-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic.mean(dim=('lat','lon')),
    label='DIC-driven pCO₂ anomaly',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal.mean(dim=('lat','lon')),
    label='Salinity-driven pCO₂',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in NWIO')
plt.xlim(time.min(), time.max())
plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
#%% nwio
dpco2_dt_cal_nwio = dpco2_dt_cal.sel(lat = slice(5,22.5),lon = slice(45,65))

T_term_nwio = T_term.sel(lat = slice(5,22.5),lon = slice(45,65))

dic_term_nwio = dic_term.sel(lat = slice(5,22.5),lon = slice(45,65))

alk_term_nwio = alk_term.sel(lat = slice(5,22.5),lon = slice(45,65))

S_term_nwio = S_term.sel(lat = slice(5,22.5),lon = slice(45,65))
#%% nas
dpco2_dt_cal_nas = dpco2_dt_cal.sel(lat = slice(22.5,28),lon = slice(56,70))

T_term_nas = T_term.sel(lat = slice(22.5,28),lon = slice(56,70))

dic_term_nas = dic_term.sel(lat = slice(22.5,28),lon = slice(56,70))

alk_term_nas = alk_term.sel(lat = slice(22.5,28),lon = slice(56,70))

S_term_nas = S_term.sel(lat = slice(22.5,28),lon = slice(56,70))
#%%
dpco2_dt_cal_esio = dpco2_dt_cal.sel(lat = slice(-6.6,8),lon = slice(92,109))

T_term_esio = T_term.sel(lat = slice(-6.6,8),lon = slice(92,109))

dic_term_esio = dic_term.sel(lat = slice(-6.6,8),lon = slice(92,109))

alk_term_esio = alk_term.sel(lat = slice(-6.6,8),lon = slice(92,109))

S_term_esio = S_term.sel(lat = slice(-6.6,8),lon = slice(92,109))
#%%
dpco2_dt_cal_eio = dpco2_dt_cal.sel(lat = slice(-6.5,5),lon = slice(49,92))

T_term_eio = T_term.sel(lat = slice(-6.5,5),lon = slice(49,92))

dic_term_eio = dic_term.sel(lat = slice(-6.5,5),lon = slice(49,92))

alk_term_eio = alk_term.sel(lat = slice(-6.5,5),lon = slice(49,92))

S_term_eio = S_term.sel(lat = slice(-6.5,5),lon = slice(49,92))

#%%  

uatm_temp_nwio = T_term_nwio.cumsum(dim='time') * dt
uatm_dic_nwio = dic_term_nwio.cumsum(dim='time') * dt
uatm_alk_nwio = alk_term_nwio.cumsum(dim='time') * dt
uatm_sal_nwio = S_term_nwio.cumsum(dim='time') * dt

#%%

pco2_recon_nwio = dpco2_dt_cal_nwio.cumsum(dim='time')*dt

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon_nwio.mean(dim=('lat','lon')),
    label='Reconstructed pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_nwio.mean(dim=('lat','lon')) ,
    label='temperature-driven pCO$_2$ ',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_nwio.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_nwio.mean(dim=('lat','lon')) ,
    label='DIC-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_nwio.mean(dim=('lat','lon')),
    label='Salinity-driven pCO$_2$',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in Northwestern Indian Ocean')

plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

#%% NAS

uatm_temp_nas = T_term_nas.cumsum(dim='time') * dt
uatm_dic_nas = dic_term_nas.cumsum(dim='time') * dt
uatm_alk_nas = alk_term_nas.cumsum(dim='time') * dt
uatm_sal_nas = S_term_nas.cumsum(dim='time') * dt

pco2_recon_nas = dpco2_dt_cal_nas.cumsum(dim='time')*dt

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon_nas.mean(dim=('lat','lon')),
    label='Reconstructed pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_nas.mean(dim=('lat','lon')) ,
    label='temperature-driven pCO$_2$ ',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_nas.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_nas.mean(dim=('lat','lon')) ,
    label='DIC-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_nas.mean(dim=('lat','lon')),
    label='Salinity-driven pCO$_2$',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in Northern Arabian Sea')

plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
#%%

uatm_temp_eio = T_term_eio.cumsum(dim='time') * dt
uatm_dic_eio = dic_term_eio.cumsum(dim='time') * dt
uatm_alk_eio = alk_term_eio.cumsum(dim='time') * dt
uatm_sal_eio = S_term_eio.cumsum(dim='time') * dt

pco2_recon_eio = dpco2_dt_cal_eio.cumsum(dim='time')*dt

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon_eio.mean(dim=('lat','lon')),
    label='Reconstructed pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_eio.mean(dim=('lat','lon')) ,
    label='temperature-driven pCO$_2$ ',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_eio.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_eio.mean(dim=('lat','lon')) ,
    label='DIC-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_eio.mean(dim=('lat','lon')),
    label='Salinity-driven pCO$_2$',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in Equatorial Indian Ocean')

plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

#%%

uatm_temp_esio = T_term_esio.cumsum(dim='time') * dt
uatm_dic_esio = dic_term_esio.cumsum(dim='time') * dt
uatm_alk_esio = alk_term_esio.cumsum(dim='time') * dt
uatm_sal_esio = S_term_esio.cumsum(dim='time') * dt

pco2_recon_esio = dpco2_dt_cal_esio.cumsum(dim='time')*dt

plt.figure(figsize=(12, 5))

plt.plot(
    time,
    pco2_recon_esio.mean(dim=('lat','lon')),
    label='Reconstructed pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_temp_esio.mean(dim=('lat','lon')) ,
    label='temperature-driven pCO$_2$ ',
    linewidth=2
)

plt.plot(
    time,
    uatm_alk_esio.mean(dim=('lat','lon')),
    label='Alkalinity-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_dic_esio.mean(dim=('lat','lon')) ,
    label='DIC-driven pCO$_2$',
    linewidth=2
)

plt.plot(
    time,
    uatm_sal_esio.mean(dim=('lat','lon')),
    label='Salinity-driven pCO$_2$',
    linewidth=2
)


plt.xlabel('Year')
plt.ylabel('pCO₂ (µatm)')
plt.title('Contribution of pCO2 by different variables in Eastern Indain Ocean')

plt.legend(frameon=False, ncol=2)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
