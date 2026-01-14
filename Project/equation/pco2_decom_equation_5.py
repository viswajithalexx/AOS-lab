#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 10:53:38 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

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

#%%

import numpy as np
import xarray as xr
import PyCO2SYS as pyco2

# allocate output in memory
pco2_all = np.empty(dic_nwio.shape, dtype=np.float32)

out = pyco2.sys(
    par1=dic_nwio.values,
    par2=alk_nwio.values,
    par1_type=2,
    par2_type=1,
    temperature=sst_nwio.values,
    salinity=sss_nwio.values,
    pressure=np.zeros_like(sst_nwio.values),
    opt_k_carbonic=10,
    opt_k_bisulfate=1
)

pco2_all[:] = out["pCO2"].astype("float32")

# wrap in xarray
ds_out = xr.Dataset(
    data_vars=dict(
        pCO2=(("time", "lat", "lon"), pco2_all)
    ),
    coords=dict(
        time=dic_nwio.time,
        lat=dic_nwio.lat,
        lon=dic_nwio.lon
    )
)

ds_out.to_netcdf("pCO2_ref_nwio.nc")