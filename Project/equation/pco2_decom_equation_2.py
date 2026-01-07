#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 19:12:44 2026

@author: bobco-08
"""

import xarray as xr
import gsw as gsw

DIC = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/DIC_2deg.nc", engine="netcdf4")
ALK = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/ALK_2deg.nc", engine="netcdf4")
SST = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SST_2deg.nc", engine="netcdf4")
SSS = xr.open_dataset("/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/SSS_2deg.nc", engine="netcdf4")

#variables

dic = DIC['DIC']
alk = ALK['ALK']
sss = SSS['SSS']
sst = SST['SST']
#%%
import xarray as xr
import gsw
import PyCO2SYS as pyco2
import numpy as np

# Broadcast lon/lat
lon2d, lat2d = xr.broadcast(sss.lon, sss.lat)

# Absolute Salinity
SA = gsw.SA_from_SP(sss, 0, lon2d, lat2d)

# Density (kg/m3)
rho = gsw.rho(SA, sst, 0)

# Convert mol/m3 → µmol/kg
dic_ukg = dic * 1e6 / rho
alk_ukg = alk * 1e6 / rho
#%%
out_ref = pyco2.sys(
    par1=dic_ukg,
    par2=alk_ukg,
    par1_type=2,      # DIC
    par2_type=1,      # ALK
    temperature=sst,  # in-situ temperature (°C)
    salinity=sss,
    pressure=xr.zeros_like(sst),
    opt_k_carbonic=10,   # Lueker et al. (2000)
    opt_k_bisulfate=1
)

pco2_ref = out_ref["pCO2"]

#%%
time_chunk = 30   # 30 days (safe)

out_ds = xr.Dataset(
    coords=dict(
        time=dic.time,
        lat=dic.lat,
        lon=dic.lon
    )
)

out_ds["pCO2"] = xr.DataArray(
    np.full(dic.shape, np.nan, dtype=np.float32),
    dims=("time","lat","lon")
)

out_ds.to_netcdf("pCO2_ref.nc", mode="w")
#%%
import PyCO2SYS as pyco2

for i in range(0, dic.sizes["time"], time_chunk):

    sl = slice(i, i + time_chunk)

    out = pyco2.sys(
        par1=dic_ukg.isel(time=sl),
        par2=alk_ukg.isel(time=sl),
        par1_type=2,
        par2_type=1,
        temperature=sst.isel(time=sl),
        salinity=sss.isel(time=sl),
        pressure=0,
        opt_k_carbonic=10,
        opt_k_bisulfate=1
    )

    pco2_chunk = out["pCO2"]

    with xr.open_dataset("pCO2_ref.nc", mode="a") as ds:
        ds["pCO2"][sl,:,:] = pco2_chunk

    del out, pco2_chunk

