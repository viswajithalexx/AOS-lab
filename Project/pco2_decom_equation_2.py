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

CT = gsw.CT_from_pt(SA,sst)
# Density (kg/m3)
rho = gsw.rho(SA,CT, 0)

# Convert mol/m3 → µmol/kg
dic_ukg = dic * 1e6 / rho
alk_ukg = alk * 1e6 / rho

dic_ukg = dic_ukg.where(dic_ukg > 0)
alk_ukg = alk_ukg.where(alk_ukg > 0)

#%%
import numpy as np
import xarray as xr
import PyCO2SYS as pyco2

# allocate output in memory
pco2_all = np.empty(dic.shape, dtype=np.float32)

out = pyco2.sys(
    par1=dic_ukg.values,
    par2=alk_ukg.values,
    par1_type=2,
    par2_type=1,
    temperature=sst.values,
    salinity=sss.values,
    pressure=np.zeros_like(sst.values),
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
        time=dic.time,
        lat=dic.lat,
        lon=dic.lon
    )
)

ds_out.to_netcdf("pCO2_ref.nc")


#%%
out_test = pyco2.sys(
    par1=float(dic_ukg.isel(time=0, lat=10, lon=10)),
    par2=float(alk_ukg.isel(time=0, lat=10, lon=10)),
    par1_type=2,
    par2_type=1,
    temperature=float(sst.isel(time=0, lat=10, lon=10)),
    salinity=float(sss.isel(time=0, lat=10, lon=10)),
    pressure=0,
    opt_k_carbonic=10,
    opt_k_bisulfate=1
)

print(out_test["pCO2"])
