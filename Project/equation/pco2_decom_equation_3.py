#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 17:42:25 2026

@author: bobco-08
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


file = '/home/bobco-08/24cl05012/CO2/codes/pCO2_codes/decomposition/pCO2_ref.nc'

data = xr.open_dataset(file)

time = data['time']
lat = data['lat']
lon = data['lon']
pco2_ref = data['pCO2']
#%%
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

import PyCO2SYS as pyco2

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
dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1
#%% DIC sensitivity
dDIC_da = xr.DataArray(dDIC, 
                       dims=dic_ukg.dims, 
                       coords=dic_ukg.coords )
dic_per = dic_ukg + dDIC
#%%
dic_pv = dic_per.values
dic_ukgv = dic_ukg.values

#%%
pco2_dic = np.empty(dic.shape, dtype=np.float32)

out = pyco2.sys(
    par1=dic_per.values,
    par2=alk_ukg.values,
    par1_type=2,
    par2_type=1,
    temperature=sst.values,
    salinity=sss.values,
    pressure=np.zeros_like(sst.values),
    opt_k_carbonic=10,
    opt_k_bisulfate=1
)
pco2_dic[:] = out["pCO2"].astype("float32")

# wrap in xarray
ds_out = xr.Dataset(
    data_vars=dict(
        pCO2_dic=(("time", "lat", "lon"), pco2_dic)
    ),
    coords=dict(
        time=dic.time,
        lat=dic.lat,
        lon=dic.lon
    )
)

ds_out.to_netcdf("pCO2_dic.nc")

#%%

dALK_da = xr.DataArray(dALK, 
                       dims=dic_ukg.dims, 
                       coords=dic_ukg.coords )
alk_per = alk_ukg + dALK

#%%
alk_pv = alk_per.values
alk_ukgv = alk_ukg.values
#%% alk sensitivity
pco2_alk = np.empty(dic.shape, dtype=np.float32)

out = pyco2.sys(
    par1=dic_ukg.values,
    par2=alk_per.values,
    par1_type=2,
    par2_type=1,
    temperature=sst.values,
    salinity=sss.values,
    pressure=np.zeros_like(sst.values),
    opt_k_carbonic=10,
    opt_k_bisulfate=1
)
pco2_alk[:] = out["pCO2"].astype("float32")

# wrap in xarray
ds_out = xr.Dataset(
    data_vars=dict(
        pCO2_alk=(("time", "lat", "lon"), pco2_alk)
    ),
    coords=dict(
        time=dic.time,
        lat=dic.lat,
        lon=dic.lon
    )
)

ds_out.to_netcdf("pCO2_alk.nc")
#%%

dT_da = xr.DataArray(dT, 
                       dims=dic_ukg.dims, 
                       coords=dic_ukg.coords )
T_per = sst + dT



#%% alk sensitivity
pco2_T = np.empty(dic.shape, dtype=np.float32)

out = pyco2.sys(  
    par1=dic_ukg.values,
    par2=alk_ukg.values,
    par1_type=2,
    par2_type=1,
    temperature=T_per.values,
    salinity=sss.values,
    pressure=np.zeros_like(sst.values),
    opt_k_carbonic=10,
    opt_k_bisulfate=1
)
pco2_T[:] = out["pCO2"].astype("float32")

# wrap in xarray
ds_out = xr.Dataset(
    data_vars=dict(
        pCO2_T=(("time", "lat", "lon"), pco2_T)
    ),
    coords=dict(
        time=dic.time,
        lat=dic.lat,
        lon=dic.lon
    )
)

ds_out.to_netcdf("pCO2_T.nc")
#%%
dS_da = xr.DataArray(dS, 
                       dims=dic_ukg.dims, 
                       coords=dic_ukg.coords )
S_per = sss + dS

#%%

#%% alk sensitivity
pco2_S = np.empty(dic.shape, dtype=np.float32)

out = pyco2.sys(
    par1=dic_ukg.values,
    par2=alk_ukg.values,
    par1_type=2,
    par2_type=1,
    temperature=sst.values,
    salinity=S_per.values,
    pressure=np.zeros_like(sst.values),
    opt_k_carbonic=10,
    opt_k_bisulfate=1
)
pco2_S[:] = out["pCO2"].astype("float32")

# wrap in xarray
ds_out = xr.Dataset(
    data_vars=dict(
        pCO2_S=(("time", "lat", "lon"), pco2_S)
    ),
    coords=dict(
        time=dic.time,
        lat=dic.lat,
        lon=dic.lon
    )
)

ds_out.to_netcdf("pCO2_S.nc")