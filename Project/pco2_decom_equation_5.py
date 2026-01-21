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

dic_nwio = dic.sel(lat = slice(-6.6,8),lon = slice(92,109))

alk_nwio = alk.sel(lat = slice(-6.6,8),lon = slice(92,109))

sst_nwio = sst.sel(lat = slice(-6.6,8),lon = slice(92,109),time = slice('1994-01-01','2024-12-01'))

sss_nwio = sss.sel(lat = slice(-6.6,8),lon = slice(92,109),time = slice('1994-01-01','2024-12-01'))

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

pco2_list = []

for t in range(dic_nwio.sizes["time"]):

    out = pyco2.sys(
        par1=dic_nwio.isel(time=t).values,
        par2=alk_nwio.isel(time=t).values,
        par1_type=2,
        par2_type=1,
        temperature=sst_nwio.isel(time=t).values,
        salinity=sss_nwio.isel(time=t).values,
        pressure= np.zeros_like(dic_nwio.isel(time=0).values),
        opt_k_carbonic=10,
        opt_k_bisulfate=1
    )

    pco2_list.append(
        xr.DataArray(
            out["pCO2"].astype("float32"),
            dims=("lat", "lon"),
            coords={
                "lat": dic_nwio.lat,
                "lon": dic_nwio.lon,
                "time": dic_nwio.time.isel(time=t)
            }
        )
    )

# combine time slices
pco2_da = xr.concat(pco2_list, dim="time")

# save
pco2_da.to_netcdf("pCO2_ref_esio.nc")

#%%
dDIC = 1.0
dALK = 1.0
dT   = 0.1
dS   = 0.1
#%% DIC sensitivity
dDIC_da = xr.DataArray(dDIC, 
                       dims=dic_nwio.dims, 
                       coords=dic_nwio.coords )
dic_per = dic_nwio + dDIC_da

#%%

pco2_dic_list = []

for t in range(dic_nwio.sizes["time"]):

    out = pyco2.sys(
        par1=dic_per.isel(time=t).values,
        par2=alk_nwio.isel(time=t).values,
        par1_type=2,
        par2_type=1,
        temperature=sst_nwio.isel(time=t).values,
        salinity=sss_nwio.isel(time=t).values,
        pressure= np.zeros_like(dic_nwio.isel(time=0).values),
        opt_k_carbonic=10,
        opt_k_bisulfate=1
    )

    pco2_dic_list.append(
        xr.DataArray(
            out["pCO2"].astype("float32"),
            dims=("lat", "lon"),
            coords={
                "lat": dic_nwio.lat,
                "lon": dic_nwio.lon,
                "time": dic_nwio.time.isel(time=t)
            }
        )
    )

# combine time slices
pco2_dic_da = xr.concat(pco2_dic_list, dim="time")

# save
pco2_dic_da.to_netcdf("pCO2_dic_esio.nc")
#%%

dALK_da = xr.DataArray(dALK, 
                       dims=dic_nwio.dims, 
                       coords=dic_nwio.coords )
alk_per = alk_nwio + dALK_da


#%% alk sensitivity

pco2_alk_list = []

for t in range(dic_nwio.sizes["time"]):

    out = pyco2.sys(
        par1=dic_nwio.isel(time=t).values,
        par2=alk_per.isel(time=t).values,
        par1_type=2,
        par2_type=1,
        temperature=sst_nwio.isel(time=t).values,
        salinity=sss_nwio.isel(time=t).values,
        pressure= np.zeros_like(dic_nwio.isel(time=0).values),
        opt_k_carbonic=10,
        opt_k_bisulfate=1
    )

    pco2_alk_list.append(
        xr.DataArray(
            out["pCO2"].astype("float32"),
            dims=("lat", "lon"),
            coords={
                "lat": dic_nwio.lat,
                "lon": dic_nwio.lon,
                "time": dic_nwio.time.isel(time=t)
            }
        )
    )

# combine time slices
pco2_alk_da = xr.concat(pco2_alk_list, dim="time")

# save
pco2_alk_da.to_netcdf("pCO2_alk_esio.nc")
#%%

dT_da = xr.DataArray(dT, 
                       dims=dic_nwio.dims, 
                       coords=dic_nwio.coords )
sst_per = sst_nwio + dT_da



#%% alk sensitivity

pco2_temp_list = []

for t in range(dic_nwio.sizes["time"]):

    out = pyco2.sys(
        par1=dic_nwio.isel(time=t).values,
        par2=alk_per.isel(time=t).values,
        par1_type=2,
        par2_type=1,
        temperature=sst_per.isel(time=t).values,
        salinity=sss_nwio.isel(time=t).values,
        pressure= np.zeros_like(dic_nwio.isel(time=0).values),
        opt_k_carbonic=10,
        opt_k_bisulfate=1
    )

    pco2_temp_list.append(
        xr.DataArray(
            out["pCO2"].astype("float32"),
            dims=("lat", "lon"),
            coords={
                "lat": dic_nwio.lat,
                "lon": dic_nwio.lon,
                "time": dic_nwio.time.isel(time=t)
            }
        )
    )

# combine time slices
pco2_temp_da = xr.concat(pco2_temp_list, dim="time")

# save
pco2_temp_da.to_netcdf("pCO2_temp_esio.nc")
#%%
dS_da = xr.DataArray(dS, 
                       dims=dic_nwio.dims, 
                       coords=dic_nwio.coords )
sss_per = sss_nwio + dS_da

#%% alk sensitivity

pco2_sss_list = []

for t in range(dic_nwio.sizes["time"]):

    out = pyco2.sys(
        par1=dic_nwio.isel(time=t).values,
        par2=alk_per.isel(time=t).values,
        par1_type=2,
        par2_type=1,
        temperature=sst_nwio.isel(time=t).values,
        salinity=sss_per.isel(time=t).values,
        pressure= np.zeros_like(dic_nwio.isel(time=0).values),
        opt_k_carbonic=10,
        opt_k_bisulfate=1
    )

    pco2_sss_list.append(
        xr.DataArray(
            out["pCO2"].astype("float32"),
            dims=("lat", "lon"),
            coords={
                "lat": dic_nwio.lat,
                "lon": dic_nwio.lon,
                "time": dic_nwio.time.isel(time=t)
            }
        )
    )

# combine time slices
pco2_sss_da = xr.concat(pco2_sss_list, dim="time")

# save
pco2_sss_da.to_netcdf("pCO2_sss_esio.nc")

#%%

