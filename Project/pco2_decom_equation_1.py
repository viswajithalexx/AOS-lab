#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 16:29:42 2026

@author: bobco-08
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import xesmf as xe

dic_025 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m_talk-dissic_35.00E-110.00E_30.00S-30.00N_0.49-11.41m_2021-11-12-2025-12-26.nc')


dic_025 = dic_025.rename({'latitude':'lat', 'longitude': 'lon'})

dic = dic_025['dissic'].isel(depth = 0)
dic = dic.sel(time = slice('2021-11-12','2024-12-31'))
#%%

pco2_2 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/oc_v2025.pCO2.nc')

lat = pco2_2['lat']
lon = pco2_2['lon']
t = pco2_2['mtime']
pco2 = pco2_2['pCO2']

pco2_io = pco2.sel(mtime = slice('2021-11-12','2024-12-31'), lat = slice(-31,31),lon = slice(35,110))


#%%
sst_025 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_phy-all_my_0.25deg_P1D-m_thetao_oras_35.00E-110.00E_30.00S-30.00N_0.51m_1994-01-01-2024-12-31.nc')
sss_025 = xr.open_dataset('/home/bobco-08/24cl05012/CO2/data/data_1/cmems_mod_glo_phy-all_my_0.25deg_P1D-m_so_oras_35.00E-110.00E_30.00S-30.00N_0.51m_1994-01-01-2024-12-31.nc')

sst_025 = sst_025.rename({'latitude':'lat', 'longitude': 'lon'})
sss_025 = sss_025.rename({'latitude':'lat', 'longitude': 'lon'})

sst = sst_025['thetao_oras'].isel(depth= 0)
sst = sst.sel(time = slice('2021-11-12','2024-12-31'))
sss = sss_025['so_oras'].isel(depth= 0)
sss = sss.sel(time = slice('2021-11-12','2024-12-31'))
alk = dic_025['talk'].isel(depth= 0)
alk = alk.sel(time = slice('2021-11-12','2024-12-31'))
#%%

source_grid = xr.Dataset(coords={"lat": dic.lat.values,"lon": dic.lon.values })

target_grid = xr.Dataset(coords = {"lat": ("lat",pco2_io.lat.values),'lon': ('lon',pco2_io.lon.values)})

#%%
regridder = xe.Regridder(source_grid,target_grid,method = 'bilinear',periodic = False,ignore_degenerate =True,filename="weights_025_to_2deg.nc",
    reuse_weights=False)

regridder = xe.Regridder(source_grid,target_grid,method = 'bilinear',periodic = False,ignore_degenerate =True,filename="weights_025_to_2deg.nc",
    reuse_weights=True)

#%%
dic_2deg = regridder(dic)
dic_2deg.name = "DIC"

dic_2deg.to_netcdf(
    "DIC_2deg.nc",
    engine="netcdf4"
)

del dic_2deg

#%%
sst_2deg = regridder(sst)
sst_2deg.name = "SST"

sst_2deg.to_netcdf(
    "SST_2deg.nc",
    engine="netcdf4"
)

del sst_2deg
#%%
sss_2deg = regridder(sss)
sss_2deg.name = "SSS"

sss_2deg.to_netcdf(
    "SSS_2deg.nc",
    engine="netcdf4"
)

del sss_2deg

#%%
alk_2deg = regridder(alk)
alk_2deg.name = "ALK"

alk_2deg.to_netcdf(
    "ALK_2deg.nc",
    engine="netcdf4"
)

del alk_2deg

#%%
dic_025.close()
sst_025.close()
sss_025.close()
