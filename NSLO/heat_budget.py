# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 16:31:24 2025

@author: user
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

file1 = 'C:/Users/user/Documents/24CL05012/nslo/data/dswrf.ntat.mon.1981-2010.ltm.nc'
file2  = 'C:/Users/user/Documents/24CL05012/nslo/data/ulwrf.ntat.mon.1981-2010.ltm.nc'
file3 = 'C:/Users/user/Documents/24CL05012/nslo/data/uswrf.ntat.mon.1981-2010.ltm.nc'

data1 = xr.open_dataset(file1)
data2 = xr.open_dataset(file2)
data3 = xr.open_dataset(file3)

lon = data1['lon']
lat = data1['lat']
dswrf = data1['dswrf']
ulwrf = data2['ulwrf']
uswrf = data3['uswrf']


rt =  dswrf - uswrf-ulwrf

rem = rt.mean(dim =('time','lon'))
swr = dswrf - uswrf
LAT = np.linspace(-90,90,94)

plt.figure(figsize = (12,8),dpi =100)
plt.plot(lat, swr.mean(dim=("time","lon")), color='green', label='Mean SWR')
plt.plot(lat, ulwrf.mean(dim=("time","lon")), color='blue', label='Mean LWR')
plt.plot(lat, rem, color='red', label='Mean Radiation')
plt.xlabel("Latitude (°)")
plt.ylabel("Net Shortwave Radiation (W/m²)")
plt.title("Latitudinal distribution of annual mean radiation at the Top of the Atmosphere (TOA)")
plt.grid(True)
plt.legend()
#%%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# File paths
uswrf = 'C:/Users/user/Documents/24CL05012/nslo/data/uswrf.ntat.mon.1981-2010.ltm.nc'
ulwrf = 'C:/Users/user/Documents/24CL05012/nslo/data/ulwrf.ntat.mon.1981-2010.ltm.nc'
dswrf = 'C:/Users/user/Documents/24CL05012/nslo/data/dswrf.ntat.mon.1981-2010.ltm.nc'

#load datasets
data1 = xr.open_dataset(uswrf)
data2 = xr.open_dataset(ulwrf)
data3 = xr.open_dataset(dswrf)

#set variables
lat = data1['lat']
lon = data1['lon']
time = data1['time']
uswr = data1['uswrf']
ulwr = data2['ulwrf']
dswr = data3['dswrf']



def compute_heat_transport(lat,net_toa):
    
    #radiation absorbed by latitudes (differ in surface area)
    lat_ = net_toa['lat'].values
    netflux = net_toa.values
    #net radiation associated with each lat
    lat_rad = np.deg2rad(lat)
    
    surface_area = np.cos(lat_rad)
    
    net_rad = netflux*surface_area
    
    #south to north transport of energy
    
    lat_diff = np.diff(lat_rad)
    
    lat_diff = np.append(lat_diff,lat_diff[-1])
    
    area =  2*np.pi*(6.371e6)**2
    
    heat_cumulative = np.cumsum(net_rad*lat_diff)
    
    heat_transport = (area*heat_cumulative)
    return heat_transport
#total mean energy remaining on earth

net_toa = (dswr - uswr - ulwr).mean(dim =('time','lon'))

heat_transport_unbalanced = compute_heat_transport(lat, net_toa)

#plot
plt.figure(figsize = (10,5),dpi =100)
plt.plot(lat,heat_transport_unbalanced / 10**15,color = 'blue')
plt.axhline(0, color='red', linestyle='--')
plt.title('Total Meridional Heat Transport in the Earth System (Unbalanced)')
plt.xlabel('Latitude')
plt.ylabel('Heat Transport (PW)')
plt.grid(True)
plt.tight_layout()
plt.show()

#%%

net_toa = (dswr - uswr - ulwr).mean(dim =('time','lon'))
lat_ = net_toa['lat'].values
netflux = net_toa.values
#net radiation associated with each lat
lat_rad = np.deg2rad(lat)

surface_area = np.cos(lat_rad)

net_rad = netflux*surface_area

average_rad = np.sum(net_rad)/np.sum(surface_area)

net_rad_ano = net_toa - average_rad

heat_transport_balanced = compute_heat_transport(lat,net_rad_ano)

plt.figure(figsize = (10,5),dpi =100)
plt.plot(lat, heat_transport_balanced / 1e15, color='darkblue', label='Heat Transport (from TOA Anomaly)')
plt.axhline(0, color='red', linestyle='--', linewidth=1)
plt.title('Meridional Heat Transport (Balanced)\nNCEP Reanalysis (1981–2010 Climatology)')
plt.xlabel('Latitude')
plt.ylabel('Heat Transport (PW)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
