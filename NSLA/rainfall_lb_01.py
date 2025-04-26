# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 09:23:13 2025

@author: user
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft
import cartopy.crs as ccrs
import cartopy.feature as cf


data =  xr.open_mfdataset("C:/Users/user/Documents/24CL05012/nsla/rainfall/*.nc")

lon = data['LONGITUDE']
lat = data['LATITUDE']
time = data['TIME']
rf = data['RAINFALL']


rf= rf.convert_calendar('noleap',dim = 'TIME')
# to remove the leap years
rf_365 = rf.groupby('TIME.dayofyear').mean(dim =('TIME'))

#%%

rf_JJAS = rf_365.isel(dayofyear = slice(151,272)).mean(dim = ('dayofyear'))

plt.figure(figsize = (12,8),dpi = 100)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(0,30,10)
im =  plt.contourf(lon,lat,rf_JJAS,levels,cmap = 'coolwarm_r',extend='both')
ax.add_feature(cf.LAND,color = 'grey')
# ax.add_feature(cf.BORDERS)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('Rainfall (mm)')
ax.set_title('JJAS average climatology (1991-2000)',pad= 30)

#%%

rft_JJAS =  rf_365.mean(dim = ('LONGITUDE','LATITUDE'))
plt.figure(figsize = (12,8),dpi = 100)
plt.plot(rft_JJAS)
plt.xlabel('TIME')
plt.ylabel('Rainfall (mm)')
plt.title('Time Series of Rainfall on JJAS (1991-2000)')
plt.legend()


#%%
rf_365_broadcasted = rf_365.isel(dayofyear=rf['TIME'].dt.dayofyear - 1)


rf_ano = rf - rf_365_broadcasted

rf_2000= rf_ano.sel(TIME = slice('2000-01-01','2000-12-31')).mean(dim = 'TIME')

plt.figure(figsize = (12,8),dpi = 100)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(-6,8,50)
im =  plt.contourf(lon,lat,rf_2000,levels,cmap = 'coolwarm_r')
ax.add_feature(cf.LAND,color = 'grey')
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('Rainfall (mm)')
plt.title('Spatial Rainfall Anomaly Plot over India in 2000',pad = 40)
#%%

rf_2000_sp = rf_ano.sel(TIME = slice('2000-01-01','2000-12-31')).mean(dim = ('LATITUDE','LONGITUDE'))

plt.figure(figsize = (12,8),dpi = 100)
plt.plot(rf_2000_sp)
plt.xlabel('TIME')
plt.ylabel('Rainfall (mm)')
plt.title('Time Series of Rainfall anomaly over India in 2000')


#%%


rf_lat = 10
rf_lon = 76.5

rf_area = rf.sel(LATITUDE = rf_lat,LONGITUDE = rf_lon,method ='nearest')


t = time[0:-3]

plt.figure(figsize = (12,8),dpi = 100)
plt.plot(t,rf_area)
plt.xlabel('TIME')
plt.ylabel('Rainfall (mm)')
plt.title('Time Series of Rainfall')
#%%

yf = np.fft.fft(rf_area)

N = np.size(rf_area)

omega = np.fft.fftfreq(N,1)

tp =1/omega

t = np.linspace(0,3650,3650)

plt.figure(figsize= (12,4),dpi =100)
plt.xlim(0,4000)
plt.plot(tp,np.abs(yf))
plt.xlabel('Time')
plt.ylabel('Power')
plt.title('Power Spectra of Rainfall')
plt.grid()

#%%
ff_ano = rf_area - rf.mean(dim =('LATITUDE','LONGITUDE'))

yf1 = np.fft.fft(ff_ano)

N1 = np.size(ff_ano)

omega1 = np.fft.fftfreq(N1,1)

tp1 =1/omega

t = np.linspace(0,3650,1826)

plt.figure(figsize= (12,4),dpi =100)
plt.xlim(0,4000)
plt.plot(tp1,np.abs(yf1))
plt.xlabel('Time')
plt.ylabel('Power')
plt.title('Power Spectra of anomaly in daily rainfall')
plt.grid()
#%% inverse fft

yf = np.fft.fft(rf_area)

# ck =abs(yf)

yf_inv = np.fft.ifft(yf)

t = np.linspace(0,3650,3650)

plt.figure(figsize= (12,4),dpi =100)
plt.plot(t,yf_inv)
plt.xlabel('Time')
plt.ylabel('rainfall')
plt.title('Time series of Rainfall(inverse method)')
plt.grid()
#%% overlapping original and reconstructed

plt.figure(figsize = (12,8),dpi =100)
plt.plot(t,rf_area)
plt.plot(t,yf_inv)

rf_diff = rf_area - yf_inv
plt.figure(figsize=(12,8),dpi =100)

plt.plot(t,rf_diff)


#%% annual cycle

omega_inv = 1/omega



sel_coeff = np.zeros_like(yf)

sel_coeff[0]= yf[0]
sel_coeff[10]= yf[10]
sel_coeff[-10]= yf[-10]

coeff_inv = np.fft.ifft(sel_coeff)
plt.plot(rf_area)
plt.plot(coeff_inv)
plt.xlabel('TIME')
plt.ylabel('RAINFALL ')
#%% semi - annual cycle

sel_coeff[20] = yf[20]
sel_coeff[-20] = yf[-20]

plt.figure(figsize =(12,8),dpi =100)
coeff_inv = np.fft.ifft(sel_coeff)
plt.plot(rf_area)
plt.plot(coeff_inv, label = 'smooth climatology')
plt.xlabel('TIME')
plt.ylabel('RAINFALL (mm)')
plt.title('Time series of Rainfall(Smooth climatology)')
plt.legend()
#%%

anomaly_smooth_clim = rf_area - coeff_inv

plt.plot(anomaly_smooth_clim)

plt.xlabel('TIME')
plt.ylabel('RAINFALL')
plt.title('Anomalies Computed by Removing Smoothed Climatology')



#%% intra seasonal 

sel_coeff2 = np.zeros_like(yf)

sel_coeff2[41:182] = yf[41:182]
sel_coeff2[-182:-41] = yf[-182:-41]

coeff_inv2 = np.fft.ifft(sel_coeff2)
plt.plot(rf_area)
plt.plot(coeff_inv2)
plt.xlabel('TIME')
plt.ylabel('RAINFALL')

#%%

rf_1991 = rf.sel(TIME = slice('1991-01-01','1991-12-31'),LATITUDE = rf_lat,LONGITUDE = rf_lon)

anomaly = rf_1991 - coeff_inv[0:365]


# plt.plot(anomaly,label = '1991 rainfall data')

plt.figure(figsize = (12,5),dpi =100)
plt.plot(anomaly_smooth_clim, label = 'Anomalies Computed by Removing Smoothed Climatology')
plt.plot(coeff_inv2, label = 'intra - seasonal variability',color = 'red')
plt.xlim(0,365)
plt.xlabel('TIME')
plt.ylabel('RAINFALL (mm)')
plt.title('Intra-seasonal oscillations overlayed above anomalies captured')
plt.grid()
plt.legend()
plt.show()
#%%
# data2 = xr.open_dataset('C:/Users/user/Documents/24CL05012/nslo/data/seasurfaceheight_1993_2020.nc')

# lon2 = data2['longitude']
# lat2 = data2['latitude']
# time2 = data2['time']
# sla = data2['sla']

# sla= sla.convert_calendar('noleap',dim = 'time')
# sla_spatial = sla.sel(time = slice('1993-01-01','2000-12-31'))
#%%

coeff_a1 = np.zeros_like(rf[0,:,:])
for i in range(len(lat)):
  for j in range(len(lon)):
        
        covariance = np.cov(rf.isel(LATITUDE = 10 ,LONGITUDE = 76),coeff_inv2)
        var_rf = np.var(rf.mean(dim =('LONGITUDE','LATITUDE')).fillna(0).values)
        coeff_a1= covariance/var_rf
           
plt.figure(figsize=(12, 8))
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(0,30,10)
im =  plt.contourf(coeff_a1,cmap = 'coolwarm_r',extend='both')
ax.add_feature(cf.LAND,color = 'grey')
# ax.add_feature(cf.BORDERS)
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)