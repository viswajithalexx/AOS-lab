# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 15:48:19 2025

@author: user
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import gsw as gsw


file = "C:/Users/user/Documents/24CL05012/nslo/data/ORAS5_1deg_SST_1990_2008.nc"
file2 = 'C:/Users/user/Documents/24CL05012/nslo/data/ORAS5_1deg_SSS_1990_2008.nc'

data1 = xr.open_dataset(file)
data2 = xr.open_dataset(file2)

lon3 = data1['LONN179_181']
lat3 = data1['LAT']
time3 = data1['TIME']
sst  = data1['SOSSTSST']
sss = data2['SOSALINE']

#%%

corr = xr.corr(sst,sss,dim = 'TIME')

plt.figure(figsize = (15,10),dpi = 100)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
levels = np.linspace(-1,1,9)
im =  plt.contourf(lon3,lat3,corr,levels,cmap = 'RdBu_r',)
ax.add_feature(cf.LAND,color = 'grey')
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('CORRELATION ')
ax.set_title('CORRELATION BETWEEN SST AND SSS ', pad = 30)

#%%

p = 0

SA = gsw.SA_from_SP(sss,p,lon3,lat3)

CT = gsw.CT_from_pt(SA,sst)

sigma0 = gsw.sigma0(SA,CT)

#%%

sst_bob = sst.sel(LAT = slice(7,30),LONN179_181 = slice(75,100))
sss_bob = sss.sel(LAT = slice(7,30),LONN179_181 = slice(75,100))


x= np.linspace(2.5,45,10)
y = np.linspace(2.5,45,10)

x,y= np.meshgrid(x,y)

sigma_grid = gsw.sigma0(x,y)
plt.figure(figsize = (15,10),dpi =100)
a = plt.contour(x,y,sigma_grid,levels =20)
plt.scatter(sss_bob,sst_bob, s = 0.6,alpha = 0.5,color = 'red')
plt.clabel(a, inline=True, fontsize=8, fmt="%.1f")

plt.xlabel("SSS (PSU)")
plt.ylabel("SST (°C)")
plt.title("Density Contours over SST-SSS Scatter of Bay of Bengal")

plt.show()

#%%

sst_as = sst.sel(LAT = slice(7,30),LONN179_181 = slice(55,70))
sss_as = sss.sel(LAT = slice(7,30),LONN179_181 = slice(55,70))


m= np.linspace(2.5,45,10)
n = np.linspace(2.5,45,10)

m,n= np.meshgrid(m,n)

sigma_as = gsw.sigma0(m,n)
plt.figure(figsize = (15,10),dpi =100)
b = plt.contour(m,n,sigma_as,levels = 11)
plt.scatter(sss_as,sst_as, s =0.4,alpha = 0.5,color = 'blue')
plt.clabel(b, inline=True, fontsize=8, fmt="%.1f")

plt.xlabel("SST (°C)")
plt.ylabel("SSS (PSU)")
plt.title("Density Contours over SST-SSS Scatter of Arabian Sea")

plt.show()
#%%
p = 0
sf = 0
sa_range = np.linspace(0,45,100)
ct_range = np.linspace(-7,45,100)

sea_freez = gsw.CT_freezing_poly(sa_range,p,sf)

max_den = gsw.CT_maxdensity(sa_range,p)



m,n = np.meshgrid(sa_range,ct_range)

sigma_as = gsw.sigma0(m,n)



plt.figure(figsize = (15,10),dpi =100)
c = plt.contour(m,n,sigma_as,colors = 'black',levels = 11)
plt.scatter(sss_as,sst_as,s = 0.4,alpha = 0.5,color ='violet')
plt.scatter(sa_range,sea_freez,s = 0.8,color = 'red')
plt.scatter(sa_range,max_den,s = 0.8,color = 'blue')
plt.clabel(c, inline=True, fontsize=8, fmt="%.1f")

plt.xlabel("SSS (PSU)")
plt.ylabel("SST (°C)")
plt.title("Density Contours over SST-SSS Scatter of Arabian Sea")

plt.show()

#%%

sf = 0
sa_range = np.linspace(0,40,100)
ct_range = np.linspace(-7,40,100)

sea_freez = gsw.CT_freezing_poly(sa_range,p,sf)

max_den = gsw.CT_maxdensity(sa_range,p)



m,n = np.meshgrid(sa_range,ct_range)

sigma_as = gsw.sigma0(m,n)



plt.figure(figsize = (15,10),dpi =100)
c = plt.contour(m,n,sigma_as,colors = 'black',levels = 20)
plt.scatter(sss_bob,sst_bob,s = 0.4,alpha = 0.5,color ='violet')
plt.scatter(sa_range,sea_freez,s = 0.6,color = 'red')
plt.scatter(sa_range,max_den,s = 0.6,color = 'blue')
plt.clabel(c, inline=True, fontsize=8, fmt="%.1f")

plt.xlabel("SSS (PSU)")
plt.ylabel("SST (°C)")
plt.title("Density Contours over SST-SSS Scatter of Bay of Bengal")

plt.show()

#%% TS diagram on the freshwater region

seasons = ['DJF','MAM','JJA','SON']

for i in range(len(seasons)):
    
        
        bob_fw_sst = sst.sel(TIME = time3.dt.season == seasons[i] ,LAT = slice(18,22),LONN179_181 = slice(87,92))
        bob_fw_sss= sss.sel(TIME = time3.dt.season == seasons[i] ,LAT = slice(18,22),LONN179_181 =slice(87,92))
        
        
        sa_r1= np.linspace(10,40,10)
        ct_r1 = np.linspace(10,40,10)
        
        x1,y1= np.meshgrid(sa_r1,ct_r1)
        
        sigma_grid1 = gsw.sigma0(x1,y1)
        plt.figure(figsize = (12,8),dpi =100)
        a = plt.contour(x1,y1,sigma_grid,levels =20,colors = 'black')
        plt.scatter(bob_fw_sss,bob_fw_sst, s = 0.6,alpha = 0.5,color = 'red')
        plt.clabel(a, inline=True, fontsize=8, fmt="%.1f")
        
        plt.xlabel("SSS (PSU)")
        plt.ylabel("SST (°C)")
        plt.title(f'Density Contours over SST-SSS Scatter of 18-22N and 87-92E during {seasons[i]}')
        plt.show()

#%%

seasons = ['DJF','MAM','JJA','SON']

for i in range(len(seasons)):
    
        
        bob_mx_sst = sst.sel(TIME = time3.dt.season == seasons[i] ,LAT = slice(0,12),LONN179_181 = slice(70,90))
        bob_mx_sss= sss.sel(TIME = time3.dt.season == seasons[i] ,LAT = slice(0,12),LONN179_181 =slice(70,90))
        
        
        sa_r1= np.linspace(10,40,10)
        ct_r1 = np.linspace(10,40,10)
        
        x1,y1= np.meshgrid(sa_r1,ct_r1)
        
        sigma_grid1 = gsw.sigma0(x1,y1)
        plt.figure(figsize = (12,8),dpi =100)
        a = plt.contour(x1,y1,sigma_grid,levels =20,colors = 'black')
        plt.scatter(bob_mx_sss,bob_mx_sst, s = 0.6,alpha = 0.5,color = 'red')
        plt.clabel(a, inline=True, fontsize=8, fmt="%.1f")
        
        plt.xlabel("SSS (PSU)")
        plt.ylabel("SST (°C)")
        plt.title(f'Density Contours over SST-SSS Scatter of 0-12N and 70-90E during {seasons[i]}')
        plt.show()

#%%
beta = gsw.alpha(SA,CT,p)*(10**4)

beta_t = beta.mean(dim =('TIME'))

plt.figure(figsize = (15,10),dpi = 100)
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
im =  plt.contourf(lon3,lat3,beta_t)
ax.add_feature(cf.LAND,color = 'grey')
ax.gridlines(visible=True,draw_labels=True)
cbar = plt.colorbar(im)
cbar.set_label('TEC')
ax.set_title('Surface Thermal Expansion Coefficient(TEC)', pad = 30)



#%%

x= np.linspace(2.5,45,10)
y = np.linspace(2.5,45,10)

x,y= np.meshgrid(x,y)

sigma_grid = gsw.sigma0(x,y)
plt.figure(figsize = (15,10),dpi =100)
a = plt.contour(x,y,sigma_grid,levels =20,colors = 'black')
plt.clabel(a, inline=True, fontsize=8, fmt="%.1f")

plt.xlabel("SSS (PSU)")
plt.ylabel("SST (°C)")
plt.title("Density Contours over SST-SSS Scatter of Bay of Bengal")

plt.show()
#%%
p1 = 4000
p = 0

m= np.linspace(2.5,35,10)
n =np.linspace(2.5,32.5,10)
m,n= np.meshgrid(m,n)

rho = gsw.rho(m,n,p)-1000
rho1 = gsw.rho(m,n,p1)-1000

plt.subplot(1,2,1)
a = plt.contour(m,n,rho,levels =40,colors ='black')
plt.scatter(20,15,color = 'red')
plt.scatter(18.66,10,color = 'blue')
plt.xlim(15,25)
plt.ylim(5,20)
plt.clabel(a, inline=True, fontsize=8, fmt="%.1f")
plt.xlabel("Absolute Salinity (g/kg)")
plt.ylabel("Conservative Temperature(°C)")
plt.title("Density Contours over SST-SSS Scatter of Bay of Bengal at 0m depth")
plt.subplot(1,2,2)
b = plt.contour(m,n,rho1,levels =40,colors ='black')
plt.scatter(20,15,color = 'red')
plt.scatter(18.66,10,color = 'blue')
plt.xlim(15,25)
plt.ylim(5,20)
plt.clabel(b, inline=True, fontsize=8, fmt="%.1f")
plt.xlabel("Absolute Salinity (g/kg)")
plt.ylabel("Conservative Temperature(°C)")
plt.title("Density Contours over SST-SSS Scatter of Bay of Bengal at 4000m depth")


plt.show()