import xarray as xr
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# =========================
# 1. Load Data
# =========================
file1 = "/home/bobco-08/24cl05012/CO2/data/data_1/cmems/bio_var/cmems_obs-mob_glo_bgc-car_my_irr-i_1772520787062.nc"

ds = xr.open_dataset(file1)
ds = ds.rename({'latitude': 'lat', 'longitude': 'lon'})

co2_flux = ds['fgco2'] * 12   # convert to gC/m2/yr

# Time selection
co2_flux = co2_flux.sel(time=slice('1994-01-01', '2024-12-31'))

# =========================
# 2. Define Regions
# =========================
regions = { 
    'NWIO':  dict(lat=slice(5,22.5), lon=slice(45,65)),         #selecting regions 
    'NAS':  dict(lat=slice(22.5,28), lon=slice(56,70)),
    'EIO': dict(lat = slice(-6.5,5),lon = slice(49,92)),
    'ESIO': dict(lat=slice(-6.6,8), lon=slice(92,109)),
}

# =========================
# 3. Seasonal Function
# =========================
def seasonal_timeseries(data, months, season_name):
    
    if season_name == 'DJF':                                   #selecting DJF
        ds_season = data.resample(time='QS-DEC').mean()
        ds_season = ds_season.sel(time=ds_season.time.dt.month == 12)
    else:
        ds_season = data.sel(time=data.time.dt.month.isin(months)) #selecting JJAS,MAM,ON
        ds_season = ds_season.resample(time='YE').mean()

    # spatial mean
    ds_mean = ds_season.mean(dim=('lat', 'lon'))

    # remove NaNs
    ds_mean = ds_mean.dropna(dim='time')

    return ds_mean

# =========================
# 4. Seasons
# =========================
seasons = {
    'DJF':  [12,1,2],
    'MAM':  [3,4,5],
    'JJAS': [6,7,8,9],
    'ON':   [10,11]
}

# =========================
# 5. Plot Function
# =========================
def plot_region(region_name, region_data):

    fig, axs = plt.subplots(2, 2, figsize=(30,14))
    axs = axs.flatten()

    for i, (season, months) in enumerate(seasons.items()):

        ax = axs[i]

        ts = seasonal_timeseries(region_data, months, season)

        years = ts.time.dt.year.values
        values = ts.values

        # trend
        coeff = np.polyfit(years, values, 1)
        trend = np.polyval(coeff, years)

        stats_out = stats.linregress(years, values)
        pval = stats_out.pvalue

        label = f"({coeff[0]:.2e}, p={pval:.3f})"

        # Plot
        ax.plot(years, values, 'o-',)
        ax.plot(years, trend, 'r--', lw=3, label=label)

        # ===== STYLE IMPROVEMENTS =====
        
        ax.set_title(season, fontsize=26, fontweight='bold', pad=15)

        ax.tick_params(axis='both', which='major', labelsize=20)
        
        # 👉 ADD THIS BLOCK HERE
        ax.set_xlabel('')
        ax.set_ylabel('')
        
        if i in [0, 2]:
            ax.set_ylabel('gC m$^{-2}$ yr$^{-1}$', fontsize=22)
        
        if i in [2, 3]:
            ax.set_xlabel('Years', fontsize=22)
        
        ax.grid(True, linestyle=':')
        ax.legend(fontsize=18)

    # ===== SUBPLOT ARRANGEMENT =====
    
    # Adjust spacing manually
    plt.subplots_adjust(
        left=0.093,
        right=0.95,
        top=0.90,
        bottom=0.089,
        wspace=0.094,   # horizontal spacing
        hspace=0.30    # vertical spacing
    )

    # Optional super title
    fig.suptitle(f"{region_name} CO$_2$ Flux Trend (1994–2024)", fontsize=26, y = 0.99)

    plt.show()

# =========================
# 6. Run for All Regions
# =========================
for name, bounds in regions.items():

    region_data = co2_flux.sel(**bounds)
    plot_region(name, region_data)
    
