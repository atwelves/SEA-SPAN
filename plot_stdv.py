# Last edited 21.10.26

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import csv
import cmocean

# List of stations
stations = ['Kemi','Oulu','Raahe',
            'Pietarsaari','Vaasa',
            'Kaskinen','Pori','Rauma',
            'Degerby',
            'Hanko','Helsinki','Porvoo','Hamina']

# This is just for use as the land-sea mask...
ds = xr.open_dataset('anomaly_maps/diff_mld_202009_mean.nc')
mld = ds.mld.values
lat = ds.nav_lat_grid_T.values
lon = ds.nav_lon_grid_T.values

lat[lat < 0] = np.nan
lon[lon < 0] = np.nan

lat = np.nanmean(lat,1)
lon = np.nanmean(lon,0)

mask = [(0.2,0.2,0.2)]

mld[mld==0] = np.nan
mld[mld>0] = 1
mld[mld<0] = 1
plt.figure(figsize=(8,10))
pcol = plt.contourf(lon[500:],lat[300:],mld[300:,500:],[-30,30],colors=mask,alpha=0.2)
plt.xticks(np.linspace(10,30,5),fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('$\degree$W',fontsize=25)
plt.ylabel('$\degree$N',fontsize=25)
i = -1

# Declare empty arrays
lat = np.zeros((13))
lon = np.zeros((13))
std = np.zeros((13))

# Now loop over stations
for st in stations:
    print(st)
    i = i+1

    # Open station time series
    ds = xr.open_dataset('data/BO_TS_TG_{}_201910.nc'.format(st))
    lat[i] = ds.LATITUDE.values
    lon[i] = ds.LONGITUDE.values
    with open('{}_hourlystats_N2000.txt'.format(st), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ',quotechar='|')
        date_time = [ row[1] for row in reader ]
    dates = np.asarray(date_time,dtype='str')
    start_date = np.argwhere(dates=='20191001T000000')
    end_date   = np.argwhere(dates=='20201001T000000')
    start_date = int(start_date)
    end_date   = int(end_date)
    with open('{}_hourlystats_N2000.txt'.format(st), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ',quotechar='|')
        ssh = [ row[2] for row in reader ]
    station_data           = np.asarray(ssh,dtype=np.float32)
    station_data = np.squeeze(station_data[start_date:end_date+1])
    data = station_data - np.nanmean(station_data)
    # Convert mm to m
    data = data/1000
    # Calculate standard deviation
    datastd = data.std(ddof=1)
    print(datastd)
    # Convert m to cm
    std[i] = 100*datastd
    print(std)
#    if st == "Degerby":
#        plt.scatter(lon,lat,s=250,marker=markers[i],label='Föglö',c='c',edgecolors='k')
#    else:
#        plt.scatter(lon,lat,s=250,marker=markers[i],label='{}'.format(st),c='c',edgecolors='k')
# Plot all stations on map
plt.scatter(lon,lat,s=250,marker='o',label='{}'.format(st),c=std,cmap='cmo.amp',edgecolors='k')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=30)
plt.savefig('stdv_map.png')
