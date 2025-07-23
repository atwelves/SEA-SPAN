import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean

stations = ['Kemi','Oulu','Raahe',
            'Pietarsaari','Vaasa',
            'Kaskinen','Pori','Rauma',
            'Degerby',
            'Hanko','Helsinki','Porvoo','Hamina']

markers = ['v','^','>',
           '<','8',
           'p','s','h',
           'o',
           '*','P','X','D']

# set colormap - different colours for Gulf of Bothnia, Gulf of Finland and Degerby
custom = [( 27/255,158/255,119/255),( 27/255,158/255,119/255),( 27/255,158/255,119/255),
          ( 27/255,158/255,119/255),( 27/255,158/255,119/255),
          ( 27/255,158/255,119/255),( 27/255,158/255,119/255),( 27/255,158/255,119/255),
          (217/255, 95/255,  2/255),
          (117/255,112/255,179/255),(117/255,112/255,179/255),(117/255,112/255,179/255),(117/255,112/255,179/250)]

plt.figure(figsize=(10,10))
i = -1
mes_obs = np.zeros((13))
syn_obs = np.zeros((13))
for st in stations:
    i = i+1
    ds = xr.open_mfdataset('{}_obs_meso*'.format(st))
    mes_obs[i] = 10000*ds.__xarray_dataarray_variable__.sum(dim='freq').mean(dim='date').values
    ds = xr.open_mfdataset('{}_obs_syno*'.format(st))
    syn_obs[i] = 10000*ds.__xarray_dataarray_variable__.sum(dim='freq').mean(dim='date').values
    plt.scatter(mes_obs[i],syn_obs[i],s=250,marker=markers[i],c=custom[i],edgecolors='k')

gob_fit = np.polyfit(mes_obs[0:8],syn_obs[0:8],1)
gob_x = np.linspace(0,1000,10)
gob_y = gob_fit[1] + gob_fit[0]*gob_x
plt.plot(gob_x,gob_y,c=custom[0],linestyle=':',label='GoB: syn={:.0f}+{:.1f}*meso'.format(gob_fit[1],gob_fit[0]),linewidth=5,alpha=0.5)

gof_fit = np.polyfit(mes_obs[9:13],syn_obs[9:13],1)
gof_x = np.linspace(0,1000,10)
gof_y = gof_fit[1] + gof_fit[0]*gof_x
plt.plot(gof_x,gof_y,c=custom[9],linestyle=':',label='GoF: syn={:.0f}+{:.1f}*meso'.format(gof_fit[1],gof_fit[0]),linewidth=5,alpha=0.5)

plt.legend(fontsize=20)

plt.ylabel('Synoptic scale variability (cm²)',fontsize=25)
plt.xlabel('Mesoscale variability (cm²)',fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid()
plt.ylim(0,3000)
plt.xlim(0,1000)
plt.savefig('ratios_obs.png')

# now plot model results

plt.figure(figsize=(10,10))
i = -1
mes_mod = np.zeros((13))
syn_mod = np.zeros((13))
for st in stations:
    i = i+1
    ds = xr.open_mfdataset('{}_mod_meso*'.format(st))
    mes_mod[i] = 10000*ds.__xarray_dataarray_variable__.sum(dim='freq').mean(dim='date').values
    ds = xr.open_mfdataset('{}_mod_syno*'.format(st))
    syn_mod[i] = 10000*ds.__xarray_dataarray_variable__.sum(dim='freq').mean(dim='date').values
    plt.scatter(mes_mod[i],syn_mod[i],s=250,marker=markers[i],c=custom[i],edgecolors='k')

gob_fit = np.polyfit(mes_mod[0:8],syn_mod[0:8],1)
gob_x = np.linspace(0,1000,10)
gob_y = gob_fit[1] + gob_fit[0]*gob_x
plt.plot(gob_x,gob_y,c=custom[0],linestyle=':',label='GoB: syn={:.0f}+{:.1f}*meso'.format(gob_fit[1],gob_fit[0]),linewidth=5,alpha=0.5)

gof_fit = np.polyfit(mes_mod[9:13],syn_mod[9:13],1)
gof_x = np.linspace(0,1000,10)
gof_y = gof_fit[1] + gof_fit[0]*gof_x
plt.plot(gof_x,gof_y,c=custom[9],linestyle=':',label='GoF: syn={:.0f}+{:.1f}*meso'.format(gof_fit[1],gof_fit[0]),linewidth=5,alpha=0.5)

plt.legend(fontsize=20)

plt.ylabel('Synoptic scale variability (cm²)',fontsize=25)
plt.xlabel('Mesoscale variability (cm²)',fontsize=25)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(0,3000)
plt.xlim(0,1000)
plt.savefig('ratios_mod.png')

# now plot anomaly (UNFILTERED)

plt.figure(figsize=(10,10))
i = -1
mes = np.zeros((13))
syn = np.zeros((13))
plt.fill_between([-200,-100,0],[-200,-100,0],[0,0,0],color=(0.9,0.9,0.9))#label='syn = meso')
for st in stations:
    i = i+1
    mes = mes_mod[i] - mes_obs[i]
    syn = syn_mod[i] - syn_obs[i]
    plt.scatter(mes,syn,s=250,marker=markers[i],c=custom[i],edgecolors='k')
plt.ylabel('Synoptic scale variability anomaly (cm²)',fontsize=25)
plt.xlabel('Mesoscale variability anomaly (cm²)',fontsize=25)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-800,0)
plt.xlim(-200,0)
plt.savefig('ratios_anom.png')

plt.figure(figsize=(10,10))
i = -1
mes = np.zeros((13))
syn = np.zeros((13))
for st in stations:
    i = i+1
    mes = (mes_mod[i] - mes_obs[i]) / mes_obs[i]
    syn = (syn_mod[i] - syn_obs[i]) / syn_obs[i]
    plt.scatter(100*mes,100*syn,s=250,marker=markers[i],c=custom[i],edgecolors='k')
plt.ylabel('Synoptic scale dispersion (m²)',fontsize=25)
plt.xlabel('Mesoscale dispersion (m²)',fontsize=25)
plt.grid()
plt.savefig('ratios_anom_rel.png')

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
for st in stations:
    i = i+1
    ds = xr.open_dataset('data/BO_TS_TG_{}_201910.nc'.format(st))
    lat = ds.LATITUDE.values
    lon = ds.LONGITUDE.values
    plt.scatter(lon,lat,s=250,marker=markers[i],label='{}'.format(st),c=custom[i],edgecolors='k')

plt.legend(fontsize=15)
