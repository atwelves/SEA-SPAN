# 17.06.2025
# Andrew G Twelves

import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean as cm

station = 'Pori'
month   = 2

# time series
ds   = xr.open_dataset('time_series/data/BO_TS_TG_{}_processed_UNFILTERED.nc'.format(station))
print(ds)
obs  = ds.ssh_obs
obs  = obs - np.nanmean(obs)
print(obs)
mod  = ds.ssh_mod
mod  = mod - np.nanmean(mod)
print(mod)

obs_spectrum = np.zeros((29,16))
mod_spectrum = np.zeros((29,16))

# transforms
meso           = xr.open_dataarray('{}_obs_meso_{}.nc'.format(station,month)).values
# synoptic
syno           = xr.open_dataarray('{}_obs_syno_{}.nc'.format(station,month)).values
# concatenate
obs_spectrum[:,0:8]  = meso
obs_spectrum[:,8:16] = syno

# mesoscale
meso           = xr.open_dataarray('{}_mod_meso_{}.nc'.format(station,month)).values
# synoptic
syno           = xr.open_dataarray('{}_mod_syno_{}.nc'.format(station,month)).values
# concatenate
mod_spectrum[:,0:8]  = meso
mod_spectrum[:,8:16] = syno

dif_spectrum = mod_spectrum - obs_spectrum

plt.figure(figsize=(30,10))
#plt.pcolormesh(10000*np.transpose(dif_spectrum),cmap='cmo.delta_r')
#plt.colorbar()
#plt.plot(np.linspace(0,29,29*24),10*obs[123*24:152*24],linewidth=3,color=(230/255,97/255,1/255))
#plt.plot(np.linspace(0,29,29*24),10*mod[123*24:152*24],linewidth=3,color=(230/255,97/255,1/255),linestyle='--')

plt.pcolormesh(np.linspace(1,30,30),np.linspace(0,16,17),10000*np.transpose(dif_spectrum),shading='flat',cmap='cmo.tarn_r')
#plt.contourf(10000*np.transpose(dif_spectrum),cmap='cmo.tarn_r')
plt.grid()
plt.xticks(np.linspace(1,29,29),fontsize=20)
plt.yticks([-12,-9.5,-7,-4.5,-2,  0,4,8,12,  18.5,21,23.5,26,28.5,31,33.5],['-50','','0','','50', '4','16','64','256',  '','0','','50','','100',''],fontsize=25)
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
plt.colorbar()
#cb=plt.colorbar(orientation='horizontal')
#cb.ax.tick_params(labelsize=30) 
plt.plot(np.linspace(1,30,29*24),10*obs[123*24:152*24]+19,linewidth=3,color='k')
plt.plot(np.linspace(1,30,29*24),10*mod[123*24:152*24]+19,linewidth=3,color='k',linestyle='--')
plt.plot(np.linspace(1,30,29*24),10*(mod[123*24:152*24]-obs[123*24:152*24])-7,linewidth=3,color='k')
plt.plot(np.linspace(1,30,29*24),np.linspace(0,0,29*24),'k')
plt.plot(np.linspace(1,30,29*24),np.linspace(16,16,29*24),'k')
plt.plot(np.linspace(1,30,29*24),np.linspace(21,21,29*24),'k')
plt.plot(np.linspace(1,30,29*24),np.linspace(-7,-7,29*24),'k')
#plt.title('Wavelet power anomaly (cm$^{2}$)')
plt.clim(-350,350)
plt.xlabel('February 2020',fontsize=30)
plt.ylabel('Period (hours)    ',fontsize=30)
plt.title('{}'.format(station),fontsize=30)
plt.ylim(-13,33)
plt.savefig('{}_transform_zoom_{}'.format(station,month))

print(10000*np.nanmax(dif_spectrum))
print(10000*np.nanmin(dif_spectrum))

#plt.yticks(np.linspace(-8,16,13))
#plt.figure(figsize=(20,10))
#plt.plot(obs[123*24:152*24])
#plt.plot(mod[123*24:152*24])
#plt.xticks(np.linspace(0,29*24,29))

#plt.grid()
#plt.savefig('{}zoom_{}'.format(station,month))
