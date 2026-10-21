# 17.06.2025
# Andrew G Twelves

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean as cm

custom = [(179/255, 88/255,  6/255),
          (224/255,130/255, 20/255),
          (253/255,184/255, 99/255),
          (254/255,224/255,182/255),
          (216/255,218/255,235/255),
          (178/255,171/255,210/255),
          (128/255,115/255,172/255),
          ( 84/255, 39/255,136/255)]

# Select station
station = 'Hamina'
# Zoom on February
month   = 2

# Open time series
ds   = xr.open_dataset('time_series/data/BO_TS_TG_{}_processed_UNFILTERED.nc'.format(station))
# Observed time series
obs  = ds.ssh_obs
# Centre relative to mean
obs  = obs - np.nanmean(obs)
# Model time series
mod  = ds.ssh_mod
# Centre relative to mean
mod  = mod - np.nanmean(mod)

# Calculate RMS and correlation
rms = np.sqrt(np.nanmean(np.square(mod - obs)))
cor = np.corrcoef(np.ravel(mod), np.ravel(obs))

# Declare empty arrays
obs_spectrum = np.zeros((31,17))
mod_spectrum = np.zeros((31,17))

# Open mesoscale part of observed spectrum
meso           = xr.open_dataarray('{}_obs_meso_{}.nc'.format(station,month)).values
# Open synoptic-scale part of observed spectrum
syno           = xr.open_dataarray('{}_obs_syno_{}.nc'.format(station,month)).values
# Concatenate
obs_spectrum[1:30,0:8]  = meso
obs_spectrum[1:30,8:16] = syno

# Open mesoscale part of model spectrum
meso           = xr.open_dataarray('{}_mod_meso_{}.nc'.format(station,month)).values
# Open synoptic-scale part of model spectrum
syno           = xr.open_dataarray('{}_mod_syno_{}.nc'.format(station,month)).values
# Concatenate
mod_spectrum[1:30,0:8]  = meso
mod_spectrum[1:30,8:16] = syno

# Find absolute and relative anomalies in wavelet amplitudes
dif_spectrum_abs = np.sqrt(mod_spectrum) - np.sqrt(obs_spectrum)
dif_spectrum_rel = np.divide(np.sqrt(mod_spectrum) - np.sqrt(obs_spectrum),np.sqrt(obs_spectrum))

# One figure with three panels
plt.figure(figsize=(30,10))

# Plot spectral amplitude anomaly
plt.contour(np.linspace(0.5,30.5,31),np.linspace(0,16,17),100*np.transpose(dif_spectrum_abs),levels=[0],linewidths=5,colors='gray')
plt.contourf(np.linspace(0.5,30.5,31),np.linspace(0,16,17),100*np.transpose(dif_spectrum_abs),levels=[-12,-8,-4,0,4,8,12],colors=custom,extend='both')

plt.grid()
plt.xticks(np.linspace(1,29,29),fontsize=20)
plt.yticks([-12,-9.5,-7,-4.5,-2,  0,4,8,12,  18.5,21,23.5,26,28.5,31,33.5],['-50','','0','','50', '4','16','64','256',  '','0','','50','','100',''],fontsize=25)
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
cb=plt.colorbar()
plt.xlim(1,30)
cb=plt.colorbar(orientation='horizontal')
cb.ax.tick_params(labelsize=30) 

top_offset = 19
# Plot observed time series in offset top panel
plt.plot(np.linspace(1,30,29*24),10*obs[123*24:152*24]+top_offset,linewidth=3,color='k',label='obs')
# Plot model time series in offset top panel
plt.plot(np.linspace(1,30,29*24),10*mod[123*24:152*24]+top_offset,linewidth=3,color='k',linestyle='--',label='mod')

bottom_offset = -7
# Plot difference between time series in offset bottom panel
plt.plot(np.linspace(1,30,29*24),10*(mod[123*24:152*24]-obs[123*24:152*24])+bottom_offset,linewidth=3,color='k')

plt.plot(np.linspace(1,30,29*24),np.linspace(0,0,29*24),'k')
plt.plot(np.linspace(1,30,29*24),np.linspace(16,16,29*24),'k')
plt.plot(np.linspace(1,30,29*24),np.linspace(21,21,29*24),'k')
plt.plot(np.linspace(1,30,29*24),np.linspace(-7,-7,29*24),'k')
plt.xlabel('February 2020',fontsize=30)
plt.ylabel('Period (hours)    ',fontsize=30)
plt.ylim(-13,33)
plt.savefig('{}_transform_zoom_{}'.format(station,month))

# Find model error in capturing extreme sea level
diff = np.squeeze(mod[123*24:152*24]) - np.squeeze(obs[123*24:152*24])
print(100*np.nanmax(diff[400:]))
print(100*np.nanmin(diff[400:]))
