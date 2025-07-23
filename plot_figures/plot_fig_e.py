# 17.06.2025
# Andrew G Twelves

import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean as cm

station = 'Rauma'
month   = 2

# time series
#obs_series   = xr.open_dataarray
#mod_series   = xr.open_dataarray

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

plt.figure()
plt.pcolormesh(10000*np.transpose(dif_spectrum),cmap='cmo.delta_r')
plt.colorbar()
plt.title('Wavelet power anomaly (cm$^{2}$)')
plt.clim(-350,350)
plt.xlabel('February 2020')
plt.ylabel('Scale')
plt.savefig('{}_transform_zoom_{}'.format(station,month))

