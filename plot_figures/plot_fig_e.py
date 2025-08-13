# 17.06.2025
# Andrew G Twelves

import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean as cm

# February 2020 record high sea level at Rauma
station = 'Rauma'
month   = 2

# time series
#obs_series   = xr.open_dataarray
#mod_series   = xr.open_dataarray

# declare array for wavelet transform of tide gauge data
obs_spectrum = np.zeros((29,16))
# declare array for wavelet transform of model output
mod_spectrum = np.zeros((29,16))

# mesoscale component of tide gauge transform 
meso           = xr.open_dataarray('{}_obs_meso_{}.nc'.format(station,month)).values
# synoptic component
syno           = xr.open_dataarray('{}_obs_syno_{}.nc'.format(station,month)).values
# concatenate mesoscale and synoptic components of tide gauge transform
obs_spectrum[:,0:8]  = meso
obs_spectrum[:,8:16] = syno

# mesoscale component of model output transform
meso           = xr.open_dataarray('{}_mod_meso_{}.nc'.format(station,month)).values
# synoptic component of model output transform
syno           = xr.open_dataarray('{}_mod_syno_{}.nc'.format(station,month)).values
# concatenate mesoscale and synoptic components
mod_spectrum[:,0:8]  = meso
mod_spectrum[:,8:16] = syno

# find anomaly in wavelet space
dif_spectrum = mod_spectrum - obs_spectrum

# plot wavelet transform anomaly for Rauma during February 2020
plt.figure()
plt.pcolormesh(10000*np.transpose(dif_spectrum),cmap='cmo.delta_r')
plt.colorbar()
plt.title('Wavelet power anomaly (cm$^{2}$)')
plt.clim(-350,350)
plt.xlabel('February 2020')
plt.ylabel('Scale')
plt.savefig('{}_transform_zoom_{}'.format(station,month))

