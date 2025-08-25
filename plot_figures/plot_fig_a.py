# 16.06.2025
# Andrew G Twelves

# Script to plot time-averaged wavelet spectra

import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

stations = ['Kemi','Oulu','Raahe',
            'Pietarsaari','Vaasa',
            'Kaskinen','Pori','Rauma',
            'Degerby',
            'Hanko','Helsinki','Porvoo','Hamina']

# Create empty arrays to hold spectra
obs_spectrum = np.zeros((16))
mod_spectrum = np.zeros((16))

# tau for normalisation
exponent = np.linspace(0,7.5,16)
tau      = np.power(2,exponent)/6

for st in stations:
    # read mesoscale part of tide gauge spectrum
    meso           = xr.open_mfdataset('{}_obs_meso_*'.format(st)).mean(dim='date').__xarray_dataarray_variable__.values
    # read synoptic scale part of tide gauge spectrum
    syno           = xr.open_mfdataset('{}_obs_syno_*'.format(st)).mean(dim='date').__xarray_dataarray_variable__.values
    # concatenate
    obs_spectrum[0:8]  = meso
    obs_spectrum[8:16] = syno
    # normalise by tau
    obs_spectrum       = obs_spectrum/tau
    # find peak, ignoring highest frequency
    obs_peak_position  = np.argmax(obs_spectrum[1:]) + 1
    obs_peak_freq      = tau[obs_peak_position]
    obs_peak_max       = obs_spectrum[obs_peak_position]

    # read mesoscale part of model spectrum
    meso           = xr.open_mfdataset('{}_mod_meso_*'.format(st)).mean(dim='date').__xarray_dataarray_variable__.values
    # read synoptic scale part of model spectrum
    syno           = xr.open_mfdataset('{}_mod_syno_*'.format(st)).mean(dim='date').__xarray_dataarray_variable__.values
    # concatenate
    mod_spectrum[0:8]  = meso
    mod_spectrum[8:16] = syno
    # normalise by tau
    mod_spectrum       = mod_spectrum/tau
    # find peak, ignoring highest frequency
    mod_peak_position  = np.argmax(mod_spectrum[1:]) + 1
    mod_peak_freq      = tau[mod_peak_position]
    mod_peak_max       = mod_spectrum[mod_peak_position]
    
    # plot figure
    plt.figure(figsize=(20,8))
    # shading for mesoscale part of spectrum
    plt.fill_betweenx([0,1.05*10000*np.nanmax(obs_spectrum)],[0,0],[2,2],color=(0.8,0.8,0.8))
    # tide gauge
    plt.plot(tau,10000*obs_spectrum,linewidth=5,color=(27/255,158/255,119/255))
    plt.scatter(tau,10000*obs_spectrum,s=1000,color=(27/255,158/255,119/255),alpha=0.5,facecolor='none',linewidth=4)
    plt.scatter(obs_peak_freq,10000*obs_peak_max,s=1000,color=(27/255,158/255,119/255),alpha=0.5,label='{:.1f} hrs    {:.1f} cm² day⁻¹'.format(obs_peak_freq*24,obs_peak_max*10000))
    # NEMO
    plt.plot(tau,10000*mod_spectrum,linewidth=5,color=(117/255,112/255,179/255))
    plt.scatter(tau,10000*mod_spectrum,s=1000,color=(117/255,112/255,179/255),alpha=0.5,facecolor='none',linewidth=4)
    plt.scatter(mod_peak_freq,10000*mod_peak_max,s=1000,color=(117/255,112/255,179/255),alpha=0.5,label='{:.1f} hrs    {:.1f} cm² day⁻¹'.format(mod_peak_freq*24,mod_peak_max*10000))
    # formatting
    plt.xticks([5,10,15,20,25,30],fontsize=40)
    plt.yticks(fontsize=40)
    plt.xlim(0,31)
    plt.xlabel('Days',fontsize=50)
    plt.ylim(0,1.05*10000*np.nanmax(obs_spectrum))
    plt.ylabel(r'$\overline{W}$² (cm² day⁻¹)', fontsize=50)
    plt.grid()
    plt.locator_params(axis='y', nbins=4)
    plt.legend(fontsize=50)
    plt.tight_layout()
    plt.savefig('{}_tavgd_spectrum.png'.format(st))
