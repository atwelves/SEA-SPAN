# Last update 20.10.2026
# Andrew G Twelves

# Script to plot spectral amplitude anomalies

# Load modules

import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean

i = -1

# List of tide gauge stations
stations = ['Kemi','Oulu','Raahe',
            'Pietarsaari','Vaasa',
            'Kaskinen','Pori','Rauma',
            'Degerby',
            'Hanko','Helsinki','Porvoo','Hamina']

letter = ['a','b','c','d','e','f','g','h','i','j','k','l','m']

for st in stations:

    # Declare empty arrays
    obs = np.zeros((16,366))
    mod = np.zeros((16,366))
    i = i+1

    # Open mesoscale part of observed spectra
    ds = xr.open_mfdataset('{}_obs_meso*'.format(st))
    time_stamp = ds.date.values
    mes = ds.__xarray_dataarray_variable__.values

    # Open synoptic-scale part of observed spectra
    ds = xr.open_mfdataset('{}_obs_syno*'.format(st))
    syn = ds.__xarray_dataarray_variable__.values
    
    # Concatenate and transpose arrays
    obs[0:8,:] = np.transpose(mes[:,:])
    obs[8:16,:] = np.transpose(syn[:,:])

    # Open mesoscale part of model spectra
    ds = xr.open_mfdataset('{}_mod_meso*'.format(st))
    mes = ds.__xarray_dataarray_variable__.values

    # Open synoptic-scale part of model spectra
    ds = xr.open_mfdataset('{}_mod_syno*'.format(st))
    syn = ds.__xarray_dataarray_variable__.values

    # Concatenate and transpose arrays
    mod[0:8,:] = np.transpose(mes[:,:])
    mod[8:16,:] = np.transpose(syn[:,:])

    # Difference in spectral amplitude (square root of spectral power)
    dif_spectrum_abs = np.sqrt(mod) - np.sqrt(obs)

    # Mask out padding at start and end of time period
    dif_spectrum_abs[dif_spectrum_abs==0] = np.nan

    plt.figure(figsize=(20,5))
    # Lines to mark seasonal division
    plt.plot([92,92],[0,15],'k',linewidth=5,linestyle='--')
    plt.plot([183,183],[0,15],'k',linewidth=5,linestyle='--')
    plt.plot([274,274],[0,15],'k',linewidth=5,linestyle='--')
    # This contour divides underestimates of spectral amplitude from overestimates
    plt.contour( 100*dif_spectrum_abs,levels=[0],linewidths=2,colors='gray') 
    
    plt.contourf(100*dif_spectrum_abs,levels=[-12,-8,-4,0,4,8,12],colors=custom,extend='both')
    plt.grid()
    plt.xlim(0,366)
    #cb=plt.colorbar()
    #cb=plt.colorbar(orientation='horizontal')
    #cb.ax.tick_params(labelsize=40)
    plt.xticks([0,31,61,92,123,152,183,213,244,274,305,336],['       O'
                                                            ,'       N'
                                                            ,'       D'
                                                            ,'       J'
                                                            ,'       F'
                                                            ,'       M'
                                                            ,'       A'
                                                            ,'       M'
                                                            ,'       J'
                                                            ,'       J'
                                                            ,'       A'
                                                            ,'       S'],fontsize=0)
    plt.yticks([0,4,8,12],['4','16','64','256'],fontsize=30)
    plt.ylabel('Period (hours)',fontsize=40)
    #plt.title('({}) {}'.format(letter[i],st),fontsize=50)
    plt.savefig('{}_amplitude'.format(st))
