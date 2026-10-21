# Andrew G Twelves
# Last edit 21.10.26

# Load modules
import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean

i = -1

custom = [(166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),
          (166/255,206/255,227/255),( 31/255,120/255,180/255),
          (166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),
          (166/255,206/255,227/255),
          (166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),( 51/255,160/255,44/250)]

season_colors = [(166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),( 51/255,160/255,44/250)]

stations = ['Kemi','Oulu','Raahe',
            'Pietarsaari','Vaasa',
            'Kaskinen','Pori','Rauma',
            'Degerby',
            'Hanko','Helsinki','Porvoo','Hamina']

letter = ['a','b','c','d','e','f','g','h','i','j','k','l','m']

markers = ['v','^','>',
           '<','8',
           'p','s','h',
           'o',
           '*','P','X','D']

for st in stations:

    i = i+1
    
    # Open mesoscale part of observed spectrum
    ds = xr.open_mfdataset('{}_obs_meso*'.format(st))
    # Sum over mesoscale frequencies
    mes = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    
    # Open synoptic-scale part of observed spectrum
    ds = xr.open_mfdataset('{}_obs_syno*'.format(st))
    # Sum over synoptic-scale frequencies
    syn = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    
    # Select date of extreme sea level
    tt = ds.date.isel(date=144).values

    # Open mesoscale part of model spectrum
    ds = xr.open_mfdataset('{}_mod_meso*'.format(st))
    # Sum over mesoscale frequencies
    mes_mod = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    # Open synoptic-scale part of model spectrum
    ds = xr.open_mfdataset('{}_mod_syno*'.format(st))
    # Sum over synoptic-scale frequencies
    syn_mod = ds.__xarray_dataarray_variable__.sum(dim='freq').values

    plt.figure(figsize=(10,9))
    # Mesoscale anomaly in spectral power
    mes_dif = mes_mod-mes
    # Synoptic-scale anomaly in spectral power
    syn_dif = syn_mod-syn
    plt.plot([-6,6],[0,0],'k',linewidth=2)
    plt.plot([0,0],[-6,6],'k',linewidth=2)
    # Plot synoptic-scale anomalies against mesoscale anomalies
    plt.scatter(mes_dif/np.std(mes),syn_dif/np.std(syn),s=250,label='{}'.format(st),c='gray',alpha=0.5,edgecolors='none')
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    # Highlight date of extreme sea level
    plt.scatter(mes_dif[144]/np.std(mes),syn_dif[144]/np.std(syn),s=1500,label='{}'.format(st),facecolors='none',edgecolors='k',linewidths=3)
    plt.xticks(fontsize=30,rotation=60)
    plt.yticks(fontsize=30,rotation=60)
    plt.ylabel('$\Delta W²_{syn}$ / $_{STDEV}(W²_{syn})$',fontsize=50)
    plt.xlabel('$\Delta W²_{meso}$ / $_{STDEV}(W²_{meso})$',fontsize=50)
    plt.grid()
    plt.locator_params(axis='y', nbins=7)
    plt.locator_params(axis='x', nbins=7)
    plt.tight_layout()
    plt.savefig('ratio_anom_seasonal_{}.png'.format(st))

    plt.figure(figsize=(10,5))
    plt.plot([92,92],[-6,6],'k',linewidth=3,linestyle='--')
    plt.plot([183,183],[-6,6],'k',linewidth=3,linestyle='--')
    plt.plot([274,274],[-6,6],'k',linewidth=3,linestyle='--')
    
    # Plot time series of mesoscale spectral power anomaly
    plt.plot(mes_dif/np.std(mes),label='meso',c=(166/255,206/255,227/255),linewidth=3)
    # Plot time series of synoptic-scale spectral power anomaly
    plt.plot(syn_dif/np.std(syn),label='syn',c=( 31/255,120/255,180/255),linestyle=':',linewidth=3)
    
    plt.xlim(0,366)
    plt.ylim(-6,6)
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
                                                            ,'       S'],fontsize=25)
    plt.yticks(fontsize=20)
    plt.ylabel('$\Delta W²$ / $_{STDEV}(W²)$',fontsize=25)
    plt.legend(fontsize=20)
    plt.grid()
    plt.savefig('mes_syn_series_{}.png'.format(st))
