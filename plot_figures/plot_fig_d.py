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

for st in stations:
    # declare array for wavelet transform of tide gauge data
    obs = np.zeros((16,366))
    # declare array for wavelet transform of model output
    mod = np.zeros((16,366))
    i = i+1

    # tide gauge
    ds = xr.open_mfdataset('{}_obs_meso*'.format(st))
    # read in times for x axis
    time_stamp = ds.date.values
    # read in mesoscale component of tide gauge transform
    mes = ds.__xarray_dataarray_variable__.values
    ds = xr.open_mfdataset('{}_obs_syno*'.format(st))
    # read in synoptic component of tide gauge transform
    syn = ds.__xarray_dataarray_variable__.values
    # concatenate mesoscale and synoptic components
    obs[0:8,:] = np.transpose(mes[:,:])
    obs[8:16,:] = np.transpose(syn[:,:])

    # model output
    ds = xr.open_mfdataset('{}_mod_meso*'.format(st))
    # read in mesoscale component of tide gauge transform
    mes = ds.__xarray_dataarray_variable__.values
    ds = xr.open_mfdataset('{}_mod_syno*'.format(st))
    # read in synoptic component of tide gauge transform
    syn = ds.__xarray_dataarray_variable__.values
    # concatenate mesoscale and synoptic components
    mod[0:8,:] = np.transpose(mes[:,:])
    mod[8:16,:] = np.transpose(syn[:,:])

    dif = mod - obs
    
    plt.figure(figsize=(20,7))
    plt.contourf(10000*obs,[0,100,100000],colors=('w',season_colors[3],season_colors[3]),alpha=0.5)
    plt.contourf(10000*mod,[0,100,100000],colors=('w',season_colors[1],season_colors[1]),alpha=0.5)
    plt.grid()
    plt.xticks([0,31,61,92,123,152,183,213,244,274,305,336],['O','N','D','J','F','M','A','M','J','J','A','S'],fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('({}) {}'.format(letter[i],st),fontsize=40)
    plt.savefig('{}_threshold'.format(st))

