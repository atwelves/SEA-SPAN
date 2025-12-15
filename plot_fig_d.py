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
            'Föglö',
            'Hanko','Helsinki','Porvoo','Hamina']

letter = ['a','b','c','d','e','f','g','h','i','j','k','l','m']

markers = ['v','^','>',
           '<','8',
           'p','s','h',
           'o',
           '*','P','X','D']

for st in stations:
    obs = np.zeros((16,366))
    mod = np.zeros((16,366))
    i = i+1
    
    ds = xr.open_mfdataset('{}/{}_obs_meso*'.format(st,st))
    time_stamp = ds.date.values
    mes = ds.__xarray_dataarray_variable__.values
    ds = xr.open_mfdataset('{}/{}_obs_syno*'.format(st,st))
    syn = ds.__xarray_dataarray_variable__.values
    obs[0:8,:] = np.transpose(mes[:,:])
    obs[8:16,:] = np.transpose(syn[:,:])

    ds = xr.open_mfdataset('{}/{}_mod_meso*'.format(st,st))
    mes = ds.__xarray_dataarray_variable__.values
    ds = xr.open_mfdataset('{}/{}_mod_syno*'.format(st,st))
    syn = ds.__xarray_dataarray_variable__.values
    mod[0:8,:] = np.transpose(mes[:,:])
    mod[8:16,:] = np.transpose(syn[:,:])

    dif = mod - obs
    
    plt.figure(figsize=(20,5))
    plt.plot([92,92],[0,15],'k',linewidth=5,linestyle='--')
    plt.plot([183,183],[0,15],'k',linewidth=5,linestyle='--')
    plt.plot([274,274],[0,15],'k',linewidth=5,linestyle='--')
    plt.contourf(10000*obs,[0,100,100000],colors=('w',season_colors[3],season_colors[3]),alpha=0.5)
    plt.contourf(10000*mod,[0,100,100000],colors=('w',season_colors[1],season_colors[1]),alpha=0.5)
    plt.grid()
    plt.xlim(0,366)
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
                                                            ,'       S'],fontsize=30)
    plt.yticks([0,4,8,12],['4','16','64','256'],fontsize=30)
    plt.ylabel('Period (hours)',fontsize=40)
    #plt.title('({}) {}'.format(letter[i],st),fontsize=50)
    plt.savefig('fig_d/{}_threshold'.format(st))

