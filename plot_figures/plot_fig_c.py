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
    plt.figure(figsize=(10,10))
    ds = xr.open_mfdataset('{}_obs_meso*'.format(st))
    tt = ds.date.isel(date=144).values
    mes = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    ds = xr.open_mfdataset('{}_obs_syno*'.format(st))
    syn = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    color_range = np.arange(366)
    color_range[0:31]=1
    color_range[31:61]=2
    color_range[61:92]=3
    color_range[92:123]=4
    color_range[123:152]=5
    color_range[152:183]=6
    color_range[183:213]=7
    color_range[213:244]=8
    color_range[244:274]=9
    color_range[274:305]=10
    color_range[305:336]=11
    color_range[336:366]=12
    
    plt.scatter(mes,syn,s=250,marker=markers[i],label='{}'.format(st),c=color_range,cmap='twilight',edgecolors='k')
    plt.scatter(mes[144],syn[144],s=500,marker=markers[i],label='{}'.format(st),c='red',edgecolors='k')
    plt.colorbar()
    plt.ylabel('Synoptic scale dispersion (m²)',fontsize=25)
    plt.xlabel('Mesoscale dispersion (m²)',fontsize=25)
    plt.xlim(-0.01,0.5)
    plt.ylim(-0.01,0.5)
    plt.grid()
    plt.savefig('ratios_seasonal_{}.png'.format(st))

    plt.figure(figsize=(10,10))
    plt.contourf([color_range,color_range],np.linspace(0.5,12.5,13),cmap='twilight')
    plt.colorbar()
    plt.savefig('discrete_colormap.png')

    plt.figure(figsize=(10,10))
    col = np.zeros((4,4))
    col[0,:]=0.5
    col[1,:]=1.5
    col[2,:]=2.5
    col[3,:]=3.5
    plt.contourf(col,np.linspace(0,4,5))
    plt.colorbar()
    plt.savefig('colorbar.png')

    # now plot model results

    plt.figure(figsize=(10,10))
    ds = xr.open_mfdataset('{}_mod_meso*'.format(st))
    mes_mod = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    ds = xr.open_mfdataset('{}_mod_syno*'.format(st))
    syn_mod = ds.__xarray_dataarray_variable__.sum(dim='freq').values
    plt.scatter(mes_mod,syn_mod,s=250,marker=markers[i],label='{}'.format(st),c=color_range,cmap='twilight',edgecolors='k')
    plt.scatter(mes_mod[144],syn_mod[144],s=500,marker=markers[i],label='{}'.format(st),c='red',edgecolors='k')
    plt.colorbar()
    plt.ylabel('Synoptic scale dispersion (m²)',fontsize=25)
    plt.xlabel('Mesoscale dispersion (m²)',fontsize=25)
    plt.grid()
    plt.xlim(-0.01,0.5)
    plt.ylim(-0.01,0.5)
    plt.savefig('ratio_mod_seasonal_{}.png'.format(st))

# now plot anomaly (UNFILTERED)

    plt.figure(figsize=(10,10))
    ds = xr.open_dataset('wavelet_transforms/data/BO_TS_TG_{}_Ricker_UNFILTERED.nc'.format(st))
    mes_dif = (mes_mod-mes)
    syn_dif = (syn_mod-syn)
    lim = np.nanmax([np.nanmax(abs(mes_dif)),np.nanmax(abs(syn_dif))])/0.9
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)
    plt.fill_between([-lim,0],[0,0],[lim,lim],color=(0.9,0.9,0.9))
    plt.fill_between([0,lim],[-lim,-lim],[0,0],color=(0.9,0.9,0.9))
    plt.scatter(mes_dif,syn_dif,s=250,marker=markers[i],label='{}'.format(st),c=color_range,cmap='twilight',edgecolors='k')
    plt.scatter(mes_dif[144],syn_dif[144],s=750,marker=markers[i],label='{}'.format(st),c=144,cmap='twilight',edgecolors='k')
    plt.colorbar()

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylabel('Synoptic scale variability anomaly (m²)',fontsize=20)
    plt.xlabel('Mesoscale variability anomaly (m²)',fontsize=20)
    plt.title('({}) {}'.format(letter[i],st),fontsize=20)
    plt.grid()
    plt.savefig('ratio_amom_seasonal_{}.png'.format(st))



