import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import cmocean

i = -1

#custom = [(166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),
#          (166/255,206/255,227/255),( 31/255,120/255,180/255),
#          (166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),
#          (166/255,206/255,227/255),
#          (166/255,206/255,227/255),( 31/255,120/255,180/255),(178/255,223/255,138/255),( 51/255,160/255,44/250)]

custom = [(179/255, 88/255,  6/255),
          (224/255,130/255, 20/255),
          (253/255,184/255, 99/255),
          (254/255,224/255,182/255),
          (216/255,218/255,235/255),
          (178/255,171/255,210/255),
          (128/255,115/255,172/255),
          ( 84/255, 39/255,136/255)]

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

total_pos = np.zeros((13))
total_neg = np.zeros((13))

major_pos = np.zeros((13))
major_neg = np.zeros((13))

five_pos = np.zeros((13))
five_neg = np.zeros((13))

ten_pos = np.zeros((13))
ten_neg = np.zeros((13))

fift_pos = np.zeros((13))
fift_neg = np.zeros((13))

denom = np.zeros((13))

threshold = 25

avg_dif = np.zeros((13))
std_dif = np.zeros((13))

top_half = np.zeros((13))
corr     = np.zeros((13))
corr_t   = np.zeros((13))

rms      = np.zeros((13))
rms_t    = np.zeros((13))

for st in stations:

    # time series
    ds   = xr.open_dataset('time_series/data/BO_TS_TG_{}_processed_UNFILTERED.nc'.format(st))
    print(ds)
    obs  = ds.ssh_obs
    obs  = obs - np.nanmean(obs)
    mod  = ds.ssh_mod
    mod  = mod - np.nanmean(mod)
    corr_matrix = np.corrcoef(mod,obs)
    corr_t[i] = corr_matrix[0,1]
    rms_t[i] = np.sqrt(np.nanmean(np.square(mod - obs))) 

#

    obs = np.zeros((16,366))
    mod = np.zeros((16,366))
    i = i+1

    ds = xr.open_mfdataset('{}_obs_meso*'.format(st))
    time_stamp = ds.date.values
    mes = ds.__xarray_dataarray_variable__.values
    ds = xr.open_mfdataset('{}_obs_syno*'.format(st))
    syn = ds.__xarray_dataarray_variable__.values
    obs[0:8,:] = np.transpose(mes[:,:])
    obs[8:16,:] = np.transpose(syn[:,:])

    ds = xr.open_mfdataset('{}_mod_meso*'.format(st))
    mes = ds.__xarray_dataarray_variable__.values
    ds = xr.open_mfdataset('{}_mod_syno*'.format(st))
    syn = ds.__xarray_dataarray_variable__.values
    mod[0:8,:] = np.transpose(mes[:,:])
    mod[8:16,:] = np.transpose(syn[:,:])

    obs_std = np.std(obs,1)
    print(np.shape(obs_std))
    obs_std = np.tile(obs_std,(366,1))
    obs_std = np.transpose(obs_std)
    print(np.shape(obs_std))

    corr_matrix = np.corrcoef(np.ravel(np.sqrt(mod)),np.ravel(np.sqrt(obs)))
    corr[i] = corr_matrix[0,1]

    dif_spectrum_abs = np.sqrt(mod) - np.sqrt(obs)

    pos = np.copy(dif_spectrum_abs)
    pos[pos<0] = 0
    pos[np.isnan(pos)] = 0
    total_pos[i] = np.count_nonzero(pos)

    neg = np.copy(dif_spectrum_abs)
    neg[neg>0] = 0
    neg[np.isnan(neg)] = 0
    total_neg[i] = np.count_nonzero(neg)

    denom[i] = total_neg[i] + total_pos[i]

    dif_spectrum_rel = np.divide(dif_spectrum_abs,np.sqrt(obs))
    
    rms[i] = np.sqrt(np.nanmean(np.square(np.ravel(mod) - np.ravel(obs))))

    std_dif[i] = np.std(dif_spectrum_abs)/np.std(np.sqrt(obs))
    avg_dif[i] = np.nanmean(dif_spectrum_abs)/np.nanmean(np.sqrt(obs))
    
    bias = np.nanmean(dif_spectrum_abs)

    std = np.std(dif_spectrum_abs)

    # mask out padding
    dif_spectrum_abs[dif_spectrum_abs==0] = np.nan


    plt.figure(figsize=(20,5))
    plt.plot([92,92],[0,15],'k',linewidth=5,linestyle='--')
    plt.plot([183,183],[0,15],'k',linewidth=5,linestyle='--')
    plt.plot([274,274],[0,15],'k',linewidth=5,linestyle='--')
    plt.contour( 100*dif_spectrum_abs,levels=[0],linewidths=2,colors='gray')
   
    # calculate proportion of energy contained in most energetic parts of spectrum
    sorted_spectrum = np.sort(np.ravel(dif_spectrum_abs))
    prop = int(np.floor(366*16/10))
    top_half[i]     = np.nansum(np.square(sorted_spectrum[0:prop])) / np.nansum(np.square(sorted_spectrum))

    # now apply mask...
    #dif_spectrum_abs[np.abs(dif_spectrum_rel)<threshold/100] = np.nan
    #dif_spectrum_abs[np.abs(dif_spectrum_abs) < 2*rms[i]] = np.nan
    #dif_spectrum_abs[np.abs(dif_spectrum_abs) < std] = np.nan
    #dif_spectrum_abs[np.abs(dif_spectrum_abs) < 10000*bias] = np.nan
    print(np.shape(dif_spectrum_abs))
    print(np.shape(obs_std))
    #dif_spectrum_abs[np.square(dif_spectrum_abs)<np.square(5*obs_std)] = np.nan
    #dif_spectrum_abs = np.divide(dif_spectrum_abs,obs_std)

    pos = np.copy(dif_spectrum_abs)
    pos[pos<0] = 0
    pos[np.isnan(pos)] = 0
    major_pos[i] = np.count_nonzero(pos)

    neg = np.copy(dif_spectrum_abs)
    neg[neg>0] = 0
    neg[np.isnan(neg)] = 0
    major_neg[i] = np.count_nonzero(neg)

    pos = np.copy(dif_spectrum_abs)
    pos[pos<0.05] = 0
    pos[np.isnan(pos)] = 0
    five_pos[i] = np.count_nonzero(pos)

    neg = np.copy(dif_spectrum_abs)
    neg[neg>-0.05] = 0
    neg[np.isnan(neg)] = 0
    five_neg[i] = np.count_nonzero(neg)

    pos = np.copy(dif_spectrum_abs)
    pos[pos<0.1] = 0
    pos[np.isnan(pos)] = 0
    ten_pos[i] = np.count_nonzero(pos)

    neg = np.copy(dif_spectrum_abs)
    neg[neg>-0.1] = 0
    neg[np.isnan(neg)] = 0
    ten_neg[i] = np.count_nonzero(neg)

    pos = np.copy(dif_spectrum_abs)
    pos[pos<0.15] = 0
    pos[np.isnan(pos)] = 0
    fift_pos[i] = np.count_nonzero(pos)

    neg = np.copy(dif_spectrum_abs)
    neg[neg>-0.15] = 0
    neg[np.isnan(neg)] = 0
    fift_neg[i] = np.count_nonzero(neg)
    
    plt.contourf(100*dif_spectrum_abs,levels=[-12,-8,-4,0,4,8,12],colors=custom,extend='both')
    #plt.contourf(10000*obs,[0,100,100000],colors=('w',season_colors[3],season_colors[3]),alpha=0.5)
    #plt.contourf(10000*mod,[0,100,100000],colors=('w',season_colors[1],season_colors[1]),alpha=0.5)
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
    plt.savefig('{}_threshold'.format(st))

lw = 3

plt.figure(figsize=(20,10))
#plt.bar(stations,100*np.divide(total_pos,denom),facecolor='none',edgecolor='gray',linewidth=lw)
plt.bar(stations,100*np.divide(major_pos,denom),facecolor=custom[4],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(five_pos,denom),facecolor=custom[5],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(ten_pos,denom),facecolor=custom[6],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(fift_pos,denom),facecolor=custom[7],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(total_pos,denom),facecolor='none',edgecolor='gray',linewidth=lw)

#plt.bar(stations,100*np.divide(-total_neg,denom),facecolor='none',edgecolor='gray',linewidth=lw)
plt.bar(stations,100*np.divide(-major_neg,denom),facecolor=custom[3],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(-five_neg,denom),facecolor=custom[2],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(-ten_neg,denom),facecolor=custom[1],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(-fift_neg,denom),facecolor=custom[1],edgecolor='none',linewidth=lw)
plt.bar(stations,100*np.divide(-total_neg,denom),facecolor='none',edgecolor='gray',linewidth=lw)

plt.yticks([-50,-25,0,25,50],['50 %','','0','','50 %'],fontsize=30)
plt.xticks(fontsize=30,rotation='90')
plt.ylim(-75,75)
plt.grid()
plt.tight_layout()
plt.savefig('areas_{}.png'.format(threshold))

plt.figure(figsize=(20,10))
plt.bar(stations,100*avg_dif)
plt.bar(stations,100*std_dif)
plt.xticks(fontsize=30,rotation='90')
plt.grid()
plt.tight_layout()
plt.savefig('metric.png')

plt.figure(figsize=(20,10))
plt.bar(stations,100*top_half)
plt.xticks(fontsize=30,rotation='90')
plt.grid()
plt.tight_layout()
plt.savefig('metric2.png') 

plt.figure(figsize=(20,10))
plt.scatter(100*rms_t,100*rms)
#plt.xticks(fontsize=30,rotation='90')
plt.grid()
plt.tight_layout()
plt.savefig('rmse.png')


plt.figure(figsize=(20,10))
plt.scatter(corr_t,corr)
#plt.xticks(fontsize=30,rotation='90')
plt.grid()
plt.tight_layout()
plt.savefig('correlation.png')

