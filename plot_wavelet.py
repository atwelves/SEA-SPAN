# new script to validate sea level from NEMO 4.2.1 tuning
# compare against tide gauge data

### --- Station indices in NEMO-nordic --- ###

# *Station               *i-index       * j-index
Aarhus             =    [      518     ,       460     ]
Degerby            =    [      883     ,       692     ]
Forsmark           =    [      806     ,       715     ]
Frederikshavn      =    [      530     ,       537     ]
GoteborgTorshamnen =    [      574     ,       552     ]
Hamina             =    [     1128     ,       724     ]
Hanko              =    [      977     ,       680     ]
Helsinki           =    [     1048     ,       699     ]
Hornbaek           =    [      598     ,       457     ]
KalixStoron        =    [      981     ,      1032     ]
Kemi               =    [     1032     ,      1031     ]
Kobenhavn          =    [      604     ,       433     ]
Kronstadt          =    [     1221     ,       689     ]
Kungsholmsfort     =    [      710     ,       457     ]
Kungsvik           =    [      551     ,       630     ]
LandsortNorra      =    [      792     ,       617     ]
OlandsNorraUdde    =    [      765     ,       532     ] 
Parnu              =    [     1031     ,       593     ]
Pietarsaari        =    [      966     ,       913     ]
Raahe              =    [     1028     ,       970     ]
Ronne              =    [      678     ,       397     ]
Sassnitz           =    [      641     ,       361     ]
Sillamae           =    [     1149     ,       656     ]
Skagen             =    [      531     ,       553     ]
StPetersburg       =    [     1236     ,       687     ]
Stockholm          =    [      801     ,       651     ]
Tallinn            =    [     1041     ,       658     ]
Viken              =    [      602     ,       459     ]
Visby              =    [      807     ,       548     ]

### ------ ###

# import python ftp libraries
import ftplib
from ftplib import FTP

import os
import xarray as xr
import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np

import scipy
from scipy.fft import fft, ifft
from scipy import stats
from scipy import signal

import calendar
from calendar import monthrange

import matplotlib as mpl
import matplotlib.pyplot as plt

import cmocean

# cycle over all stations

# specify station
station       = 'Kemi'
station_locat =  Kemi

# months of interest
months = ['202001','202002','202003','202004','202005','202006','202007','202008','202009','202010','202011','202012','202101','202102','202103','202104','202105','202106','202107','202108','202109','202110','202111','202112','202201','202202','202203','202204','202205','202206','202207','202208','202209','202210','202211','202212','202301','202302','202303','202304','202305','202306','202307','202308','202309','202310','202311','202312']
#months = ['202105','202106','202107','202108']#,'202109','202110','202111','202112']
### --- Tide gauges --- ###

# tide gauge prefix: 'BO' or 'NO'
nobo = 'BO'
# check if folder exists
newpath = 'station_data_{}'.format(station)
# if not, download from ftp
if not os.path.exists(newpath):
    print('got here')
    os.makedirs(newpath)
    os.chdir(newpath)

    ### --- CMEMS FTP --- ###

    ftp = FTP('nrt.cmems-du.eu')
    ftp.login(user='atwelves1',passwd='Shamzi69!')
    # navigate to tide gauge data
    ftp.cwd('Core/INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032')
    ftp.cwd('cmems_obs-ins_bal_phybgcwav_mynrt_na_irr/monthly/TG')
    for month_stamp in months:
        ftp.cwd(month_stamp)
        fname = '{}_TS_TG_{}_{}.nc'.format(nobo,station,month_stamp)
        # download month of data
        with open(fname,'wb') as fp:
            ftp.retrbinary('RETR {}'.format(fname),fp.write)
        ftp.cwd('..')
    # quit ftp server
    ftp.quit()
    os.chdir('..')    
### ------ ###

print('GOT HERE')

print(os.getcwd())
os.chdir('station_data_{}'.format(station))
print(os.getcwd())
slv = np.zeros((24))
tms = np.arange('2020-10-01','2020-10-02', dtype='datetime64[h]')
#tms = np.arange('2021-03-31','2021-04-01', dtype='datetime64[h]')
print(tms)
days=-1
for month_stamp in months:
    year_str  = month_stamp[0:4]
    year_int  = int(year_str)
    month_str = month_stamp[4:6]
    month_int = int(month_str)
    m_len = monthrange(year_int,month_int)
    days  = days + m_len[1]
    fname = '{}_TS_TG_{}_{}.nc'.format(nobo,station,month_stamp)
    print(fname)
    ds = xr.open_dataset(fname)
    slv_in = np.squeeze(ds.SLEV.values)
    print(np.shape(slv_in))
    slv_qc = np.squeeze(ds.SLEV_QC.values)
    print(np.shape(slv_qc))
    slv_in[slv_qc>1] = np.nan
    print(np.shape(slv_qc))
    tms_in = np.squeeze(ds.TIME.values)
    # quality control
    tms_in = tms_in[~np.isnan(slv_in)]
    slv_in = slv_in[~np.isnan(slv_in)]
    slv = np.concatenate((slv,slv_in))
    tms = np.concatenate((tms,tms_in))
    print(days)
os.chdir('..')

print('number of days')
print(days)

slv = slv - np.nanmean(slv)

### --- --- ###


# declare array to hold model outputs
count= 0
expt_list = ['00']#,'00r','00s','10','sm2','glp2','sm3','map','ctrl']

# declare array to hold model outputs
ssh_locat = np.zeros((2,3648))
#timer     = np.zeros((2,3648))


#start = 181 # Start in April

#for expt in expt_list:
#    print(expt)
#    expt_name          = 'TEST_{}'.format(expt)
    # check if file already exists
#    extracted_name     = '{}/{}_{}.nc'.format(expt_name,expt_name,station)
#    if not os.path.exists(extracted_name):
#        ds                 = xr.open_dataset('{}/{}_1h_20201002_20211231_ssh_T.nc'.format(expt_name,expt_name))
#        ds_locat           = ds.isel(y=station_locat[1],x=station_locat[0])
#        ssh_in             = ds_locat.zos.values
#        ssh_locat[count,:] = ssh_in[start*24:(start+days)*24]
        #timer_in           = ds.time_instant.values
        #timer[count,:]     = timer_in[:364*24]
        #print(timer[0,0])
#        ncfile = Dataset(extracted_name,mode='w')
#        t                  = ncfile.createDimension('t',np.size(ssh_locat[count,:]))
#        sim_ssh            = ncfile.createVariable('sim_ssh',np.float32,('t'))
#        sim_ssh[:]         = ssh_locat[count,:]
#        print(ncfile)
#        ncfile.close()
#    else:
#        ds_locat           = xr.open_dataset(extracted_name)
#        print(count)
#        if (count<30):
#            ssh_locat[count,:] = ds_locat.sim_ssh.values
#        else:
#            ssh_locat[count,:] = ds_locat.ssh.values[:days*24]
#        ssh_locat[count,:] = ssh_locat[count,:] - np.nanmean(ssh_locat[count,:])
#        count              = count+1

#del ssh_locat

ssh_locat = np.zeros((len(slv)))
ds = xr.open_dataset('calval24_{}/calval24_{}.nc'.format(station,station))
ssh_locat[8760-61*24:2*8760-144] = np.squeeze(ds.SSH_inst.values)
ssh_locat = ssh_locat - np.nanmean(ssh_locat)

ssh_locat[ssh_locat==0] = np.nan

slv = np.squeeze(slv)
#slv = np.squeeze(slv[:10080])
ssh_locat = np.squeeze(ssh_locat)
#ssh_locat = np.squeeze(ssh_locat[:10080])

#plt.figure()
#plt.plot(ssh_locat,label='model')
#plt.plot(slv,label='obs')
#plt.legend()
#plt.show()

ssh_tr = signal.cwt(ssh_locat,signal.ricker,np.arange(1,720))
print('shape:{}'.format(np.shape(ssh_tr)))
slv_tr = signal.cwt(slv,signal.ricker,np.arange(1,720))
print(np.shape(slv_tr))

#ssh_anom = ssh_locat - slv

slv_tr1 = np.nansum(np.abs(slv_tr[:48,:]),0)
ssh_tr1 = np.nansum(np.abs(ssh_tr[:48,:]),0)
slv_tr2 = np.nansum(np.abs(slv_tr[48:,:]),0)
ssh_tr2 = np.nansum(np.abs(ssh_tr[48:,:]),0)

average_slv1 = []
average_ssh1 = []
average_slv2 = []
average_ssh2 = []
average_ssh = []
average_slv = []
window = 720
for ind in range(len(slv_tr1)-window+1):
#    average_anom.append(np.mean(ssh_anom[ind:ind+window]))
    average_ssh.append(np.mean(ssh_locat[ind:ind+window]))
    average_slv.append(np.mean(slv[ind:ind+window]))
    average_slv1.append(np.mean(slv_tr1[ind:ind+window]))
    average_ssh1.append(np.mean(ssh_tr1[ind:ind+window]))
    average_slv2.append(np.mean(slv_tr2[ind:ind+window]))
    average_ssh2.append(np.mean(ssh_tr2[ind:ind+window]))

#plt.figure()
#plt.scatter(average_slv2,np.abs(average_anom))
#plt.show()

plt.figure(figsize=(40,15))
#plt.plot(average_slv1/sum(average_slv1),color='b')
#plt.plot(average_ssh1/sum(average_slv1),color='b',linestyle=':')
plt.plot(average_slv1,color='g',linewidth=5,label='tide gauge')
plt.plot(np.linspace(8760-61*24,2*8760-144-720,8760+61*24-144-720),average_ssh1[8760-61*24:2*8760-144-720],color='b',linewidth=5,linestyle=':',label='nemo')
plt.xticks(np.linspace(0,8760*4,5)-window/2,['                                         2020',
                                             '                                         2021',
                                             '                                         2022',
                                             '                                         2023'
                                             ,''],fontsize=40)
plt.yticks(fontsize=40)
plt.ylabel('Sum of wavelet components, monthly running average',fontsize=30)
plt.grid()
plt.legend(fontsize=40)
plt.title('Mesoscale (<2 days) sea level variability, {}'.format(station),fontsize=40)
plt.savefig('ssh_mesoscale_{}.png'.format(station))

plt.figure(figsize=(40,15))
#plt.plot(average_slv1/sum(average_slv1),color='b')
#plt.plot(average_ssh1/sum(average_slv1),color='b',linestyle=':')
plt.plot(average_slv2,color='g',linewidth=5,label='tide gauge')
plt.plot(np.linspace(8760-61*24,2*8760-144-1440,8760+61*24-144-1440),average_ssh2[8760-61*24:2*8760-144-1440],color='b',linewidth=5,linestyle=':',label='nemo')
plt.xticks(np.linspace(0,8760*4,5)-window/2,['                                         2020',
                                             '                                         2021',
                                             '                                         2022',
                                             '                                         2023'
                                             ,''],fontsize=40)
plt.yticks(fontsize=40)
plt.ylabel('Sum of wavelet components, monthly running average',fontsize=30)
plt.grid()
plt.legend(fontsize=40)
plt.title('Synoptic (>2 days, <1 month) sea level variability, {}'.format(station),fontsize=40)
plt.savefig('ssh_synoptic_{}.png'.format(station))

plt.figure(figsize=(40,15))
#plt.plot(average_slv1/sum(average_slv1),color='b')
#plt.plot(average_ssh1/sum(average_slv1),color='b',linestyle=':')
plt.plot(slv,color='g',linewidth=5)
plt.plot(np.linspace(8760-61*24,2*8760-144-720,8760+61*24-144-720),ssh_locat[8760-61*24:2*8760-144-720],color='b',linewidth=5,linestyle=':')
plt.xticks(np.linspace(0,8760*4,5),['                                         2020',
                                             '                                         2021',
                                             '                                         2022',
                                             '                                         2023'
                                             ,''],fontsize=40)
plt.yticks(fontsize=40)
plt.ylabel('Sea level (m)',fontsize=40)
plt.xlim(24*365,24*730)
plt.grid()
plt.title('{} sea level time series'.format(station),fontsize=40)
plt.savefig('ssh_comparison_{}.png'.format(station))

### ------ ###

#ssh_transform = signal.cwt(np.squeeze(ssh_locat[0,:]),signal.ricker,np.arange(1,240))
ssh_transform = signal.cwt(ssh_locat,signal.ricker,np.arange(1,240))

ssh_trcm = np.cumsum(ssh_transform,0)

print(np.shape(slv))
print(np.shape(np.nansum(ssh_transform[1:48,:],0)))

ssh_transform = np.nansum(np.abs(ssh_transform[0:24,:]),0)
#mv_avg = ssh_transform[0:-9]+ssh_transform[1:-8]+ssh_transform[2:-7]+ssh_transform[3:-6]+ssh_transform[4:-5]+ssh_transform[5:-4]+ssh_transform[6:-3]+ssh_transform[7:-2]+ssh_transform[8:-1]+ssh_transform[9:]
#mv_avg = mv_avg/10
lo_lim
up_lim=720
#mv_avg = np.zeros((np.size(np.squeeze(ssh_locat[0,:]))-up_lim))
mv_avg = np.zeros((np.size(ssh_locat)-up_lim-lo_lim))
for it in range(0,up_lim):
    mv_avg = mv_avg + ssh_transform[it:-up_lim+it]

mod1_avg=mv_avg/up_lim



#ssh_transform = signal.cwt(np.squeeze(ssh_locat[1,:]),signal.ricker,np.arange(1,240))
ssh_transform = signal.cwt(ssh_locat,signal.ricker,np.arange(1,240))

ssh_trcm = np.cumsum(ssh_transform,0)

print(np.shape(slv))
print(np.shape(np.nansum(ssh_transform[1:48,:],0)))

ssh_transform = np.nansum(np.abs(ssh_transform[0:24,:]),0)
#mv_avg = ssh_transform[0:-9]+ssh_transform[1:-8]+ssh_transform[2:-7]+ssh_transform[3:-6]+ssh_transform[4:-5]+ssh_transform[5:-4]+ssh_transform[6:-3]+ssh_transform[7:-2]+ssh_transform[8:-1]+ssh_transform[9:]
#mv_avg = mv_avg/10
up_lim=720
#mv_avg = np.zeros((np.size(np.squeeze(ssh_locat[0,:]))-up_lim))
mv_avg = np.zeros((np.size(ssh_locat)-up_lim))
for it in range(0,up_lim):
    mv_avg = mv_avg + ssh_transform[it:-up_lim+it]

mod2_avg=mv_avg/up_lim

ssh_transform = signal.cwt(slv,signal.ricker,np.arange(1,240))

ssh_trcm = np.cumsum(ssh_transform,0)

#plt.figure()
#pcol = plt.contourf(np.log10(np.abs(ssh_transform)),np.linspace(-3,1,11),cmap='cmo.solar')
#pcol = plt.pcolormesh(np.abs(ssh_trcm),cmap='cmo.solar')
#plt.xticks(np.cumsum(np.array([0,31,28,31,30,31,30,31,31,30,31,30,31])*31),['J','F','M','A','M','J','J','A','S','O','N','D',''])
#plt.xlabel('Date')
#plt.yticks(np.linspace(0,240,6),['0','2','4','6','8','10'])
#plt.ylabel('Days')
#plt.ylim(0,24)
#cbar=plt.colorbar(pcol)
#plt.clim(0,200)
#plt.show()

print(np.shape(slv))
print(np.shape(np.nansum(ssh_transform[1:48,:],0)))

ssh_transform = np.nansum(np.abs(ssh_transform[0:24,:]),0)
#mv_avg = ssh_transform[0:-9]+ssh_transform[1:-8]+ssh_transform[2:-7]+ssh_transform[3:-6]+ssh_transform[4:-5]+ssh_transform[5:-4]+ssh_transform[6:-3]+ssh_transform[7:-2]+ssh_transform[8:-1]+ssh_transform[9:]
#mv_avg = mv_avg/10
up_lim=720
mv_avg = np.zeros((np.size(slv)-up_lim))
ts_avg = np.zeros((np.size(slv)-up_lim))
for it in range(0,up_lim):
    mv_avg = mv_avg + ssh_transform[it:-up_lim+it]
    ts_avg = ts_avg + slv[it:-up_lim+it]

obs_avg=mv_avg/up_lim
ts_avg = ts_avg/up_lim

plt.figure(figsize=(60,10))
#plin = plt.plot(slv)
plin = plt.plot(np.arange(24*90,24*90+np.size(mod1_avg)),np.log10(mod1_avg),linewidth=1)
plin = plt.plot(np.arange(24*90,24*90+np.size(mod2_avg)),np.log10(mod2_avg),linewidth=1)
plin = plt.plot(np.arange(0,np.size(obs_avg)),np.log10(obs_avg),linewidth=1)
plt.xticks(np.cumsum(np.array([0,31,28,31,30,31,30,31,31,30,31,30,31])*31),['J','F','M','A','M','J','J','A','S','O','N','D',''])
plt.xlabel('Date')
plt.savefig('wavelet_{}.png'.format(station))

plt.figure(figsize=(60,10))
plin = plt.plot(np.arange(0,np.size(ts_avg)),ts_avg,linewidth=1)
plin = plt.plot(np.arange(24*90,24*90+np.size(np.squeeze(ssh_locat[0,:]))),ssh_locat[0,:]-slv[24*90:24*90+np.size(np.squeeze(ssh_locat[0,:]))],linewidth=1)
plin = plt.plot(np.arange(24*90,24*90+np.size(np.squeeze(ssh_locat[0,:]))),ssh_locat[0,:],linewidth=1)
plin = plt.plot(np.arange(24*90,24*90+np.size(np.squeeze(ssh_locat[1,:]))),ssh_locat[1,:],linewidth=1)
plin = plt.plot(np.arange(0,np.size(slv)),slv,linewidth=1)
plt.xticks(np.cumsum(np.array([0,31,28,31,30,31,30,31,31,30,31,30,31])*31),['J','F','M','A','M','J','J','A','S','O','N','D',''])
plt.xlabel('Date')
plt.savefig('tseries_{}.png'.format(station))


