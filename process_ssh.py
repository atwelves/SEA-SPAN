### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SEA-SPAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### --- DESCRIPTION --------------------------------------------------- ###
# Script to process sea level from tide gauges and model outputs.
# 18.09.2024 | 
# A G Twelves
### ------------------------------------------------------------------- ###

### --- GLOBAL VARIABLES ---------------------------------------------- ###
#           |
#obs/mod
setting     = 'mod'
#
nemo_name   = "UNFILTERED"    # Name of NEMO model output file
nemo_freq   = "1h"
#nemo_freq   = "1d" 
#           |
#tide_name   =       # Tide gauge file naming convention
#           |
#output_name =       # Output netcdf naming convention
#           |
lat_min     = 65.6      # Southern edge of region of interest     
lat_max     = 80      # Northern edge of region of interest
lon_min     = 0      # Western edge of region of interest
lon_max     = 50      # Eastern edge of region of interest
#           |
year_min    = 2019      # First year of interest
year_max    = 2020      # Last year of interest
month_min   = 10      # First month of interest
month_max   = 9      # Last month of interest
#           |   
freq_mode   =  'meso'
if (freq_mode == 'meso'):
    scale_min   =  0.0      # Shortest timescale of interest
    scale_max   =  3.5      # Longest timescale of interest
elif (freq_mode == 'syno'): 
    scale_min   =  4.0      # Shortest timescale of interest
    scale_max   =  7.5      # Longest timescale of interest
#           |
time_res    = 24      # Time resolution of transform
freq_res    = 8     # Number of scales to consider
#           |
qc_pass     = 1      # Quality control flags deemed good
qc_tol      = 2500      # Number of bad flagged values allowed
#           |
#basis_fn    = "Cosine"      # Basis function for transform
#basis_fn    = 'Gaussian'
basis_fn   = 'Ricker'
#           |
centre_ssh  = 1      # Subtract mean from time series (1) or not (0)
#           { 1
#             0 }
#           |
### ------------------------------------------------------------------ ###

### --- Load modules ------------------- ###

import xarray as xr
import numpy as np
from pathlib import Path
import hashlib
import pandas as pd
import netCDF4 as nc
from netCDF4 import Dataset
import csv

### ------------------------------------ ###

### --- Top level functions --------------------------------------------- ###

### --- Set bounds in latitude, longitude, time, and frequency domain --- ###
# GLOBAL VARIABLES:
# lat_min
# lat_max
# lon_min
# lon_max
# year_min
# month_min
def within_bounds(source_dir):
    path = Path(source_dir)                              # path of file to find tide gauge data
    count = 0                                            # count number of files
    extension = "{}{:02d}.nc".format(year_min,month_min) # generate file names
    # loop over all matching filenames
    for filename in path.glob(f"*{extension}"):
        ds                 = xr.open_dataset(filename)
        tg_lat             = ds.LATITUDE.values          # read station latitude 
        tg_lon             = ds.LONGITUDE.values         # read station longitude
        # check if station is within bounded area
        if(lat_min < tg_lat < lat_max and lon_min < tg_lon < lon_max):
            station_name = str(filename)                 # extract station name from filename
            station_name = station_name[:-10]            # remove time extension from name
            # check if station has good quality data
            if (quality_pass(station_name) == 1):
                if(count>0):
                    station_list.append(station_name)    # add station to list
                else:
                    station_list = [station_name]        # initialise list
                count = count + 1
    return (station_list)
### --------------------------------------------------------------------- ###

### --- Read in tide gauge data ----------------------------------------- ###
# GLOBAL VARIABLES:
# tide_name
def read_the_tides(station_id):
#    ds                    = xr.open_mfdataset('{}*'.format(station_id))#,year_min))
#    station_data          = ds.SLEV.values           # read sea level from tide gauge
    with open('{}_hourlystats_N2000.txt'.format(station_id), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ',quotechar='|')
        date_time = [ row[1] for row in reader ]
    dates = np.asarray(date_time,dtype='str')
    start_date = np.argwhere(dates=='20191001T000000')
    end_date   = np.argwhere(dates=='20201001T000000')
    start_date = int(start_date)
    end_date   = int(end_date)
    with open('{}_hourlystats_N2000.txt'.format(station_id), newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ',quotechar='|')
        ssh = [ row[2] for row in reader ]
    station_data           = np.asarray(ssh,dtype=np.float32)
    station_data = np.squeeze(station_data[start_date:end_date+1])
    #station_data          = sub_marine(station_data) # subtract mean sea level
    #convert units!
    station_data = station_data/1000
    return station_data
### --------------------------------------------------------------------- ###

### --- Extract time series from point in NEMO output closest to station  ###
# GLOBAL VARIABLES:
# nemo_name
# nemo_freq
def sea_extract(station_id):
    ds        = xr.open_mfdataset('{}/full_results/{}_{}_*'.format(nemo_name,nemo_name,nemo_freq))
    mod_ind   = match_point(station_id)  # find point in model domain closest to tide gauge
    sl        = ds.isel(y=mod_ind[0], x=mod_ind[1]) # extract this point from the model output
    ssh_mod   = sl.zos                       # read in sea level from model output
    ssh_mod_mean = ssh_mod.mean(dim='time_counter')
    ssh_mod   = ssh_mod - ssh_mod_mean
    ssh_mod   = ssh_mod#.values
    #ssh_mod   = sub_marine(ssh_mod)                 # subtract mean level
    plat      = ds.nav_lat.isel(y=mod_ind[0], x=mod_ind[1]).values
    plon      = ds.nav_lon.isel(y=mod_ind[0], x=mod_ind[1]).values
    print(plat)
    print(plon)
    print(mod_ind[0])
    print(mod_ind[1])
    return ssh_mod
### --------------------------------------------------------------------- ###

### --- Sum convolution to complete transform --------------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def do_transform(time_series):
    # for Fourier transform...
    if(basis_fn=="Cosine"): 
        raw_transform
    # for wavelet transform...
    else:
        print('apply Ricker transform')
        convolution = convolute_it(time_series)                            # convolute time series with wavelet basis
        convolv_array = convolution[0]                                     # convolved array
        normal_array  = convolution[1]                                     # normalisation factor
        raw_transform = np.nansum(convolv_array,2)                         # integrate over length of time series
        #normal_tot    = np.nanmin(np.nansum(np.nansum(normal_array,2),0))
        #raw_transform = np.divide(raw_transform,normal_tot)                # normalise wavelet transform
        raw_transform = np.square(raw_transform)                           # spectrum of wavelet transform
        raw_transform = np.transpose(raw_transform)                        # transpose
    transform = sort_bins(raw_transform)                                   # reduce time resolution of output
    return transform
### --------------------------------------------------------------------- ###

### --- Sum convolution to complete transform --------------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def do_transform_mod(time_series):
    # for Fourier transform...
    if(basis_fn=="Cosine"):
        raw_transform
    # for wavelet transform...
    else:
        print('apply Ricker transform')
        raw_transform = convolute_mod(time_series)                            # convolute time series with wavelet basis
        #raw_transform = convolv_array.sum(dim='time_counter')                         # integrate over length of time series
        #normal_tot    = np.nanmin(np.nansum(np.nansum(normal_array,2),0))
        #raw_transform = np.divide(raw_transform,normal_tot)                # normalise wavelet transform
        #raw_transform = np.square(raw_transform)                           # spectrum of wavelet transform
        #raw_transform = raw_transform.transpose
        #print(raw_transform)
        #raw_transform.to_netcdf('transform_Rauma.nc')
        # transpose
    #transform = sort_bins(raw_transform)                                   # reduce time resolution of output
    return
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### --- Secondary functions --------------------------------------------- ###

### --- Read tide gauge quality control flags --------------------------- ###
# GLOBAL VARIABLES:
# tide_name
# qc_pass
# qc_tol
def quality_pass(station_id):
    ds                 = xr.open_mfdataset('{}*'.format(station_id))#,year_min))
    qc_flags           = ds.SLEV_QC.values                      # read quality control values
    qc_bad             = np.count_nonzero(qc_flags !=  qc_pass) # count bad quality data
    # Set threshold for number of 'bad' values
    if(qc_bad > qc_tol):
        station_filter = 0 
    else:
        station_filter = 1
    xr.Dataset.close(ds)
    return station_filter
### --------------------------------------------------------------------- ### 

### --- Match NEMO grid points to tide gauge locations ------------------ ###
def match_point(station_id):
    mod_lat = get_grid("lat")                                # extract latitude axis from model domain
    mod_lon = get_grid("lon")                                # extract longitude axis from model domain
    stat_lat = get_point(station_id,"lat")                   # extract latitude coordinate for station
    stat_lon = get_point(station_id,"lon")                   # extract longitude coordinate for station
    # calculate distance from each model point to station
    diff_latlon = np.square(mod_lat-stat_lat)+np.square((mod_lon-stat_lon)*np.cos(stat_lat*np.pi/180))
    ds        = xr.open_mfdataset('{}/{}_{}_*'.format(nemo_name,nemo_name,nemo_freq))
    sl        = ds.isel(time_counter=0)                                 # extract initial time slice
    ssh_mod   = sl.zos.values                                           # extract initial sea level field
    diff_latlon[ssh_mod==0] = np.nan                                    # disregard land pointd
    diff_latlon[np.isnan(diff_latlon)]=1e20                             # convert NaNs
    mod_ind = np.unravel_index(diff_latlon.argmin(), diff_latlon.shape) # find model index closest to station
    return mod_ind
### --------------------------------------------------------------------- ###

### --- Convolute basis functions with time series ---------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def convolute_it(time_series):
    sp = make_space(time_series)         
    # set grid for transform
    t_space = sp[0]       
    s_space = sp[1]
    tiled_series = np.tile(time_series,(np.size(s_space),1))                           # duplicate time series over s axis
#    tiled_series = np.transpose(tiled_series)                                      
    convolv_array = np.zeros((np.size(t_space),np.size(s_space),np.size(t_space))) # declare array for convolution
    norm_factor   = np.zeros((np.size(t_space),np.size(s_space),np.size(t_space))) # declare array for normalisation
    # Loop over time dimension for wavelet centre
    for tau in range(0,np.size(t_space)):
        #print(tau)
        wavelets = make_waves(t_space,s_space,tau)                                 # generate wavelets
        convolv_array[tau,:,:] = np.multiply(tiled_series,wavelets)                # convolute time series with wavelets
        # normalisation - is this correct??  
        norm_factor[tau,:,:] = wavelets                                            # normalisation factor
        #convolv_array = np.divide(convolve_array,norm_factor)
    return convolv_array, norm_factor
### --------------------------------------------------------------------- ###

### --- Convolute basis functions with time series ---------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def convolute_mod(time_series):
    time_coords = np.arange('2019-10-01T00','2020-10-01T00',dtype='datetime64[h]')
    sp = make_space(time_series)
    t_space = sp[0]
    s_space = sp[1]
    #tiled_series = time_series * xr.DataArray(np.ones([np.size(t_space),np.size(s_space)]), dims=("tau_dim", "scale_dim"))
    #tiled_series.to_netcdf('tiled_series.nc')
    # make wavelets
    
    if (setting=='mod'):
        series = xr.open_dataset('time_series/data/BO_TS_TG_{}_processed_UNFILTERED.nc'.format(station))
        series = series.ssh_mod.values
    elif (setting=='obs'):
        series = np.copy(time_series)

    timer        = xr.DataArray(t_space,[("time_counter",t_space)])
    tiled_timer  = timer * xr.DataArray(np.ones([np.size(t_space),np.size(s_space)]), dims=("tau_dim","scale_dim"))

    tau          = xr.DataArray(t_space,[("tau_dim", t_space)])
    tiled_tau    = tau * xr.DataArray(np.ones([np.size(t_space),np.size(s_space)]), dims=("time_counter","scale_dim"))
    
    scale        = xr.DataArray(np.power(2,s_space),[("scale_dim", s_space)])
    tiled_scale  = scale * xr.DataArray(np.ones([np.size(t_space),np.size(t_space)]), dims=("time_counter","tau_dim"))

    # now combine
    gauss = np.exp(-np.square(tiled_timer-tiled_tau)/np.square(tiled_scale)/2)
    gauss.to_netcdf('gauss_test.nc')

    wavelets = (1 - np.square((tiled_timer-tiled_tau)/tiled_scale) ) * gauss
    # normalise
    wavelets = wavelets / tiled_scale
    wavelets = wavelets
        # convert coordinate
    #wavelets['time_counter'] = pd.DatetimeIndex(wavelets['time_counter'].values)
    wavelets['time_counter']  = time_coords
    wavelets.to_netcdf('wavelets.nc')
    # convolute
    scale = scale.values
    del wavelets
    
    wavelets = xr.open_dataarray('wavelets.nc')
    wt = np.zeros((np.size(t_space),np.size(s_space)))

    scale = np.ceil(scale).astype(int)
    t_space = t_space.astype(int)
    # try for loops
    for scale_index in range(0,np.size(s_space)):
        tau_min = t_space[4*scale[scale_index]]
        tau_max = t_space[-1-4*scale[scale_index]]
        #print(tau_min)
        #print(tau_max)
        tau_range = np.squeeze(t_space[tau_min:tau_max])
        for tau_index in tau_range:
            #print(tau_index)
            t_min = tau_index-4*scale[scale_index]
            #print(t_min)
            t_max = tau_index+4*scale[scale_index]
            #print(t_max)
            time_range = np.squeeze(t_space[t_min:t_max])
            wavelets_local = wavelets.isel(time_counter=slice(time_range[0],time_range[-1])).isel(tau_dim=tau_index).isel(scale_dim=scale_index).values
            #print(wavelets_local)
            tiled_series_local = series[time_range[0]:time_range[-1]]
            #print(tiled_series_local)
            transform_local = np.multiply(wavelets_local,tiled_series_local)
            #print(transform_local)
            transform_local = np.nansum(transform_local)
            #transform_local = np.square(transform_local)
            #print(transform_local)
            wt[tau_index,scale_index] = transform_local

    wt = np.square(wt)
    
    # convert scale to frequency for output
    scale = 4*scale

    da = xr.DataArray(data=wt,dims=["time","freq"],coords=dict(time=time_coords,freq=scale))
    da = da.groupby("time.date").mean()
    da['date'] = np.arange('2019-10-01','2020-10-01',dtype='datetime64[D]')
    for m in [10,11,12,1,2,3,4,5,6,7,8,9]:
        dm = da.groupby('date.month')[m]
        dm.to_netcdf('{}_{}_{}_{}.nc'.format(station,setting,freq_mode,m))
    return 
### --------------------------------------------------------------------- ###


### --- Coarsen transform output by binning ----------------------------- ###
# GLOBAL VARIABLES:
# time_res
# freq_res
def sort_bins(unbinned):
    # Change from two dimensions (t_old x f_old) to four dimensions (t_new x t_res x f_new x f_res)
    t_old    = np.size(unbinned,1)                                  # length of old time dimension
    f_old    = np.size(unbinned,0)                                  # length of old scale dimension
    t_new    = np.int(t_old/time_res)                               # length of new time dimension
    #f_new    = np.int(f_old/freq_res)                               # length of new scale dimension
    unbinned = np.reshape(unbinned,(f_old,1,t_new,time_res)) # Convert 2d array to 4d 
    # Sum within each bin to return to two dimensions (t_new x f_new)
    binned   = np.nansum(unbinned,(1,3))
    return binned
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### --- Tertiary functions ---------------------------------------------- ###

### --- Get NEMO grid coordinates --------------------------------------- ###
# GLOBAL VARIABLES:
# nemo_name
# nemo_freq
def get_grid(mod_axis):
    ds        = xr.open_mfdataset('{}/full_results/{}_{}_*'.format(nemo_name,nemo_name,nemo_freq))
    if (mod_axis=="lat"):
        mod_crd = ds.nav_lat.values # extract latitude values from model
    if (mod_axis=="lon"):
        mod_crd = ds.nav_lon.values # extract longitude values from model
    mod_crd[mod_crd==-1] = np.nan   # set land to NaN
    return mod_crd
### --------------------------------------------------------------------- ### 

### --- Get NEMO grid coordinates --------------------------------------- ###
# GLOBAL VARIABLES:
# year_min
# month_min
def get_point(station_id,stat_axis):
    ds      = xr.open_dataset('data/BO_TS_TG_{}_{}{:02d}.nc'.format(station_id,year_min,month_min))
    if (stat_axis=="lat"):
        stat_crd = ds.LATITUDE.values  # extract latitude value for station
    if (stat_axis=="lon"):
        stat_crd = ds.LONGITUDE.values # extract longitude value for station
    return stat_crd
### --------------------------------------------------------------------- ### 

### --- Make spaces for time and scale ---------------------------------- ###
# GLOBAL VARIABLES:
# scale_min
# scale_max
def make_space(time_series):
    t_space = np.linspace(0,np.size(time_series)-1,np.size(time_series)) # generate time axis values
    s_space = np.linspace(scale_min,scale_max,freq_res)     # generate scale axis values
    return(t_space,s_space)
### --------------------------------------------------------------------- ###

### --- Construct wavelet transform basis functions --------------------- ###
# GLOBAL VARIABLES:
# basis_fn
# scale_min
# scale_max
def make_waves(t_space,s_space,tau):                  
    wavelets = np.zeros((np.size(s_space),np.size(t_space)))                          # declare array for wavelets
#   if (basis_fn=="Cosine"):
#        wavelets[:]   =  
    if (basis_fn=="Gaussian"):
        # Loop over different scales
        for s in range(0,np.size(s_space)):                  
            scale         = np.power(2,s_space[s])                                  # use dyadic convention
            wavelets[s,:] = gaussian(t_space,tau,scale)
            #norm_factor   = np.sqrt(np.nansum(np.square(wavelets[s,:])))
            norm_factor = np.sqrt(scale) #wrong!
            wavelets[s,:] = wavelets[s,:] / norm_factor                             # normalisation
    if (basis_fn=="Ricker"):
        # Loop over different scales        
        for s in range(0,np.size(s_space)):
            scale = np.power(2,s_space[s])                                          # use dyadic convention
            wavelets[s,:] = np.multiply((1 - np.square((t_space[:]-tau)/scale)),gaussian(t_space,tau,scale))
            #norm_factor   = np.sqrt(np.nansum(np.square(wavelets[s,:])))
            norm_factor = np.sqrt(scale) #wrong!
            wavelets[s,:] = wavelets[s,:] / norm_factor                             # normalisation
    return wavelets
### --------------------------------------------------------------------- ###

### --- Construct wavelet transform basis functions --------------------- ###
# GLOBAL VARIABLES:
# basis_fn
# scale_min
# scale_max
def make_waves_mod(t_space,s,tau):
    wavelets = np.zeros((np.size(t_space)))                          # declare array for wavelets
#   if (basis_fn=="Cosine"):
#        wavelets[:]   =  
    if (basis_fn=="Gaussian"):
        # Loop over different scales
        for s in range(0,np.size(s_space)):
            scale         = np.power(2,s_space[s])                                  # use dyadic convention
            wavelets[s,:] = gaussian(t_space,tau,scale)
            #norm_factor   = np.sqrt(np.nansum(np.square(wavelets[s,:])))
            norm_factor = np.sqrt(scale)
            wavelets[s,:] = wavelets[s,:] / norm_factor                             # normalisation
    if (basis_fn=="Ricker"):
        # Loop over different scales        
        #for s in range(0,np.size(s_space)):
        scale = np.power(2,s)                                          # use dyadic convention
        wavelets[:] = np.multiply((1 - np.square((t_space[:]-tau)/scale)),gaussian(t_space,tau,scale))
        #norm_factor   = np.sqrt(np.nansum(np.square(wavelets[s,:])))
        norm_factor = np.sqrt(scale)
        wavelets[:] = wavelets[:] / norm_factor                             # normalisation
    return wavelets
### --------------------------------------------------------------------- ###


### --------------------------------------------------------------------- ###

### --- Supplementary functions ----------------------------------------- ###

### --- Construct gaussian function ------------------------------------- ###
def gaussian(t_space,tau,scale):
    gauss = np.exp(-np.square(t_space[:]-tau)/np.square(scale)/2)
    return gauss
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### --- Main --------------------------- ###
fmi_dir ="data" # directory of tide gauge data
print(fmi_dir)
#station_list = within_bounds(fmi_dir) # this extracts a full list of stations within bounds
station_list = ['Rauma']
print("Processing {} tide gauge stations".format(np.size(station_list)))
print(station_list)
for station in station_list:
    # read in tide gauge data
    if (setting == 'obs'):
        tg            = read_the_tides(station)
        tg_transform  = do_transform_mod(tg) # transform obs
    # read in model output
    elif (setting == 'mod'):
        mod           = sea_extract(station)
        mod_transform = do_transform_mod(mod) # transform mod

### ------------------------------------ ###

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SEA-SPAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
