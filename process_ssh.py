### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SEA-SPAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### --- DESCRIPTION --------------------------------------------------- ###
# Script to process sea level from tide gauges and model outputs.
# 18.09.2024 | 
# A G Twelves
### ------------------------------------------------------------------- ###

### --- GLOBAL VARIABLES ---------------------------------------------- ###
#           |
nemo_name   = "WW_ORIGINAL"    # Name of NEMO model output file
nemo_freq   = "1h"
#nemo_freq   = "1d" 
#           |
#tide_name   =       # Tide gauge file naming convention
#           |
#output_name =       # Output netcdf naming convention
#           |
lat_min     = 50      # Southern edge of region of interest     
lat_max     = 80      # Northern edge of region of interest
lon_min     = 0      # Western edge of region of interest
lon_max     = 50      # Eastern edge of region of interest
#           |
year_min    = 2019      # First year of interest
year_max    = 2019      # Last year of interest
month_min   = 1      # First month of interest
month_max   = 4      # Last month of interest
#           |   
scale_min   = 0      # Shortest timescale of interest
scale_max   = 8      # Longest timescale of interest
#           |
time_res    = 24      # Time resolution of transform
freq_res    = 1      # Frequency resolution of transform
#           |
qc_pass     = 1      # Quality control flags deemed good
qc_tol      = 25      # Number of bad flagged values allowed
#           |
#basis_fn    = "Cosine"      # Basis function for transform
#basis_fn   = 'Ricker'
basis_fn   = 'Gaussian'
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
            if (quality_pass(source_dir,station_name) == 1):
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
def read_the_tides(source_dir,station_id):
    ds                    = xr.open_mfdataset('{}/{}*2019*'.format(source_dir,station_id))
    station_data          = ds.SLEV.values           # read sea level from tide gauge
    station_data          = sub_marine(station_data) # subtract mean sea level
    return station_data
### --------------------------------------------------------------------- ###

### --- Extract time series from point in NEMO output closest to station  ###
# GLOBAL VARIABLES:
# nemo_name
# nemo_freq
def sea_extract(source_dir,station_id):
    ds        = xr.open_mfdataset('{}/{}/{}_{}_*'.format(source_dir,nemo_name,nemo_name,nemo_freq))
    mod_ind   = match_point(source_dir,station_id)  # find point in model domain closest to tide gauge
    sl        = ds.isel(y=mod_ind[0], x=mod_ind[1]) # extract this point from the model output
    ssh_mod   = sl.zos.values                       # read in sea level from model output
    ssh_mod   = sub_marine(ssh_mod)                 # subtract mean level
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
        convolution = convolute_it(time_series)                            # convolute time series with wavelet basis
        convolv_array = convolution[0]                                     # convolved array
        normal_array  = convolution[1]                                     # normalisation factor
        raw_transform = np.nansum(convolv_array,2)                         # integrate over length of time series
        raw_transform = np.divide(raw_transform,np.nansum(normal_array,2)) # normalise wavelet transform
        raw_transform = np.square(raw_transform)                           # spectrum of wavelet transform
        raw_transform = np.transpose(raw_transform)                        # transpose
    transform = sort_bins(raw_transform)                                   # reduce time resolution of output
    return transform
### --------------------------------------------------------------------- ###

### --- Write arrays to netcdf ------------------------------------------ ###
# GLOBAL VARIABLES:
# output_name
def write_all(station, time_series, transformed_series, time_series2, transformed_series2):#, mod_series):
    sdiff             = np.squeeze(time_series2) - np.squeeze(time_series) # model sea level anomaly
    # 
    outfile           = '{}_proc.nc'.format(station)
    ncfile            = Dataset(outfile,mode='w')
    t                 = ncfile.createDimension('t',np.size(time_series))
    slev_obs          = ncfile.createVariable('ssh_obs',np.float32,('t'))
    slev_mod          = ncfile.createVariable('ssh_mod',np.float32,('t'))
    slev_diff         = ncfile.createVariable('ssh_diff',np.float32,('t'))
    slev_obs[:]       = time_series[:]
    slev_mod[:]       = time_series2[:]
    slev_diff[:]      = sdiff[:]
    ncfile.close(); print('Dataset is closed')
    # 
    tdiff               = transformed_series2 - transformed_series         # model transform anomaly  
    outfile             = '{}_transform.nc'.format(station)
    ncfile              = Dataset(outfile,mode='w')
    t                   = ncfile.createDimension('t',np.size(transformed_series,1))
    s                   = ncfile.createDimension('s',np.size(transformed_series,0))
    transform_obs       = ncfile.createVariable('obs_transform',np.float32,('s','t'))
    transform_mod       = ncfile.createVariable('mod_transform',np.float32,('s','t'))
    transform_diff      = ncfile.createVariable('diff_transform',np.float32,('s','t'))
    transform_obs[:,:]  = transformed_series[:,:]
    transform_mod[:,:]  = transformed_series2[:,:]
    transform_diff[:,:] = tdiff[:,:]
    ncfile.close(); print('Dataset is closed')
    return
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### --- Secondary functions --------------------------------------------- ###

### --- Read tide gauge quality control flags --------------------------- ###
# GLOBAL VARIABLES:
# tide_name
# qc_pass
# qc_tol
def quality_pass(source_dir,station_id):
    ds                 = xr.open_mfdataset('{}/{}*2019*'.format(source_dir,station_id))
    qc_flags           = ds.SLEV_QC.values
    qc_bad             = np.count_nonzero(qc_flags !=  qc_pass)
    if(qc_bad > qc_tol):
        station_filter = 0
    else:
        station_filter = 1
    xr.Dataset.close(ds)
    return station_filter
### --------------------------------------------------------------------- ### 

### --- Match NEMO grid points to tide gauge locations ------------------ ###
def match_point(source_dir,station_id):
    mod_lat = get_grid(source_dir,"lat")
    mod_lon = get_grid(source_dir,"lon")
    # flatten to 1d
    #i#mod_lat = np.nanmean(mod_lat,1)
    #mod_lon = np.nanmean(mod_lon,0)
    stat_lat = get_point(source_dir,station_id,"lat")
    stat_lon = get_point(source_dir,station_id,"lon")
    diff_latlon = np.square(mod_lat-stat_lat)+np.square((mod_lon-stat_lon)*np.cos(stat_lat*np.pi/180))
    # test for land points
    ds        = xr.open_mfdataset('{}/{}/{}_{}_*'.format(source_dir,nemo_name,nemo_name,nemo_freq))
    sl        = ds.isel(time_counter=0)
    ssh_mod   = sl.zos.values
    diff_latlon[ssh_mod==0] = np.nan
    #diff_lat = np.abs(mod_lat - stat_lat)
    #diff_lon = np.abs(mod_lon - stat_lon)
    #diff_lat[np.isnan(diff_lat)] = 999
    #diff_lon[np.isnan(diff_lon)] = 999
    diff_latlon[np.isnan(diff_latlon)]=1e20
    mod_ind = np.unravel_index(diff_latlon.argmin(), diff_latlon.shape)
    #mod_y = np.argmin(diff_lat)
    #mod_x = np.argmin(diff_lon)
    return mod_ind
### --------------------------------------------------------------------- ###

### --- Subtract mean sea level from time series, removes NaNs ---------- ###
# GLOBAL VARIABLES:
# centre_ssh
def sub_marine(time_series):
    if(centre_ssh==1):
        time_series = time_series - np.nanmean(time_series)
    time_series[np.isnan(time_series)] = 0
    return time_series
### --------------------------------------------------------------------- ###

### --- Convolute basis functions with time series ---------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def convolute_it(time_series):
    print(np.shape(time_series))
    # DECLARE ARRAY:
    sp = make_space(time_series)
    t_space = sp[0]
    s_space = sp[1]
    tiled_series = np.tile(time_series,np.size(s_space))
    tiled_series = np.transpose(tiled_series)
    convolv_array = np.zeros((np.size(t_space),np.size(s_space),np.size(t_space)))
    norm_factor   = np.zeros((np.size(t_space),np.size(s_space),np.size(t_space)))
    for tau in range(0,np.size(t_space)):
        wavelets = make_waves(t_space,s_space,tau)
        print(np.shape(wavelets))
        convolv_array[tau,:,:] = np.multiply(tiled_series,wavelets)
        # normalisation
        norm_factor[tau,:,:] = np.square(wavelets)
        #convolv_array = np.divide(convolve_array,norm_factor)
    return convolv_array, norm_factor
### --------------------------------------------------------------------- ###

### --- Coarsen transform output by binning ----------------------------- ###
# GLOBAL VARIABLES:
# time_res
# freq_res
def sort_bins(unbinned):
    # unbinned has two dimensions, t x f
    # change to four dimensions, tnew x tres x fnew x fres
    t_old    = np.size(unbinned,1)
    f_old    = np.size(unbinned,0)
    t_new    = np.int(t_old/time_res)
    f_new    = np.int(f_old/freq_res)
    # Convert 2d to 4d  
    unbinned = np.reshape(unbinned,(f_new,freq_res,t_new,time_res))
    # Sum within each bin
    binned   = np.nansum(unbinned,(1,3))
    return binned
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### --- Tertiary functions ---------------------------------------------- ###

### --- Get NEMO grid coordinates --------------------------------------- ###
# GLOBAL VARIABLES:
# nemo_name
def get_grid(source_dir,mod_axis):
    ds        = xr.open_mfdataset('{}/{}/{}_{}_*'.format(source_dir,nemo_name,nemo_name,nemo_freq))
    if (mod_axis=="lat"):
        mod_crd = ds.nav_lat.values
    if (mod_axis=="lon"):
        mod_crd = ds.nav_lon.values
    mod_crd[mod_crd==-1] = np.nan
    return mod_crd
### --------------------------------------------------------------------- ### 

### --- Get NEMO grid coordinates --------------------------------------- ###
# GLOBAL VARIABLES:
# year_min
# month_min
def get_point(source_dir,station_id,stat_axis):
    ds      = xr.open_dataset('{}/{}_{}{:02d}.nc'.format(source_dir,station_id,year_min,month_min))
    if (stat_axis=="lat"):
        stat_crd = ds.LATITUDE.values
    if (stat_axis=="lon"):
        stat_crd = ds.LONGITUDE.values
    return stat_crd
### --------------------------------------------------------------------- ### 

### --- Make spaces for time and scale ---------------------------------- ###

def make_space(time_series):
    t_space = np.linspace(0,np.size(time_series),np.size(time_series))
    s_space = np.linspace(scale_min,scale_max,scale_max-scale_min)
    return(t_space,s_space)

### --------------------------------------------------------------------- ###

### --- Construct wavelet transform basis functions --------------------- ###
# GLOBAL VARIABLES:
# basis_fn
# scale_min
# scale_max
def make_waves(t_space,s_space,tau):
    # DECLARE ARRAY
    wavelets = np.zeros((np.size(s_space),np.size(t_space)))
#!!!!!!!!!TEST
#    wavelets = np.random.rand(np.size(s_space),np.size(t_space))
#    if (basis_fn=="Cosine"):
#        wavelets[:]   = 
#    if (basis_fn=="Ricker"):
#        wavelets[:,:] = 
    if (basis_fn=="Gaussian"):
        for s in range(0,np.size(s_space)):
            # use dyadic convention
            scale = np.power(2,s_space[s])
            wavelets[s,:] = np.exp(-np.square(t_space[:]-tau)/np.square(scale))
    return wavelets
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### --- Main --------------------------- ###
fmi_dir ="." # directory of tide gauge data
station_list = within_bounds(fmi_dir) # this extracts a full list of stations within bounds
print("Processing {} tide gauge stations".format(np.size(station_list)))
print(station_list)
for station in station_list:
    # read in tide gauge data
    tg            = read_the_tides(fmi_dir,station)
    tg_transform  = do_transform(tg) # transform obs
    # read in model output
    mod           = sea_extract(fmi_dir,station)
    mod_transform = do_transform(np.expand_dims(mod,axis=1)) # transform mod
    # write out data to netcdf
    write_all(station,tg,tg_transform,mod,mod_transform)
    #write_both("transformed")
    print("{} processed".format(station))

### ------------------------------------ ###

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SEA-SPAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
