### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SEA-SPAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### --- DESCRIPTION --------------------------------------------------- ###
# Script to process sea level from tide gauges and model outputs.
# 18.09.2024 | 
# A G Twelves
### ------------------------------------------------------------------- ###

### --- GLOBAL VARIABLES ---------------------------------------------- ###
#           |
nemo_name   = "NORDIC"    # Name of NEMO model output file
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
basis_fn    = "Cosine"      # Basis function for transform
#           { 'Cosine'
#             'Ricker'
#             'Gaussian' }
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
# year_min
# year_max
# month_min
# month_max
# tide_name
def within_bounds(source_dir):
    path = Path(source_dir)
    count = 0 # count number of files
    #station_list = np.empty((1000),dtype=str)
    # read all unique stations
    extension = "{}{:02d}.nc".format(year_min,month_min)
    for filename in path.glob(f"*{extension}"):
        #ds                 = xr.open_mfdataset('{}/*_{}{:02d}.nc'.format(source_dir,year_min,month_min))    
        ds                 = xr.open_dataset(filename)
        tg_lat             = ds.LATITUDE.values
        tg_lon             = ds.LONGITUDE.values
        tg_time            = ds.TIME.values
        print(tg_time)
        if(lat_min < tg_lat < lat_max and lon_min < tg_lon < lon_max):
            #station_byte        = ds.STATION.values
            #print(station_byte)
            #station_name        = station_byte.astype('UTF-8') 
            #print(type(station_name))
            station_name = str(filename)
            # remove year+month
            station_name = station_name[:-10]
            if (quality_pass(source_dir,station_name) == 1):
                print(station_name)
                if(count>0):
                    station_list.append(station_name)
                else:
                    station_list = [station_name]
                count = count + 1
    #station_list = np.squeeze(station_list[0:count])        
    print(station_list)
    return (station_list)
### --------------------------------------------------------------------- ###

### --- Read in tide gauge data ----------------------------------------- ###
# GLOBAL VARIABLES:
# tide_name
def read_the_tides(source_dir,station_id):
    ds                    = xr.open_mfdataset('{}/{}*2019*'.format(source_dir,station_id))
    print(ds)
    station_data          = ds.SLEV.values # account for qc flag here too?
    station_data          = sub_marine(station_data) 
    return station_data
### --------------------------------------------------------------------- ###

### --- Extract time series from point in NEMO output closest to station  ###
# GLOBAL VARIABLES:
# nemo_name
def sea_extract(station):
    ds        = xr.open_dataset('{}'.format(nemo_name))
    mod_ind   = match_point(stat_lat, stat_lon)
    sl        = ds.isel(j=mod_ind[0], i=mod_ind[1])
    ssh_mod   = sl.ssh.values
    ssh_mod   = sub_marine(ssh_mod)
    return ssh_mod
### --------------------------------------------------------------------- ###

### --- Sum convolution to complete transform --------------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def do_transform(time_series):
    if(basis_fn=="cosine"):
        raw_transform
    else:
        convolv_array = convolute_it(time_series)
        raw_transform = np.nansum(convolv_array,0)
    transform = sort_bins(raw_transform)
    # Normalisation?
    return transform
### --------------------------------------------------------------------- ###

### --- Write arrays to netcdf ------------------------------------------ ###
# GLOBAL VARIABLES:
# output_name
def write_all(station, time_series, transformed_series):#, mod_series):
    #transform_out = np.concat(obs_transform,mod_transform)
    #time_series.to_netcdf('{}'.format(station))
    ##! Need to do this the long way!!
    outfile           = '{}_proc.nc'.format(station)
    ncfile            = Dataset(outfile,mode='w')
    t                 = ncfile.createDimension('t',np.size(time_series))
    slev              = ncfile.createVariable('ssh',np.float32,('t'))
    slev[:]           = time_series[:]
    ncfile.close(); print('Dataset is closed')

    outfile           = '{}_transform.nc'.format(station)
    ncfile            = Dataset(outfile,mode='w')
    t                 = ncfile.createDimension('t',np.size(transformed_series,0))
    s                 = ncfile.createDimension('s',np.size(transformed_series,1))
    slev_transform    = ncfile.createVariable('ssh_transform',np.float32,('t','s'))
    slev_transform[:] = transformed_series[:,:]
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
    print(station_id)
    ds                 = xr.open_mfdataset('{}/{}*2019*'.format(source_dir,station_id))
    print(ds)
    qc_flags           = ds.SLEV_QC.values
    qc_bad             = np.count_nonzero(qc_flags !=  qc_pass)
    print("qc_bad:{}".format(qc_bad))
    if(qc_bad > qc_tol):
        station_filter = 0
    else:
        station_filter = 1
    xr.Dataset.close(ds)
    return station_filter
### --------------------------------------------------------------------- ### 

### --- Match NEMO grid points to tide gauge locations ------------------ ###
def match_point(station):
    mod_lat = get_grid("lat")
    mod_lon = get_grid("lon")
    stat_lat = get_point("lat")
    stat_lon = get_point("lon")
    mod_y = np.argmin(mod_lat - stat_lat)
    mod_x = np.argmin(mod_lon - stat_lon)
    return mod_y, mod_x
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
    # DECLARE ARRAY:
    sp = make_space(time_series)
    t_space = sp[0]
    s_space = sp[1]
    tiled_series = np.tile(time_series,np.size(s_space))
    tiled_series = np.transpose(tiled_series)
    convolv_array = np.zeros((np.size(t_space),np.size(s_space),np.size(t_space)))
    for tau in range(0,np.size(t_space)):
        wavelets = make_waves(t_space,s_space)
        convolv_array[tau,:,:] = np.multiply(tiled_series,wavelets)
    return convolv_array
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
def get_grid(mod_axis):
    ds      = xr.open_dataset(nemo_name)
    if (mod_axis=="lat"):
        mod_crd = ds.nav_lat.values
    if (mod_axis=="lon"):
        mod_crd = ds.nav_lon.values
    return mod_crd
### --------------------------------------------------------------------- ### 

### --- Get NEMO grid coordinates --------------------------------------- ###
# GLOBAL VARIABLES:
def get_point(stat_axis):
    ds      = xr.open_dataset('{}/*_{}{}.nc'.format(source_dir,year_min,month_min))
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
def make_waves(t_space,s_space):
    # DECLARE ARRAY
#    wavelets = np.zeros((np.size(s_space),np.size(t_space)))
#!!!!!!!!!TEST
    wavelets = np.random.rand(np.size(s_space),np.size(t_space))
#    if (basis_fn=="Cosine"):
#        wavelets[:]   = 
#    if (basis_fn=="Ricker"):
#        wavelets[:,:] = 
    if (basis_fn=="Gaussian"):
        for s in range(0,np.size(s_space)):
            # use dyadic convention
            scale = np.power(2,s_space)
            wavelets[s,:] = np.exp(-np.square(tspace[:])/np.square(scale))
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
    #mod           = sea_extract(station)
    #mod_transform = do_transform(mod) # transform mod
    # write out data to netcdf
    write_all(station,tg,tg_transform)
    #write_both("transformed")
    print("{} processed".format(station))

### ------------------------------------ ###

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SEA-SPAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
