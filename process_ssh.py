### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SEA-SPAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### --- DESCRIPTION --------------------------------------------------- ###
# Script to process sea level from tide gauges and model outputs.
# 18.09.2024 | 
# A G Twelves
### ------------------------------------------------------------------- ###

### --- GLOBAL VARIABLES ---------------------------------------------- ###
#           |
nemo_name   =       # Name of NEMO model output file
#           |
tide_name   =       # Tide gauge file naming convention
#           |
output_name =       # Output netcdf naming convention
#           |
lat_min     =       # Southern edge of region of interest     
lat_max     =       # Northern edge of region of interest
lon_min     =       # Western edge of region of interest
lon_max     =       # Eastern edge of region of interest
#           |
year_min    =       # First year of interest
year_max    =       # Last year of interest
month_min   =       # First month of interest
month_max   =       # Last month of interest
#           |   
scale_min   =       # Shortest timescale of interest
scale_max   =       # Longest timescale of interest
#           |
time_res    =       # Time resolution of transform
freq_res    =       # Frequency resolution of transform
#           |
qc_pass     =       # Quality control flags deemed good
qc_tol      =       # Number of bad flagged values allowed
#           |
basis_fn    =       # Basis function for transform
#           { 'Cosine'
#             'Ricker'
#             'Gaussian' }
#           |
centre_ssh  =       # Subtract mean from time series (1) or not (0)
#           { 1
#             0 }
#           |
### ------------------------------------------------------------------ ###

### --- Load modules ------------------- ###

import xarray as xr
import numpy as np

### ------------------------------------ ###

### --- Main --------------------------- ###
fmi_dir ="." # directory of tide gauge data
station_list = within_bounds(fmi_dir) # this extracts a full list of stations within bounds
print("Processing {} tide gauge stations".format(np.size(station_list))
for station in station_list:
    # read in tide gauge data
    tg            = read_the_tides(station)
    tg_transform  = do_transform(tg) # transform obs
    # read in model output
    mod           = sea_extract(station)
    mod_transform = do_transform(mod) # transform mod
    # write out data to netcdf
    write_both("original")
    write_both("transformed")
    print("{} processed".format(station))

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
    count = 0 # count number of files
    station_list = np.empty(dtype='string')
    # read all unique stations
    ds                 = xr.open_mfdataset('{}/*_{}{}.nc'.format(source_dir,*,year_min,month_min))    
    tg_lat             = ds.LATITUDE.values
    tg_lon             = ds.LONGITUDE.values
    if(lat_min < tg_lat < lat_max .and. lon_min < tg_lon < lon_max):
        station_name        = ds.STATION.values
        station_list[count] = '{}'.format(station_name)# extract station name
        count = count + 1
        station_list = np.squeeze(station_list[0:count])        
    return (station_list)
### --------------------------------------------------------------------- ###

### --- Read in tide gauge data ----------------------------------------- ###
# GLOBAL VARIABLES:
# tide_name
def read_the_tides(station_id):
    if (quality_pass(station) = 1)
        ds                    = xr.open_mfdataset('{}/*_{}*'.format(source_dir,station_id))
        station_data          = ds.slev.values # account for qc flag here too?
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
def do_transform(time_series)
    if(basis=="cosine"):
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
def write_both(station, obs_transform, mod_transform)
    transform_out = np.concat(obs_transform,mod_transform)
    xr.to_netcdf('{}'.format(station))
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
    ds                 = xr.open_mfdataset('{}/{}'.format(source_dir,tide_name))
    qc_flags           = ds.qc.values
    qc_bad             = np.count_nonzero(qc_flags not in  qc_pass)
    if(qc_bad > qc_tol)
        station_filter = 0
    else
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
def convolute_it(time_series)
    # DECLARE ARRAY:
    make_space(time_series)
    tiled_series = np.tile(time_series,np.size(s_space))
    for tau in range(0,np.size(t_space))
        make_waves(tau)
        convolv_array(tau,:,:) = tiled_series(:,:) * wavelets(:,:)
    return convolv_array
### --------------------------------------------------------------------- ###

### --- Coarsen transform output by binning ----------------------------- ###
# GLOBAL VARIABLES:
# time_res
# freq_res
def sort_bins(unbinned)
    # unbinned has two dimensions, t x f
    # change to four dimensions, tnew x tres x fnew x fres
    t_old    = np.size(unbinned,1)
    f_old    = np.size(unbinned,0)
    t_new    = t_old/t_res
    f_new    = f_old/t_res
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
    ds      = xr.open_dataset('{}/*_{}{}.nc'.format(source_dir,*,year_min,month_min))
    if (stat_axis=="lat"):
        stat_crd = ds.LATITUDE.values
    if (stat_axis=="lon"):
        stat_crd = ds.LONGITUDE.values
    return stat_crd
### --------------------------------------------------------------------- ### 

### --- Make spaces for time and scale ---------------------------------- ###

def make_space(time_series)
    tspace = np.linspace(0,np.size(time_series),np.size(time_series))
    sspace = np.linspace(scale_min,scale_max,scale_max-scale_min)
    return(t_space,s_space)

### --------------------------------------------------------------------- ###

### --- Construct wavelet transform basis functions --------------------- ###
# GLOBAL VARIABLES:
# basis_fn
# scale_min
# scale_max
def make_waves(t_space,s_space)
    # DECLARE ARRAY
    wavelets = np.zeros((np.size(s_space),np.size(t_space))
#    if (basis_fn=="Cosine"):
#        wavelets[:]   = 
#    if (basis_fn=="Ricker"):
#        wavelets[:,:] = 
    if (basis_fn=="Gaussian"):
        for s in range(0,np.size(s_space))
        # use dyadic convention
        scale = np.power(2,s_space)
        wavelets[s,:] = np.exp(-np.square(tspace[:])/np.square(scale))
    return wavelets
### --------------------------------------------------------------------- ###

### --------------------------------------------------------------------- ###

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SEA-SPAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
