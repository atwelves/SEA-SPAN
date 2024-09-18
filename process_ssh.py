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
freq_bins   =       # Number of frequency bins in transform
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

### --- Set bounds in latitude, longitude, time, and frequency domain --- ###
# GLOBAL VARIABLES:
# year_min
# year_max
# month_min
# month_max
def set_bounds(t_init, t_final):
    for (x in nx):
    return (t_init, t_final)
### --------------------------------------------------------------------- ###

### --- Read tide gauge quality control flags --------------------------- ###
# GLOBAL VARIABLES:
# tide_name
# qc_pass
# qc_tol
def quality_pass(station):
    ds                 = xr.open_dataset('{}'.format(tide_name))
    qc_flags           = ds.qc.values
    qc_bad             = np.count_nonzero(qc_flags not in  qc_pass)
    if(qc_bad > qc_tol)
        station_filter = 0
    else
        station_filter = 1
    xr.Dataset.close(ds)
    return station_filter
### --------------------------------------------------------------------- ###

### --- Read in tide gauge data ----------------------------------------- ###
# GLOBAL VARIABLES:
# tide_name
def read_the_tides(station):
    if (quality_pass(station) = 1)
        ds                    = xr.open_dataset('{}'.format(tide_name))
        station_data          = ds.slev.values
        station_data          = sub_marine(station_data)
    return station_data
### --------------------------------------------------------------------- ###

### --- Get NEMO grid coordinates --------------------------------------- ###
# GLOBAL VARIABLES:
# nemo_name
def get_grid():
    ds      = xr.open_dataset(nemo_name)
    mod_lat = ds.nav_lat.values
    mod_lon = da.nav_lon.values
    return (mod_lat, mod_lon)
### --------------------------------------------------------------------- ### 

### --- Match NEMO grid points to tide gauge locations ------------------ ###
def match_points(mod_lat, mod_lon, stat_lat, stat_lon):
    mod_y = np.argmin(mod_lat - stat_lat)
    mod_x = np.argmin(mod_lon - stat_lon)
    return (mod_y, mod_x)
### --------------------------------------------------------------------- ###

### --- Extract time series from point in NEMO output closest to station  ###
# GLOBAL VARIABLES:
# nemo_name
def sea_extract(mod_y, mod_x):
    ds        = xr.open_dataset('{}'.format(nemo_name))
    sl        = ds.isel(j=mod_lat, i=mod_lon)
    ssh_mod   = sl.ssh.values
    ssh_mod   = sub_marine(ssh_mod)
    return ssh_mod
### --------------------------------------------------------------------- ###

### --- Subtract mean sea level from time series, removes NaNs ---------- ###
# GLOBAL VARIABLES:
# centre_ssh
def sub_marine(time_series):
    if(centre_ssh==1):
        time_series = time_series - np.nanmean(time_series)
    time_series[np.isnan(time_series)] = 0
    return time_series

### --- Construct wavelet transform basis functions --------------------- ###
# GLOBAL VARIABLES:
# basis_fn
# scale_min
# scale_max
def make_waves(t_init, t_final)
    for t in range(t_init, t_final)
    return
### --------------------------------------------------------------------- ###

### --- Convolute basis functions with time series ---------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def convolute_it(time_series, wavelets)
    return convolv_array
### --------------------------------------------------------------------- ###

### --- Sum convolution to complete transform --------------------------- ###
# GLOBAL VARIABLES:
# basis_fn
def do_transform(convolv_array)
    return transform
### --------------------------------------------------------------------- ###

### --- Do binning?
def collect_bins(transform)
    return binned_transform
### --------------------------------------------------------------------- ###

### --- Write arrays to netcdf ------------------------------------------ ###
# GLOBAL VARIABLES:
# output_name
def write_both(station, )
    xr.to_netcdf('{}'.format(station))
    return 
### --------------------------------------------------------------------- ###

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SEA-SPAN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
