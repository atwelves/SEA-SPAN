### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SEA-SPAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

### --- DESCRIPTION ---------------------------------------------------------------------------- ###
# Script to calculate model-obs RMSE, as a function of scale and magnitude of sea level variability
# 04.11.2024 |
# A G Twelves
### -------------------------------------------------------------------------------------------- ###

### --- Load modules ----------------------- ###

import xarray as xr
import numpy as np
from pathlib import Path
import hashlib
import pandas as pd
import netCDF4 as nc
from netCDF4 import Dataset

### ---------------------------------------- ###

#
station = "BO_TS_TG_Helsinki"
Nq      = 4

##! This example uses 12 scales, 181 time bins of 24 hours each, and 4 quant(quart)iles
# Could be expanded to 10 quantiles??
# Need to test on real data...

### --- Functions -------------------------------------------------------------------- ###

def binary_choice(quant,scale_exp,raw_transform):
    binary       = np.zeros((np.size(diff2)))
    threshold_lo = ds_quantiles.isel(quantile=quant,s=scale_exp).obs_transform.values
    threshold_hi = ds_quantiles.isel(quantile=quant+1,s=scale_exp).obs_transform.values
    print(quant,scale_exp,threshold_hi)
    binary[raw_transform>threshold_lo] = 1
    binary[raw_transform>threshold_hi] = 0
    return binary

def rmse(difference):
    output = np.sqrt(np.nanmean(difference))
    return output 

def write_out(error_matrix):
    print(np.shape(error_matrix))
    outfile             = '{}_rmse_.nc'.format(station)
    ncfile              = Dataset(outfile,mode='w')
    q                   = ncfile.createDimension('q',np.size(error_matrix,1))    # time dimension
    s                   = ncfile.createDimension('s',np.size(error_matrix,0))    # scale dimension
    ssh_rmse            = ncfile.createVariable('ssh_rmse',np.float32,('s','q'))
    ssh_rmse[:,:]       = error_matrix[:,:]
    ncfile.close(); print('Dataset is closed')
    return 

### ---------------------------------------------------------------------------------- ###

### --- Main -------------------------------------------------------- ###

# Time series
ds_series = xr.open_dataset('{}_proc.nc'.format(station))
diff      = ds_series.ssh_diff.values  # Read in model-obs error
diff2     = np.square(diff)            # Square it...
diff2     = np.reshape(diff2,(181,24)) # reshape...
diff2     = np.nansum(diff2,1)         # and do daily sum

# Transform
ds_transform = xr.open_dataset('{}_Ricker.nc'.format(station))

# Find quantiles for each scale
ds_quantiles = ds_transform.quantile(np.linspace(0,1,Nq+1),dim="t")

error_matrix = np.zeros((12,Nq)) # this will be the output

# loop over scales
for scale_exp in range(0,12):
    raw_transform = ds_transform.isel(s=scale_exp).obs_transform.values # Read in transform for given scale
    for quant in range(0,Nq): 
        filtr = binary_choice(quant,scale_exp,raw_transform)            # Select only times in this quantile
        error_matrix[scale_exp,quant] = rmse(np.multiply(filtr,diff2))  # Calculate array of rmse values

write_out(error_matrix)

### ----------------------------------------------------------------- ###

# This script will enable questions like:

# Does rmse error increase monotonically with variability?
# Does the relation between error and variability depend on the scale at which variability is assessed?
# Is there a characteristic scale at which the variability most strongly affects model skill?
