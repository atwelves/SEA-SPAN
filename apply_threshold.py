# script to apply threshold

import numpy as np
import xarray as xr
import matplotlib as mpl
from matplotlib import pyplot as plt

##! This example uses 12 scales, 181 time bins of 24 hours each, and 4 quant(quart)iles
# Could be expanded to 10 quantiles??
# Need to test on real data...

# read in time series
ds_series = xr.open_dataset('BO_TS_TG_Helsinki_proc.nc')
diff      = ds_series.ssh_diff.values
diff2     = np.square(diff)
print(np.shape(diff2))
diff2     = np.reshape(diff2,(181,24))
print(np.shape(diff2))
diff2 = np.nansum(diff2,1)

# read in transform
ds_transform = xr.open_dataset('BO_TS_TG_Helsinki_Ricker.nc')

# find quantiles for each scale
ds_quantiles = ds_transform.quantile(np.linspace(0,1,5),dim="t")

rmse = np.zeros((4,12))

# loop over scales
for scale in range(0,12):
    raw_transform = ds_transform.isel(s=scale).obs_transform.values
    threshold = np.zeros((5))
    for q in range(0,4):
        binary = np.zeros((181))
        threshold[q+1] = ds_quantiles.isel(quantile=q+1,s=scale).obs_transform.values
        scale_hrs = np.power(2,scale)
#        print("{}th quantile at scale {} hours is {}".format(q,scale_hrs,threshold[q+1]))
        binary[raw_transform>threshold[q-1]] = 1
        binary[raw_transform>threshold[q]]   = 0
 #       print(binary)
        rmse[q,scale] = np.sqrt(np.nanmean(np.multiply(binary,diff2)))

plt.pcolormesh(np.transpose(rmse)); plt.colorbar(); plt.show()

# This script will enable questions like:

# Does rmse error increase monotonically with variability?
# Does the relation between error and variability depend on the scale at which variability is assessed?
# Is there a characteristic scale at which the variability most strongly affects model skill
