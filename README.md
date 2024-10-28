**SEA-SPAN** (**SEA** level **SP**ectral **AN**alysis) is a set of Python tools that can be used to extract SSH at tide gauge stations from NEMO model output, and compare against tide gauge data using Fourier and (Ricker) wavelet analysis. 

- 1 script to download (?)
- 1 script to process
  
Top level functions
  
      - "within_bounds" takes a given directory and searches for data belonging to tide gauges within specified latitude and longitude limits, then checks if these stations have a satisfactory number of data values that pass quality control. It returns the matching stations as a list.
           
      - "read_the_tides" takes a given directory and station name, reads in sea level data, and converts to a mean-centred time series without gaps.  It returns this series as output.
           
      - "sea_extract" takes a given directory and station name and extracts, from the gridded model output, a mean-centred time series without gaps.  It returns this series as output.
          
      - "do_transform" takes a given time series, performs a Fourier or wavelet transform, and returns the normalised transform as an output.
           
      - "write_all" takes a given station name, plus original and transformed series from both tide gauge data and model outputs, and writes out netCDF files.

Second level functions

      - "quality_pass" takes a given directory and station name, loops over all the files belonging to that station in the directory, and counts the number of data points with bad quality control values.  Depending on whether the number of these points exceeds a user defined threshold, a logical one/zero is returned.
      
      - "match_point" takes a given directory and station name, and finds the model indices describing the point closest to the location of the tide gauge station, returning these indices as output.

      - "sub_ marine" takes a given time series, subtracts the mean, then sets any NaN values to zero and returns the corrected series as output.

      - "convolute_it" takes a given time series, convolutes it with a wavelet basis, and returns a 3 dimensional array. 

      - "sort_bins" takes a given two dimensional array, and sorts it into bins along each axis, effectively coarsening the resolution, and returns this coarsened array as output.
          
- 1 script to plot outputs

Processing script produces a NetCDF file with dimensions 

**2** x **N_{s}** x (**L_{t}**/**N_{t}**) x (**L_{f}**/**N_{f}**), 

where **N_{s}** is the number of stations selected, **L_{t}** is the length of the time series, **N_{t}** is the number of time bins within the transform (= 1 for Fourier), **L_{f}** is the range of frequencies that the transform is conducted over, and **N_{f}** is the number of frequency bins within the transform, and the factor of **2** accounts for the concurrence of model outputs and observational data.

*Plotting script does several operations, including differencing over the first dimension (model-obs error), averaging over the third dimension (to generate a spectrum), summing over the fourth dimension (to produce a time series), thresholding (to produce a binary that can be convoluted with the original time series to seperate out stable and unstable conditions.*  
