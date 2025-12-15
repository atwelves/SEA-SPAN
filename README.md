**SEA-SPAN** (**SEA** level **SP**ectral **AN**alysis) is a set of Python tools that can be used to extract SSH at tide gauge stations from NEMO model output, and compare against tide gauge data using Fourier and (Ricker) wavelet analysis. 

- 1 script to process sea level data, match with model output locations and transform both model and observations into time-frequency coordinates.
- 5 scripts to plot figures 

**process_ssh.py**:

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
          
Third level functions

      - "get_grid" takes a given directory and model axis (lat/lon), and extracts the model grid from files in that directory, returning latitude or longitude grid as output.

      - "get_point" takes a given directory, station name and axis (lat/lon), extracts the coordinates of the station, and returns either the station latitude or longitude as output.

      - "make_waves" takes a given of values in time and scale space, plus a wavelet centering time (tau), and constructs a series of waveforms centred at tau, with logarithmically varying scale.  It returns this as output in the form of a two dimensional array.

Fourth level functions

      - "gaussian" takes a given time range, centre (tau) and scale, then constructs a one-dimensional gaussian distribution centred at tau, which it returns as output. 
