**SEA-SPAN** (**SEA** level **SP**ectral **A**nalysis for **N**emo) is a set of Python tools that can be used to extract SSH at tide gauge stations from NEMO model output, and compare against tide gauge data using Fourier and (Ricker) wavelet analysis. 

- 1 script to download (?)
- 1 script to process
    - 1 function to select bounds in lat, lon, time, and frequency domain
    - 1 function to read tide gauge qc flags
    - 1 function to read NEMO grid
    - 1 function to read 
- 1 script to plot outputs

Processing script produces a NetCDF file with dimensions 

**2** x **N_{s}** x (**L_{t}**/**N_{t}**) x (**L_{f}**/**N_{f}**), 

where **N_{s}** is the number of stations selected, **L_{t}** is the length of the time series, **N_{t}** is the number of time bins within the transform (= 1 for Fourier), **L_{f}** is the range of frequencies that the transform is conducted over, and **N_{f}** is the number of frequency bins within the transform, and the factor of **2** accounts for the concurrence of model outputs and observational data.

*Plotting script does several operations, including differencing over the first dimension (model-obs error), averaging over the third dimension (to generate a spectrum), summing over the fourth dimension (to produce a time series), thresholding (to produce a binary that can be convoluted with the original time series to seperate out stable and unstable conditions.*  
