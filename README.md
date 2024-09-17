**SEA-SPAN** (**SEA** level **SP**ectral **A**nalysis for **N**emo) is a set of Python tools that can be used to extract SSH at tide gauge stations from NEMO model output, and compare against tide gauge data using Fourier and (Ricker) wavelet analysis. 

- 1 script to download (?)
- 1 script to process
    - 1 function to select bounds in lat, lon, time, and frequency domain
    - 1 function to read tide gauge qc flags
    - 1 function to read NEMO grid
    - 1 function to read 
- 1 script to plot outputs

Processing script produces an array with dimensions **2** x **N_{s}** x (**L_{t}**/**R_{t}**) x (**L_{f}**/**R_{f}**)
