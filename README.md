# nimrod-toolbox
Python tools for manipulating Met Office NIMROD data, creating timeseries from it, converting to ascii, creating hydro index for CAESAR-Lisflood and LSDCatchmentModel.

**Author's note**: *I wrote these tools during my PhD to regrid rainfall radardata from the UK and Irish radar network - there are a few rough edges here and there but hopefully they may be of use to other users of Met Office NIMROD data. A few years after writing them, a newer set of tools appeared on the CEDA/BADC website that were not available to me at the time. I would recommend using the newer tools on CEDA if you have access to that service (I.e. you're a UK academic/scientist/researcher/student). (Presumably if you have access to NIMROD data outside the Met Office you have a  CEDA/BADC account)*

## Workflow for creating Rainfall Timeseries and Hydroindex File

### NIMROD_convert.py

This converts the MetOffice binary format radar .dat files into .ascii format.
It also converts the coordinates to UTM from British Grid.

OUTPUT: A bunch of ASCII files. (Make sure you have space if processing a lot of data!)

### NIMROD_timeseries.py

This extracts the rainfall data from the region of interest. You must supply a DEM file, and the script will extract the radar data from the UK image based on the header extents of the DEM. 

The cropped radar data is then extracted and compiled into a timeseries file. The number of columns in the file corespons to the number of rainfall 'zones' or cells in the radar data.

OuTPUT: timeseries.txt, sample_radar_extent.asc

### create_hydroindex.py

This takes the sample_radar_extent.asc file and creates an appropriate hydroindex.asc file for use in CAESAR-Lisflood and LSDCatchmentModel.


