# -*- coding: utf-8 -*-
"""
Created on Sat Apr 3 2015

This tool takes the nimrod radar data (already converted from native format
to ascii, using the NIMpy tool), extracts a predefined region from it (note the column and row
indexing convetions below), and converts it into an hourly radar rainfall
timeseries in the format for the CAESAR-Lisflood landscape evoltuion model, 

This script is also designed to create a hydroindex file from a clipped radar
file in ascii format. Hydroindex files are require for use in the LSDCatchmentModel
and CAESAR-Lisflood models when investigating spatially variable rainfall.

Run this script in the same directory as your clipped (and correctly georeffed)
radar ascii data file. You will need to set the upsampling resolution. i.e. the 
grid cell resolution of your corresponding DEMs that will be used in the model
simulation.

So if your DEMs are at 10m resolution, the upscaling resultion needs to be set at 10.

This script does not yet clip the resulting hydroindex to the catchment extent.

DAV, March 2015

To Do:
1) Error handling if Cellsize is not an integer
2) Writing to different formats (default is ascii grid .asc)

@author: Declan Valters
"""

## File for Converting the 5 minute NIMROD data (already converted to ascii format)
## into hourly rainfall rates

from __future__ import division
# because we want to get the answer 5/2 = 2.5 (Python 3), not 5/2 = 2 (Python 2.x)
from linecache import getline

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.ndimage
import math
import glob
import os
import pyproj

#import create_hydroindex.py

# File name of ASCII digital elevation model
#path = "/home/dvalters/asciifiles/"
# As long as you are sure all your radar files have been converted similarly, any sample radar (in ascii format) will work.
#radarpath = "D:\\DATASETS\\NIMROD\\NIMROD\\1km-Composite\\2005\\19June\\datfiles\\"
radarpath = "/home/dav/DATADRIVE/DATASETS/NIMROD/NIMROD/BOSCASTLE/Boscastleflood/asciifiles/"
radarsource = radarpath + "metoffice-c-band-rain-radar_uk_200408160000_1km-composite.txt"

#basinpath = 'D:\\CODE_DEV\\PyToolsPhD\\Radardata_tools\\multiple_radar_test\\'
basinpath = "/home/dav/DATADRIVE/CODE_DEV/PyToolsPhD/Radardata_tools/"
basinsource = basinpath + 'boscastle5m_bedrock_fill_outlet.asc'

CROPPED_RADAR_DEM = "cropped_test_radar.asc"
TERRAIN_DEM = "ryedale20m_fillcrop.asc"

"""
Reads in header information from an ASCII DEM
"""
def read_ascii_header(ascii_raster_file):
    
    with open(ascii_raster_file) as f:
        header_data = [float(f.next().split()[1]) for x in xrange(6)] #read the first 6 lines
    return header_data
    
"""
Parse the header using a loop and
the built-in linecache module
"""
def parse_header(filesource):
    hdr = [getline(filesource, i) for i in range(1,7)]
    values = [float(h.split(" ")[-1].strip()) \
     for h in hdr]
    cols,rows,lx,ly,cell,nd = values
    return cols,rows,lx,ly,cell,nd
    

def set_header_values():
    pass

"""
The np.savetxt function takes a string header as an argumnet, use this to
convert the array into a long string
"""
def convert_array_hdr_to_str_hdr(array_header):
    str_header = 'NCols ' + str(array_header[0]) + '\n' + 'NRows ' + str(array_header[1]) + '\n' + 'xllcorner ' + str(array_header[2]) + '\n' + 'yllcorner ' + str(array_header[3]) + '\n' + 'cellsize ' + str(array_header[4]) + '\n' + 'NODATA_value ' + str(array_header[5])
    return str_header

"""
This is to transform the osgb (Great Britain Ordnance Survey) coords into UTM ones.
Often, DEM data and RADAR data is supplied in OSGB if you get it from a UK source.
Check both sources though!
"""
def convert_OSGB36_to_UTM30(xcoord, ycoord):
    # Set up the coord systems
    OSGB36 = pyproj.Proj("+init=EPSG:27700")
    UTM30N = pyproj.Proj("+init=EPSG:32630")
    utm_x, utm_y = pyproj.transform(OSGB36, UTM30N, xcoord, ycoord)
    return utm_x, utm_y

"""
Creates the base hydroindex grid from a cropped radar image tile
"""
def create_base_indexgrid(CROPPED_RADAR_DEM):
    ncols = read_ascii_header(CROPPED_RADAR_DEM)[0]
    nrows = read_ascii_header(CROPPED_RADAR_DEM)[1]
    no_of_cells = ncols*nrows
    # The base grid is based on the coarse radar data clipped to the right area.
    # So it should only have a low number of cells really
    # +1 at the end because we want the numbering to start at 1 instead of zero 
    # for later purposes.
    base_indexgrid = np.arange(no_of_cells).reshape(nrows,ncols) + 1   
    
    return base_indexgrid
  
"""
Resamples the base indexgrid to the resolution of you terrain DEM.
Returns the upscaled hydroindex.
"""  
def create_upscaled_hydroindex(base_hydroindex, CROPPED_RADAR_DEM):
    radar_cellsize = read_ascii_header(CROPPED_RADAR_DEM)[4]
    terrain_cellsize = read_ascii_header(TERRAIN_DEM)[4]
    # We want to scale up the base grid to the resolution of the terrain DEM
    # that will be used in the model run by a mult_factor
    #mult_factor = radar_cellsize/terrain_cellsize
    #print mult_factor
    radarNcols = read_ascii_header(CROPPED_RADAR_DEM)[0]
    radarNrows = read_ascii_header(CROPPED_RADAR_DEM)[1]
    terrainNcols = read_ascii_header(TERRAIN_DEM)[0]
    terrainNrows = read_ascii_header(TERRAIN_DEM)[1]
    
    #mult_factor =  math.sqrt( (terrainNcols*terrainNrows) )/ math.sqrt((radarNcols*radarNrows) )
    #print mult_factor
    
    # Yes! this us the way to do it!
    mult_factor_cols = terrainNcols / radarNcols
    mult_factor_rows = terrainNrows / radarNrows
    
    print mult_factor_cols, mult_factor_rows
    
    # Now resample. Using nearest neighbour (order=0) as we don't want any fancy interpolation going on
    # (Although this could be a later feature) (order!=0)
    # The mult_factor can have different values for each axis, and they can be float values as I understand
    upscaled_hydroindex = scipy.ndimage.zoom(base_hydroindex, (mult_factor_rows, mult_factor_cols), order=0)
    print upscaled_hydroindex.shape
    
    return upscaled_hydroindex

"""
Crops the upscaled hydroindex grid to the outline of your terrain basin
Only needed if you have a clipped catchment DEM
"""    
def crop_upscaled_hydroindex_to_basin(upscaled_hydroindex, TERRAIN_DEM):
    #load the terrain DEM
    terrain_array = np.loadtxt(TERRAIN_DEM, skiprows=6)
    terrain_NODATA = read_ascii_header(TERRAIN_DEM)[5]
    # Create a mask from the terrain DEM
    terrain_nodata_mask = np.ma.masked_where(terrain_array<-1, terrain_array)
    #plt.imshow(terrain_nodata_mask)
    
    # Apply the mask to the hydroindex grid
    cropped_hydroindex = np.ma.masked_array(upscaled_hydroindex, terrain_nodata_mask.mask)
    
    
    crop_with_nodataval_hydroindex = np.ma.filled(cropped_hydroindex, terrain_NODATA)
    plt.imshow(crop_with_nodataval_hydroindex)
    
    return crop_with_nodataval_hydroindex

"""
Writes the hydroindex grid to an ascii DEM with the same header
information as the terrain DEM
"""
def write_hydroindex_to_file(hydroindex_filename, hydroindex_masked):
    #hydroindex_header = np.zeros(6)
    header = "ncols " + str(hydroindex_masked.shape[1]) + "\n"
    header += "nrows " + str(hydroindex_masked.shape[0]) + "\n"
    header += "xllcorner " + str(read_ascii_header(TERRAIN_DEM)[2]) + "\n"
    header += "yllcorner " + str(read_ascii_header(TERRAIN_DEM)[3]) + "\n"
    header += "cellsize " + str(read_ascii_header(TERRAIN_DEM)[4]) + "\n"
    header += "NODATA_value " + str(read_ascii_header(TERRAIN_DEM)[5])
    # Remember Hydroindex must be integer
    np.savetxt(hydroindex_filename, hydroindex_masked, fmt="%1i", header=header, comments='')

"""
This function calculates what the x and y coordinates
should be for the cropped radar data, based on the header
data in the radar raster and terrain DEMs
"""
def calculate_crop_coords(basin_header, radar_header):
    # set values for easier calcs
    y0_rad = round(radar_header[3])
    x0_rad = round(radar_header[2])
    
    y0_bas = round(basin_header[3])
    x0_bas = round(basin_header[2])
    
    x0_bas, y0_bas = convert_OSGB36_to_UTM30(x0_bas, y0_bas)
    
    nrows_rad = radar_header[1]
    ncols_rad = radar_header[0]
    
    nrows_bas = basin_header[1]
    ncols_bas = basin_header[0]

    cellres_rad = radar_header[4]
    cellres_bas = basin_header[4]
    
    xp = x0_bas - x0_rad
    yp = y0_bas - y0_rad
    
    xpp = ncols_bas * cellres_bas
    ypp = nrows_bas * cellres_bas
    
    start_col = xp / cellres_rad
    end_col = (xpp + xp) / cellres_rad
    
    start_row = nrows_rad - ( (yp + ypp)/cellres_rad )
    end_row = nrows_rad - (yp/cellres_rad)
    
    return start_col, start_row, end_col, end_row

"""
This loops over the radar rasters and extracts the relevant sub-section
of the radar data area of interest.
It then creates the rainfall timeseries as it goes along (maybe split this into separate function??)

It outputs:
1) A cumulative rainfall raster for the area of extraction. (Total mm/hr rainfall)
2) A 5min rainfall timeseries
3) An hourly rainfall timeseries
"""
def extract_cropped_rain_data():
    rainfile = []
    cum_rain_totals = np.zeros((1,1)) # shape is arbitrary here, it gets changed later
    for f in glob.iglob('/home/dav/DATADRIVE/DATASETS/NIMROD/NIMROD/BOSCASTLE/Boscastleflood/asciifiles/*composite.txt'):
        #print f
        basin_header = read_ascii_header(basinsource)  # this does not change as there is only one file
        radar_header = read_ascii_header(f)   # need to check the header each time for radar
        
        # since grid size can change (Thanks, Met Office!), these values have to be recaclulated for every iteration...ugh
        start_col, start_row, end_col, end_row = calculate_crop_coords(basin_header, radar_header)
        
        print start_col, start_row, end_col, end_row
        start_col = int(round(start_col))
        start_row = int(round(start_row)) 
        end_col = int(round(end_col))
        end_row = int(round(end_row))
        
        ##### WARNING MAsSIVE FUDGE CODE WITH MAGIC NUMBERS
        if radar_header[0] < 1000:
            end_col = end_col + 1
        ##### WARNING THIS IS A FUDGE!

        #print start_col, start_row, end_col, end_row        
        
        # Load in the entire rainfall radar grid for this timestep
        cur_rawgrid = np.loadtxt(f, skiprows=6)
        # perhaps make sure that -1 values and NODATA values are masked or made = 0 
        # (should not be issue over land but better safe than sorry)
        
        # Crop the rainfall radar grid of the whole country to the subset area
        cur_croppedrain = cur_rawgrid[start_row:end_row, start_col:end_col]###/32   division done in NimPy now 
        
        
        # Create the first row of the rainfile by stacking the rainfall rates side by side with hstack
        cur_rainrow = np.hstack(cur_croppedrain)
        
        # Add this to the rainfile list (i.e. add the next row for the next timestep)
        rainfile.append(cur_rainrow)
        
        if cum_rain_totals.shape != cur_croppedrain.shape:
            new_shape = cur_croppedrain.shape
            cum_rain_totals = np.zeros(new_shape)
        
        cum_rain_totals += cur_croppedrain
        
    #print total_rainfall
    np.savetxt('rainfall_totals_boscastle.asc', cum_rain_totals, delimiter=' ', fmt='%1.1f')
    
    #print rainfile in current format
    rainfile_arr = np.vstack(rainfile)
    np.savetxt('test_rainfile_5min.txt', rainfile_arr, delimiter=' ', fmt='%1.1f')
    
    ## Optional, if your array(rows) are not divisible by an integer
    ## (This would happen if you had an odd number of time steps from the radar)
    rainfile_arr_rows, rainfile_arr_cols = rainfile_arr.shape
    if rainfile_arr_rows % 2 != 0:            # check if odd
        extra = np.zeros(rainfile_arr_cols)   # number of hydroindex cells
        rainfile_arr = np.vstack([rainfile_arr, extra])  # append extra row to make even
        print rainfile_arr.shape #just to check, should be even now
        rainfile_arr_rows, rainfile_arr_cols = rainfile_arr.shape # double check we have even no. or rows
    
    ## Now we are going to calculate hourly rainfall from the 5 min data by taking binned averages
    # First, reshape the array into 3d bins:
    reshaped_rain = np.reshape(rainfile_arr, ((rainfile_arr_rows/12),12,rainfile_arr_cols)) # assuming 
    # shape = (number of timesteps in new rainfall file, radar interval to new timestep interval [i.e 5mins to 1hour steps would be 60/5 =12], no of cols
    
    hourly_mean = reshaped_rain.mean(axis=1)  # we are taking the mean down the columns. Remember that we have already binned into hourly bins previously.
    print hourly_mean.shape # just to check (it should print the number of hours (rows) by grid cells (columns))
    
    rounded_hourly_rain = np.around(hourly_mean, decimals=1)
    
    ## Save to file in the CAESAR required format
    np.savetxt('test_rainfile_hourly.txt', rounded_hourly_rain, delimiter=' ', fmt='%1.1f')  # the fmt bit gives it to 1 decimal place. Oddly, the rounding step above is non-permanent??

"""
This creates a mean rainfall timeseries (uniform rainfall), by simply
averaging the spatial file. It does NOT take into account raincells that only
partially cover the catchment
"""    
def create_catchment_mean_rainfall(rainfile):
    # For the Uniform rainfall test cases in CAESAR
    # Takes the rainfile txt file and takes the average of all raincells over the catchment, returning a single hourly timeseries for
    # the whole catchment.
    rainfile_arr = np.loadtxt(rainfile)
    average_rain_arr = np.mean(rainfile_arr, axis=1)
    np.savetxt("RYEDALE_rainfile_uniform72hr_5min.txt", average_rain_arr, delimiter=' ', fmt='%1.1f')

"""
For mean rainfall, cells not entirely within the catchmnet boundaries should be weighted appropriately.  
I.e. If you include radar raincells that don't cover the entire catchment in the average, you are
giving them disproportionate weigthing.

This function examines how much of a radar raincell covers the catchment, and calculates a
weghting factor to reduce its contribution to the catchment-wide mean.

(Thanks to Chris Skinner for pointing this out!)
"""  
def create_catchment_weighted_rainfall(rainfile, terrain_dem):
    rainfile_arr = np.loadtxt(rainfile)
    weighting_grid = calculate_weighting_array(terrain_dem)
    # stack the weighting grid
    # we are going to broadcast a horizontal row to all the columns in the next step,
    # so we need to hstack our weighting_grid.
    weighting_flattened = np.hstack(weighting_grid)
    # open the variable rainfile, multiply columns by weighting amounts
    weighted_rainfall_arr = rainfile_arr * weighting_flattened
    # average along axis 1.
    average_weighted_rain_array = np.mean(weighted_rainfall_arr, axis=1)
    np.savetxt("WEIGHTED_UNIFORM_RAINFALL.txt", average_weighted_rain_array, delimiter=' ', fmt='%1.1f')

"""
Creates an array of weighting factors based on radar rain cells that 
do not fully cover the catchment domain. This is used in conjuction with create_catchment_weighted_rainfall()
"""
def calculate_weighting_array(terrain_dem):
    terrain_array = np.loadtxt(terrain_dem, skiprows=6)
    terrain_NODATA = read_ascii_header(terrain_dem)[5]
    
    # no need to reload sample image, use base index grid
    baseindexgrid = create_hydroindex.create_base_indexgrid(CROPPED_RADAR_DEM)
    weighting_array = np.ones(np.shape(baseindexgrid))
    
    upscaled_baseindexgrid = create_hydroindex.create_upscaled_hydroindex(baseindexgrid, CROPPED_RADAR_DEM)
    
    # Iterate row-wise over the entire arrat, allow writing values to array
    # Specify that we want to iterate in C-style order, i.e. row by row.
    for i in np.nditer(weighting_array, op_flags=['readwrite'], order=['C']):
        # Check that the arrays are the same shape
        print terrain_array.shape == upscaled_baseindexgrid.shape
        # Create a boolean array where the terrain data is elevations (not nodata) 
        ## and where we are in the current radar grid cell. 
        cells_in_catchment_boolean = (terrain_array != terrain_NODATA) & (upscaled_baseindexgrid == i)
        # Sum the boolean array to get number of hydroindex cells in catchmnet
        number_of_catchment_cells = cells_in_catchment_boolean.sum()
        # Will be 1 if all cells are in the catchment
        weighting_ratio = number_of_catchment_cells / cells_in_catchment_boolean.size
        # set weighting ratio in array
        # Note: http://docs.scipy.org/doc/numpy/reference/arrays.nditer.html
        i[...] = weighting_ratio
    
    print weighting_array
    return weighting_array

"""
This writes a sample cropped radar raster for checking purposes.
"""
def write_sample_radar_img():
    basin_header = read_ascii_header(basinsource)  # this does not change as there is only one file
    radar_header = read_ascii_header(radarsource)   # need to check the header each time for radar
    
    # since grid size can change, these values have to be recaclulated for every iter.
    start_col, start_row, end_col, end_row = calculate_crop_coords(basin_header, radar_header)
    rawgrid = np.loadtxt(radarsource, skiprows=6)
    
    # perhaps make sure that -1 values and NODATA values are masked or made = 0 
    #(should not be issue over land but better safe than sorry)
    croppedrain = rawgrid[start_row:end_row ,start_col:end_col]
    mpl.pylab.imshow(croppedrain, interpolation='none')
    
    nrows, ncols = croppedrain.shape # get rows and cols from the shape of the cropped radar array
    sampleradar_header = basin_header # bring in the same georef data as the basin DEM
    sampleradar_header[4] = radar_header[4] # set the resolution to be the same as the radar
    sampleradar_header[5] = radar_header[5] # set no data to be the same as well
    sampleradar_header[0] = ncols
    sampleradar_header[1] = nrows
    print "SAMPLE RADAR HEADER", sampleradar_header
    sampleradarhdr = convert_array_hdr_to_str_hdr(sampleradar_header)    
    
    np.savetxt("cropped_test_radar.asc", croppedrain, header=sampleradarhdr, comments='', delimiter=' ', fmt='%1.1f')

#-=-=-=-=-=-=#
# MAIN -=-=-=#
#-=-=-=-=-=-=#

# no. Met office data changes format throughout the day for some datasets, you need a more
# robust solution that calculates header and start cols each time!!!
# EXTRACT THE TIME SERIES FROM THE ASCII RADAR FILES
extract_cropped_rain_data()   

# This is ok if you just want a sample, the hydroindex is not affected
# As long as resolution is the same in each file! (this should be more reliable)
# WRITE A SAMPLE CROPPED RADAR FILE
##write_sample_radar_img()

# Create an average rainfall file
#create_catchment_mean_rainfall("C:\\DATA\\PROJECTS\\HYDRO-GEOMORPHIC_STORM_RESPONSE\\CASE_STUDIES\\RYEDALE\\SPATIAL_72hr\\RYEDALE_spatial_72hr_5min.txt")

##### TO DO #####
# Write sample raster cropped ascii with correct header info
# georef stuff can come from the basin_header - is the same spatial extent by definition
# plot overlay of radar onto terrain dem
# create total radar accumultion plots
    
#=-=-=-=-=-=-=-=-=-=
# CREATE HYDRPINDEX
#=-=-=-=-=-=-=-=-=-=

base_hydroindex = create_base_indexgrid(CROPPED_RADAR_DEM)
upscaled_hydroindex = create_upscaled_hydroindex(base_hydroindex, CROPPED_RADAR_DEM)

hydroindex_masked = crop_upscaled_hydroindex_to_basin(upscaled_hydroindex, TERRAIN_DEM)


#plt.imshow(terrain_mask_test, interpolation='none')
#plt.imshow(upscaled_hydroindex, interpolation='none')

#print base_hydroindex
#print upscaled_hydroindex

write_hydroindex_to_file("ryedale_hydroindex_test2.asc", hydroindex_masked)