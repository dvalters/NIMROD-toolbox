from __future__ import division, print_function
from linecache import getline

# -*- coding: utf-8 -*-

'''
This script takes GPM or other precipitation radar data and produces an hourly rainfall timeseries and hydroindex file for use as spatially variable rainfall data in the HAIL-CAESAR hydrodynamic model. The script:
1. Reads in .asc rainfall radar files in their designated timesteps (e.g. UK NIMROD = 5 mins, GPM = 30 mins)
2. Extracts a predefined region from it
3. Converts the predefined row and column region to an hourly rainfall timeseries file

The output from this radar cropping and timeseries tool is to be used to create a hydroindex file (grid fit across a river catchment DEM derived from the grid cells of the rainfall data) for use in spatially variable rainfall modelling (LSDCatchmentModel and HAIL-CAESAR).

The script and its functions are derived from the original script created by D. Valters (April 2015).

F. Muir, May 2018
'''

'''
--- FUNCTION AND PACKAGE IMPORTS ---
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.ndimage
import math
import glob
import os
import pyproj


'''
--- INPUT DATA NAMES ---
'''

# Assigns radar file name and path to its directory to the following variables.
# Any sample radar file (in ascii format) will work.
radar_sample_file = "3B-HHR-GIS.asc.asc.ascRG.asc0913-S000000-E002959.asc.asc.asc"
#radar_sample_file = "3B-HHR-GIS.0913-S160000-E162959_con.asc"
radarpath = "/Users/freyamuir/Documents/Masters/MDiss/primary/max/asc/"
#radarpath = "/Users/freyamuir/Documents/Masters/MDiss/primary/max/asc_scaled/"
radarsource = radarpath + radar_sample_file

# These variables are for accessing the rest of your radar files in the series.
radar_wildcard_file = "*.asc"
radar_mult_source = radarpath + radar_wildcard_file

# Assigns the file name and directory for the DEM that the radar_sample_file is to be clipped to.
basinfile_name = "comonwsh_2.asc" 
basinpath = "/Users/freyamuir/Documents/Masters/MDiss/primary/max/extra/"
basinsource = basinpath + basinfile_name

'''
--- OUTPUT DATA NAMES ---
'''

# Name to write the summed rainfall levels for the entire period to
cumulative_rainfall_raster_name = 'hurricanemax_small_rainfall_totals.asc'

# Name to write rainfall timeseries file names to
original_rainfall_spatial_timeseries_name = 'hurricanemax_small_rainfile_30min.txt'
hourly_spatial_rainfall_timeseries_name = 'hurricanemax_small_rainfile_1hour.txt'

# Name to write the cropped radar file to (for use in the hydroindex.py tool)
cropped_test_radar_name = "small_cropped_radar_test.asc"


'''
--- DATA PREPROCESSING FUNCTIONS ---
'''

# Reads in header information from DEM
def read_ascii_header(ascii_raster_file):
    print("FUNCTION: read_ascii_header")
    with open(ascii_raster_file) as f:
        header_data = [float(f.next().split()[1]) for x in range(6)] #reads the first 6 lines
    return header_data
    
# Parses the header using a loop and the built-in linecache (getline()) module
def parse_header(filesource):
    print("FUNCTION: parse_header")
    hdr = [getline(filesource, i) for i in range(1,7)]
    values = [float(h.split(" ")[-1].strip()) \
     for h in hdr]
    cols,rows,lx,ly,cell,nd = values
    return cols,rows,lx,ly,cell,nd
    
def set_header_values():
    pass

# Converts header array into string (np.savetxt function takes a string header as an argument)
def convert_array_hdr_to_str_hdr(array_header):
    print("FUNCTION: convert_array_hdr_to_str_hdr")
    str_header = 'NCols ' + str(array_header[0]) + '\n' + 'NRows ' + str(array_header[1]) + '\n' + 'xllcorner ' + str(array_header[2]) + '\n' + 'yllcorner ' + str(array_header[3]) + '\n' + 'cellsize ' + str(array_header[4]) + '\n' + 'NODATA_value ' + str(array_header[5]) # long string used for conversion
    return str_header

# Transforms the OSGB coordinates into UTM if required (check DEM and radar sources)
def convert_OSGB36_to_UTM30(xcoord, ycoord):
    print("FUNCTION: convert_OSGB36_to_UTM30")
    # Set up the coord systems
    OSGB36 = pyproj.Proj("+init=EPSG:27700")
    UTM30N = pyproj.Proj("+init=EPSG:32630")
    utm_x, utm_y = pyproj.transform(OSGB36, UTM30N, xcoord, ycoord)
    return utm_x, utm_y

'''
--- CROPPED RADAR FILE CREATION ---
'''

# Calculates what the x and y coordinates should be for the cropped radar data, 
# based on the header data in the radar raster and terrain DEM
def calculate_crop_coords(basin_header, radar_header):
    print("FUNCTION: calculate_crop_coords")
    # Sets values for easier calculations
    y0_radar = radar_header[3]
    x0_radar = radar_header[2]
    print("radar (x,y):", x0_radar, y0_radar)
    
    y0_basin = basin_header[3]
    x0_basin = basin_header[2]
    
    x0_basin_UTM, y0_basin_UTM = convert_OSGB36_to_UTM30(x0_basin, y0_basin)
    print("UTM (x,y):", x0_basin_UTM, y0_basin_UTM)
    
    x0_radar_UTM, y0_radar_UTM = convert_OSGB36_to_UTM30(x0_radar, y0_radar)
    print("radar UTM (x,y):", x0_radar_UTM, y0_radar_UTM)
    
    nrows_radar = radar_header[1]
    ncols_radar = radar_header[0]
    
    nrows_basin = basin_header[1]
    ncols_basin = basin_header[0]

    cellres_radar = radar_header[4]
    cellres_basin = basin_header[4]
    
    xp = x0_basin_UTM - x0_radar_UTM
    yp = y0_basin_UTM - y0_radar_UTM
    print("xp:", xp, yp)

    xpp = ncols_basin * cellres_basin
    ypp = nrows_basin * cellres_basin
    print("ypp:", xpp, ypp) 
    
    # Floor and Ceiling are used to ensure all rainfall area covering basin is captured
    start_col = np.floor( xp / cellres_radar )   # Should be -1 as indexing starts at 0
    end_col = np.ceil( (xpp + xp) / cellres_radar )
    
    start_row = np.floor(nrows_radar - ( (yp + ypp)/cellres_radar ))
    end_row = np.ceil(nrows_radar - (yp/cellres_radar))
    
    print("floor-ceiling:", start_col, start_row, end_col, end_row)
    return int(start_col), int(start_row), int(end_col), int(end_row)


# Do I need this? 
def calculate_crop_coords2(basin_header, radar_header):
    # set values for easier calcs
    print("FUNCTION: calculate_crop_coords2")
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

'''
--- MAIN RADAR CROPPING FUNCTION ---
'''

# Loops over radar rasters and extracts the relevant area of interest.
# Then creates the rainfall timeseries as it goes along, producing:
# 1) Cumulative rainfall raster for the area of extraction. (Total mm/hr rainfall)
# 2) Defined minute interval rainfall timeseries (e.g. 5 mins, 30 mins)
# 3) Hourly rainfall timeseries
def extract_cropped_rain_data():
    print("FUNCTION: extract_cropped_rain_data")
    rainfile = []
    cum_rain_totals = np.zeros((1,1)) # Shape is arbitrary here
    for f in glob.iglob(radar_mult_source):
        print("Filename:", f)
        basin_header = read_ascii_header(basinsource)  # Only one file so this does not change
        radar_header = read_ascii_header(f)   # Checks header for each radar file
        
        # Grid size can change for NIMROD, so these values have to be recalculated for each iteration
        start_col, start_row, end_col, end_row = calculate_crop_coords(basin_header, radar_header)
        
        print(start_col, start_row, end_col, end_row)
        start_col = int(round(start_col))
        start_row = int(round(start_row) )
        end_col = int(round(end_col) )
        end_row = int(round(end_row) )
        print("START COL, START ROW, END COL, END ROW:")
        print(start_col, start_row, end_col, end_row)        
        
        # Loads in the entire rainfall radar grid for this timestep
        # Makes sure that -1 values and NODATA values are masked or made = 0
        cur_rawgrid = np.genfromtxt(f, skip_header=6, filling_values=0.0, loose=True, invalid_raise=False)
        print(cur_rawgrid)
        print(cur_rawgrid.shape)

        # Crops the rainfall radar grid to the subset area
        cur_croppedrain = cur_rawgrid[start_row:end_row, start_col:end_col]#/32 # division by 32 done in NIMROD_convert.py 
        print(cur_croppedrain)
        print(cur_croppedrain.shape)
        
        # Creates the first row of the rainfile by stacking the rainfall rates side by side
        cur_rainrow = np.hstack(cur_croppedrain)
        
        # Adds first row to the rainfile list (i.e. adds the next row for the next timestep)
        rainfile.append(cur_rainrow)
        
        if cum_rain_totals.shape != cur_croppedrain.shape:
            new_shape = cur_croppedrain.shape
            cum_rain_totals = np.zeros(new_shape)
        
        # NOTE: Rain total depends on the temporal resolution of the radar data. 
        # Divide by 12 if 5 min data (12*5 mins = 1 hour)
        # Divide by 2 if 30 min data (2*30 mins = 1 hour)
        cum_rain_totals += cur_croppedrain / 2
        
    # Prints total_rainfall
    np.savetxt(cumulative_rainfall_raster_name, cum_rain_totals, delimiter=' ', fmt='%1.1f')
    
    # Prints rainfile in current format
    rainfile_arr = np.vstack(rainfile)
    np.savetxt(original_rainfall_spatial_timeseries_name, rainfile_arr, delimiter=' ', fmt='%1.1f')
    
    # Optional: if your array(rows) are not divisible by an integer
    # (e.g. if you had an odd number of time steps from the radar)
    rainfile_arr_rows, rainfile_arr_cols = rainfile_arr.shape
    if rainfile_arr_rows % 2 != 0:            # Check if odd number of hydroindex cells
        extra = np.zeros(rainfile_arr_cols)
        rainfile_arr = np.vstack([rainfile_arr, extra])  # Appends extra row to make even
        print(rainfile_arr.shape) # Just to check, should be even number of rows
        rainfile_arr_rows, rainfile_arr_cols = rainfile_arr.shape
    
    # Calculates hourly rainfall from the x min data by taking binned averages (reshapes the array into 3D bins)
    # Shape = number of timesteps in new rainfall file, or radar interval to new timestep interval
    # i.e 5 mins to 1 hour steps would be 60/5 = 12, no of cols
    # When using GPM, radar is in 30 min data. So hourly rainfall from 30 min data is 60/30 = 2
    reshaped_rain = np.reshape(rainfile_arr, ((rainfile_arr_rows//2),2,rainfile_arr_cols)) 
    # [Error of array division not resulting in integer solved with // instead of /]
    
    # Takes the mean down the columns
    hourly_mean = reshaped_rain.mean(axis=1)  
    print(hourly_mean.shape) # Checks by printing number of hours (rows) by grid cells (columns)
    
    rounded_hourly_rain = np.around(hourly_mean, decimals=1)
    
    # Saves to file in the HAIL-CAESAR required format of .asc
    np.savetxt(hourly_spatial_rainfall_timeseries_name, rounded_hourly_rain, delimiter=' ', fmt='%1.1f')  # Gives to 1 decimal place

'''
--- MEAN RAINFALL TIMESERIES CREATION ---
'''

# Creates a mean rainfall timeseries (uniform rainfall), by averaging the spatial file. 
# Does not take into account raincells that only partially cover the catchment
def create_catchment_mean_rainfall(rainfile):
    # Takes the rainfile text file and takes the average of all raincells over the catchment
    # Returns a single hourly timeseries for the whole catchment.
    print("FUNCTION: create_catchment_mean_rainfall")
    rainfile_arr = np.loadtxt(rainfile)
    average_rain_arr = np.mean(rainfile_arr, axis=1)
    np.savetxt(uniform_hourly_rainfall_name, average_rain_arr, delimiter=' ', fmt='%1.1f')

# Examines how much of a radar raincell covers the catchment and calculates a weighting factor to reduce its overall contribution
# If you include radar raincells that don't cover the entire catchment in the average, you are giving them disproportionate weigthing.
def create_catchment_weighted_rainfall(rainfile, terrain_dem):
    print("FUNCTION: create_catchment_weighted_rainfall")
    rainfile_arr = np.loadtxt(rainfile)
    # Mult ratio = ratio of total radar cells to catchment cells to DEM/RADAR raster size
    mult_ratio, weighting_grid = calculate_weighting_array(terrain_dem)
    # A horizontal row will be broadcast to all the columns in the next step, so need to hstack the weighting_grid
    weighting_flattened = np.hstack(weighting_grid)
    # Opens the variable rainfile, multiplies columns by weighting amounts
    weighted_rainfall_arr = rainfile_arr * weighting_flattened
  
    # Prints weighted_rainfall_arr average along axis 1
    print(rainfile_arr)
    average_weighted_rain_array = np.mean(weighted_rainfall_arr, axis=1)*mult_ratio
    np.savetxt(weighted_uniform_hourly_rainfall_name, average_weighted_rain_array, delimiter=' ', fmt='%1.1f')

# Creates an array of weighting factors based on radar rain cells that do not fully cover the catchment domain. 
# This is used in conjuction with create_catchment_weighted_rainfall()
def calculate_weighting_array(terrain_dem):
    print("FUNCTION: calculate_weighting_array")
    terrain_array = np.loadtxt(terrain_dem, skiprows=6)
    terrain_NODATA = read_ascii_header(terrain_dem)[5]
    print(terrain_array)
    
    # No need to reload sample image, just use base indexgrid
    baseindexgrid = create_base_indexgrid(CROPPED_RADAR_DEM)
    weighting_array = np.ones(np.shape(baseindexgrid))
    upscaled_baseindexgrid = create_upscaled_hydroindex(baseindexgrid, CROPPED_RADAR_DEM)
    
    # Number of grid cells in radar image
    no_radar_cells = upscaled_baseindexgrid.size
    # Number of actual catchment cells (total)
    no_catchment_cells = (terrain_array != terrain_NODATA).sum()
    
    # Iterates row-wise over the entire arrat, allow writing values to array
    # Want to iterate in C-style order (row by row)
    for i in np.nditer(baseindexgrid, op_flags=['readwrite'], order='C'):
        print(i)
        
        # Creates a boolean array for where the terrain data is elevations (not NODATA) 
        # and the location in the current radar grid cell. 
        cells_in_catchment_boolean = (terrain_array != terrain_NODATA) & (upscaled_baseindexgrid == i)
        
        # Counts the number of cells in this hydroindex zone
        total_cells_this_hydrozone = np.count_nonzero(upscaled_baseindexgrid == i)
        print("Total number of cells in hydrozone ", i.__int__(), ": ", total_cells_this_hydrozone)

        # Sums the boolean array to get number of hydroindex cells in catchment
        # Gives the number of cells in the current hydroindex zone that are not over NODATA cells
        number_of_catchment_cells = cells_in_catchment_boolean.sum()
        print("No. of catchment cells within hydrozone: ", number_of_catchment_cells)
        # Should be 1 if all the cells are in the catchment
        weighting_ratio = number_of_catchment_cells / total_cells_this_hydrozone
        # Sets weighting ratio in array
        # NOTE: http://docs.scipy.org/doc/numpy/reference/arrays.nditer.html
        print ("Weighting ratio: ", weighting_ratio)
        
        weighting_array.flat[int(i)-1] = weighting_ratio
    
    # Corrects for the size of catchment relative to radar grid rectangular size
    mult_ratio = no_radar_cells / no_catchment_cells
    
    # Prints weighting_array
    return mult_ratio, weighting_array

'''
--- CROPPED RADAR FILE WRITING ---
'''
# Writes a sample cropped radar for use in the hydroindex section
def write_sample_radar_img():
    print("FUNCTION: write_sample_radar_img")
    basin_header = read_ascii_header(basinsource)  # This does not change as there is only one file
    radar_header = read_ascii_header(radarsource)   # Need to check the header each time for radar
    
    # Since grid size can change, these values have to be recalculated for every iteration
    start_col, start_row, end_col, end_row = calculate_crop_coords(basin_header, radar_header)
    print("START COL {0}, START_ROW {1}, END_COL {2}, END_ROW {3}".format(start_col, start_row, end_col, end_row))
    rawgrid = np.loadtxt(radarsource, skiprows=6)
    
    # Make sure that -1 values and NODATA values are masked or made = 0
    croppedrain = rawgrid[start_row:end_row ,start_col:end_col]
    plt.imshow(croppedrain, interpolation='none')
    
    nrows, ncols = croppedrain.shape # Gets rows and cols from the shape of the cropped radar array
    sampleradar_header = basin_header # Brings in the same georeferenced data as the basin DEM
    sampleradar_header[4] = radar_header[4] # Sets the resolution to be the same as the radar
    sampleradar_header[5] = radar_header[5] # Sets NODATA to be the same as well
    sampleradar_header[0] = ncols
    sampleradar_header[1] = nrows
    print("SAMPLE RADAR HEADER", sampleradar_header)
    sampleradarhdr = convert_array_hdr_to_str_hdr(sampleradar_header)    
    
    np.savetxt(cropped_test_radar_name, croppedrain, header=sampleradarhdr, comments='', delimiter=' ', fmt='%1.1f')

'''
--- EXECUTE FUNCTIONS ---
'''
# EXTRACTS THE TIMESERIES FROM THE ASCII RADAR FILES
extract_cropped_rain_data()   

# If just a sample is required, the hydroindex is not affected
# WRITES A SAMPLE CROPPED RADAR FILE
write_sample_radar_img()

# CREATES AN AVERAGE RAINFALL FILE
#create_catchment_mean_rainfall("C:/INSERT FILE PATH/and name.txt")