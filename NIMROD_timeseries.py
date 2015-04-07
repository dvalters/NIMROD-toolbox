# -*- coding: utf-8 -*-
"""
Created on Sat Apr 3 2015

This tool takes the nimrod radar data (already converted from native format
to ascii, using the NIMpy tool), extracts a predefined region from it (note the column and row
indexing convetions below), and converts it into an hourly radar rainfall
timeseries in the format for the CAESAR-Lisflood landscape evoltuion model, 


@author: Declan Valters
"""

## File for Converting the 5 minute NIMROD data (already converted to ascii format)
## into hourly rainfall rates

import glob
import os
import matplotlib as mpl

from linecache import getline
import numpy as np

import pyproj


# File name of ASCII digital elevation model
#path = "/home/dvalters/asciifiles/"
# As long as you are sure all your radar files have been converted similarly, any sample radar (in ascii format) will work.
radarpath = "D:\\DATASETS\\NIMROD\\NIMROD\\1km-Composite\\2005\\19June\\datfiles\\"
radarsource = radarpath + "metoffice-c-band-rain-radar_uk_200506191500_1km-composite.dat.asc"

basinpath = 'D:\\CODE_DEV\\PyToolsPhD\\Radardata_tools\\multiple_radar_test\\'
basinsource = basinpath + 'ryedale_20m_basin.asc'

# This is to transform the osgb coords into UTM ones.
def convert_OSGB36_to_UTM30(xcoord, ycoord):
    # Set up the coord systems
    OSGB36 = pyproj.Proj("+init=EPSG:27700")
    UTM30N = pyproj.Proj("+init=EPSG:32630")
    utm_x, utm_y = pyproj.transform(OSGB36, UTM30N, xcoord, ycoord)
    return utm_x, utm_y

# Parse the header using a loop and
# the built-in linecache module
def parse_header(filesource):
    hdr = [getline(filesource, i) for i in range(1,7)]
    values = [float(h.split(" ")[-1].strip()) \
     for h in hdr]
    cols,rows,lx,ly,cell,nd = values
    return cols,rows,lx,ly,cell,nd
    
# Another header reading method
def read_ascii_header(filesource):
    with open(filesource) as f:
        header_data = [float(f.next().split()[1]) for x in xrange(6)] 
        return header_data

def set_header_values():
    pass

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
    

rainfile = []

# Now do this in a file-by-file loop
def extract_cropped_rain_data():
    for f in glob.iglob('*.asc'):
        print f
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

        print start_col, start_row, end_col, end_row        
        
        cur_rawgrid = np.loadtxt(f, skiprows=6)
        print cur_rawgrid.shape
        # perhaps make sure that -1 values and NODATA values are masked or made = 0 
        #(should not be issue over land but better safe than sorry)
        cur_croppedrain = cur_rawgrid[start_row:end_row, start_col:end_col]###/32   division done in NimPy now 
        print cur_croppedrain.shape
        cur_rainrow = np.hstack(cur_croppedrain)
        #print cur_rainrow
        rainfile.append(cur_rainrow)
    
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

# The np.savetxt function takes a string header as an argumnet, use this to
# convert the array into a long string
def convert_array_hdr_to_str_hdr(array_header):
    str_header = 'NCols ' + str(array_header[0]) + '\n' + 'NRows ' + str(array_header[1]) + '\n' + 'xllcorner ' + str(array_header[2]) + '\n' + 'yllcorner ' + str(array_header[3]) + '\n' + 'cellsize ' + str(array_header[4]) + '\n' + 'NODATA_value ' + str(array_header[5])
    return str_header


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
write_sample_radar_img()

##### TO DO #####
# Write sample raster cropped ascii with correct header info
# georef stuff can come from the basin_header - is the same spatial extent by definition
# plot overlay of radar onto terrain dem
# create total radar accumultion plots
    
