#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 20:17:45 2014

This tool creates ascii rasters of Met Office 'NIMROD' rainfall radar
data. It can be used to convert a single NIMROD .dat file (These are binary,
proprietary format files) or an entire directory of NIMROD files. The default
output format is ASCII raster (.asc) but this could be adapted to suit the user's
purpose.


@author: Declan Valters
@author: Adapted from Charles Kilburn's script 2008
"""
# Read a NIMROD file
# Look up all the meanings of the header values in the NIMORD format document available at the BADC
# Charles Kilburn Aug 2008
# Extended: Declan Valters Apr 2015

import os, stat, re, sys, time
import struct
import array
import numpy as np
import matplotlib as mpl
import pyproj 
from glob import glob
    
# This is to transform the osgb coords into UTM ones.
def convert_OSGB36_to_UTM30(xcoord, ycoord):
    # Set up the coord systems
    OSGB36 = pyproj.Proj("+init=EPSG:27700")
    UTM30N = pyproj.Proj("+init=EPSG:32630")
    utm_x, utm_y = pyproj.transform(OSGB36, UTM30N, xcoord, ycoord)
    return utm_x, utm_y
    
# check to see if gride cell resolution is same in both dimensions
def check_horiz_vert_resolution(row_interval_res, col_interval_res):
    if row_interval_res != col_interval_res:
        print("Warning: the row resolution ", row_interval_res, "is not the \
        same resolution as the column resolution",  col_interval_res, " in this \
        dataset. Row resolution has been used to set the ASCII 'cellsize'. Check your data.")
        
#write the array and header to an ascii file
def write_NIMROD_toASCII(asciiname, radararray, header):
    np.savetxt(asciiname, radararray, header=header, comments='', fmt='%1.1f')

#Horizontal grid types
grid_types = {'0':"British National Grid", \
'1':"Lat/Long", \
'2':"Space View", \
'3':"Polar Stereo", \
'4':"UTM32", \
'5':"Rotated Lat/Lon", \
'6':"Other"}

def ingest_NIMROD_file(full_filepath):
    file_id = open(full_filepath,"rb")
    record_length, = struct.unpack(">l", file_id.read(4))
    if record_length != 512: 
        raise( "Unexpected record length", record_length)
    
    # Set up the arrays, with the correct type
    gen_ints = array.array("h")
    gen_reals = array.array("f")
    spec_reals = array.array("f")
    characters = array.array("c")
    spec_ints = array.array("h")
    
    # read in the data from the open file
    gen_ints.read(file_id, 31) 
    gen_reals.read(file_id, 28)
    spec_reals.read(file_id, 45)
    characters.read(file_id, 56)
    spec_ints.read(file_id, 51)
    
    gen_ints.byteswap()
    gen_reals.byteswap()
    spec_reals.byteswap()
    spec_ints.byteswap()
    
    row_interval_res = gen_ints[3]
    col_interval_res = gen_ints[5]
    
    date_yy = gen_ints[0]
    date_mm = gen_ints[1]
    date_dd = gen_ints[2]
    time_hh = gen_ints[3]
    time_mm = gen_ints[4]
    grid_rows = gen_ints[15]
    grid_cols = gen_ints[16]
    
    grid_typeID = gen_ints[14]
    grid_type = grid_types[str(grid_typeID)] #(Dictionaries need a string to look up)
    print(grid_type)
    
    record_length, = struct.unpack(">l", file_id.read(4))
    if record_length != 512: 
        raise( "Unexpected record length", record_length)
        
    chars = characters.tostring()
    
    units_rain = chars[0:8]
    radar_source = chars[8:32]
    radar_data_type = chars[32:55]
    
    #Read the Data
    array_size = gen_ints[15] * gen_ints[16]
    nrows = gen_ints[15]
    ncols = gen_ints[16]
    #xllcornerNG = spec_reals[7]
    #yllcornerNG = spec_reals[6]  # 'NG' = British National Grid, or 'OSGB36' to be precise
    cellsize = gen_reals[3]
    nodata_value = gen_reals[6]
    
    ## Special Case when the extended header has NODATA values. ##
    # In this case the header does not contain the xll/yll-corner values for ASCII files
    # Instead, you are given the top left corner, and corresponding x-value for this
    # In effect ytlcorner and xtlcorner.
    # So to get yllcorner = ytlcorner - nrows*resolution
    ytlcorner = gen_reals[2] # get the top left corner y co-ord
    xtlcorner = gen_reals[4] # get the top left corner x co-ord
    # calculate the lowerleft corner from the top left corner
    yllcornerNG = ytlcorner - (nrows*cellsize)
    xllcornerNG = xtlcorner # they are effectively the same thing...left-most x coordinate in the grid
    
    xllcornerUTM, yllcornerUTM = convert_OSGB36_to_UTM30(xllcornerNG, yllcornerNG)
    
    #Note if you use the data in spec_reals, the co-ordnates are 500m apart...probably not big enough to worry about    
    record_length, = struct.unpack(">l", file_id.read(4))
    if record_length != array_size * 2: 
        raise( "Unexpected record length", record_length)
    
    data = array.array("h")
    try:
        data.read(file_id, array_size) # read() is deprecated. And it only reads linearly through a file
        record_length, = struct.unpack(">l", file_id.read(4))
        if record_length != array_size * 2: raise( "Unexpected record length", record_length)
        data.byteswap()
        #print "First 100 values are", data[:100]
        #print "Last 100 values are", data[-100:]
    except:
        print( "Read failed")
        
    radararray = np.reshape(data, (nrows,ncols)) 
    radararray = radararray / 32.0 # This is due to a strange NIMROD convention where everything is *32
    
    # Make the header
    header = 'NCols ' + str(ncols) + '\n' + 'NRows ' + str(nrows) + '\n' + 'xllcorner ' + str(xllcornerUTM) + '\n' + 'yllcorner ' + str(yllcornerUTM) + '\n' + 'cellsize ' + str(cellsize) + '\n' + 'NODATA_value ' + str(nodata_value)
    
    file_id.close()
    return radararray, header
    
def convert_multiple_files(path_to_dir, basename):
    # The base name will should be appended to a wildcard to narrow down the search.
    for filename in glob(path_to_dir + basename):
        print( filename)
        thisradararray, thisheader = ingest_NIMROD_file(filename)
        asciiname = filename + '.asc'
        write_NIMROD_toASCII(asciiname, thisradararray, thisheader)
    
def convert_single_file(path, fname, ext, asciiname):
    full_fname = path + fname + '.' + ext
    radararray, header = ingest_NIMROD_file(full_fname)
    write_NIMROD_toASCII(asciiname, radararray, header)
    print("NIMROD file converted to ASCII: ", asciiname)
    
def show_radar_image(full_filepath):
    #full_fname = path + fname + '.' + ext
    radararray, header = ingest_NIMROD_file(full_filepath)
    mpl.pyplot.imshow(radararray)
    
########
# MAIN #
########  

# File name and extension
#fname = '201201010305_nimrod_ng_radar_rainrate_composite_1km_UK'  # Testing the current files
fpath = ''
fname = "201201160000_nimrod_ng_radar_rainrate_composite_1km_UK"
file_ext = ''
asciiname = "NIMtest.asc"

# for multiple files. You do not need to specify a basename if you want to conver ALL files in dir.
basename = "*composite.dat"

convert_single_file(fpath, fname, file_ext, asciiname)
#convert_multiple_files(fpath, basename)


# To do:
# Improve the plotting function (base map of UK, mask the nodata values)
# In multiple file converter, strip the file extension before appending ".asc"










