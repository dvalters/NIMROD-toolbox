# -*- coding: utf-8 -*-
"""
This script is designed to create a hydroindex file from a clipped radar
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
"""
from __future__ import division
# because we want to get the answer 5/2 = 2.5 (Python 3), not 5/2 = 2 (Python 2.x)


import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
import math



"""
There's some existing raster tools in the LSD python files, you would have to
add the folder to your PYTHONPATH and then import as below if you wanted to 
do it this way:

import LSDVisualisation.trunk.raster_tools as rastertool
rastertool.read_ascii_header

you need __init__.py files in each directory.

"""

RADAR_DEM = "cropped_test_radar.asc"
TERRAIN_DEM = "ryedale20m_fillcrop.asc"


def read_ascii_header(ascii_raster_file):
    
    with open(ascii_raster_file) as f:
        header_data = [float(f.next().split()[1]) for x in xrange(6)] #read the first 6 lines
         
    #raster_data = np.genfromtxt(ascii_raster_file, delimiter=' ', skip_header=6)
    #raster_data = raster_data.reshape(header_data[1], header_data[0]) #rows, columns
    
    return header_data

def create_base_indexgrid():
    ncols = read_ascii_header(RADAR_DEM)[0]
    nrows = read_ascii_header(RADAR_DEM)[1]
    no_of_cells = ncols*nrows
    # The base grid is based on the coarse radar data clipped to the right area.
    # So it should only have a low number of cells really
    # +1 at the end because we want the numbering to start at 1 instead of zero 
    # for later purposes.
    base_indexgrid = np.arange(no_of_cells).reshape(nrows,ncols) + 1   
    
    return base_indexgrid
    
def create_upscaled_hydroindex(base_hydroindex):
    radar_cellsize = read_ascii_header(RADAR_DEM)[4]
    terrain_cellsize = read_ascii_header(TERRAIN_DEM)[4]
    # We want to scale up the base grid to the resolution of the terrain DEM
    # that will be used in the model run by a mult_factor
    #mult_factor = radar_cellsize/terrain_cellsize
    #print mult_factor
    radarNcols = read_ascii_header(RADAR_DEM)[0]
    radarNrows = read_ascii_header(RADAR_DEM)[1]
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

#=-=-=-=-=-
# MAIN
#=-=-=-=-=-

base_hydroindex = create_base_indexgrid()
upscaled_hydroindex = create_upscaled_hydroindex(base_hydroindex)

hydroindex_masked = crop_upscaled_hydroindex_to_basin(upscaled_hydroindex, TERRAIN_DEM)


#plt.imshow(terrain_mask_test, interpolation='none')
#plt.imshow(upscaled_hydroindex, interpolation='none')

#print base_hydroindex
#print upscaled_hydroindex

write_hydroindex_to_file("ryedale_hydroindex_test2.asc", hydroindex_masked)






