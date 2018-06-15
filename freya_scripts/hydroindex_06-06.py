
from __future__ import division, print_function
# because we want to get the answer 5/2 = 2.5 (Python 3), not 5/2 = 2 (Python 2.x)

# -*- coding: utf-8 -*-
'''
This script is designed to create a hydroindex file from a clipped radar
file in ascii format. Hydroindex files are required for use in the LSDCatchmentModel
and HAIL-CAESAR models when investigating spatially variable rainfall in flooding and landscape evolution.

Run this script in the same directory as a clipped (and correctly georeferenced)
radar ascii data file (obtained from timeseries.py). 

The upsampling resolution must be set 
(the grid cell resolution of the corresponding DEM that will be used in the model simulation).
I.e. if DEM is at 10m resolution, the upscaling resolution = 10

The script has been altered from work by D. Valters (March 2015)

F. Muir, May 2018
'''

'''
--- FUNCTION AND PACKAGE IMPORTS ---
'''



import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
import math

'''
--- INPUT DATA NAMES ---
'''

CROPPED_RADAR_DEM = "/Users/freyamuir/Documents/Masters/MDiss/primary/max/small_cropped_radar_test.asc"
TERRAIN_DEM = "/Users/freyamuir/Documents/Masters/MDiss/primary/max/extra/comonwsh_2.asc"

'''
--- DATA PREPROCESSING FUNCTIONS ---
'''

def read_ascii_header(ascii_raster_file):
    
    with open(ascii_raster_file) as f:
        header_data = [float(f.next().split()[1]) for x in range(6)] #read the first 6 lines
         
    #raster_data = np.genfromtxt(ascii_raster_file, delimiter=' ', skip_header=6)
    #raster_data = raster_data.reshape(header_data[1], header_data[0]) #rows, columns
    
    return header_data

''' 
--- BASE INDEXGRID CREATION ---
'''
# Creates the base hydroindex grid from a cropped radar image tile
def create_base_indexgrid(CROPPED_RADAR_DEM):
    ncols = read_ascii_header(CROPPED_RADAR_DEM)[0]
    nrows = read_ascii_header(CROPPED_RADAR_DEM)[1]
    no_of_cells = ncols*nrows
    # The base grid is based on coarse radar data clipped to the right area.
    # Should have a low number of cells at the end (need the numbering to start at 1 instead of zero).
    base_indexgrid = np.arange(int(no_of_cells)).reshape(int(nrows),int(ncols)) + 1   
    
    return base_indexgrid

'''
--- HYDROINDEX GRID CREATION ---
'''

# Resamples the base indexgrid to the resolution of the DEM, and returns the upscaled hydroindex
def create_upscaled_hydroindex(base_hydroindex, CROPPED_RADAR_DEM, TERRAIN_DEM):
    radar_cellsize = read_ascii_header(CROPPED_RADAR_DEM)[4]
    terrain_cellsize = read_ascii_header(TERRAIN_DEM)[4]
    # Want to scale up the base indexgrid to the resolution of the terrain DEM by a mult_factor
    radarNcols = read_ascii_header(CROPPED_RADAR_DEM)[0]
    radarNrows = read_ascii_header(CROPPED_RADAR_DEM)[1]
    terrainNcols = read_ascii_header(TERRAIN_DEM)[0]
    terrainNrows = read_ascii_header(TERRAIN_DEM)[1] 
    mult_factor_cols = terrainNcols / radarNcols
    mult_factor_rows = terrainNrows / radarNrows
    
    print(mult_factor_cols, mult_factor_rows)
    
    # Resampling using nearest neighbour (order = 0)
    # The mult_factor can have different values for each axis, and can be float values
    upscaled_hydroindex = scipy.ndimage.zoom(base_hydroindex, (mult_factor_rows, mult_factor_cols), order=0)
    print(upscaled_hydroindex.shape)
    
    return upscaled_hydroindex

# Crops the upscaled hydroindex grid to terrain basin (only needed if using a clipped DEM)
def crop_upscaled_hydroindex_to_basin(upscaled_hydroindex, TERRAIN_DEM): 
    # Loads the terrain DEM
    terrain_array = np.loadtxt(TERRAIN_DEM, skiprows=6)
    terrain_NODATA = read_ascii_header(TERRAIN_DEM)[5]
    # Creates a mask from the terrain DEM
    terrain_nodata_mask = np.ma.masked_where(terrain_array<-1, terrain_array)
    #plt.imshow(terrain_nodata_mask)
    
    # Applies the mask to the hydroindex grid
    cropped_hydroindex = np.ma.masked_array(upscaled_hydroindex, terrain_nodata_mask.mask)

    crop_with_nodataval_hydroindex = np.ma.filled(cropped_hydroindex, terrain_NODATA)
    plt.imshow(crop_with_nodataval_hydroindex)
    
    return crop_with_nodataval_hydroindex
    
# Writes the hydroindex grid to an ascii DEM with the same header as the terrain DEM
def write_hydroindex_to_file(hydroindex_filename, hydroindex_masked):
    #hydroindex_header = np.zeros(6)
    header = "ncols " + str(hydroindex_masked.shape[1]) + "\n"
    header += "nrows " + str(hydroindex_masked.shape[0]) + "\n"
    header += "xllcorner " + str(read_ascii_header(TERRAIN_DEM)[2]) + "\n"
    header += "yllcorner " + str(read_ascii_header(TERRAIN_DEM)[3]) + "\n"
    header += "cellsize " + str(read_ascii_header(TERRAIN_DEM)[4]) + "\n"
    header += "NODATA_value " + str(read_ascii_header(TERRAIN_DEM)[5])
    # Hydroindex must be integer
    np.savetxt(hydroindex_filename, hydroindex_masked, fmt="%1i", header=header, comments='')

'''
--- EXECUTE FUNCTIONS ---
'''

# ASSIGN BASE HYDROINDEX FUNCTION RESULTS
base_hydroindex = create_base_indexgrid(CROPPED_RADAR_DEM)

# ASSIGN UPSCALED HYDROINDEX FILE
upscaled_hydroindex = create_upscaled_hydroindex(base_hydroindex, CROPPED_RADAR_DEM, TERRAIN_DEM)

# ASSIGN MASKED HYDROINDEX FILE
hydroindex_masked = crop_upscaled_hydroindex_to_basin(upscaled_hydroindex, TERRAIN_DEM)

# CREATE HYDROINDEX FILE
write_hydroindex_to_file("hurricanemax_small_hydroindex.asc", hydroindex_masked)

