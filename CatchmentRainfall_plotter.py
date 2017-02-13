# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 09:16:08 2015

@author: dav
"""

import numpy as np
import matplotlib.pyplot as plt

rainfall_totals = np.loadtxt("/mnt/WORK/Dev/PyToolsPhD/Radardata_tools/rainfall_totals_boscastle_downscaled.asc")
#rainfall_totals = rainfall_totals/12
# Have you accounted for 5 minute data? If not, divide above by 12.

def read_ascii_raster(ascii_raster_file):
    import numpy as np
    
    with open(ascii_raster_file) as f:
        header_data = [float(f.next().split()[1]) for x in xrange(6)] #read the first 6 lines
         
    raster_data = np.genfromtxt(ascii_raster_file, delimiter=' ', skip_header=6)
    raster_data = raster_data.reshape(header_data[1], header_data[0]) #rows, columns
    
    return raster_data, header_data

def simple_radar_totals_plot():
    plt.imshow(rainfall_totals, interpolation="None", alpha=1, cmap="cool")
    plt.colorbar()
    plt.ylabel("Distance (m)")
    plt.xlabel("Distance (m)")
    plt.title("Cumulative rainfall (mm/72hrs)")


def cumulative_rainfall_catchment(hillshade_file, radar_data_totals):
    """
    Plots the catchment hillshade and overlays the total rainfalls accumulated
    during the model run.
    """
    label_size = 20
    #title_size = 30
    axis_size = 28

    import matplotlib.pyplot as pp
    import numpy as np
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    from matplotlib import rcParams
    import matplotlib.lines as mpllines
    
    #get data
    #hillshade, hillshade_header = read_flt(hillshade_file)
    
    hillshade, hillshade_header = read_ascii_raster(hillshade_file)
    rainfall_totals = np.loadtxt(radar_data_totals)
    
    #ignore nodata values    
    hillshade = np.ma.masked_where(hillshade == -9999, hillshade)    
    
    #fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size      
    
    fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
    ax = fig.add_subplot(1,1,1)
    
    plt.imshow(hillshade, vmin=0, vmax=255, cmap=cmx.gray)
    plt.imshow(rainfall_totals, interpolation="none", alpha=0.2)
