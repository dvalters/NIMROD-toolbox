# -*- coding: utf-8 -*-
"""
Created on Thu May 12 08:09:02 2016

@author: dav

Class contents:

Convert NIMROD binary to ascii (single file, multiple files)
Create arrays from rainfall files (single, multiple)

Plot radar image for inspection
Extract subset of radar data (single file, multiple files)

Write a sample (cropped) radar image (ascii format)
Write multiple (cropped) radar images (ascii format)

Save as rainfall timeseries (spatially variable sones (caesar))
Save as rainfall timeseries (uniform, uniform weighted by pixel contribution to catchment)

Create a hydroindex file (CAESAR format)

Later:

Create artificial rainfall data from spatial and tempora distribution parameters.
Use NOAA return periods.


"""

import numpy as np

class Raingrid:
    def __init__(self):
        pass
    
    def from_file(self, rainfile):
        pass
    
    def from_NIMROD(self, nimrod_files):
        pass
    
    