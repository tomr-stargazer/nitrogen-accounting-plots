"""
This will load and display the ALMA science verification (SV) data of IRAS 16293 in Band 9, showing the H13CN 8-7 line.
I intend for the source size to have been estimated by-hand, in CASA using a polygon select tool...
https://casaguides.nrao.edu/index.php?title=IRAS16293_Band9_-_Imaging_for_CASA_4.2#Create_a_velocity_cube_for_H13CN_.28J.3D8-7.29

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from config import alma_data_path

def show_alma_h13cn_image():

