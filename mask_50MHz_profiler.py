#!/usr/bin/env python
#
# File: mask_50MHz_profiler.py
#
# Author: D. Adriaansen
#
# Date: 28 Apr 2016
#
# Purpose: Read in 50 MHz profiler data and mask it using 920 MHz profiler data
#
# Notes:
#________________________________________________________________________________________

# Python libs
import netCDF4 as nc
import numpy as np
import glob
import matplotlib.pyplot as plt

################################## User Config ##################################

# Path to profiler data files
data = "/d1/dadriaan/paper/20141008_Christopher_Williams/netcdf"

# TEST 1: z > 4 km && valid dbz
t1 = True
t1z = 4000.0

# TEST 2: abs(dopvel) > 1.0 m/s
t2 = False
t2v = 10.0

# TEST 3: dbz > 20.0 dbz
t3 = False
t3r = 20.0

# TEST 4: abs(omega) > 1.0 m/s and Z < 8 km
t4 = False
t4v = 5.0
t4z = 8000.0

#################################################################################

# Get a list of the files available
flist = glob.glob('%s/*.nc' % data)

# Open a single file to read variables that don't change between files 
sf = nc.Dataset(flist[0])
mf = nc.MFDataset(flist, aggdim='time')

# Store data from the single file
bf = sf.variables['bad_flag'][:]
palt = sf.variables['prof_site_altitude'][:]
plat = sf.variables['prof_site_latitude'][:]
plon = sf.variables['prof_site_longitude'][:]
pagl = sf.variables['prof_height_100m_AGL'][:]

# Get data from all files
y = mf.variables['prof_time_year'][:]
m = mf.variables['prof_time_month'][:]
jd = mf.variables['prof_time_dayofyear'][:]
d = mf.variables['prof_time_dayofmonth'][:]
h = mf.variables['prof_time_hour'][:]
mn = mf.variables['prof_time_hour'][:]

# Read in 50 MHz data from all files
omean = mf.variables['prof_dar50_Omega_mean'][:]
zmean = mf.variables['prof_dar50_Zonal_mean'][:]
mmean = mf.variables['prof_dar50_Merid_mean'][:]

# Read in 920 MHz data from all files
raindbz = mf.variables['prof_dar920_vert_zdb'][:]

# Create arrays with bad data turned off (masked, using the value of "bf")
omean2 = np.ma.masked_values(omean,bf)
raindbz2 = np.ma.masked_values(raindbz,bf)

# Create a new numpy array to hold mask information, initialize to 1
mask = np.empty_like(omean)
mask[::] = 1

# Store the height as a 3D variable instead of a single vector
heights = np.empty_like(omean)
heights[::] = pagl

# Mark the bad data locations with a 2
mask[omean == bf] = 2

# Perform test 1 using 920 data
t1mask1 = (heights > t1z)
t1mask2 = (~raindbz2.mask)
combo = np.logical_and(t1mask1,t1mask2)
mask[combo] = 3

# Now, use the mask array to set any data to "bf" where it's not 1
good_vals = (mask == 1)

# Why are these a vector and not a 2D thing?
goodw = omean[good_vals]
goodp = raindbz[good_vals]
print(goodw.shape)

