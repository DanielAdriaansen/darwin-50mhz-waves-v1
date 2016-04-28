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

#################################################################################

# Get a list of the files available
flist = glob.glob('%s/*.nc' % data)

# Open a single file to read variables that don't change between files 
sf = nc.Dataset(flist[0])
mf = nc.MFDataset(flist, aggdim='time')

# Store data from the single file
bf = sf.variables['bad_flag']
palt = sf.variables['prof_site_altitude']
plat = sf.variables['prof_site_latitude']
plon = sf.variables['prof_site_longitude']
pagl = sf.variables['prof_height_100m_AGL']

# Get data from all files
y = mf.variables['prof_time_year'][:]
m = mf.variables['prof_time_month'][:]
jd = mf.variables['prof_time_dayofyear'][:]
d = mf.variables['prof_time_dayofmonth'][:]
h = mf.variables['prof_time_hour'][:]
mn = mf.variables['prof_time_hour'][:]

# Read in 50 MHz data from all files
omean = mf.variables['prof_dar50_Omega_mean'][:]
ouncert = mf.variables['prof_dar50_Omega_mean_uncert'][:]
osnr = mf.variables['prof_dar50_Omega_snr'][:]
owid = mf.variables['prof_dar50_Omega_wid'][:]
ozdb = mf.variables['prof_dar50_Omega_zdb'][:]

omean2 = np.ma.masked_values(omean,bf)

