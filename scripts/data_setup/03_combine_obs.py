# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:48:58 2021

@author: Sarah

Combine precip and temp netCDFs into single  netCDF
Output will be used for debiasing and downscaling 
"""

import os
import pandas as pd
import requests
import xarray as xr
import scipy


# precipitation
ds = xr.open_mfdataset('D:/Research/SWAT/swat-omo/data/rainfall/CHIRPS/0.05-degree/data/nc-files/*.nc',combine = 'by_coords')
da = ds.precip.sel(latitude=slice(4,10), longitude=slice(34,39)) # units are mm

da.to_netcdf('D:/Research/Debias/obs/pcp_CHIRPS_1981_2020.nc')


# tasmin 
tmin = xr.open_mfdataset('D:/Research/SWAT/swat-omo/data/temperature/nc-files/resampled/min/*.nc', combine='by_coords')
# tmin_da = tmin.t2m_C.sel(latitude=slice(4,10), longitude=slice(34,39)) # units are degrees celsius 

tmin.t2m_C.to_netcdf('D:/Research/Debias/obs/tasmin_ERA5_1981_2019.nc')


# tas max
tmax = xr.open_mfdataset('D:/Research/SWAT/swat-omo/data/temperature/nc-files/resampled/max/*.nc', combine='by_coords')
# tmax_da = tmax.t2m_C.sel(latitude=slice(4,10), longitude=slice(34,39)) # units are degrees celsius 

tmax.t2m_C.to_netcdf('D:/Research/Debias/obs/tasmax_ERA5_1981_2019.nc')

# combine climate 
access = xr.open_mfdataset('D:/Research/Debias/cmipProjs/test/*.nc', combine='by_coords')
access.pr.to_netcdf('D:/Research/Debias/cmipProjs/test/pr_ssp245_combo.nc')
