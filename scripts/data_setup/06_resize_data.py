# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:25:03 2021

@author: Sarah

Resize the resampled netCDF files so that they cover the
same extent (problems because of pre-downloaded pr and temp data)
"""

from esd.util import buffer_cdo_lonlat_grid_cfg, mkdir_p
import esd
import glob
import os
import subprocess
import xarray as xr 
import pandas as pd

if __name__ == '__main__':
    # start with cmip data
    path_in = os.path.join(esd.cfg.path_cmip6_resample, 'remapbic') #remapbil or remapbic
    paths_rcp = sorted(glob.glob(os.path.join(path_in, "historical*_merged")))
    path_out = esd.cfg.path_cmip6_resize
    
    for path_rcp in paths_rcp:
        
        fpaths_tasmax = sorted(glob.glob(os.path.join(path_rcp, 'tasmax', 'tasmax.day*')))
        fpaths_tasmin = sorted(glob.glob(os.path.join(path_rcp, 'tasmin', 'tasmin.day*')))
        fpaths_pr = sorted(glob.glob(os.path.join(path_rcp, 'pr', 'pr.day*')))

        fpaths_all = fpaths_tasmax + fpaths_tasmin + fpaths_pr
        
        path_out_rcp = os.path.join(path_out, os.path.basename(path_rcp))
        mkdir_p(path_out_rcp)
        
        for fpath_in in fpaths_all:
            temp = xr.open_dataset(fpath_in)
            var = os.path.basename(fpath_in).split('.')[0]
            temp = temp[var].sel(lat=slice(3,10), lon=slice(34,39), time=slice("1981-01-01", "2099-12-31"))

            fpath_out = os.path.join(path_out_rcp, var, os.path.basename(fpath_in))
            mkdir_p (os.path.join(path_out_rcp, var))
            
            temp.to_netcdf(fpath_out)
    
    # now obs 
    path_in_obs = os.path.join(esd.cfg.path_obs_resample)
    obs_files = sorted(glob.glob(os.path.join(path_in_obs, '*.nc')))
    
    for f in obs_files:
        o = xr.open_dataset(f)
        var = os.path.basename(f).split('_')[0]
        if var == 'pcp':
            var = 'precip'
        if var[0:3] == 'tas':
            var = 't2m_C'
        o = o[var].sel(lat=slice(3,10), lon=slice(34,39))
        
        fpath_out_obs = os.path.join(esd.cfg.path_obs_resize, os.path.basename(f))
        o.to_netcdf(fpath_out_obs)
