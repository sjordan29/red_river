'''
Script for and prototyping and testing temperature downscaling approaches.
'''
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/downscale/')
from esd.downscale.tair import setup_data_for_analog, downscale_analog_anoms,TairDownscale
import esd
import numpy as np
import os
import pandas as pd
import xarray as xr
from esd.downscale.common import _convert_times

if __name__ == '__main__':
    
    fpath_tair_obs = esd.cfg.fpath_obs_tasmax
    fpath_tair_obsc = os.path.join(esd.cfg.path_obs_resample,
                                   'tasmax_ERA5_1981_2019_p10deg_remapbic.nc')
    fpath_tair_mod = os.path.join('/scratch/smj5vup/CMIP6/Debias//cmip6_debiased/resampled/tasmax/tasmax.day.ACCESS1-0.rcp45.r1i1p1.19810101-20991231.historical_merged.nc')
    
    base_start_year = '1981'
    base_end_year = '2018'
    downscale_start_year = '1981'
    downscale_end_year = '2030'
        
#     base_start_year = '1978'
#     base_end_year = '2007'
#     downscale_start_year = '1961'
#     downscale_end_year = '1977'
         
    
    tair_d = TairDownscale(fpath_tair_obs, fpath_tair_obsc, base_start_year,
                           base_end_year, base_start_year, base_end_year,
                           [(downscale_start_year,downscale_end_year)])


    ds = xr.open_dataset(fpath_tair_mod, decode_cf=False)
    vname = 'tasmax' # variable name
    ds[vname].attrs.pop('missing_value')
    ds = xr.decode_cf(ds)
    latlon_dict_ds = {'latitude': ds.latitude.round(3), 'longitude': ds.longitude.round(3)}
    ds = ds.assign_coords(latlon_dict_ds)
    da = ds[vname].load()
    da_ds = tair_d.downscale(da)
    
    # Add metadata and create dataset
    da_ds.attrs = da.attrs
    ds_out = da_ds.to_dataset(name=vname)
    ds_out.attrs = ds.attrs

    fname = 'test.nc'
    path_out = '/sfs/lustre/bahamut/scratch/smj5vup/CMIP6/omoCMIP/scripts/downscale/'
    fpath_out = os.path.join(path_out, fname)
    ds_out.to_netcdf(fpath_out)
