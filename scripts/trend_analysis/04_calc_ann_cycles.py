'''
Script to calculate historical annual cycles of tair and prcp for 
APHRODITE and CMIP5 realizations.
'''
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
sys.path.append('/home/smj5vup/.local/bin')

import esd
import glob
import itertools
import numpy as np
import os
import pandas as pd
import xarray as xr

if __name__ == '__main__':
    
    # Historical annual cycles for aphrodite 1deg and CMIP5 realizations that
    # have NOT been biased corrected and downscaled
    path_in_cmip5 = esd.cfg.path_cmip6_cleaned
    fpath_aphrodite_tmin = os.path.join(esd.cfg.fpath_obs_tasmin)
    fpath_aphrodite_tmax = os.path.join(esd.cfg.fpath_obs_tasmax)
    fpath_aphrodite_prcp = os.path.join(esd.cfg.fpath_obs_pr)

    fpath_out_tmin = os.path.join(esd.cfg.path_cmip6_trends, 'hist_ann_cycle_tmin.csv')
    fpath_out_tmax = os.path.join(esd.cfg.path_cmip6_trends, 'hist_ann_cycle_tmax.csv')
    fpath_out_prcp = os.path.join(esd.cfg.path_cmip6_trends, 'hist_ann_cycle_prcp.csv')
    
    # Historical annual cycles for original aphrodite .25deg and CMIP5 realizations that
    # have been biased corrected and downscaled
    # path_in_cmip5 = esd.cfg.path_cmip5_downscaled
    # fpath_aphrodite_tair = os.path.join(esd.cfg.path_aphrodite_resample,
    #                                     'aprhodite_redriver_sat_1961_2007_p25deg.nc')
    # fpath_aphrodite_prcp = os.path.join(esd.cfg.path_aphrodite_resample,
    #                                     'aprhodite_redriver_pcp_1961_2007_p25deg.nc')
    # fpath_out_tair = os.path.join(esd.cfg.path_cmip5_trends, 'hist_ann_cycle_tair_p25deg.csv')
    # fpath_out_prcp = os.path.join(esd.cfg.path_cmip5_trends, 'hist_ann_cycle_prcp_p25deg.csv')
    
    
    ############################################################################
    
    paths_in = sorted(glob.glob(os.path.join(path_in_cmip5, 'historical*_merged', '*')))
        
    fpaths_tmin = sorted(list(itertools.chain.
                             from_iterable([glob.glob(os.path.join(apath,'tasmin.day*'))
                                            for apath in paths_in])))
    fpaths_tmax = sorted(list(itertools.chain.
                             from_iterable([glob.glob(os.path.join(apath,'tasmax.day*'))
                                            for apath in paths_in])))
    fpaths_pr = sorted(list(itertools.chain.
                            from_iterable([glob.glob(os.path.join(apath,'pr.day*'))
                                           for apath in paths_in])))

    
    ds = xr.open_dataset(fpaths_tmin[0])
    mask_lat = np.logical_and(ds.lat.values >= esd.cfg.bbox[1],
                              ds.lat.values <= esd.cfg.bbox[-1])
    mask_lon = np.logical_and(ds.lon.values >= esd.cfg.bbox[0],
                              ds.lon.values <= esd.cfg.bbox[2]) 
    ds.close()
    
    def mth_tair_norms(fpath, vname='tas', apply_mask=True,name=None, obs=False):
        if obs == False:
            ds = xr.open_dataset(fpath)
        else:
            ds = fpath

        if apply_mask:
            da = ds[vname].loc['1981':'2005',mask_lat,mask_lon].load()
        else:
            da = ds[vname].loc['1981':'2005',:,:].load()

        if obs == False:
            da =  da.sel(lat=slice(4,10),lon=slice(32,39))
        else:
            da = da.sel(latitude=slice(10,4), longitude=slice(32,39))
            
        da_mthly = da.resample(time='MS').mean()

        if obs == False: 
            norms = da_mthly.groupby('time.month').mean(dim='time').mean(dim=['lat','lon'])
        else:
            norms = da_mthly.groupby('time.month').mean(dim='time').mean(dim=['latitude','longitude']) # observed data has different names

        norms = norms.to_pandas()
        
        if name is None:
            name = os.path.basename(fpath)
        norms.name = name
        ds.close()
        
        print(norms.name)
        
        return norms
    
    
    norms_tmin = [mth_tair_norms(fpath, 'tasmin', False, None) for fpath in fpaths_tmin]
    norms_tmin = pd.concat(norms_tmin, axis=1)

    norms_tmax = [mth_tair_norms(fpath, 'tasmax', False, None) for fpath in fpaths_tmax]
    norms_tmax = pd.concat(norms_tmax, axis=1)
    
    def mth_prcp_norms(fpath, vname='pr', apply_mask=False,name=None, obs=False):
        if obs == False:
            ds = xr.open_dataset(fpath)
        else:
            ds = fpath


        if apply_mask:
            da = ds[vname].loc['1981':'2005',mask_lat,mask_lon].load()
        else:
            da = ds[vname].loc['1981':'2005',:,:].load()


        if obs == False:
            da =  da.sel(lat=slice(4,10),lon=slice(32,39))
        else:
            da = da.sel(latitude=slice(4,10), longitude=slice(32,39))

        da_mthly = da.resample(time='MS').sum()

        if obs == False: 
            norms = da_mthly.groupby('time.month').mean(dim='time').mean(dim=['lat','lon'])
        else:
            norms = da_mthly.groupby('time.month').mean(dim='time').mean(dim=['latitude','longitude']) 
        norms = norms.to_pandas()
        if name is None:
            name = os.path.basename(fpath)
        norms.name = name
        ds.close()
        
        print(norms.name)
        
        return norms
    
    norms_prcp = [mth_prcp_norms(fpath) for fpath in fpaths_pr]
    norms_prcp = pd.concat(norms_prcp, axis=1)
        
    ds_tmin = xr.open_dataset(fpath_aphrodite_tmin)
    ds_tmax = xr.open_dataset(fpath_aphrodite_tmax)
    ds_prcp = xr.open_dataset(fpath_aphrodite_prcp)
    
    def convert_times(ds):
        
        if ds.time.values.dtype.name != 'datetime64[ns]':
            # Convert times to datetime format
            times = pd.to_datetime(ds.time.to_pandas().astype(np.str),
                                   format='%Y%m%d.0', errors='coerce')        
            ds['time'] = times.values
            
        return ds
    
    ds_tmin = convert_times(ds_tmin)
    ds_tmax = convert_times(ds_tmax)
    ds_prcp = convert_times(ds_prcp)
    
    aphro_norms_tmin = mth_tair_norms(ds_tmin, 't2m', False, 'era5', True)
    aphro_norms_tmax = mth_tair_norms(ds_tmax, 't2m', False, 'era5', True)
    aphro_norms_prcp = mth_prcp_norms(ds_prcp, 'precip', False, 'chirps05', True)
    
    norms_prcp = pd.concat([norms_prcp, aphro_norms_prcp],axis=1)
    norms_tmax = pd.concat([norms_tmax, aphro_norms_tmax],axis=1)
    norms_tmin = pd.concat([norms_tmin, aphro_norms_tmin],axis=1)
    
    norms_prcp.to_csv(fpath_out_prcp)
    norms_tmin.to_csv(fpath_out_tmin)
    norms_tmax.to_csv(fpath_out_tmax)