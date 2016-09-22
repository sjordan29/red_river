'''
Functions for downscaling temperature and validation of the downscaled results.
'''

import numpy as np
import xarray as xr
from .common import _convert_times, _window_masks

def downscale_analog_anoms(ds, win_masks):
    
    da_mod_anoms = ds.mod_anoms.loc[ds.downscale_start_year:ds.downscale_end_year]
    da_obsc_anoms = ds.obsc_anoms.loc[ds.base_start_year:ds.base_end_year]
    da_obs_anoms = ds.obs_anoms
    da_obs_clim = ds.obs_clim
    
    # DataArray to store downscaled results
    da_mod_anoms_d = da_mod_anoms.copy()
    
    for a_date in da_mod_anoms.time.to_pandas().index:
        
        print a_date
        
        vals_mod = da_mod_anoms.loc[a_date]
        analog_pool = da_obsc_anoms[win_masks.loc[a_date.strftime('%m-%d')].values]        
        rmse_analogs = np.sqrt((np.square(vals_mod - analog_pool)).mean(dim=('lon', 'lat')))
        vals_analog = analog_pool[int(rmse_analogs.argmin())]
        
        s = vals_mod - vals_analog
        vals_obs_anoms = da_obs_anoms.loc[vals_analog.time.values]
        vals_d = vals_obs_anoms + s
                   
        da_mod_anoms_d.loc[a_date] = vals_d.values
    
    da_mod_d = da_mod_anoms_d.groupby('time.month') + da_obs_clim
    
    return da_mod_d

def to_anomalies(da_tair, base_startend_yrs=None):
    
    da_tair_mthly = da_tair.resample('MS', dim='time', how='mean', skipna=False)
    
    if base_startend_yrs is None:
        da_tair_clim = da_tair_mthly.groupby('time.month').mean(dim='time') 
    else:
        da_tair_clim = da_tair_mthly.loc[base_startend_yrs[0]:base_startend_yrs[1]].groupby('time.month').mean(dim='time')    
    
    da_tair_anoms = da_tair.groupby('time.month') - da_tair_clim
    da_tair_anoms = da_tair_anoms.drop('month')

    return da_tair_anoms,da_tair_clim

def setup_data_for_analog(fpath_tair_obs, fpath_tair_obsc, fpath_tair_mod,
                          base_start_year, base_end_year, downscale_start_year,
                          downscale_end_year):
    
    da_obs = xr.open_dataset(fpath_tair_obs).SAT
    da_obsc = xr.open_dataset(fpath_tair_obsc).SAT.load()
    da_mod = xr.open_dataset(fpath_tair_mod).SAT.load()

    da_obs = da_obs.loc[:, da_mod.lat.values, da_mod.lon.values].load()

    da_obs['time'] = _convert_times(da_obs)
    da_obsc['time'] = _convert_times(da_obsc)
    da_mod['time'] = _convert_times(da_mod)

    da_obs_anoms,da_obs_clim = to_anomalies(da_obs, (base_start_year,base_end_year))
    da_obsc_anoms = to_anomalies(da_obsc, (base_start_year,base_end_year))[0]
    da_mod_anoms = to_anomalies(da_mod, (base_start_year,base_end_year))[0]
        
    da_obsc_anoms = da_obsc_anoms.loc[base_start_year:base_end_year]
    da_obsc = da_obsc.loc[base_start_year:base_end_year]
    win_masks = _window_masks(da_obsc_anoms.time.to_pandas().index)
    
    da_obs.name = 'obs'
    da_obs_anoms.name = 'obs_anoms'
    da_obs_clim.name = 'obs_clim'
    da_obsc.name = 'obsc'
    da_obsc_anoms.name = 'obsc_anoms'
    da_mod.name = 'mod'
    da_mod_anoms.name = 'mod_anoms'
    
    ds = xr.merge([da_obs, da_obs_anoms, da_obs_clim, da_obsc, da_obsc_anoms, da_mod, da_mod_anoms])
    
    ds.attrs['base_start_year'] = base_start_year
    ds.attrs['base_end_year'] = base_end_year
    ds.attrs['downscale_start_year'] = downscale_start_year
    ds.attrs['downscale_end_year'] = downscale_end_year
    
    return ds, win_masks

def validate_avg_tair(mod_d, obs, months=None):
    
    mod_d_mthly = mod_d.resample('1MS', dim='time', how='mean', skipna=False)
    obs_mthly = obs.resample('1MS', dim='time', how='mean', skipna=False)
    
    if months is None:
        
        # Calculate annual avg 
        mod_d_avg = mod_d_mthly.resample('AS', dim='time', how='mean',skipna=False).mean(dim='time')
        obs_avg = obs_mthly.resample('AS', dim='time', how='mean',skipna=False).mean(dim='time')
        
    else:
                
        mask_season = np.in1d(mod_d_mthly['time.month'].values, months)
        mask_season = xr.DataArray(mask_season, coords=[mod_d_mthly.time])
        
        mod_d_mthly = mod_d_mthly.where(mask_season)
        obs_mthly = obs_mthly.where(mask_season)
        
        nmths = len(months)
        mod_d_roller = mod_d_mthly.rolling(min_periods=nmths, center=False, time=nmths)
        obs_roller = obs_mthly.rolling(min_periods=nmths, center=False, time=nmths)
        
        mod_d_avg = mod_d_roller.mean(skipna=False).mean(dim='time')
        obs_avg = obs_roller.mean(skipna=False).mean(dim='time')
        
    # Calculate error
    err = mod_d_avg-obs_avg
    obs_avg.name = 'obs'
    mod_d_avg.name = 'mod'
    err.name = 'err'
    
    return xr.merge([err, obs_avg, mod_d_avg])

def validate_daily_tair(mod_d, obs, months=None):
    
    if months is not None:
        
        mask_season = mod_d['time.month'].to_pandas().isin(months).values
        mod_d = mod_d[mask_season,:,:]
        obs = obs[mask_season,:,:]
        
    err = mod_d-obs
    mae = np.abs(err).mean(dim='time', skipna=False)
    bias = err.mean(dim='time', skipna=False)
    mae.name = 'mae'
    bias.name = 'bias'
        
    return xr.merge([mae, bias])

