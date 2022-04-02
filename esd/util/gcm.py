import numpy as np
import pandas as pd
import xarray as xr
# import datetime as dt


def unit_convert_pr(da):
    # Convert from kg m-2 s-1 to mm per day
    assert da.units == 'kg m-2 s-1'
    # Check for possible invalid negative values from bicubic interpolation
    da.values[da.values < 0] = 0
    da = da*86400
    da.attrs['units'] = 'mm day-1'
    return da

def unit_convert_tas(da):
    
    # Convert from Kelvin to Celsius
    assert da.units == 'K'
    da = da - 273.15
    da.attrs['units'] = 'C'
    return da

def unit_convert_none(da):
    return da

def clean_gcm_da(da, unit_convert_func=unit_convert_none, startend_dates=None):
    da = da.sel(lat=slice(1,13), lon=slice(29,42))
    # Get rid of duplicated days
    da = da.isel(time=np.nonzero((~da.time.to_pandas().duplicated()).values)[0])

    # # convert dates
    # datetimeindex = da.indexes['time'].to_datetimeindex()
    # da['time'] = datetimeindex

    # Get rid of days > 2101 due to pandas Timestamp range limitation
    # da = da.sel(time=da.time[da.time <= pd.to_datetime("2101-01-01")])
    
    # Perform unit conversion
    da = unit_convert_func(da)
    
    # Convert times to datetime format
    # times = pd.to_datetime(da['time'].dt.strftime("%Y%m%d").to_pandas(),   #.astype(np.str),
    #                         format='%Y%m%d.5', errors='coerce')   
    # print(times)     
    # da['time'] = times.values
    # Get rid of invalid dates (occurs with 360-day calendar)
    datetimeindex = da.indexes['time'].strftime("%Y%m%d").values.tolist()
    datetimeindex = pd.to_datetime(datetimeindex,
                           format='%Y%m%d', errors='coerce')
    da['time'] = datetimeindex


    # Get rid of invalid dates (occurs with 360-day calendar)
    da = da.isel(time=np.nonzero(da.time.notnull().values)[0])
    times = da.time.to_pandas()
    
    # Full sequence of valid dates
    # Todo: get 12-31
    first_day = times.iloc[0]
    last_day = times.iloc[-1]
    
    if (last_day.day == 30) and (last_day.month==12):
        last_day = last_day + pd.Timedelta(days=1)
    
    times_full = pd.date_range(first_day, last_day)
    
    # Reindex to full sequence and place NA for missing values
    da = da.reindex(time=times_full)
    # da = da.to_array()
    s = da.to_series() #index: time, lat, lon
    s = s.reorder_levels(['lat','lon','time'])
    s = s.sort_index(0, sort_remaining=True)
    s = s.unstack(['lat','lon'])
    s = s.interpolate(method='time')
    da_cleaned = xr.DataArray.from_series(s.stack(['lat','lon']))
    da_cleaned.attrs = da.attrs
    
    if startend_dates is not None:
        da_cleaned = da_cleaned.loc[startend_dates[0]:startend_dates[-1]]
    
    return da_cleaned