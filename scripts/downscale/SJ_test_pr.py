import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/downscale/')


from esd.downscale.prcp import PrcpDownscale
from esd.util import Unbuffered, StatusCheck, mkdir_p
import esd
import numpy as np
import os
import pandas as pd
import xarray as xr
import esd
import glob
import itertools
import numpy as np
import os
import sys
import xarray as xr
from esd.downscale.common import _convert_times

sys.stdout = Unbuffered(sys.stdout)

print("hello!")

print("in main")



def get_cmip5_fpaths(a_path):
    
    paths_in = sorted(glob.glob(os.path.join(a_path, '*')))
        
    # fpaths_tasmin = sorted(list(itertools.chain.
    #                          from_iterable([glob.glob(os.path.join(apath,'tasmin.day*'))
    #                                         for apath in paths_in])))

    # fpaths_tasmax = sorted(list(itertools.chain.
    #                          from_iterable([glob.glob(os.path.join(apath,'tasmax.day*'))
    #                                         for apath in paths_in])))

    fpaths_pr = sorted(list(itertools.chain.
                            from_iterable([glob.glob(os.path.join(apath,'pr.day*'))
                                           for apath in paths_in])))
    
    # fpaths_all = np.concatenate((fpaths_pr, fpaths_tasmin, fpaths_tasmax))  #  
    
    return fpaths_pr

def downscale_stuff(n):
    path_in_cmip5 = os.path.join(esd.cfg.path_cmip6_debiased, 'resampled')
    print("01",path_in_cmip5)
    path_out = esd.cfg.path_cmip6_downscaled
    print("02", path_out)
    fpaths_all = get_cmip5_fpaths(path_in_cmip5)
    print("total files", len(fpaths_all))


    
    fpath_prcp_obs = esd.cfg.fpath_obs_pr
    fpath_pr_obsc = os.path.join(esd.cfg.path_obs_resample,
                                   'pr_CHIRPS05_1981_2019_p05deg_remapbic.nc')
    
    base_start_yr = str(esd.cfg.start_date_baseline.year)
    base_end_yr = str(esd.cfg.end_date_baseline.year)
    train_start_yr = str(esd.cfg.start_date_train_downscale.year)
    train_end_yr = str(esd.cfg.end_date_train_downscale.year)
    ds_start_yr = str(esd.cfg.start_date_downscale.year)
    ds_end_yr = str(esd.cfg.end_date_downscale.year)

    prcp_d = PrcpDownscale(esd.cfg.fpath_obs_pr, fpath_pr_obsc,
                           base_start_yr, base_end_yr, train_start_yr, train_end_yr,
                           ds_start_yr, ds_end_yr)

    print("prcp_d")

    downscalers = {'pr': prcp_d}

    fpath_cmip5 = fpaths_all[n]
    print("Starting", fpath_cmip5)
    
    
    ds = xr.open_dataset(fpath_cmip5, decode_cf=False)
    vname = list(ds.data_vars.keys())[0] # variable name
    ds[vname].attrs.pop('missing_value')
    ds = xr.decode_cf(ds)
    latlon_dict_ds = {'latitude': ds.latitude.round(3), 'longitude': ds.longitude.round(3)}
    ds = ds.assign_coords(latlon_dict_ds)
    da = ds[vname].load()
    print("about to downscale", fpath_cmip5)
    da_ds = downscalers[vname].downscale(da)
    print("finished downscaling")

    # Add metadata and create dataset
    da_ds.attrs = da.attrs
    ds_out = da_ds.to_dataset(name=vname)
    ds_out.attrs = ds.attrs
    print("add metadata")

    NC_COMMENT_ATTR = ("Downscaled version of CMIP5 model realization "
                   "(0.05 or 0.1 deg resolution) for the Omo River Basin, Ethiopia. "
                   "All CMIP5 ensemble members were resampled to a standard "
                   "1.0deg resolution grid and then bias corrected using a variation "
                   "of equidistant quantile matching "
                   "(EDCDFm; Li et al. 2010 10.1029/2009JD012882; "
                   "Pierce et al 2015 10.1175/JHM-D-14-0236.1). A constructed "
                   "analog method similar to that of Pierce et al. (2014) "
                   "10.1175/JHM-D-14-0082.1 was then used to downscale the bias "
                   "ensemble members using the CHIRPS05 or ERA5 gridded product as the "
                   "observation calibration dataset.")


    ds_out.attrs['comment'] = NC_COMMENT_ATTR

    subdir = os.path.split(os.path.split(fpath_cmip5)[0])[-1]
    fname = os.path.basename(fpath_cmip5)
    path_out = os.path.join(esd.cfg.path_cmip6_downscaled, subdir)
    mkdir_p(path_out)
            
    fname_splt = fname.split('.')
    fname_splt.insert(5, 'omo_downscaled_0p05deg')
    fname = ".".join(fname_splt)
    fpath_out = os.path.join(path_out, fname)
    print("fpath_out")
    ds_out.to_netcdf(fpath_out)

if __name__ == '__main__':
    i = sys.argv[1]
    n = int(i) - 1
    downscale_stuff(n)
   