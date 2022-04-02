# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:18:38 2021

@author: Sarah
"""
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')

from esd.downscale.prcp import PrcpDownscale
from esd.downscale.tair import TairDownscale
from esd.util import Unbuffered, StatusCheck, mkdir_p
import esd
import glob
import itertools
import numpy as np
import os
import sys
import xarray as xr


def proc_work(rank):
    # load in years of interest 
    base_start_yr = str(esd.cfg.start_date_baseline.year)
    base_end_yr = str(esd.cfg.end_date_baseline.year)
    train_start_yr = str(esd.cfg.start_date_train_downscale.year)
    train_end_yr = str(esd.cfg.end_date_train_downscale.year)
    ds_start_yr = str(esd.cfg.start_date_downscale.year)
    ds_end_yr = str(esd.cfg.end_date_downscale.year)
    
    # define downscale windows 
    downscale_wins =  [('1981','2018'), ('2019','2039'), ('2040','2069'), ('2070','2099')]
    
    # data that has been upsampled and then downsampled
    # from original resolution to 1 degree back to original 
    fpath_tasmax_obsc = os.path.join(esd.cfg.path_obs_resize,
                                   'tasmax_ERA5_1981_2019_p10deg_remapbic.nc')
    fpath_tasmin_obsc = os.path.join(esd.cfg.path_obs_resize,
                                     'tasmin_ERA5_1981_2019_p10deg_remapbic.nc')
    fpath_pr_obsc = os.path.join(esd.cfg.path_obs_resize,
                                   'pcp_CHIRPS_1981_2020_p05deg_remapbic.nc')

    # 
    tasmax_d = TairDownscale(esd.cfg.fpath_obs_tasmax, fpath_tasmax_obsc,
                           base_start_yr, base_end_yr,
                           train_start_yr, train_end_yr, downscale_wins)
    tasmin_d = TairDownscale(esd.cfg.fpath_obs_tasmin, fpath_tasmin_obsc,
                           base_start_yr, base_end_yr,
                           train_start_yr, train_end_yr, downscale_wins)
    prcp_d = PrcpDownscale(esd.cfg.fpath_obs_pr, fpath_pr_obsc,
                           base_start_yr, base_end_yr, train_start_yr, train_end_yr,
                           ds_start_yr, ds_end_yr)
    
    downscalers = {'tasmax': tasmax_d, 'tasmin':tasmin_d, 'pr': prcp_d}

    # cmip6 
    fpath_cmip6 = 'D:/Research/Debias/test_debiasedreample_t.nc'
    ds = xr.open_dataset(fpath_cmip6, decode_cf=False)
    vname = list(ds.data_vars.keys())[0] # variable name
    ds[vname].attrs.pop('missing_value')
    ds = xr.decode_cf(ds)
    da = ds[vname].load()
    da_ds = downscalers[vname].downscale(da)
    
    # Add metadata and create dataset
    da_ds.attrs = da.attrs
    ds_out = da_ds.to_dataset(name=vname)
    ds_out.attrs = ds.attrs
    ds_out.attrs['comment'] = NC_COMMENT_ATTR

    subdir = os.path.split(os.path.split(fpath_cmip6)[0])[-1]
    fname = os.path.basename(fpath_cmip6)
                
    MPI.COMM_WORLD.send((subdir, fname, ds_out), dest=RANK_WRITE, tag=TAG_DOWORK)
    MPI.COMM_WORLD.send(rank, dest=RANK_COORD, tag=TAG_DOWORK)

                
def proc_coord(nwrkers):
    
    def get_cmip5_fpaths(a_path):
        
        paths_in = sorted(glob.glob(os.path.join(a_path, '*')))
            
        fpaths_tas = sorted(list(itertools.chain.
                                 from_iterable([glob.glob(os.path.join(apath,'tas.day*'))
                                                for apath in paths_in])))
        fpaths_pr = sorted(list(itertools.chain.
                                from_iterable([glob.glob(os.path.join(apath,'pr.day*'))
                                               for apath in paths_in])))
        
        fpaths_all = np.concatenate((fpaths_pr, fpaths_tas))
        
        return fpaths_all
        
    
    path_in_cmip5 = os.path.join(esd.cfg.path_cmip5_debiased, 'resampled')
    path_out = esd.cfg.path_cmip5_downscaled
    fpaths_all = get_cmip5_fpaths(path_in_cmip5)
    fpaths_done =  get_cmip5_fpaths(path_out)
    fnames_all = np.array([os.path.basename(fpath) for fpath in fpaths_all])
    fnames_done = np.array([os.path.basename(fpath) for fpath in fpaths_done])
    mask_done = np.in1d(fnames_all, fnames_done)
    fnames_all = fnames_all[~mask_done]
    fpaths_all = fpaths_all[~mask_done]
    
    print "COORD: %d datasets already downscaled. %d to go."%(mask_done.sum(), fnames_all.size)
    
    MPI.COMM_WORLD.bcast(fpaths_all.size, root=RANK_COORD)
        
    print "COORD: Done initialization. Starting to send work."
    
    cnt = 0
    nrec = 0
    
    for a_fpath in fpaths_all:
        
        if cnt < nwrkers:
            dest = cnt + N_NON_WRKRS
        else:
            dest = MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
            nrec += 1

        MPI.COMM_WORLD.send(a_fpath, dest=dest, tag=TAG_DOWORK)
        cnt += 1
    
    for w in np.arange(nwrkers):
        MPI.COMM_WORLD.send(None, dest=w + N_NON_WRKRS, tag=TAG_STOPWORK)
        
    print "COORD: done"

def proc_write(nwrkers):

    status = MPI.Status()
    
    bcast_msg = None
    bcast_msg = MPI.COMM_WORLD.bcast(bcast_msg, root=RANK_COORD)
    n_datasets = bcast_msg  
    
    nwrkrs_done = 0
                        
    stat_chk = StatusCheck(n_datasets, 5)
    
    while 1:

        subdir, fname, ds_out = MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE,
                                                    tag=MPI.ANY_TAG,
                                                    status=status)
        
        if status.tag == TAG_STOPWORK:
            
            nwrkrs_done += 1
            if nwrkrs_done == nwrkers:
                print "WRITER: Finished"
                return 0
        else:
            path_out = os.path.join(esd.cfg.path_cmip5_downscaled, subdir)
            mkdir_p(path_out)
            
            fname_splt = fname.split('.')
            fname_splt.insert(5, 'red_river_downscaled_0p25deg')
            fname = ".".join(fname_splt)
            
            fpath_out = os.path.join(path_out, fname)
            print "WRITER: Writing %s..."%fpath_out
            ds_out.to_netcdf(fpath_out)
            
            stat_chk.increment()

if __name__ == '__main__':
            
    rank = MPI.COMM_WORLD.Get_rank()
    nsize = MPI.COMM_WORLD.Get_size()
        
    if rank == RANK_COORD:
        proc_coord(nsize - N_NON_WRKRS)
    elif rank == RANK_WRITE:
        proc_write(nsize - N_NON_WRKRS)
    else:
        proc_work(rank)

    MPI.COMM_WORLD.Barrier()