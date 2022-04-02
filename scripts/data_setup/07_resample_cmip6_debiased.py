'''
Script to smooth debiased CMIP6 model data from 1deg resolution to that of the 0.05 degree CHIPRS05 degree data (pr)
or the 0.1 degree ERA5 ddegree data (tasmin, tasmax). These smoothed versions of the debiased data are used
in the analog downscaling procedure.
'''
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
sys.path.append('/home/smj5vup/.local/bin')

from esd.util import mkdir_p
import esd
import glob
import os
import subprocess

if __name__ == '__main__':
    
    # Upsample 1deg debiased CMIP6 datasets to Red River 0.05 deg CHIRPS grid with
    # bicubic or bilinear interpolation (remapbil or remapbic) or 0.25 deg ERA5 temperature grid
    # These are used in the analog downscaling procedure
    
    path_in = os.path.join(esd.cfg.path_cmip6_debiased)
    resample_method = 'remapbic' #remapbil or remapbic
    path_out = os.path.join(esd.cfg.path_cmip6_debiased, 'resampled')
    mkdir_p(path_out)
    
    fpath_cdo_grid_def_pr = os.path.join(esd.cfg.path_data_bbox, 'cdo_grid_005deg_omo') # different resolution for precip and temp
    fpath_cdo_grid_def_tas = os.path.join(esd.cfg.path_data_bbox, 'cdo_grid_01deg_omo')

    paths_rcp = sorted(glob.glob(os.path.join(path_in, "*")))
    
    for path_rcp in paths_rcp:
        
        fpaths_tasmax = sorted(glob.glob(os.path.join(path_rcp, 'tasmax.day*')))
        fpaths_tasmin = sorted(glob.glob(os.path.join(path_rcp, 'tasmin.day*')))
        fpaths_pr = sorted(glob.glob(os.path.join(path_rcp, 'pr.day*')))
        fpaths_tas  = fpaths_tasmin + fpaths_tasmax
        
        path_out_rcp = os.path.join(path_out, os.path.basename(path_rcp))
        mkdir_p(path_out_rcp)
        
        for fpath_in in fpaths_tas:
            
            fpath_out = os.path.join(path_out_rcp, os.path.basename(fpath_in))
            cmd = "cdo %s,%s %s %s"%(resample_method, fpath_cdo_grid_def_tas, fpath_in, fpath_out)
            print(cmd)
            subprocess.call(cmd, shell=True)

        for fpath_in in fpaths_pr:
            fpath_out = os.path.join(path_out_rcp, os.path.basename(fpath_in))
            cmd = "cdo %s,%s %s %s"%(resample_method, fpath_cdo_grid_def_pr, fpath_in, fpath_out)
            print(cmd)
            subprocess.call(cmd, shell=True)
            
        
        