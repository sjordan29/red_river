# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 18:52:09 2021

@author: Sarah

Search and load CMIP6 data to xarray
Convert xarray to netCDF
Combine historical and ssp into single netCDF
Adapted from Jared Oyler 
"""

### PACKAGES ####################
from pyesgf.search import SearchConnection
import os
import pandas as pd
import requests
import xarray as xr
import scipy
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
import esd
from esd.util import clean_gcm_da, unit_convert_pr, unit_convert_tas, mkdir_p, StatusCheck
from pyesgf.logon import LogonManager # if getting keyerror need to change 'HOME' to 'USERPROFILE' if using windows in \lib\site-packages\pyesgf\logon.py; line 58, 62
import os.path
home_folder = os.path.expanduser('~')
lm = LogonManager()
lm.logoff()
lm.logon_with_openid('https://esgf-data.dkrz.de/esgf-idp/openid/smj5vup', password='76Crooked!', bootstrap=True)
# lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
print(lm.is_logged_on())


## search all files ###########################################################
data_node = 'https://esgf-data.dkrz.de/esg-search'
conn = SearchConnection(data_node, distrib=True)


def parse_fpath_info(fpath):
           try:
           
               bname = os.path.basename(fpath)
               
               if not bname.endswith('.nc'):
                   bname = bname + ".nc"
                
               # s
               vals = bname.split('_')
               
               elem = vals[0]   # pr, tasmin, or tasmax  
               model_name = vals[2] 
               rcp_name = vals[3] # historical, ssp126, ssp245, ssp370, or ssp585-- had to be changed for cmip6 format
               realization_name = vals[4] # ensemble 
               
               start_yr = pd.to_datetime(vals[-1][0:8]) # had to change these for cmip6 format
               end_yr = pd.to_datetime(vals [-1][9:-3])
           
           except Exception:
               return tuple([fpath,bname]+[None]*5)
           
           return fpath, bname, elem, model_name, rcp_name, realization_name, start_yr, end_yr



def ESGF_to_netCDF(proj, exp_id, vb, var_lab):
    ''' 
    Function to search ESGF for climate projection, subset, and convert to netCDF
    proj = projection name
    exp_id = experiment_id (historical, ssp126, ssp245, ssp370, ssp585)
    vb = variable (pr, tasmin, tasmax)
    var_lab = ensemble 
    '''

    # ctx = conn.new_context(
    #     project='CMIP6',
    #     source_id=proj,
    #     experiment_id=exp_id, 
    #     variable=vb,
    #     frequency='day',
    #     variant_label=var_lab,
    #     data_node='esgf-data3.ceda.ac.uk')
    
    
    ctx = conn.new_context(
        project='CMIP5',
        experiment=exp_id,
        model=proj,
        ensemble=var_lab,
        realm = 'atmos',
        time_frequency='day',
        data_node='esgf-data1.ceda.ac.uk'
        )
    
    
    # get all files in search 
    result = ctx.search()[0]
    files = result.file_context().search()
    

    
    # if len([f for f in files]) > 0:
        # combine all datasets by coordinates
    start_end = (esd.cfg.start_date_baseline, pd.Timestamp('2099-12-31'))

    
    # put all GCMs in a dataframe 
    fpaths_df = pd.DataFrame([parse_fpath_info(f.opendap_url) for f in files],
                         columns=['fpath','fname','elem','model_name','rcp_name',
                                  'realization_name','start_yr','end_yr'])
    
    fpaths_df = fpaths_df[fpaths_df.elem == vb]
    
    if len(fpaths_df.end_yr) != len(fpaths_df.end_yr.unique()):
        print("There are probably duplicates! Projection:", proj, vb)
    

    
    # Remove GCM files with no fpath info or incomplete records of the specified
    # period of interest
    mask_keep = ((fpaths_df.start_yr <= start_end[-1]) &
             (fpaths_df.end_yr >= start_end[0]))
    fpaths_df = fpaths_df[mask_keep]
    
    if len(fpaths_df.end_yr) > 0:
    
        ds_agg = xr.open_mfdataset([f for f in fpaths_df.fpath], chunks={'time': 120}, combine='by_coords')
        
        # subset over omo and timeframe
        da = ds_agg[vb].sel(lat=slice(2,12), lon=slice(32, 42))
        
        # output file name and path name
        fname_out = "%s.day.%s.%s.%s.%s-%s.nc"%(vb,
                                                           proj,
                                                           exp_id,
                                                           var_lab,
                                                           min(fpaths_df.start_yr).strftime('%Y%m%d'),
                                                           min([max(fpaths_df.end_yr), start_end[-1]]).strftime('%Y%m%d'))
       
        mkdir_p(os.path.join(esd.cfg.path_cmip6_archive,vb))
        fpath_out = os.path.join(esd.cfg.path_cmip6_archive,
                                    vb,
                                    fname_out)
        # save to netCDF
        da.load().to_netcdf(fpath_out)
        
        return fpath_out


def combine_xarray(hist_path, rcp_path):
    # combine by 
    full = xr.open_mfdataset([hist_path, rcp_path], combine='by_coords')
    b = os.path.basename(rcp_path).split('.')[0:-2] 
    name_str = ""
    for item in b:
        name_str = name_str+item+'.'
    name_str = name_str + hist_path.split('.')[-2].split('-')[0] + '-' + os.path.basename(rcp_path).split('.')[-2].split('-')[1]
    fname_out = name_str + '.historical_merged.nc'
    
    # convert units 
    unit_convert_funcs = {'pr':unit_convert_pr,
                          'tasmin':unit_convert_tas,
                          'tasmax':unit_convert_tas}
    var = fname_out.split('.')[0]
    # func = unit_convert_funcs[var]
    # full = func(full, var) 
    
    # clean dates
    da = clean_gcm_da(full[var], unit_convert_funcs[var])
    ds_out = da.to_dataset(name=var)
    ds_out.attrs.update(full.attrs)
    
    # make directories -- cleaned
    mkdir_p(os.path.join(esd.cfg.path_cmip6_cleaned,
                                 'historical+%s_merged'%rcp_path.split('.')[4]))
    mkdir_p(os.path.join(esd.cfg.path_cmip6_cleaned,
                                 'historical+%s_merged'%rcp_path.split('.')[4],
                                 var))
    # file name
    fpath_out = os.path.join(esd.cfg.path_cmip6_cleaned,
                                 'historical+%s_merged'%rcp_path.split('.')[4],
                                 var,
                                 fname_out)
    # save to netCDF
    ds_out.to_netcdf(fpath_out)
    
    
        
        
def full_merge(projection_name, ensemble_name):
    # precip
    print("Starting", projection_name)
    try: 
        hist_pr_path = ESGF_to_netCDF(projection_name, 'historical', 'pr', ensemble_name)
        print("Historical pr data done.")
        hist_precip_exists = True
    except IndexError:
        print("No historical pr data.")
        hist_precip_exists = False

            
    # tasmin
    try:
        hist_tasmin_path = ESGF_to_netCDF(projection_name, 'historical', 'tasmin', ensemble_name)
        print("Historical tasmin data done.")
        hist_tasmin_exists = True
    except IndexError:
        print("No historical tasmin data.")
        hist_tasmin_exists = False
    
    # tasmax 
    try:
        hist_tasmax_path = ESGF_to_netCDF(projection_name, 'historical', 'tasmax', ensemble_name)
        print("Historical Data Done!")
        hist_tasmax_exists = True 
    except IndexError:
        print("No historical tasmax data.")
        hist_tasmax_exists = False 
        
    # dictionary to define historical comparison 
    xr_hist_select = {}
    if hist_precip_exists == True:
        xr_hist_select['pr'] = hist_pr_path
    if hist_tasmin_exists == True:
        xr_hist_select['tasmin'] = hist_tasmin_path
    if hist_tasmax_exists == True:
        xr_hist_select['tasmax'] = hist_tasmax_path
     
    
    for r in ['rcp26', 'rcp45', 'rcp60', 'rcp85']: # ssp126','ssp126', 'ssp245', 
        for v in ['pr', 'tasmin', 'tasmax']:
        
            try:
                f = ESGF_to_netCDF(projection_name, r, v, ensemble_name)
                try:
                    
                    h = xr_hist_select[v]
                    combine_xarray(h,f)
                    print(projection_name, r, v, "done!")
                except:
                    print("Can't combine - no historical data for", projection_name, v)
                
                
            except (AttributeError, IndexError) as e:
                print("No projections found for " + projection_name + ' ' + r + ' ' + ensemble_name + ' ' + v)

    
 


proj_ls = ["BCC-CSM1.1(m)", "BNU-ESM", "CCSM4", "CESM1(BGC)", "CESM1(CAM5)",
           "CESM1(FASTCHEM)", "CESM1(WACCM)", "CFSv2-2011", "CMCC-CESM", "CMCC-CM", "CMCC-CMS", "CNRM-CM5",
           "CNRM-CM5-2", "CSIRO-Mk3.6.0", "CSIRO-Mk3L-1-2", "EC-EARTH", "GEOS-5", "FGOALS-g2", 
           "FGOALS-gl", "FGOALS-s2", "GFDL-CM2.1", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H",
           "GISS-E2-H-CC", "GISS-E2-R", "GISS-E2-R-CC", "HadCM3", "HadGEM2-A", "HadGEM2-AO",
           "HadGEM2-CC", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",
           "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC4h", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MPI-ESM-P",
           "MRI-AGCM3-2H", "MRI-AGCM3-2S", "MRI-CGCM3", "MRI-ESM1", "NorESM1-M", "NorESM1-ME"] #"ACCESS1.0", "ACCESS1.0", "ACCESS1.3", 


# import sys 

i = sys.argv[1]
print(i)
n = int(i) - 1
full_merge(proj_ls[n], 'r1i1p1')            
    

