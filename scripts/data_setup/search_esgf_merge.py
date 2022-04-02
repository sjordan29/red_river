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
sys.path.append('/home/smj5vup/.local/bin')
import esd 
from esd.util import clean_gcm_da, unit_convert_pr, unit_convert_tas, mkdir_p, StatusCheck


# from pyesgf.logon import LogonManager # if getting keyerror need to change 'HOME' to 'USERPROFILE' if using windows in \lib\site-packages\pyesgf\logon.py; line 58, 62
# lm = LogonManager()
# lm.logoff()
# lm.logon_with_openid('https://esgf-data.dkrz.de/esgf-idp/openid/XXX', password='XXXX')# , bootstrap=True)
# # lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
# print(lm.is_logged_on())

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
        for i,fpath_row in fpaths_df.iterrows():


            print(fpath_row.fname)
            ds = xr.open_dataset(fpath_row.fpath)
            
            unit_convert_funcs = {'pr':unit_convert_pr,
                          'tasmin':unit_convert_tas,
                          'tasmax':unit_convert_tas}
            
            da = clean_gcm_da(ds[fpath_row.elem], unit_convert_funcs[fpath_row.elem])
            ds_out = da.to_dataset(name=fpath_row.elem)
            # ds_out.attrs.update(ds.attrs)
            
        
            # output file name and path name
            fname_out = "%s.day.%s.%s.%s.%s-%s.nc"%(fpath_row.elem,
                                                    fpath_row.model_name,
                                                    fpath_row.rcp_name,
                                                    fpath_row.realization_name,
                                                    fpath_row.start_yr.strftime('%Y%m%d'),
                                                    fpath_row.end_yr.strftime('%Y%m%d')) 

            mkdir_p(os.path.join(esd.cfg.path_cmip6_archive,vb))
            fpath_out = os.path.join(esd.cfg.path_cmip6_archive,
                                        vb,
                                        fname_out)
            # save to netCDF
            ds_out.load().to_netcdf(fpath_out)


def combine_xarray(hist_path, rcp_path, rcp):
    # combine and open 
    hist_path.extend(rcp_path)
    full = xr.open_mfdataset(hist_path, combine='by_coords')

    # define start and end dates 
    # start_date = str(min([int(f.split('.')[-2].split('-')[0]) for f in hist_path]))
    # end_date = str(max([int(f.split('.')[-2].split('-')[1]) for f in rcp_path]))

    start_date = "19810101"
    end_date = "20991231"

    # determine name 
    b = os.path.basename(rcp_path[-1]).split('.')[0:-2] 
    name_str = ""
    for item in b:
        name_str = name_str+item+'.'
    name_str = name_str + start_date + '-' + end_date
    fname_out = name_str + '.historical_merged.nc'
    full.lat.attrs['units'] = 'degrees_north'
    full.lon.attrs['units'] = 'degrees_east'

    full = full.sel(time=slice("1981-01-01", "2099-12-31"))
    

    var = fname_out.split('.')[0]

    # make directories -- cleaned
    mkdir_p(os.path.join(esd.cfg.path_cmip6_cleaned,
                                 'historical+%s_merged'%rcp))
    mkdir_p(os.path.join(esd.cfg.path_cmip6_cleaned,
                                 'historical+%s_merged'%rcp,
                                 var))
    # file name
    fpath_out = os.path.join(esd.cfg.path_cmip6_cleaned,
                                 'historical+%s_merged'%rcp,
                                 var,
                                 fname_out)
    # save to netCDF
    full.to_netcdf(fpath_out)
    
    
        
        
def full_merge(projection_name, ensemble_name):
    print("Starting", projection_name)
    for r in ['historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85']: # ssp126','ssp126', 'ssp245', 
        for v in ['pr', 'tasmin', 'tasmax']:
        
            try:
                ESGF_to_netCDF(projection_name, r, v, ensemble_name)
                print(projection_name, r, v, "data downloaded!")                
            except (AttributeError, IndexError) as e:
                print("No projections found for " + projection_name + ' ' + r + ' ' + ensemble_name + ' ' + v)

    print("Downloads are done!")

    for r in ['rcp26', 'rcp45', 'rcp60', 'rcp85']: # ssp126','ssp126', 'ssp245', 
        for v in ['pr', 'tasmin', 'tasmax']:
            projection_name_replaced = projection_name.replace(".", "-")
            projection_name_replaced = projection_name_replaced.replace("(", "-")
            projection_name_replaced = projection_name_replaced.replace(")", "")
            files = os.listdir(os.path.join(esd.cfg.path_cmip6_archive, v))
            proj_files = [os.path.join(esd.cfg.path_cmip6_archive, v, f) for f in files if projection_name_replaced in f]
            
            h = [f for f in proj_files if 'historical' in f]
            f_ls = [f for f in proj_files if r in f]

            if len(h) == 0:
                print("No historical data for", projection_name, v)
            elif len(f_ls) == 0:
                print("No projections for", projection_name, r, v)
            else:
                combine_xarray(h, f_ls, r)
                print("Finished", projection_name, r, v)



    
 

# proj_ls = ["CSIRO-Mk3.6.0", "CSIRO-Mk3L-1-2", "EC-EARTH", "GEOS-5", "FGOALS-g2", 
#            "FGOALS-gl", "FGOALS-s2", "GFDL-CM2.1", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H",
#            "GISS-E2-H-CC", "GISS-E2-R", "GISS-E2-R-CC", "HadCM3","HadGEM2-A", "HadGEM2-AO",
#            "HadGEM2-CC", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",
#            "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC4h", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MPI-ESM-P",
#            "MRI-AGCM3-2H", "MRI-AGCM3-2S", "MRI-CGCM3", "MRI-ESM1", "NorESM1-M", "NorESM1-ME"] # , "BNU-ESM",

proj_ls = ["ACCESS1.0", "BCC-CSM1.1(m)", "CCSM4", "CESM1(BGC)", 
           "CESM1(FASTCHEM)", "CESM1(WACCM)", "CFSv2-2011", "CMCC-CESM","CMCC-CM", "CMCC-CMS", "CNRM-CM5", "CNRM-CM5-2",
           "CSIRO-Mk3.6.0", "CSIRO-Mk3L-1-2", "EC-EARTH", "GEOS-5",
           "FGOALS-gl", "FGOALS-s2", "GFDL-CM2.1", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H",
           "GISS-E2-H-CC", "GISS-E2-R", "GISS-E2-R-CC", "HadCM3","HadGEM2-A", "HadGEM2-AO",
           "HadGEM2-CC", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",
           "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC4h", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MPI-ESM-P",
           "MRI-AGCM3-2H", "MRI-AGCM3-2S", "MRI-CGCM3", "MRI-ESM1", "NorESM1-M", "NorESM1-ME"]



# now merge data
i = sys.argv[1]
print(i)
n = int(i) - 1
print(n)
projection_name = proj_ls[n]

# full_merge(projection_name, "r1i1p1")
# ESGF_to_netCDF(projection_name, "historical", "pr","r1i1p1")
for r in ['rcp26', 'rcp45', 'rcp60', 'rcp85']: # ssp126','ssp126', 'ssp245', 
    for v in ['pr', 'tasmin', 'tasmax']:
        print("Starting", projection_name, r, v)

        # all files with given projection name 
        files = os.listdir(os.path.join(esd.cfg.path_cmip6_archive, v))

        projection_name_replaced = projection_name.replace(".", "-")
        projection_name_replaced = projection_name_replaced.replace("(", "-")
        projection_name_replaced = projection_name_replaced.replace(")", "")
        proj_files = [os.path.join(esd.cfg.path_cmip6_archive, v, f) for f in files if projection_name_replaced in f]

        # historical 
        h = [f for f in proj_files if 'historical' in f]
        h.sort()

        # projections
        f_ls = [f for f in proj_files if r in f]
        f_ls.sort()

       	# combine
        if len(h) == 0:
            print("No historical data for", projection_name, v)
        elif len(f_ls) == 0:
            print("No projections for", projection_name, r, v)
        else:
            combine_xarray(h, f_ls, r)

