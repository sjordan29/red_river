'''
Script to find different metrics for mid-century and late-century compared with the baseline time period. 
'''
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
sys.path.append('/home/smj5vup/.local/bin')

import esd
import glob
import numpy as np
import os
import pandas as pd
import statsmodels.formula.api as smf
import xarray as xr
from esd.util import mkdir_p



class Unbuffered:
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)



def fpaths_runnames(path_cmip5, fname_pattern, min_year, max_year):
    # define rcp paths 
    paths_rcps = sorted(glob.glob(os.path.join(path_cmip5, fname_pattern.split('.')[0])))
     
    fpaths = []
     
    for path_rcp in paths_rcps:
        fpaths.extend(sorted(glob.glob(os.path.join(path_rcp, fname_pattern))))
        
    fpaths = np.array(fpaths)
     
    def parse_run_name(fpath):
         
        try:
            vals = os.path.basename(fpath).split('.') # had to remove "." from projection names - automate this piece in the future 
             
            model_name = vals[2]
            rcp_name = vals[3]
            realization_name = vals[4]

            return "_".join([model_name,rcp_name,realization_name])
        except Exception as e:
            return None
        
    def parse_start_end_yr(fpath):
        
        return [int(i[0:4])
                for i in os.path.basename(fpath).split('.')[-3].split('-')]
         
    run_names = pd.Series([parse_run_name(fpath) for fpath in fpaths])
     
    yrs = pd.DataFrame([parse_start_end_yr(fpath) for fpath in fpaths],
                       columns=['start_yr','end_yr'])
    
    mask_incmplt = (~((yrs.start_yr <= min_year) & (yrs.end_yr >= max_year)))
    mask_rm = ((run_names.isnull()) | (mask_incmplt)).values
    
    print("Removing the following fpaths: "+str(fpaths[mask_rm]))
     
    fpaths = fpaths[~mask_rm]
    run_names = run_names[~mask_rm]

    return fpaths, run_names
    
def calc_trend(anoms, vname):
    
    anoms = pd.DataFrame(anoms)
    anoms['i'] = np.arange(len(anoms))
    
    lm = smf.ols(formula='%s ~ i'%vname, data=anoms).fit()
    
    return lm.params.i



def output_ann_trends(path_cmip5, fname_pattern, vname, unit_convert_func,
                      resample_func, anom_func, start_base,
                      end_base, start_trend, end_trend, path_out, observ):
    
    

    fpaths, run_names = fpaths_runnames(path_cmip5, fname_pattern,
                                        pd.Timestamp(start_base).year,
                                        pd.Timestamp(end_trend).year)
    # anoms_all = []
    anoms_all = pd.Series(np.nan,index=run_names)

    for run_name,fpath in zip(run_names, fpaths):
        #print run_name


        
        ds = xr.open_dataset(fpath)
        mask_lat = np.logical_and(ds.latitude.values >= 4.5,
                                  ds.latitude.values <= 9.5)
        mask_lon = np.logical_and(ds.longitude.values >= 35.5,
                                  ds.longitude.values <= 38.5) 
        da = ds[vname][:, mask_lat, mask_lon].load()
        ds.close()

        if resample_func == 'mean':
            ann = da.resample(time='AS').mean()
        elif resample_func == 'sum':
            ann = da.resample(time='AS').sum()
        
        # norms from obs
        if vname == "pr":
            vname_obs = "precip"
        elif vname == "tasmin":
            vname_obs = "t2m"
        elif vname == "tasmax":
            vname_obs = "t2m"

        ds_obs = xr.open_dataset(observ)
        mask_lat = np.logical_and(ds_obs.latitude.values >= 4.5,
                                  ds_obs.latitude.values <= 9.5)
        mask_lon = np.logical_and(ds_obs.longitude.values >= 35.5,
                                  ds_obs.longitude.values <= 38.5) 
        da_obs = ds_obs[vname_obs][:, mask_lat, mask_lon].load()
        ds_obs.close()

        if resample_func == 'mean':
            ann_obs = da_obs.resample(time='AS').mean()
        elif resample_func == 'sum':
            ann_obs = da_obs.resample(time='AS').sum()

        full_tp = ann.loc[start_trend:end_trend].mean(dim=('longitude', 'latitude','time')) # average over future time period
        norms = ann_obs.loc[start_base:end_base].mean(dim=('longitude', 'latitude','time')) # average over past
        anoms = anom_func(full_tp, norms) # anomalies 

        # anoms.name = vname
        # anoms = anoms.mean(['latitude', 'longitude']).to_dataframe()

        # anoms_all.append(pd.DataFrame({run_name:anoms[vname]}))
        anoms_all.loc[run_name] = anoms


    anoms_all.index.name = 'model_run'
    fname_out = "%s_%s_%s_ann.csv"%(vname,start_trend,end_trend)
    anoms_all.to_csv(os.path.join(path_out,fname_out))


def output_season_trends(path_cmip5, fname_pattern, season_mths, vname,
                         unit_convert_func, resample_func, anom_func,
                         start_base, end_base, start_trend, end_trend, path_out, observ,
                         rolling=False):
    
    fpaths, run_names = fpaths_runnames(path_cmip5, fname_pattern, pd.Timestamp(start_base).year, pd.Timestamp(end_trend).year)
    
    anoms_all = pd.Series(np.nan,index=run_names)
    
    for run_name,fpath in zip(run_names, fpaths):
        print(run_name)
        #print run_name
        
        ds = xr.open_dataset(fpath)
        mask_lat = np.logical_and(ds.latitude.values >= 4.5,
                                  ds.latitude.values <= 9.5)
        mask_lon = np.logical_and(ds.longitude.values >= 35.5,
                                  ds.longitude.values <= 38.5) 
        da = ds[vname][:, mask_lat, mask_lon].load()
        ds.close()

        if resample_func == 'mean':
            mthly = da.resample(time = '1MS').mean() 
        elif resample_func == 'sum':
            mthly = da.resample(time = '1MS').sum() 

        mask_season = np.in1d(mthly['time.month'].values, season_mths)
        
        if rolling:
        
            mthly = mthly.where(xr.DataArray(mask_season, coords=[mthly.time]))
            roller = mthly.rolling(min_periods=len(season_mths), center=True, time=len(season_mths))
            mthly = getattr(roller, resample_func)()
        
        else:
            
            mthly = mthly[mask_season]
        
        if resample_func == 'mean':
            ann_season = mthly.resample(time='AS').mean() 
        elif resample_func == "sum":
            ann_season = mthly.resample(time='AS').sum() 


        # norms from obs
        # norms from obs
        if vname == "pr":
            vname_obs = "precip"
        elif vname == "tasmin":
            vname_obs = "t2m"
        elif vname == "tasmax":
            vname_obs = "t2m"

        ds_obs = xr.open_dataset(observ)
        mask_lat = np.logical_and(ds_obs.latitude.values >= 4.5,
                                  ds_obs.latitude.values <= 9.5)
        mask_lon = np.logical_and(ds_obs.longitude.values >= 35.5,
                                  ds_obs.longitude.values <= 38.5) 
        da_obs = ds_obs[vname_obs][:, mask_lat, mask_lon].load()
        ds_obs.close()

        if resample_func == 'mean':
            mthly_obs = da_obs.resample(time = '1MS').mean() 
        elif resample_func == 'sum':
            mthly_obs = da_obs.resample(time = '1MS').sum() 

        mask_season_obs = np.in1d(mthly_obs['time.month'].values, season_mths)
        
        if rolling:
        
            mthly_obs = mthly_obs.where(xr.DataArray(mask_season, coords=[mthly_obs.time]))
            roller = mthly_obs.rolling(min_periods=len(season_mths), center=True, time=len(season_mths))
            mthly_obs = getattr(roller, resample_func)()
        
        else:
            
            mthly_obs = mthly_obs[mask_season_obs]


        if resample_func == 'mean':
            ann_obs = mthly_obs.resample(time='AS').mean()
        elif resample_func == 'sum':
            ann_obs = mthly_obs.resample(time='AS').sum()


        season_tp = ann_season.loc[start_trend:end_trend].mean(dim=('longitude', 'latitude','time'), skipna=True) # average over future time period    
        norms = ann_obs.loc[start_base:end_base].mean(dim=('longitude', 'latitude','time'), skipna = True)
        anoms = anom_func(season_tp, norms)
        # anoms.name = vname
        
        # anoms = anoms.mean(['lat','lon']).to_dataframe()

        # anoms_all.append(pd.DataFrame({run_name:anoms[vname]}))
        anoms_all.loc[run_name] = anoms

    anoms_all.index.name = 'model_run'
    fname_out = "%s_%s_%s_mths%s.csv"%(vname,start_trend,end_trend,
                                       "-".join([str(x) for x in season_mths]))
    anoms_all.to_csv(os.path.join(path_out,fname_out))
    
    

if __name__ == '__main__':
    
    # Trends for CMIP5 realizations that have been debiased and downscaled 
    path_cmip5_orig = esd.cfg.path_cmip6_downscaled

    path_out_orig = os.path.join(esd.cfg.path_cmip6_trends, 'MVN')
    mkdir_p(path_out_orig)
    
    paths_all = zip([path_cmip5_orig],
                    [path_out_orig])
    
    start_base = esd.cfg.start_date_baseline.strftime('%Y-%m-%d')
    end_base = esd.cfg.end_date_baseline.strftime('%Y-%m-%d')
    trend_periods = [('2040', '2069'), ('2070', '2099')]
    
    # Mean annual precipitation
    for path_cmip5,path_out in paths_all:
        
        print("#######################")
        print("Processing CMIP5 realizations for: " + path_cmip5)
        print("#######################")
        
        for a_start, a_end in trend_periods:
            
            print("Processing Prcp ann over %s to %s"%(a_start, a_end))
            
            output_ann_trends(path_cmip5,
                              fname_pattern='pr.day*.nc',
                              vname='pr',
                              unit_convert_func=None,
                              resample_func='sum',
                              anom_func= lambda mthly_grp,norms: mthly_grp / norms, 
                              start_base=start_base,
                              end_base=end_base,
                              start_trend=a_start,
                              end_trend=a_end,
                              path_out=path_out,
                              observ=esd.cfg.fpath_obs_pr)
        
        
        # Prcp seasonal
        for a_start, a_end in trend_periods:
            
            for season,use_roll in [([2,3,4,5], False), # belg: feb to may (short rains season)
                                    ([6,7,8,9], False), # kiremt: june and sept (logn rains season)
                                    ([10,11,12,1], False), # bega: dry weather 
                                    ([8,9], False),]: # recession agriculture
                
                print("Processing Prcp seasonal %s over %s to %s"%(str(season), a_start, a_end))
                
                output_season_trends(path_cmip5,
                                     fname_pattern='pr.day*.nc',
                                     season_mths=season,
                                     vname='pr',
                                     unit_convert_func=None,
                                     resample_func='sum',
                                     anom_func=lambda mthly_grp,norms: mthly_grp / norms,
                                     start_base=start_base,
                                     end_base=end_base,
                                     start_trend=a_start,
                                     end_trend=a_end,
                                     path_out=path_out,
                                     observ=esd.cfg.fpath_obs_pr,
                                    rolling=use_roll)
                
        for a_start, a_end in trend_periods:
             
            print("Processing Tasmin ann over %s to %s"%(a_start, a_end))
             
            output_ann_trends(path_cmip5,
                              fname_pattern='tasmin.day*.nc',
                              vname='tasmin',
                              unit_convert_func=None,
                              resample_func='mean',
                              anom_func= lambda mthly_grp,norms: mthly_grp - norms,
                              start_base=start_base,
                              end_base=end_base,
                              start_trend=a_start,
                              end_trend=a_end,
                              path_out=path_out,
                              observ=esd.cfg.fpath_obs_tasmin)
        
        
        # Tair seasonal
        for a_start, a_end in trend_periods:
            
            for season,use_roll in [([2,3,4,5], False),
                                    ([6,7,8,9], False),
                                    ([10,11,12,1], False),
                                    ([8,9], False),]:
                
                print("Processing tasmin seasonal %s over %s to %s"%(str(season), a_start, a_end))
                
                output_season_trends(path_cmip5,
                                     fname_pattern='tasmin.day*.nc',
                                     season_mths=season,
                                     vname='tasmin',
                                     unit_convert_func=None,
                                     resample_func='mean',
                                     anom_func=lambda mthly_grp,norms: mthly_grp - norms,
                                     start_base=start_base,
                                     end_base=end_base,
                                     start_trend=a_start,
                                     end_trend=a_end,
                                     path_out=path_out, 
                                     observ=esd.cfg.fpath_obs_tasmin,
                                     rolling=use_roll)

        # tasmax
        for a_start, a_end in trend_periods:
             
            print("Processing tasmax ann over %s to %s"%(a_start, a_end))
             
            output_ann_trends(path_cmip5,
                              fname_pattern='tasmax.day*.nc',
                              vname='tasmax',
                              unit_convert_func=None,
                              resample_func='mean',
                              anom_func= lambda mthly_grp,norms: mthly_grp - norms,
                              start_base=start_base,
                              end_base=end_base,
                              start_trend=a_start,
                              end_trend=a_end,
                              path_out=path_out,
                              observ=esd.cfg.fpath_obs_tasmax)
        
        
        # Tair seasonal
        for a_start, a_end in trend_periods:
            
            for season,use_roll in [([2,3,4,5], False),
                                    ([6,7,8,9], False),
                                    ([10,11,12,1], False),
                                    ([8,9], False),]:
                
                print("Processing tasmax seasonal %s over %s to %s"%(str(season), a_start, a_end))
                
                output_season_trends(path_cmip5,
                                     fname_pattern='tasmax.day*.nc',
                                     season_mths=season,
                                     vname='tasmax',
                                     unit_convert_func=None,
                                     resample_func='mean',
                                     anom_func=lambda mthly_grp,norms: mthly_grp - norms,
                                     start_base=start_base,
                                     end_base=end_base,
                                     start_trend=a_start,
                                     end_trend=a_end,
                                     path_out=path_out, 
                                     observ=esd.cfg.fpath_obs_tasmax,
                                     rolling=use_roll)