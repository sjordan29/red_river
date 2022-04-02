from misc import mkdir_p
import configparser # modified from ConfigParser import ConfigParser
from copy import copy
import numpy as np
import os
import pandas as pd

class OmoConfig():
    
    def __init__(self, fpath_ini):
    
        cfg = configparser.ConfigParser()
        cfg.read(fpath_ini)
        
        self.data_root = cfg.get('OMO_CONFIG', 'data_root')
        self.code_root = cfg.get('OMO_CONFIG', 'code_root')
        
        bbox_str = cfg.get('OMO_CONFIG', 'domain_bounds_omo')
        self.bbox = tuple([np.float(i) for i in bbox_str.split(',')])
        self.path_data_bbox = os.path.join(self.data_root, 'domain_bbox')
        mkdir_p(self.path_data_bbox)
        
        self.fpath_obs_pr = cfg.get('OMO_CONFIG', 'FPATH_OBS_PR')
        self.fpath_obs_tasmin = cfg.get('OMO_CONFIG', 'FPATH_OBS_TASMIN')
        self.fpath_obs_tasmax = cfg.get('OMO_CONFIG', 'FPATH_OBS_TASMAX')
        
        # self.path_cmip6_archive = cfg.get('OMO_CONFIG', 'PATH_CMIP6_ARCHIVE')
        self.path_cmip6_archive = os.path.join(self.data_root, 'cmip6_archive')
        mkdir_p(self.path_cmip6_archive)
        
        self.path_cmip6_resample = os.path.join(self.data_root, 'cmip6_resample')
        mkdir_p(self.path_cmip6_resample)
        
        self.path_cmip6_resize = os.path.join(self.data_root, 'cmip6_resize')
        mkdir_p(self.path_cmip6_resize)
        
        self.path_cmip6_trends = os.path.join(self.data_root, 'cmip6_trends')
        mkdir_p(self.path_cmip6_trends)
        
        self.path_cmip6_wetbias = os.path.join(self.data_root, 'cmip6_wetbias')
        mkdir_p(self.path_cmip6_wetbias)
        
        self.path_cmip6_cleaned = os.path.join(self.data_root, 'cmip6_cleaned')
        mkdir_p(self.path_cmip6_cleaned)
        
        self.path_cmip6_debiased = os.path.join(self.data_root, 'cmip6_debiased')
        mkdir_p(self.path_cmip6_debiased)
        
        self.path_cmip6_downscaled = os.path.join(self.data_root, 'cmip6_downscaled')
        mkdir_p(self.path_cmip6_downscaled)
        
        self.path_obs_resample = os.path.join(self.data_root, 'obs_resample')
        mkdir_p(self.path_obs_resample)
        
        self.path_obs_resize = os.path.join(self.data_root, 'obs_resize')
        mkdir_p(self.path_obs_resize)
        
        self.start_date_baseline = pd.Timestamp(cfg.get('OMO_CONFIG', 'START_DATE_BASELINE'))
        self.end_date_baseline = pd.Timestamp(cfg.get('OMO_CONFIG', 'END_DATE_BASELINE'))
        
        self.start_date_train_downscale = pd.Timestamp(cfg.get('OMO_CONFIG', 'START_DATE_TRAIN_DOWNSCALE'))
        self.end_date_train_downscale = pd.Timestamp(cfg.get('OMO_CONFIG', 'END_DATE_TRAIN_DOWNSCALE'))
        
        self.start_date_downscale = pd.Timestamp(cfg.get('OMO_CONFIG', 'START_DATE_DOWNSCALE'))
        self.end_date_downscale = pd.Timestamp(cfg.get('OMO_CONFIG', 'END_DATE_DOWNSCALE'))
        
        self.path_logs = os.path.join(self.data_root, 'logs')
        mkdir_p(self.path_logs)
        
        
    def to_str_dict(self):
        
        a_dict = copy(self.__dict__)
        for a_key in a_dict.keys():
            a_dict[a_key] = str(a_dict[a_key])
        return a_dict
            
        