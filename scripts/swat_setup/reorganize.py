# -*- coding: utf-8 -*-
"""
Sarah Jordan
02/23/2021

Reorganize data from zonal statistics for pr, tasmin, tasmax
"""

import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')

from esd.util import mkdir_p
import pandas as pd
import glob
import numpy as np
import os


def rearrange(vb):
    os.chdir('/scratch/smj5vup/CMIP6/Debias/cmip6_swat/')

    # define output path 
    output_path = '/scratch/smj5vup/CMIP6/Debias/cmip6_swat_cleaned/'
    mkdir_p(output_path)
    output_path2 = '/scratch/smj5vup/CMIP6/Debias/cmip6_swat_cleaned/' + vb
    mkdir_p(output_path2)

    # all files
    files = glob.glob(vb + '/*.csv')


    for file in files:
        print(file)
        # read data
        df = pd.read_csv(file)
        fname = str(file[:-4])
        print(fname)
        names = df.ID
        sub_names = ['subbasin-' + str(f) for f in names]

        # transpose
        df_T = df.transpose()
        df_T = df_T.iloc[2:]

        # change column names and fix index 
        df_T.columns = sub_names
        df_T.index = [f[1:] for f in df_T.index] # remove leading X
        df_T.index = pd.to_datetime(df_T.index) # to datetime

        # save
        df_T.to_csv(output_path + fname + '.csv')




rearrange('pr')
rearrange('tasmin')
rearrange('tasmax')