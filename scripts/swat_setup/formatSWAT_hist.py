# -*- coding: utf-8 -*-
"""
Sarah Jordan
02/24/2021

Organize dataframes from CMIP5 projections in to SWAT format
pcp1.pcp
tmp1.tmp

Historical data from existing pcp and tmp files for 1981-2019
Projection data fter that 

csvs reading in have 27 columns - 1 for each subbasin - and a datetime index
"""
import sys
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/')
sys.path.append('/scratch/smj5vup/CMIP6/omoCMIP/esd/util/')
from esd.util import mkdir_p
import glob
import pandas as pd
import numpy as np
import os


### FUNCTION ################################################################
# convert dates to format for swat files 
def convert_dt(df):
    '''
    Parameters
    ----------
    df : dataframe of interest, index is a datetime index 

    Returns
    -------
    combo : list of dates in the index in YYYYJJJ format -- type = integer

    '''
    yrs = df.index.year
    yr_str = [str(f) for f in yrs]
    doys = df.index.dayofyear
    doy_str = [str(f).zfill(3) for f in doys] 
    combo = [int(y + d) for y,d in zip(yr_str, doy_str)]
    return combo



### PRECIP FILES #############################################################
# get header (subbasin names, latitude, longitude, elevation)

with open('scripts/swat_setup/swat_base/precip.pcp', 'r') as f:
    lines = f.readlines()
    header = lines[0:4]
    

# csv files
pr_csvs = glob.glob('../Debias/cmip6_swat_cleaned/pr/*.csv')


for pr_csv in pr_csvs:
    # get name
    fname = os.path.basename(pr_csv)[3:-4]
    # define output path 
    fpath_out = '../Debias/cmip6_swat_input/'
    mkdir_p(fpath_out)
    fpath_out = fpath_out + fname + '/'
    mkdir_p(fpath_out)
    print(fname)

    # read as dataframe 
    df = pd.read_csv(pr_csv,index_col=0)
    df.index = pd.to_datetime(df.index) # index to datetime

    pr_combo = convert_dt(df)
    df['dates'] = pr_combo

    filename = fpath_out + 'pcp1.pcp'

    with open(filename, 'w') as p:
        for f in lines:
            p.write(f)
        for a in df.itertuples():
            if int(a[28]) >= 2020001:
                p.write("{:3d}".format(a[28]))
                for i in range(1,28):
                    p.write("{:05.1f}".format(a[i]))
                p.write('\n')
        
        

print("-------------------------------------------------------------------")
### TMP1.TMP #################################################################
with open('scripts/swat_setup/swat_base/temp.tmp', 'r') as f:
    t_lines = f.readlines()
    t_header = t_lines[0:4]

tasmax_csvs = glob.glob('../Debias/cmip6_swat_cleaned/tasmax/*.csv')
 

for tasmax_csv in tasmax_csvs: 
    # csv files
    df_tasmax = pd.read_csv(tasmax_csv, index_col=0)
    tasmin_csv = '../Debias/cmip6_swat_cleaned/tasmin/tasmin' + os.path.basename(tasmax_csv)[6:]
    print(tasmin_csv, tasmax_csv)
    df_tasmin = pd.read_csv(tasmin_csv, index_col=0)
    df_tasmax.index = pd.to_datetime(df_tasmax.index) # index to datetime
    df_tasmin.index = pd.to_datetime(df_tasmin.index) # index to datetime

    # define output path 
    fname = os.path.basename(tasmax_csv)[7:-4]
    fpath_out = '../Debias/cmip6_swat_input/' + fname + '/'
    mkdir_p(fpath_out)
    print(fname)


    # create a DT column of year and julian date (integer)
    tas_combo = convert_dt(df_tasmax)
    df_tasmax['dates'] = tas_combo
    df_tasmin['dates'] = tas_combo

    filename = fpath_out + 'tmp1.tmp'
    with open(filename, 'w') as p:
        for f in t_lines:
            p.write(f) # same header as pcp file
        for t1, t2 in zip(df_tasmin.itertuples(), df_tasmax.itertuples()):
            if int(t1[28]) >= 2020001:
                p.write("{:3d}".format(t1[28]))
                for i in range(1,28):
                    p.write("{:05.1f}{:05.1f}".format(t2[i], t1[i]))
                p.write('\n')