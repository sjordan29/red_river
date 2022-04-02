# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:43:55 2021

@author: Sarah
"""

from pyesgf.search import SearchConnection
conn = SearchConnection('https://esgf-data.dkrz.de/esg-search', distrib=True)

ctx = conn.new_context(
project='CMIP5',
experiment='rcp45',
model='HadCM3',
ensemble='r1i1p1',
time_frequency='mon',
realm='atmos',
data_node='esgf-data1.ceda.ac.uk',
)
ctx.hit_count

result = ctx.search()[0]
result.dataset_id

files = result.file_context().search()
for file in files:
    if 'tasmax' in file.opendap_url:
        tasmax_url = file.opendap_url
        print(tasmax_url)


from pyesgf.logon import LogonManager
lm = LogonManager()
lm.logoff()
lm.is_logged_on()


lm.logon(hostname='esgf-data.dkrz.de', interactive=True, bootstrap=True)
lm.is_logged_on()


import xarray as xr
ds = xr.open_dataset(tasmax_url, chunks={'time': 120})
print(ds)
