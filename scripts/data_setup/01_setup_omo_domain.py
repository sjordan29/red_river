'''
Script to create rasters representing the Omo spatial domain grid
'''

from osgeo import gdalconst
import esd
import numpy as np
import os
import xarray as xr

if __name__ == '__main__':
    # observed precipitation - mask to desired longitude and latitude 
    bnds = esd.cfg.bbox
    ds = xr.open_dataset(esd.cfg.fpath_obs_pr)
    
    mask_lat = (ds.latitude >= bnds[1]) & (ds.latitude <= bnds[3]) 
    mask_lon = (ds.longitude >= bnds[0]) & (ds.longitude <= bnds[2])
    
    da = ds.precip[:,mask_lat,mask_lon].load()
    
    da = da.mean(dim='time')
    a = da.values
    a = np.flipud(a)
    a = np.ma.masked_invalid(a)
    
    lon = da.longitude.values
    lat = da.latitude.values
    lat = lat[::-1]
    
    fpath_bnds_gtif = os.path.join(esd.cfg.path_data_bbox, 'omo_bnds.tif')
    
    esd.util.grid_wgs84_to_raster(fpath_bnds_gtif, a, lon, lat, gdalconst.GDT_Float32, ndata=-9999)
        
    ds = esd.util.RasterDataset(fpath_bnds_gtif)
    
    # Get the upper left and right point x coordinates in the raster's projection
    ulX = ds.geo_t[0] + (ds.geo_t[1] / 2.0)
    urX = ds.get_coord(0, ds.gdal_ds.RasterXSize - 1)[1]

    ulx = ds.geo_t[0]
    uly = ds.geo_t[3]
    
    ptx = ulx + .5
    pty = uly - .5
    
    # 1 degree resample
    lat = np.linspace(pty, pty - 5, num=6)
    lon = np.linspace(ptx, ptx + 6, num=7)
    
    a = np.random.rand(lat.size, lon.size)
    
    fpath_1degbnds_gtif = os.path.join(esd.cfg.path_data_bbox, 'omo_1degbnds.tif')
    
    esd.util.grid_wgs84_to_raster(fpath_1degbnds_gtif, a, lon, lat,
                                  gdalconst.GDT_Float32, ndata=-9999)
