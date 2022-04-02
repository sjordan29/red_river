# Script to bias correct CMIP5 model realizations using Aphrodite observations

Sys.setenv(TZ="UTC")
library(ncdf4)
library(RNetCDF)
library(xts)
library(doParallel)
library(foreach)
# source('esd_r/config.R') # rPython is no longer available - SJ just defines file paths in work
source('esd_r/debias.R')

# setting up paths 
path_in_cmip6 <- '/scratch/smj5vup/CMIP6/Debias/cmip6_resample'
path_out_cmip6 <- '/scratch/smj5vup/CMIP6/Debias/cmip6_debiased'
paths_in = Sys.glob(file.path(path_in_cmip6,"*",'*'))
paths_out = sapply(basename(paths_in),function(x){file.path(path_out_cmip6,x)}, USE.NAMES=FALSE)
for (apath in paths_out) dir.create(apath, showWarnings=FALSE, recursive=TRUE)
fpaths_in <- Sys.glob(file.path(paths_in,'*')) # all .nc files 


# Remove files already completed 
fpaths_out_exist <- Sys.glob(file.path(paths_out,'*'))
mask_incmplt <- !(basename(fpaths_in) %in% basename(fpaths_out_exist))
fpaths_in <- fpaths_in[mask_incmplt]

### OBSERVED PR, TASMAX, TASMIN  ########################################################
# read in observed files
ds_obs_tasmin <- nc_open(file.path('/scratch/smj5vup/CMIP6/Debias/obs_resample', 'tasmin_ERA5_1981_2019_3degbuf.nc'))
ds_obs_tasmax <- nc_open(file.path('/scratch/smj5vup/CMIP6/Debias/obs_resample', 'tasmax_ERA5_1981_2019_3degbuf.nc'))
ds_obs_pr <- nc_open(file.path('/scratch/smj5vup/CMIP6/Debias/obs_resample', 'pr_CHIRPS05_1981_2019_3degbuf.nc'))

a_obs_tasmin <-  ncvar_get(ds_obs_tasmin)
a_obs_tasmax <-  ncvar_get(ds_obs_tasmax)
a_obs_pr <-  ncvar_get(ds_obs_pr)

# times
times_obs <- utcal.nc(ds_obs_tasmin$dim$time$units, ds_obs_tasmin$dim$time$vals, type="c")
# times_obs_pr <- utcal.nc(ds_obs_pr$dim$time$units, ds_obs_pr$dim$time$vals, type="c")

nc_close(ds_obs_tasmin)
nc_close(ds_obs_tasmax)
nc_close(ds_obs_pr)

start_yr_base <- '1981'
end_yr_base <- '2018'
idx_train <- paste(start_yr_base,end_yr_base,sep='/')
idx_fut <- c(idx_train,'2019/2039','2040/2069','2070/2099')
idx_all <- unique(c(idx_train,idx_fut))
winsize = 31

## PROJECTION #####################################################################
# read in data 
ds_ptype <- nc_open(fpaths_in[1])
a_ptype <- ncvar_get(ds_ptype)
times_ptype <- utcal.nc(ds_ptype$dim$time$units, ds_ptype$dim$time$vals, type="c")
times_ptype <- as.POSIXct(times_ptype,format='%Y-%d-%m %H:%M:%S', tz="UTC") # format the same
times_ptype <- as.Date(times_ptype, format = '%Y-%d-%m')
times_ptype <- as.POSIXct(times_ptype, format=times_ptype.fmt, tz="UTC")
nc_close(ds_ptype)



mod_ptype <- xts(a_ptype[1,1,],times_ptype)
obs_ptype <- xts(a_obs_tasmax[1,1,],times_obs)  
merge_ptype <- merge.xts(obs_ptype, mod_ptype)
win_masks <- build_window_masks(merge_ptype,idx_all,winsize=winsize)
win_masks1 <- build_window_masks(merge_ptype,idx_all,winsize=1)

bias_funcs <- list('tasmax'=edqmap_tair, 'tasmin'=edqmap_tair, 'pr'=edqmap_prcp)
a_obs <- list('tasmax'=a_obs_tasmax,'tasmin'=a_obs_tasmin,'pr'=a_obs_pr)

## MERGE ###########################################################################

fpath_log <- file.path('/scratch/smj5vup/CMIP6/Debias/logs', 'bias_correct.log')
writeLines(c(""),fpath_log)
cmdArgs <- commandArgs(trailingOnly=TRUE)
numCores <- as.integer(cmdArgs[1]) - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

results <- foreach(fpath=fpaths_in) %dopar% {
  library(ncdf4)
  library(xts)
  source(file.path('/scratch/smj5vup/CMIP6/omoCMIP','esd_r','debias.R'))
  sink(fpath_log, append=TRUE)
  cat(sprintf("Started processing %s \n", fpath), file=fpath_log, append=TRUE)
  
  # open a resampled netcdf file
  ds <- nc_open(fpath) 
  elem <- names(ds$va) # pr, tasmin, or tasmax
  bias_func <- bias_funcs[[elem]] # bias correction is different for precip vs temp
  a <- ncvar_get(ds)
  
  a_debias <- array(NA,dim(a)) # create an empty array
  
  # Loop through each grid cell, debias it, and add it to the new, empty array
  for (r in 1:dim(a_debias)[1]) {
    
    for (c in 1:dim(a_debias)[2]) {
            
      mod_rc <- xts(a[r,c,],times_ptype)
      obs_rc <- xts(a_obs[[elem]][r,c,],times_obs)

      
      # Only perform bias correct if observations (or model observations) are available for the grid cell
      if (! all(is.na(obs_rc))) {
        if (! all(is.nan(mod_rc))){
        merge_rc <- merge.xts(obs_rc, mod_rc)
        xtsAttributes(merge_rc) <- list(modname=sprintf('%s_r%dc%d',basename(fpath),r,c))
        a_debias[r,c,] <- mw_bias_correct(merge_rc[,1], merge_rc[,2], idx_train,
                                          idx_fut, bias_func, win_masks, win_masks1)
        

        
        }}
    }
    
  }
  
  # Output debiased data
  
  # Create output dimensions
  dim_lon <- ncdim_def('lon',units='degrees_east',vals=ds$dim[['lon']]$vals)
  dim_lat <- ncdim_def('lat',units='degrees_north',vals=ds$dim[['lat']]$vals)
  dim_time <- ncdim_def('time', units=ds$dim[['time']]$units, ds$dim[['time']]$vals,
                        calendar=ds$dim[['time']]$calendar)

  # Create output variable
  var_out <- ncvar_def(elem, units=ds$var[[elem]]$units,
                       dim=list(dim_lon,dim_lat,dim_time), missval=NA, prec='double')
  
  fpath_split <- strsplit(fpath,.Platform$file.sep)[[1]]
  path_rcp <- fpath_split[length(fpath_split)-1]
  fname_out <- fpath_split[length(fpath_split)]
  fpath_out <- file.path(path_out_cmip6,path_rcp,fname_out)

  # Create output dataset
  ds_out <- nc_create(fpath_out, var_out, force_v4=TRUE)
  
  
  # Write data
  ncvar_put(ds_out, vals=a_debias)
  
  # Write global attributes
  attrs_global <- ncatt_get(ds,0)
  for (attname in names(attrs_global)) {
    ncatt_put(ds_out,0,attname,attrs_global[[attname]])
  }
  # Add lon/lat attributes
  ncatt_put(ds_out,'lon','standard_name','longitude')
  ncatt_put(ds_out,'lon','long_name','longitude')
  ncatt_put(ds_out,'lon','axis','X')
  ncatt_put(ds_out,'lat','standard_name','latitude')
  ncatt_put(ds_out,'lat','long_name','latitude')
  ncatt_put(ds_out,'lat','axis','Y')

  # Close
  nc_close(ds_out)
  nc_close(ds)
  
  cat(sprintf("Finished processing %s \n", fpath), file=fpath_log, append=TRUE)
  
  sink()
  
  return(1)
}
stopCluster(cl)
