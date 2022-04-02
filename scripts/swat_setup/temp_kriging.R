library(rgdal)
library(atakrig)
library(raster)
library(MASS)

######## functions ##############################
source("scripts/swat_setup/rivanna-basefiles.R")

####### file paths #############################
fpath_in_tasmin <- '/scratch/smj5vup/CMIP6/Debias/cmip6_downscaled/tasmin'
fpath_in_tasmax <- '/scratch/smj5vup/CMIP6/Debias/cmip6_downscaled/tasmax'

fpath_out <- '/scratch/smj5vup/CMIP6/Debias/cmip6_swat/'
dir.create(fpath_out)

dir.create(file.path(fpath_out,'tasmin'))
fpath_out_tasmin <- file.path(fpath_out,'tasmin')

dir.create(file.path(fpath_out,'tasmax'))
fpath_out_tasmax <- file.path(fpath_out,'tasmax')

###### projections from array number ############
tasmin_files <- list.files(fpath_in_tasmin)
tasmax_files <- list.files(fpath_in_tasmax)

cmdArgs <- commandArgs(trailingOnly=TRUE)
arrayNo <- as.integer(cmdArgs[1])
tasminPath <- tasmin_files[arrayNo]
tasmin.base <- tools::file_path_sans_ext(basename(tasminPath))
dir.create(file.path(fpath_out_tasmin, tasmin.base))
fpath_out_tasmin <- file.path(fpath_out_tasmin, tasmin.base)
tasmin.stack <- stack(file.path(fpath_in_tasmin,tasminPath))

tasmaxPath <- tasmax_files[arrayNo]
tasmax.base <- tools::file_path_sans_ext(basename(tasmaxPath))
dir.create(file.path(fpath_out_tasmax, tasmax.base))
fpath_out_tasmax <- file.path(fpath_out_tasmax, tasmax.base)
tasmax.stack <- stack(file.path(fpath_in_tasmax,tasmaxPath))

####### set paths ##############################
shpPath <- "scripts/swat_setup/shapefile"
shp <- "subs1"
maskfile <- rgdal::readOGR(dsn = shpPath, layer = shp)


####### dates #################################
startDate <- paste("1981/01/01",sep="")
endDate <- paste("2099/12/31", sep="")
days <- seq(as.Date(startDate), as.Date(endDate), "days")


##### execute function #######################
i <- 1
for (day in days){
  val<-as.Date(day, origin="1970-01-01")
  valStr <- toString(val)
  print(valStr)
  tasmin.day <- tasmin.stack[[i]]
  tasmax.day <- tasmax.stack[[i]]
  predictKrig(valStr, 10000, maskfile, 10000, tasmin.day, tasmax.day, fpath_out_tasmin, fpath_out_tasmax)
  i <- i+1
}