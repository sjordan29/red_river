####### packages ################################)
library(raster)
library(ncdf4)
library(rgdal)

####### define files #############################
fpath_in <- '/scratch/smj5vup/CMIP6/Debias/cmip6_downscaled/tasmax'
fpath_out <- '/scratch/smj5vup/CMIP6/Debias/cmip6_swat'
dir.create(fpath_out)
dir.create(file.path(fpath_out,'tasmax'))
fpath_out <- file.path(fpath_out,'tasmax')


##### translate array no to projection #########
tasmax_files <- list.files(fpath_in)
cmdArgs <- commandArgs(trailingOnly=TRUE)
arrayNo <- as.integer(cmdArgs[1])
tasmaxPath <- tasmax_files[arrayNo]
base <- tools::file_path_sans_ext(basename(tasmaxPath))
print(base)
# dir.create(file.path(fpath_out, base))
#fpath_out <- file.path(fpath_out, base)


####### set paths ##############################
shpPath <- "scripts/swat_setup/shapefile"
shp <- "subs1"
mask <- shapefile(paste(shpPath,"/",shp,".shp",sep=""))
print("mask setup!")


###### function ###############################
zonStatNC = function(poly, bn, f_out){
  tasmax.stack <- stack(file.path(fpath_in,tasmaxPath))
  print("stack")
  ex <- raster::extract(tasmax.stack, poly, fun=mean, na.rm=TRUE, df=TRUE, weights=TRUE)
  print("extract")
  write.csv(ex, file.path(f_out, paste(bn,".csv",sep="")))
}


##### execute function #######################
zonStatNC(mask, base, fpath_out)
