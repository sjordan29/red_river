####### packages ################################)
library(raster)
library(ncdf4)
library(rgdal)

####### define files #############################
fpath_in <- '/scratch/smj5vup/CMIP6/Debias/cmip6_downscaled/tasmin'
fpath_out <- '/scratch/smj5vup/CMIP6/Debias/cmip6_swat'
dir.create(fpath_out)
dir.create(file.path(fpath_out,'tasmin'))
fpath_out <- file.path(fpath_out,'tasmin')


##### translate array no to projection #########
tasmin_files <- list.files(fpath_in)
cmdArgs <- commandArgs(trailingOnly=TRUE)
arrayNo <- as.integer(cmdArgs[1])
tasminPath <- tasmin_files[arrayNo]
base <- tools::file_path_sans_ext(basename(tasminPath))
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
  tasmin.stack <- stack(file.path(fpath_in,tasminPath))
  print("stack")
  ex <- raster::extract(tasmin.stack, poly, fun=mean, na.rm=TRUE, df=TRUE, weights=TRUE)
  print("extract")
  write.csv(ex, file.path(f_out, paste(bn,".csv",sep="")))
}


##### execute function #######################
zonStatNC(mask, base, fpath_out)
