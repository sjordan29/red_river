discretizePcp = function(pcp_file, pcp_res,msk){
  # precip data
  precip <- raster(file.path(precipPath, pcp_file))
  e <- as(extent(35, 39, 4, 10), 'SpatialPolygons')
  crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  precip.crop <- crop(precip, e)
  sr <-"+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  precip.proj <- projectRaster(precip.crop, crs=sr)
  precip.mask <- mask(precip.proj, mask=msk)
  precip.d <- discretizeRaster(precip.mask,pcp_res)
  # plot(precip.mask)
  # plot(msk, add=TRUE)
  # points(precip.d$areaValues$centx, precip.d$areaValues$centy)
  # print(precip.d)
  
}

discretizeTmin = function(tmin_file, temp_res, msk) {
  tmin <- tmin_file
  crs(tmin) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  e <- as(extent(35, 39, 4, 10), 'SpatialPolygons')
  crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  tmin.crop <-  crop(tmin, e)
  sr <-"+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  tmin.proj <- projectRaster(tmin.crop, crs=sr)
  tmin.mask <- mask(tmin.proj, mask=msk)
  tmin.d <- discretizeRaster(tmin.mask,temp_res)
  #print(tmin.d)
  
}

discretizeTmax = function(tmax_file, temp_res, msk) {
  #tmax data
  print('started')
  tmax <- tmax_file
  print('tmax')
  crs(tmax) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  e <- as(extent(35, 39, 4, 10), 'SpatialPolygons')
  crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  tmax.crop <-  crop(tmax, e)
  sr <-"+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  tmax.proj <- projectRaster(tmax.crop, crs=sr)
  tmax.mask <- mask(tmax.proj, mask=msk)
  tmax.d <- discretizeRaster(tmax.mask,temp_res)
  # print(tmax.d)
}


discretizeShp = function(msk, msk_res){
  # shapefile
  subs.d <- discretizePolygon(msk, cellsize=msk_res, id="OBJECTID")
}

transformPcp = function(precip.d){
  # boxcox transformation
  
  vals <- precip.d$areaValues$value + 0.1
  # print(vals)
  Box = boxcox(vals~1, lambda=seq(-5,5,0.1))
  Cox = data.frame(Box$x, Box$y)
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]
  precip.d$areaValues$value <- (vals ^ lambda - 1)/lambda
  T_box = (vals ^ lambda - 1)/lambda
  output <- c(lambda, precip.d)
  
  #     vals = as.vector(precip.d$areaValues$value + 0.1)
  
  #     bc_obj <- boxcox(vals, standardize = TRUE)
  #     # bc_obj <- predict(bc_obj)
  #     return(bc_obj)
}

normalize = function(discreteRaster){
  mu = mean(discreteRaster$areaValues$value)
  sigma = sd(discreteRaster$areaValues$value)
  discreteRaster$areaValues$value = (discreteRaster$areaValues$value - mu) / sigma
  information <- c(mu, sigma, discreteRaster)
  return(information)
}

kriging = function(pt.list, sv.ck, vble,disArea){
  pred <- ataCoKriging(pt.list, unknownVarId=vble, unknown=disArea,
                       ptVgms=sv.ck, oneCondition=TRUE, auxRatioAdj=TRUE, showProgress = TRUE)    
}


denorm = function(prediction_df, mean_val, sd_val){
  prediction_df$pred <- prediction_df$pred * sd_val + mean_val
  return(prediction_df)
}

detrans = function(prediction_df, lambda){
  prediction_df$pred <- exp(log(lambda * prediction_df$pred  + 1) / lambda) - 0.1
  return(prediction_df)
}



prep = function(tmin_file, tmax_file, temp_res,msk,msk_res){
  ######## DISCRETIZE ######################################
  tmin.d <- discretizeTmin(tmin_file, temp_res,msk)
  print('tmin.d')
  tmax.d <- discretizeTmax(tmax_file, temp_res,msk)
  print('tmax.d')
  
  
  ####### NORMALIZE AND TRANSFORM #########################
  
  # transform and normalize tmin
  tmin_N <- normalize(tmin.d)
  tmin_mu <- tmin_N[[1]]
  tmin_sigma <- tmin_N[[2]]
  tmin_N <- tmin_N[3:4]
  tmin_fitted <- tmin_N
  
  
  # transform and normalize tmax
  tmax_N <- normalize(tmax.d)
  tmax_mu <- tmax_N[[1]]
  tmax_sigma <- tmax_N[[2]]
  tmax_N <- tmax_N[3:4]
  tmax_fitted <- tmax_N
  
  output <- c(tmin_fitted, tmax_fitted, tmin_mu, tmin_sigma, tmax_mu, tmax_sigma)
}


predictKrig = function(date, temp_res, msk, msk_res, tmin_file, tmax_file, tasmin_out, tasmax_out){
  ########## DATE TO FILES ############################
  # tmin_file  <- paste(date,"_min.tif",sep="")
  # print(tmin_file)
  # tmax_file <- paste(date, "_max.tif", sep="")
  # print(tmax_file)
  
  ########## NORMALIZE ################################
  norm <- prep(tmin_file, tmax_file, temp_res,msk,msk_res)
  
  t.d <-norm[1:2]
  T.d <- norm[3:4]
  
  tmin_mu <- norm[5][[1]]
  tmin_sigma <- norm[6][[1]]
  tmax_mu <- norm[7][[1]]
  tmax_sigma <- norm[8][[1]]
  
  subs.d <- discretizeShp(msk, msk_res)
  
  ######### KRIGING ###################################
  #  pt.list <- list(tmin = t.d, tmax=T.d)
  #  sv.ck <- deconvPointVgmForCoKriging(pt.list, model="Exp",ngroup=15, rd=0.4, fixed.range=9e4)
  #  tmin <- kriging(pt.list, sv.ck, "tmin",subs.d)
  #  tmax <- kriging(pt.list, sv.ck, "tmax",subs.d)
  vgm.tmin <- deconvPointVgm(t.d, model="Exp", ngroup=15, rd=0.4, fig=TRUE)
  tmin <- ataKriging(t.d, subs.d, vgm.tmin, showProgress = TRUE)
  vgm.tmax <- deconvPointVgm(T.d, model="Exp", ngroup=15, rd=0.4, fig=TRUE)
  tmax <- ataKriging(T.d, subs.d, vgm.tmax, showProgress = TRUE)
  
  
  ###### DENORM ######################################
  # denorm
  tmin.denorm <- denorm(tmin, tmin_mu, tmin_sigma)
  tmax.denorm <- denorm(tmax, tmax_mu, tmax_sigma)
  
  ##### TO CSV ########################################
  write.csv(tmin.denorm, file.path(tasmin_out, paste(date,".csv",sep="")))
  write.csv(tmax.denorm, file.path(tasmax_out, paste(date,".csv",sep="")))
  
}