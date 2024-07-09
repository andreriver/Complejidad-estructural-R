#################### Text_S1. R code for processing DEMs ####################
#
# Integrating Three-Dimensional Benthic Habitat Characterization Techniques
# into Ecological Monitoring of coral reefs
# A. Fukunaga, J.H.R. Burns, B.K. Craig and R.K. Kosaki 


library(raster)  # need v. 2.5-8 or ealier or uncomment the two lines to convert NaN to NA
library(rgeos)
library(ggplot2)

source("Text_S2.R")  # or run scripts in Text_S2.R
# install.pacakges("rgdal")   # need to have the rgdal pacakge installed

########## analysis prep #########
resolution <- 0.01  # enter DEM resolution in meter

### initialize output data frame
hab_data <- data.frame(file_name = character(0), 
                       orig_area = numeric(0), 
                       surface64 = numeric(0), planerS64 = numeric(0), 
                       surface128 = numeric(0), planerS128 = numeric(0), 
                       max_h = numeric(0), min_h = numeric(0), 
                       fd64 = numeric(0), fd128 = numeric(0), 
                       surface_complexity = numeric(0), 
                       mean_slope = numeric(0), 
                       mean_profile_curvature = numeric(0), 
                       mean_plan_curvature = numeric(0))

########## process files #########

file_path <- "~/R_wd/DEMs/"  # enter an appropriate path to the folder containing DEMs
file_names <- c("FFS_4297", "LAY_5046", "LIS_4199", "KUR_4059")  # enter file (site) names
file_suf <- "_DEM_1cm.tif"  # enter any suffix common in DEM file names

files <- paste(file_path, file_names, file_suf, sep = "")  # vector of pathes for the DEM files

for (k in 1:length(files)) {
  
  ras <- raster(files[k])
  ras_g <- as(ras, "SpatialGridDataFrame")
  orig_area <- sum(!is.na(ras_g@data)) * resolution ^ 2
  
  print(file_names[k])
  
  extent64 <- aggregate(ras, fac = 64, fun = mean, expand = FALSE, na.rm = FALSE)
  extent64_uniform <- extent64 > -Inf
  extent64_polygon <- rasterToPolygons(extent64_uniform, dissolve = TRUE)
  
  extent128 <- aggregate(ras, fac = 128, fun = mean, expand = FALSE, na.rm = FALSE)
  extent128_uniform <- extent128 > -Inf
  extent128_polygon <- rasterToPolygons(extent128_uniform, dissolve = TRUE)
  
  ras_clip_64 <- mask(crop(ras, extent(extent64_polygon)), extent64_polygon)
  ras_clip_128 <- mask(crop(ras, extent(extent128_polygon)), extent128_polygon)
  
  # # uncomment to plot raster
  # plot(ras)
  # plot(ras_clip_64)
  # plot(ras_clip_128)
  
  dat64 <- data.frame(fac = c(1, 2, 4, 8, 16, 32, 64), 
                      s_area = NA, area = NA, 
                      max_height = NA, min_height = NA, 
                      mean_slope = NA,  
                      mean_profile_curvature = NA, mean_plan_curvature = NA)
  dat128 <- data.frame(fac = c(1, 2, 4, 8, 16, 32, 64, 128), 
                       s_area = NA, area = NA, 
                       max_height = NA, min_height = NA, 
                       mean_slope = NA, 
                       mean_profile_curvature = NA, mean_plan_curvature = NA)
  
  dat <- list(dat64, dat128)
  
  for (j in 1:length(dat)) {
    
    if (j == 1) {
      for (i in 1:nrow(dat[[j]])) {
        
        print(dat[[j]]$fac[i])
        
        if (dat[[j]]$fac[i] == 1) {
          reef <- ras_clip_64
        } else {
          reef <- aggregate(ras_clip_64, fac = dat[[j]]$fac[i], fun = mean, 
                            expand = FALSE, na.rm = FALSE)
        }
        reef_g <- as(reef, "SpatialGridDataFrame")
        #reef_g@data[, 1][is.nan(reef_g@data[, 1])] <- NA   # uncomment if using raster 2.6-7 or later
        terrain_res <- terrain_fun(reef, resolution)
        dat[[j]]$s_area[i] <- surfaceArea(reef_g)
        dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * resolution) ^ 2
        dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
        dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
        dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature
      }
    } else if (j == 2) {
      for (i in 1:nrow(dat[[j]])) {
        
        print(dat[[j]]$fac[i])
        
        if (dat[[j]]$fac[i] == 1) {
          reef <- ras_clip_128
        } else {
          reef <- aggregate(ras_clip_128, fac = dat[[j]]$fac[i], fun = mean, 
                            expand = FALSE, na.rm = FALSE)
        }
        reef_g <- as(reef, "SpatialGridDataFrame")
        #reef_g@data[, 1][is.nan(reef_g@data[, 1])] <- NA   # uncomment if using raster 2.6-7 or later
        terrain_res <- terrain_fun(reef, resolution)
        dat[[j]]$s_area[i] <- surfaceArea(reef_g)
        dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * resolution) ^ 2
        dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
        dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
        dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature      }
    } else {
      break
    }
  }
  
  p64 <- ggplot(dat[[1]], aes(x = log(fac * 0.01), y = log(s_area))) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    xlim(c(-5, 1)) + ylim(c(4, 6)) + 
    labs(x = expression(paste("log(", delta, ")")),
         y = expression(paste("logS(", delta, ")"))) +
    ggtitle(file_names[k])
  theme_bw()
  
  p128 <- ggplot(dat[[2]], aes(x = log(fac * 0.01), y = log(s_area))) +
    geom_point() +
    geom_smooth(method = "lm", se = F) + 
    xlim(c(-5, 1)) + ylim(c(4, 6)) + 
    labs(x = expression(paste("log(", delta, ")")),
         y = expression(paste("logS(", delta, ")"))) +
    ggtitle(file_names[k])
  theme_bw()
  
  plot(p64)
  plot(p128)
  
  d64 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[1]])
  slope64 <- coef(d64)[[2]]
  fd64 <- 2 - slope64
  
  d128 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[2]])
  slope128 <- coef(d128)[[2]]
  fd128 <- 2 - slope128
  
  surface64 <- dat[[1]]$s_area[1]
  surface128 <- dat[[2]]$s_area[1]
  
  planerS64 <- dat[[1]]$area[1]
  planerS128 <- dat[[2]]$area[1]
  
  surface_complexity <- dat[[1]]$s_area[1]/dat[[1]]$area[1]
  
  max_h <- dat[[1]]$max_height[1]
  min_h <- dat[[1]]$min_height[1]
  
  mean_slope <- dat[[1]]$mean_slope[1]

  mean_profile_curvature <- dat[[1]]$mean_profile_curvature[1]
  mean_plan_curvature <- dat[[1]]$mean_plan_curvature[1]
  
  temp <- data.frame(file_name = file_names[k], 
                     orig_area = orig_area, 
                     surface64 = surface64, planerS64 = planerS64,
                     surface128 = surface128, planerS128 = planerS128, 
                     max_h = max_h, min_h = min_h, 
                     fd64 = fd64, fd128 = fd128, 
                     surface_complexity = surface_complexity, 
                     mean_slope = mean_slope, 
                     mean_profile_curvature = mean_profile_curvature, 
                     mean_plan_curvature = mean_plan_curvature)
  
  hab_data <- rbind(hab_data, temp)
  
}

hab_data

########## export output data frame #########
write.csv(hab_data, "habitat_complexity.csv", row.names = FALSE)
