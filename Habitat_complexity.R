#################### Text_S1. R code for processing DEMs ####################
#
# Integrating Three-Dimensional Benthic Habitat Characterization Techniques
# into Ecological Monitoring of coral reefs
# A. Fukunaga, J.H.R. Burns, B.K. Craig and R.K. Kosaki 


library(raster)  # need v. 2.5-8 or ealier or uncomment the two lines to convert NaN to NA
library(rgeos)
library(ggplot2)
library(RColorBrewer)

source("Text_S2_v5.R")  # or run scripts in Text_S2.R
# install.pacakges("rgdal")   # need to have the rgdal pacakge installed

########## analysis prep #########
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

files <- paste("DEMs_moc_mor_chi/", list.files("DEMs_moc_mor_chi/", pattern = ".tif"), sep = "")  # vector of pathes for the DEM files
file_names <- gsub("_DEM.tif", "", gsub("DEMs_moc_mor_chi/", "\\1", files))

for (k in 1:length(files)) {
    
    ras <- raster(files[k]) #Lee el DEM k
    
    #calibrate values to 0
    values(ras) <- values(ras) - min(values(ras), na.rm = T)
    
    resolution <- res(ras)[1] # Determina la resolucion del raster
    ras_g <- as(ras, "SpatialGridDataFrame")
    orig_area <- sum(!is.na(ras_g@data)) * resolution ^ 2 # Calculo del area exacta que ocupa el raster
    
    print(file_names[k]) # imprime el nombre del raster que esta leyendo
    
    extent64 <- aggregate(ras, fac = ceiling(64/res(ras)[1]/100), fun = mean, expand = FALSE, na.rm = FALSE)
    extent64_uniform <- extent64 > -Inf
    extent64_polygon <- rasterToPolygons(extent64_uniform, dissolve = TRUE)
    
    extent128 <- aggregate(ras, fac = ceiling(128/res(ras)[1]/100), fun = mean, expand = FALSE, na.rm = FALSE)
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
                    reef <- aggregate(ras_clip_64, fac = ceiling(dat[[j]]$fac[i]/res(ras)[1]/100), fun = mean, 
                                      expand = FALSE, na.rm = FALSE)
                } else {
                    reef <- aggregate(ras_clip_64, fac = ceiling(dat[[j]]$fac[i]/res(ras)[1]/100), fun = mean, 
                                      expand = FALSE, na.rm = FALSE)
                }
                reef_g <- as(reef, "SpatialGridDataFrame")
                reef_g@data[, 1][is.nan(reef_g@data[, 1])] <- NA   # uncomment if using raster 2.6-7 or later
                #plot(reef_g)
                terrain_res <- terrain_fun(reef, resolution)
                dat[[j]]$s_area[i] <- surfaceArea(reef_g)
                dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (res(reef)[1])^2
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
                    reef <- aggregate(ras_clip_128, fac = ceiling(dat[[j]]$fac[i]/res(ras)/100), fun = mean, 
                                      expand = FALSE, na.rm = FALSE)
                } else {
                    reef <- aggregate(ras_clip_128, fac = ceiling(dat[[j]]$fac[i]/res(ras)/100), fun = mean, 
                                      expand = FALSE, na.rm = FALSE)
                }
                reef_g <- as(reef, "SpatialGridDataFrame")
                reef_g@data[, 1][is.nan(reef_g@data[, 1])] <- NA   # uncomment if using raster 2.6-7 or later
                terrain_res <- terrain_fun(reef, resolution)
                dat[[j]]$s_area[i] <- surfaceArea(reef_g)
                dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (res(reef)[1])^2
                dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
                dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
                dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
                dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
                dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature      }
        } else {
            break
        }
    }
    
    if(is.infinite(min(log(dat[[1]]$s_area)))) { #### Agregado para estandarizar el eje Y a que siempre tenga un rango de 2 unidades
        miny <- floor(min(log(dat[[1]]$s_area[-7])))
    } else {
        miny <- floor(min(log(dat[[1]]$s_area)))
    }
    #Produce el grafico para evaluar dimension fractal a 64
    p64 <- ggplot(dat[[1]], aes(x = log(fac * 0.01), y = log(s_area))) +
        geom_point() +
        geom_smooth(method = "lm", se = F) +
        xlim(c(-5, 1)) + ylim(c(miny, miny + 2)) + 
        labs(x = expression(paste("log(", delta, ")")),
             y = expression(paste("logS(", delta, ")"))) +
        ggtitle(file_names[k]) + 
        theme_bw()
    
    
    if(is.infinite(min(log(dat[[2]]$s_area)))) { #### Agregado para estandarizar el eje Y a que siempre tenga un rango de 2 unidades
        miny <- floor(min(log(dat[[2]]$s_area[-8])))
    } else {
        miny <- floor(min(log(dat[[2]]$s_area)))
    }
    #Produce el grafico para evaluar dimension fractal a 128
    p128 <- ggplot(dat[[2]], aes(x = log(fac * 0.01), y = log(s_area))) +
        geom_point() +
        geom_smooth(method = "lm", se = F) + 
        xlim(c(-5, 1)) + ylim(c(miny, miny + 2)) + 
        labs(x = expression(paste("log(", delta, ")")),
             y = expression(paste("logS(", delta, ")"))) +
        ggtitle(file_names[k]) +
        theme_bw() 
    
    # Los siguientes dos bloques guardan los graficos en la carpeta FD_Graficos
    png(paste("FD_Graficos/", file_names[k], "64FD.png", sep = ""))
    plot(p64)
    dev.off()
    
    png(paste("FD_Graficos/", file_names[k], "64FD_LM.png", sep = ""))
    plot(log(dat[[1]]$s_area/dat[[1]]$area) ~ log(dat[[1]]$fac))
    abline(lm(log(dat[[1]]$s_area/dat[[1]]$area) ~ log(dat[[1]]$fac)))
    dev.off()
    
    png(paste("FD_Graficos/", file_names[k], "128FD.png", sep = ""))
    plot(p128)
    dev.off()
    
    png(paste("FD_Graficos/", file_names[k], "128FD_LM.png", sep = ""))
    plot(log(dat[[2]]$s_area/dat[[2]]$area) ~ log(dat[[2]]$fac))
    abline(lm(log(dat[[2]]$s_area/dat[[2]]$area) ~ log(dat[[2]]$fac)))
    dev.off()
    
    # Los siguientes bloques de codigo generan los valores que seran agregados a hab_data
    d64 <- lm(log(s_area/area) ~ log(fac * resolution), data = dat[[1]])
    slope64 <- coef(d64)[[2]]
    fd64 <- 2 - slope64
    
    d128 <- lm(log(s_area/area) ~ log(fac * resolution), data = dat[[2]])
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

hab_data2 <- hab_data
hab_data2$Site <- factor(gsub("_T[0-9]", "", as.character(hab_data2$file_name)))
hab_data2$Locality <- rep(c("Mochima", "Morrocoy", "Mochima", "Mochima", "Chichirivichi", "Morrocoy", "Morrocoy", "Chichirivichi", "Chichirivichi", "Chichirivichi", "Mochima", "Morrocoy"), each = 4)

write.csv(hab_data2, "habitat_complexity_wFactors.csv", row.names = FALSE)

ggplot(data = hab_data2, aes(x = surface_complexity, y = fd64, color = Site, shape = Locality)) +
    geom_point(size = 4) +
    scale_color_manual(values = brewer.pal(12, "Paired"))


ggplot(data = hab_data2, aes(x = surface_complexity, y = fd64, color = Locality)) +
    geom_point(size = 4) + 
    scale_color_manual(values = brewer.pal(3, "Dark2"))
