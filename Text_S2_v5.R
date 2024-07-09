##### Text_S2. R code for a custom function to obtain slope and curvature #####
#
# Integrating Three-Dimensional Benthic Habitat Characterization Techniques
# into Ecological Monitoring of Coral Reefs
# A. Fukunaga, J.H.R. Burns, B.K. Craig and R.K. Kosaki 


terrain_fun <- function(data, cell_size) {
  
  d0 <- as.matrix(data)
  if (ncol(d0) == 1 | nrow(d0) == 1) {
    terrain_list <- list(NA, NA, NA)
    names(terrain_list) <- c("mean_slope", "mean_profile_curvature", "mean_plan_curvature")
  } else {
    da <-  cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                 rbind(matrix(NA, nrow = 1, ncol = ncol(d0) - 1), 
                       matrix(d0[1 : (nrow(d0) - 1), 1 : (ncol(d0) - 1)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1)))
    db <- rbind(matrix(NA, nrow = 1, ncol = ncol(d0)), 
                matrix(d0[1 : (nrow(d0) - 1), ], nrow = nrow(d0) - 1, ncol = ncol(d0)))
    dc <- cbind(rbind(matrix(NA, nrow = 1, ncol = ncol(d0) - 1), 
                      matrix(d0[1 : (nrow(d0) - 1), 2 : ncol(d0)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1)),
                matrix(NA, nrow = nrow(d0), ncol = 1))
    dd <- cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                matrix(d0[, 1 : (ncol(d0) - 1)], nrow = nrow(d0), ncol = ncol(d0) -1))
    df <- cbind(matrix(d0[, 2 : ncol(d0)], nrow = nrow(d0), ncol = ncol(d0) - 1), 
                matrix(NA, nrow = nrow(d0), ncol = 1))
    dg <- cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                rbind(matrix(d0[2 : nrow(d0), 1 : (ncol(d0) - 1)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1), 
                      matrix(NA, nrow = 1, ncol = ncol(d0) - 1)))
    dh <- rbind(matrix(d0[2 : nrow(d0), ], nrow = nrow(d0) - 1, ncol = ncol(d0)), 
                matrix(NA, nrow = 1, ncol = ncol(d0)))
    di <- cbind(rbind(matrix(d0[2 : nrow(d0), 2 : ncol(d0)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1), 
                      matrix(NA, nrow = 1, ncol = (ncol(d0) - 1))), 
                matrix(NA, nrow = nrow(d0), ncol = 1))
    
    ## slope
    
    x_rate <- ((dc + (2 * df) + di) - (da + (2 * dd) + dg)) / (8 * cell_size)
    y_rate <- ((dg + (2 * dh) + di) - (da + (2 * db) + dc)) / (8 * cell_size)
    
    slope_degrees <- atan(sqrt(x_rate ^ 2 + y_rate ^ 2)) * 57.29578
    mean_slope <- mean(slope_degrees, na.rm = TRUE)

    ## curvature
    
    coefD <- (((dd + df) / 2) - d0) / (cell_size ^ 2)
    coefE <- (((db + dh) / 2) - d0) / (cell_size ^ 2)
    coefF <- (dc + dg -da - di) / (4 * (cell_size ^ 2))
    coefG <- (df - dd) / 2 * cell_size
    coefH <- (db - dh) / 2 * cell_size
    
    prof_curv <- -2 * ((coefD * coefG ^ 2 + coefE * coefH ^ 2 + coefF * coefG * coefH) / 
                         (coefG ^ 2 + coefH ^ 2))
    mean_prof_curv <- mean(prof_curv, na.rm = TRUE)
    
    plan_curv <- 2 * ((coefD * coefH ^ 2 + coefE * coefG ^ 2 - coefF * coefG * coefH) / 
                        (coefG ^ 2 + coefH ^ 2))
    mean_plan_curv <- mean(plan_curv, na.rm = TRUE)
    
    ## output list
    
    terrain_list <- list(mean_slope, mean_prof_curv, mean_plan_curv)
    names(terrain_list) <- c("mean_slope", "mean_profile_curvature", "mean_plan_curvature")
  }
  
  return(terrain_list)
  
}
