files <- list.files(pattern = "cloud.txt")

dir.create("Summary10cm")

SummaryDF <- data.frame(file = "", n = rep(0, length(files)), 
           Max_Z = 0, Mean_Z = 0, Var_Z = 0, Kurtosis_Z = 0, Skewness_Z = 0,
           Max_Ang = 0, Mean_Ang = 0, Var_Ang = 0, Kurtosis_Ang = 0, Skewness_Ang = 0, 
           stringsAsFactors = F)

system_time <- Sys.time()
for(j in 1:length(files)) {
    
    numrows <- 100000 #Filas a leer
    skiprows <- 0 #Filas a saltar - se actualiza en cada iteracion
    
    DF <- NULL #Matriz a ser llenada con la data resumida
    
    i <- 1 #numero de iteraciones - se actualiza en cada iteracion
    
    repeat {
        print(paste("File:", files[j], "| Iteration =", i))
        A <- read.csv(files[j], header = F, nrows = numrows, sep = " ", skip = skiprows) #lectura de la matriz
        
        D <- data.frame(table(X = round(A[, 1], 1),
                              Y = round(A[, 2], 1), 
                              Z = round(A[, 3], 1), 
                              Angle = round(A[, 9]*180/pi)))
        
        
        DF <- rbind(DF, D[D$Freq != 0, ])
        
        rm(D); gc()
        
        DF <- aggregate(Freq ~ X + Y + Z + Angle, data = DF, FUN = sum)
        
        DF <- DF[order(DF$Z), ]
        
        skiprows = numrows + skiprows #Actualiza las filas a saltar
        i <- i + 1 #Actualiza la iteracion en la que se va
        if(nrow(A) < numrows) {
            print(paste("Finished with:", files[j]))
            rm(A); gc()
            break() #Condicion para romper el loop  
        } 
        rm(A); gc()
        #if(i == 9) break()
    }
    
    
    DF$X <- factor(DF$X, levels = levels(DF$X)[order(as.numeric(levels(DF$X)))])
    DF$Y <- factor(DF$Y, levels = levels(DF$Y)[order(as.numeric(levels(DF$Y)))])
    DF$Z <- factor(DF$Z, levels = levels(DF$Z)[order(as.numeric(levels(DF$Z)))])
    DF$Angle <- factor(DF$Angle, levels = levels(DF$Angle)[order(as.numeric(levels(DF$Angle)))])
    
    DF$X <- as.numeric(as.character(DF$X))
    DF$Y <- as.numeric(as.character(DF$Y))
    DF$Z <- as.numeric(as.character(DF$Z))
    DF$Angle <- as.numeric(as.character(DF$Angle))
    
    #### Estandarizacion de Z por su minimo ####
    DF$Z <- DF$Z - min(DF$Z)
    
    #### Agregar archivo y numero de puntos
    SummaryDF$file[j] <- files[j]
    SummaryDF$n[j] <- sum(DF$Freq)
    
    #### Agregar estimadores de Z ####
    SummaryDF$Max_Z[j] <- max(DF$Z)
    SummaryDF$Mean_Z[j] <- sum(DF$Z*DF$Freq)/sum(DF$Freq)
    SummaryDF$Var_Z[j] <- sum(((DF$Z - SummaryDF$Mean_Z[j])^2)*(DF$Freq))/(SummaryDF$n[j])
    SummaryDF$Kurtosis_Z[j] <- (sum(((DF$Z - SummaryDF$Mean_Z[j])^4)*(DF$Freq))/((SummaryDF$n[j])*(sqrt(SummaryDF$Var_Z[j]))^4))-3
    SummaryDF$Skewness_Z[j] <- ((1/SummaryDF$n[j])*sum((((DF$Z - SummaryDF$Mean_Z[j])^3)*DF$Freq)))/(sqrt(SummaryDF$Var_Z[j])^3)
    
    #### Agregar estimadores de angulos ####
    SummaryDF$Max_Ang[j] <- max(abs(DF$Angle))
    SummaryDF$Mean_Ang[j] <- sum(abs(DF$Angle)*DF$Freq)/sum(DF$Freq)
    SummaryDF$Var_Ang[j] <- sum(((abs(DF$Angle) - SummaryDF$Mean_Ang[j])^2)*(DF$Freq))/(SummaryDF$n[j])
    SummaryDF$Kurtosis_Ang[j] <- (sum(((abs(DF$Angle) - SummaryDF$Mean_Ang[j])^4)*(DF$Freq))/((SummaryDF$n[j])*(sqrt(SummaryDF$Var_Ang[j]))^4))-3
    SummaryDF$Skewness_Ang[j] <- ((1/SummaryDF$n[j])*sum(((abs(DF$Angle) - SummaryDF$Mean_Ang[j])^3)*(DF$Freq)))/(sqrt(SummaryDF$Var_Ang[j])^3)
    
    
    write.table(DF, paste("Summary10cm", files[j], sep = "/"), row.names = F, sep = " ")
}

write.table(SummaryDF, paste("Summary10cm", "DataSummary.csv", sep = "/"), row.names = F, sep = " ")