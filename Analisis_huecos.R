#### Carga datos de huecos ####
require(readxl)
Data <- read_xlsx("Huecos.xlsx")

#### Revisar que no haya mediciones impares ####
for(i in 1:length(unique(Data$Transecta))) {
    print(paste(unique(Data$Transecta)[i], nrow(Data[Data$Transecta == unique(Data$Transecta)[i], ]), nrow(Data[Data$Transecta == unique(Data$Transecta)[i], ]) %% 2))
    if(nrow(Data[Data$Transecta == unique(Data$Transecta)[i], ]) %% 2 != 0) {
        stop("La ultima transecta no tenia las medidas pares")
    }
}

#### Crear matriz de huecos y sus respectivas areas ####
NewData <- data.frame(Localidad = rep("NA", nrow(Data)/2),
                      Sitio = rep("NA", nrow(Data)/2),
                      Transecta = rep("NA", nrow(Data)/2),
                      AreaHueco = as.numeric(NA), stringsAsFactors = F)


for(i in seq(1, nrow(Data), 2)) {
    NewData$Localidad[(i+1)/2] <- Data$Localidad[i]
    NewData$Sitio[(i+1)/2] <- gsub("_t[0-9]", "", Data$Transecta[i])
    NewData$Transecta[(i+1)/2] <- Data$Transecta[i]
    NewData$AreaHueco[(i+1)/2] <- as.numeric(Data[i, 3])*as.numeric(Data[i+1, 3])
}

#### Agregar Categoria de Tamaño de Hueco a la matriz anterior ####
NewData$CatHueco <- "NA"

CatStats <- boxplot.stats(NewData$AreaHueco)$stats

#NewData$CatHueco[NewData$AreaHueco < CatStats[1]] <- "0" #Si quieres agregar una categoria extra al principio
NewData$CatHueco[NewData$AreaHueco >= CatStats[1] & NewData$AreaHueco < CatStats[2]] <- "1"
NewData$CatHueco[NewData$AreaHueco >= CatStats[2] & NewData$AreaHueco < CatStats[3]] <- "2"
NewData$CatHueco[NewData$AreaHueco >= CatStats[3] & NewData$AreaHueco < CatStats[4]] <- "3"
NewData$CatHueco[NewData$AreaHueco >= CatStats[4] & NewData$AreaHueco < CatStats[5]] <- "4"
NewData$CatHueco[NewData$AreaHueco >= CatStats[5]] <- "5"

#### Añade el conteo de huecos ####
NewData$Conteo <- 1

#### Resume la tabla anterior por transecta sitio, y localidad ####
SummaryData <- cbind(aggregate(Conteo ~ Localidad + Sitio + Transecta, NewData, sum), 
AreaHuecos = aggregate(AreaHueco ~ Localidad + Sitio + Transecta, NewData, sum)[, 4])

#### Agrega las categorias de huecos como columnas de conteo a la tabla resumen ####
require(dplyr)
require(tidyr)
A <- as.data.frame(table(NewData[c(3, 5)]))
A$CatHueco <- paste("Size_", A$CatHueco, sep = "")
A <- spread(data = A, key = "CatHueco", value = "Freq")

SummaryData <- left_join(SummaryData, A, by = "Transecta")

#### Guarda la tabla en un csv ####
write.csv(SummaryData, "SummaryHuecos.csv", row.names = F)
#write.csv(NewData, "HuecosCompletas.csv", row.names = F)