#######################################################################################
## CODIGOS MEMORIA 
## 1 - SELECCION DE PUNTOS A MEDIR EN TERRENO
## Lucas Rivero Iribarne
#######################################################################################


################################################################################################
## Este codigo permite generar zonas para muestrear en terreno, de la siguiente forma:
## 
##    1) carga imagenes sentinel 2 
##    2) Recorta estas imagenes en el area de estudio
##    3) Calcula NDVI
##    4) Con el algoritmo Kmeans clasifica en K zonas (en este caso 3)
##      > es importante asignar una 4 zonas si no se hubieran eliminado caminos previamente 
##
################################################################################################

rm(list = ls()) # limpia "Global Enviroment"
graphics.off() # limpia graficos

# ---- DEFINIR DIRECTORIOS ----
setwd('D:/MEMORIA_LRI') # directorio
dir.in = paste(getwd(),"/Materiales/seleccion_pixeles_terreno", sep = "") # directorio entrada con imgs. Sentinel 2 .SAFE
dir.out = paste(getwd(),"/Resultados/seleccion_pixeles_terreno", sep = "") # directorio salida

ROI = readOGR('D:/MEMORIA_LRI/Materiales/seleccion_pixeles_terreno/shapes/ROI.shp') # Area de estudio

imgs_L2A = grep('.SAFE',list.dirs(path = paste0(dir.in,'/imgs'),full.names = T,recursive = FALSE),ignore.case=TRUE, value=TRUE) # imgs. Sentinel 2 de interes 
fechas = substr(x = imgs_L2A,start = nchar(paste0(dir.in,'/imgs'))+13 ,stop = nchar(paste0(dir.in,'/imgs'))+20) # extrae la fecha de cada img.

NDVI_temp = list() # lista vacia para llenar con indices espectrales
for (i in 1:length(imgs_L2A)) {
  wd0 = (paste0(dir.in,'/imgs'))
  setwd(imgs_L2A[i])
  setwd(paste0(getwd(),'/GRANULE'))
  setwd(paste0(getwd(),'/',list.files()[1]))
  setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
  
  red = mask(crop(raster(list.files()[4]),ROI),ROI) # banda 'red' Sentinel 2
  nir = mask(crop(raster(list.files()[10]),ROI),ROI) # banda '8A' Sentinel 2
  
  ndvi = ((nir/10000 - red/10000)/(nir/10000 + red/10000))
  
  NDVI_temp[[i]] = ndvi
  names(NDVI_temp)[i] = fechas[i]
  print(i)
}
NDVI_stack = stack(NDVI_temp)

Los_Molinos = ROI[c(1:4),] # subset ROI en area huerto 'Los Molinos'
Sofruco = ROI[c(5:10),] # subset ROI en area huerto 'Sofruco'

NDVI_molinos = mask(crop(NDVI_stack[[1]],Los_Molinos),Los_Molinos)
NDVI_sofruco = mask(crop(NDVI_stack[[2]],Sofruco),Sofruco)

a = mask(crop(NDVI_stack[[1]],Los_Molinos),Los_Molinos)
b = mask(crop(NDVI_stack[[2]],Sofruco),Sofruco)
NDVI_AOI = merge(a,b)

# ---- K MEANS PARA DEFINIR ZONAS DE MEDICION ----

### K-Means Los Molinos

{v = getValues(NDVI_molinos)
i = which(!is.na(v))
v = na.omit(v)

E = kmeans(v, 3, iter.max = 100, nstart = 10)
kmeans_raster = raster(NDVI_molinos)
kmeans_raster[i] = E$cluster
plot(kmeans_raster)

setwd(paste0(dir.out,'/Kmeans'))
writeRaster(kmeans_raster, filename="NDVI_kmeans_Los_Molinos.tif", format="GTiff", overwrite=TRUE)
}

### K-Means Sofruco

{v = getValues(NDVI_sofruco)
  i = which(!is.na(v))
  v = na.omit(v)
  
  E = kmeans(v, 3, iter.max = 100, nstart = 10)
  kmeans_raster = raster(NDVI_sofruco)
  kmeans_raster[i] = E$cluster
  plot(kmeans_raster)
  
  setwd(paste0(dir.out,'/Kmeans'))
  writeRaster(kmeans_raster, filename="NDVI_kmeans_sofruco.tif", format="GTiff", overwrite=TRUE)
}


#### FIN ####