#######################################################################################
## CODIGOS MEMORIA 
## 3 - Espacializacion y temporalizacion de variables
## Lucas Rivero Iribarne
#######################################################################################


############################################################################################################
## Este codigo espacializa y temporaliza las variables de interes 
## de la siguiente forma:
##
##  1) carga la serie temporal de sentinel 2 L2A
##  2) Calcula WDVI, TSAVI, LAI, fPAR, NDVI y TSAVI para cada escena (y guarda como .tiff en directorio)
##  3) Extrae los valores para el pixel donde hay sonda
##  4) Interpola las variables entre las fechas con escenas
##  5) carga la serie temporal de variables biofisicas de SNAP
##  6) Extrae los valores de estas variables en el pixel sonda
##  7) Interpola las variables entre las fechas con escenas
##  8) Guarda variables en .tiff en directorio
##
############################################################################################################

rm(list = ls()) # limpia "Global Enviroment"
graphics.off() # limpia graficos
t0 = Sys.time()

library(rgdal)
library(raster)
library(landsat)
library(MLmetrics)
library(hydroGOF)
library(ggplot2)
library(dplyr)
library(lubridate)

# ---- DEFINIR DIRECTORIOS ----
setwd('D:/MEMORIA_LRI') # directorio
dir.in = paste(getwd(),"/Materiales/espacializacion_y_temporalizacion_IAF_y_fPAR", sep = "") # directorio entrada con imgs. Sentinel 2 .SAFE
dir.out = paste(getwd(),"/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR", sep = "") # directorio salida
ROI = readOGR('D:/MEMORIA_LRI/Materiales/seleccion_pixeles_terreno/shapes/ROI.shp') # Area de estudio
sonda = readOGR('D:/MEMORIA_LRI/Materiales/IVs_y_modelos_regresion_IAF_fPAR/shapes/sonda.shp') # Pixel sonda
## espacializacion IAF y FPAR

## LAI = 4.38 * SAVI - 0.77
## FPAR = 1.41 * NDVI - 0.35


setwd("D:/MEMORIA_LRI/Materiales/espacializacion_y_temporalizacion_IAF_y_fPAR/temporalizacion_variables/imgs/S2_L2a") 
imgs_L2A = grep('.SAFE',list.dirs(path = getwd(), full.names = T,recursive = FALSE),ignore.case=TRUE, value=TRUE)  

nam2 = substr(x = imgs_L2A,start = nchar(getwd())+2 ,stop = nchar(getwd())+9)
# day = substr(x = nam2,start = 7, stop = 8)
# month = substr(x = nam2,start = 5, stop = 6)
# year = substr(x = nam2,start = 1, stop = 4)
# date = paste0(day, '-', month, '-' , year)
# write.csv(date,paste0('D:/MEMORIA_LRI/Materiales/meteo/fechas_imgs.csv'))



lai_temp = list()
fpar_temp = list()
#WDVI_temp = list()
#TSAVI_temp = list()
SAVI_temp = list()
NDVI_temp = list()

for (i in 1:length(imgs_L2A)) {
  
  wd0 = (getwd())
  setwd(imgs_L2A[i])
  setwd(paste0(getwd(),'/GRANULE'))
  setwd(paste0(getwd(),'/',list.files()[1]))
  setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
  
  red = crop(raster(list.files()[4]),ROI) # banda 'red' Sentinel 2
  nir = crop(raster(list.files()[10]),ROI) # banda '8A' Sentinel 2
  
  # a = data.frame(red = red@data@values/10000, nir = nir@data@values/10000)
  # soil_line = BSL(band3 = a$red, band4 = a$nir, method = "quantile", ulimit = 0.99, llimit = 0.005 )
  # alpha = soil_line$BSL[2]
  # 
  # wdvi = (nir/10000 - (alpha * red/10000))
  # 
  # beta = soil_line$BSL[1]
  # tsavi = (alpha * (nir/10000 - alpha*red/10000 - beta)/(red/10000 + alpha*nir/10000 - alpha*beta))
  # 
  # LAI = 7.43 * wdvi - 0.78
  # FPAR = 1.42 * tsavi - 0.34
  
  L = 0.5
  savi = ((nir/10000 - red/10000)/((nir/10000 + red/10000 + L))) * (1 + L)
  
  ndvi = ((nir/10000 - red/10000)/((nir/10000 + red/10000)))
  
  LAI = 4.38 * savi - 0.77
  FPAR = 1.41 * ndvi - 0.35
  
  lai_temp[[i]] = LAI
  names(lai_temp)[i] = nam2[i]
  fpar_temp[[i]] = FPAR
  names(fpar_temp)[i] = nam2[i]
  
  #WDVI_temp[[i]] = wdvi
  #names(WDVI_temp)[i] = nam2[i]
  #TSAVI_temp[[i]] = tsavi
  #names(TSAVI_temp)[i] = nam2[i]
  SAVI_temp[[i]] = savi
  names(SAVI_temp)[i] = nam2[i]
  NDVI_temp[[i]] = ndvi
  names(NDVI_temp)[i] = nam2[i]
  
  setwd(wd0)
  print(i)
}

setwd(paste0(dir.out, '/tiff/LAI'))
lai_stack = stack(lai_temp)
lai_stack = stack(mask(crop(lai_stack,ROI),ROI))
#writeRaster(lai_stack, filename=paste0(names(lai_stack),'_lai'), bylayer=TRUE,format="GTiff", overwrite = TRUE)


setwd(paste0(dir.out, '/tiff/fPAR'))
fpar_stack = stack(fpar_temp)
fpar_stack = stack(mask(crop(fpar_stack,ROI),ROI))
#writeRaster(fpar_stack, filename=paste0(names(fpar_stack),'_fPAR'), bylayer=TRUE,format="GTiff", overwrite = TRUE)

setwd(paste0(dir.out, '/tiff/WDVI'))
wdvi_stack = stack(WDVI_temp)
wdvi_stack = stack(mask(crop(wdvi_stack,ROI),ROI))
#writeRaster(wdvi_stack, filename=paste0(names(wdvi_stack),'_WDVI'), bylayer=TRUE,format="GTiff", overwrite = TRUE)


setwd(paste0(dir.out, '/tiff/TSAVI'))
tsavi_stack = stack(TSAVI_temp)
tsavi_stack = stack(mask(crop(tsavi_stack,ROI),ROI))
#writeRaster(tsavi_stack, filename=paste0(names(tsavi_stack),'_TSAVI'), bylayer=TRUE,format="GTiff", overwrite = TRUE)


setwd(paste0(dir.out, '/tiff/SAVI'))
savi_stack = stack(SAVI_temp)
savi_stack_ROI = stack(mask(crop(savi_stack,ROI),ROI))
#writeRaster(savi_stack_ROI, filename=paste0(names(savi_stack_ROI),'_SAVI'), bylayer=TRUE,format="GTiff", overwrite = TRUE)


setwd(paste0(dir.out, '/tiff/NDVI'))
ndvi_stack = stack(NDVI_temp)
ndvi_stack_ROI = stack(mask(crop(ndvi_stack,ROI),ROI))
#writeRaster(ndvi_stack_ROI, filename=paste0(names(ndvi_stack_ROI),'_NDVI'), bylayer=TRUE,format="GTiff", overwrite = TRUE)

# ---- EXTRACT SONDA ----

lai_sonda = extract(lai_stack, sonda)
fpar_sonda = extract(fpar_stack, sonda)

wdvi_sonda = extract(wdvi_stack, sonda)
tsavi_sonda = extract(tsavi_stack, sonda)


day = substr(x = nam2,start = 7 ,stop = 8)
month = substr(x = nam2,start = 5 ,stop = 6) 
year = substr(x = nam2,start = 1 ,stop = 4)
fecha = paste0(year,'-',month,'-',day)
fecha = as.Date(fecha, format = '%Y-%m-%d')
sonda_df = data.frame(fecha = fecha,
                      lai = c(lai_sonda),
                      fpar = c(fpar_sonda))#,
                      #wdvi = c(wdvi_sonda),
                      #tsavi = c(tsavi_sonda))


#write.csv(sonda_df,'D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/SONDA/sonda_df.csv')

# ---- INTERPOLACION DE VARIABLES ENTRE FECHAS ----

ApproxLAI = approxfun(x = sonda_df$fecha, y = sonda_df$lai)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitLAI = ApproxLAI(Dates)

ApproxFPAR = approxfun(x = sonda_df$fecha, y = sonda_df$fpar)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitFPAR = ApproxFPAR(Dates)

sonda_df_temporal = data.frame(fecha = Dates, 
                               lai = LinearFitLAI, 
                               fpar = LinearFitFPAR)

#write.csv(sonda_df_temporal,'D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/SONDA/sonda_df_temporal.csv')

# ---- LAI & FPAR Biophysical ----

setwd('D:/MEMORIA_LRI/Materiales/espacializacion_y_temporalizacion_IAF_y_fPAR/temporalizacion_variables/imgs/biophyscal_SNAP') #Directorio  con imgs Sentinel 2 L2A .SAFE
imgs_biophysical = grep('.tif',list.files(path = getwd(),full.names = T,recursive = FALSE),ignore.case=TRUE, value=TRUE)  

nam2 = substr(x = imgs_biophysical,start = nchar(getwd())+2 ,stop = nchar(getwd())+9)

lai_S2_temp = list()
lai_flags_S2_temp = list()
fpar_S2_temp = list()
fpar_flags_S2_temp = list()
fcover_S2_temp = list()
fcover_flags_S2_temp = list()
extent=extent(ROI)

for (i in 1:length(imgs_biophysical)) {
  
  b=brick(imgs_biophysical[i])
  b=crop(b,extent)
  
  lai_S2_temp[[i]] = b[[1]]
  names(lai_S2_temp)[i] = nam2[i]
  
  lai_flags_S2_temp[[i]] = b[[2]]
  names(lai_flags_S2_temp)[i] = nam2[i]
  
  fpar_S2_temp[[i]] = b[[3]]
  names(fpar_S2_temp)[i] = nam2[i]
  
  fpar_flags_S2_temp[[i]] = b[[4]]
  names(fpar_flags_S2_temp)[i] = nam2[i]
  
  fcover_S2_temp[[i]] = b[[5]]
  names(fcover_S2_temp)[i] = nam2[i]
  
  fcover_flags_S2_temp[[i]] = b[[6]]
  names(fcover_flags_S2_temp)[i] = nam2[i]
  
  print(i)
}

lai_S2_stack = stack(lai_S2_temp)
lai_flags_S2_stack = stack(lai_flags_S2_temp)
lai_S2_stack[lai_flags_S2_stack != 0] <- NA
lai_S2_stack = approxNA(lai_S2_stack)

fpar_S2_stack = stack(fpar_S2_temp)
fpar_flags_S2_stack = stack(fpar_flags_S2_temp)
fpar_S2_stack[fpar_flags_S2_stack != 0] <- NA
fpar_S2_stack = approxNA(fpar_S2_stack)

fcover_S2_stack = stack(fcover_S2_temp)
fcover_flags_S2_stack = stack(fcover_flags_S2_temp)
fcover_S2_stack[fcover_flags_S2_stack != 0] <- NA
fcover_S2_stack = approxNA(fcover_S2_stack)

# ---- EXTRACT VARIABLES BIOPHYSICAL SNAP ----

lai_S2_sonda = extract(lai_S2_stack, sonda)
lai_flags_S2_sonda = extract(lai_flags_S2_stack, sonda)

fpar_S2_sonda = extract(fpar_S2_stack, sonda)
fpar_flags_S2_sonda = extract(fpar_flags_S2_stack, sonda)

fcover_S2_sonda = extract(fcover_S2_stack, sonda)
fcover_flags_S2_sonda = extract(fcover_flags_S2_stack, sonda)

day = substr(x = nam2,start = 7 ,stop = 8)
month = substr(x = nam2,start = 5 ,stop = 6) 
year = substr(x = nam2,start = 1 ,stop = 4)
fecha = paste0(year,'-',month,'-',day)
fecha = as.Date(fecha, format = '%Y-%m-%d')

sonda_S2_df = data.frame(fecha = fecha,
                         lai_s2 = c(lai_S2_sonda),
                         fpar_s2 = c(fpar_S2_sonda),
                         fcover_s2 = c(fcover_S2_sonda))

#write.csv(sonda_S2_df,'D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/SONDA/sonda_S2_df.csv')

ApproxLAI_S2 = approxfun(x = sonda_S2_df$fecha, y = sonda_S2_df$lai_s2)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitLAI_S2 = ApproxLAI_S2(Dates)

ApproxFPAR_S2 = approxfun(x = sonda_S2_df$fecha, y = sonda_S2_df$fpar_s2)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitFPAR_S2 = ApproxFPAR_S2(Dates)

ApproxFCOVER_S2 = approxfun(x = sonda_S2_df$fecha, y = sonda_S2_df$fcover_s2)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitFCOVER_S2 = ApproxFCOVER_S2(Dates)

sonda_S2_df_temporal = data.frame(fecha = Dates,
                                  lai_s2 = LinearFitLAI_S2,
                                  fpar_s2 = LinearFitFPAR_S2,
                                  fcover_s2 = LinearFitFCOVER_S2)

#write.csv(sonda_S2_df_temporal,'D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/SONDA/sonda_S2_df_temporal.csv')

#### GUARDAR TIFF variables biophysical

setwd(paste0(dir.out, '/tiff/BIOPHYSICAL/LAI'))
lai_biophysical_stack = mask(crop(lai_S2_stack,ROI),ROI)
#writeRaster(lai_biophysical_stack, paste0(filename=names(lai_biophysical_stack),'_LAI_SNAP'), bylayer=TRUE,format="GTiff", overwrite = TRUE)

setwd(paste0(dir.out, '/tiff/BIOPHYSICAL/fPAR'))
fpar_biophyscal_stack = mask(crop(fpar_S2_stack,ROI),ROI)
#writeRaster(fpar_biophyscal_stack, paste0(filename=names(fpar_biophyscal_stack),'_fPAR_SNAP'), bylayer=TRUE,format="GTiff", overwrite = TRUE)

setwd(paste0(dir.out, '/tiff/BIOPHYSICAL/FCOVER'))
fcover_biophysical_stack = mask(crop(fcover_S2_stack,ROI),ROI)
#writeRaster(fcover_biophysical_stack, filename=paste0(names(fcover_biophysical_stack),'_FCOVER_SNAP'), bylayer=TRUE,format="GTiff", overwrite = TRUE)

tf = Sys.time()
print(tf - t0)
#### FIN ####