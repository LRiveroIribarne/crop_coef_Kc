#######################################################################################
## CODIGOS MEMORIA 
## 2 - IVs y modelos regresion IAF y fPAR
## Lucas Rivero Iribarne
#######################################################################################


################################################################################################
## Este codigo calcula varios indices de vegetacion y genera modelos de regresion,
## tanto lineales como logaritmicos entre los IV y variables medidas en terreno,
## luego lleva a cabo un bootstrapping y genera "beanplots" con coeficientes (R2, RMSE y PBIAS),
## de la siguiente forma:
##
################################################################################################

rm(list = ls()) # limpia "Global Enviroment"
graphics.off() # limpia graficos

library(caret)
library(rgdal)
library(raster)

# ---- DEFINIR DIRECTORIOS ----
setwd('D:/MEMORIA_LRI') # directorio
dir.in = paste(getwd(),"/Materiales/IVs_y_modelos_regresion_IAF_fPAR", sep = "") # directorio entrada con imgs. Sentinel 2 .SAFE
dir.out = paste(getwd(),"/Resultados/IVs_y_modelos_regresion_IAF_fPAR", sep = "") # directorio salida


coords = readOGR(paste0(dir.in,'/shapes'),'hass') # Puntos terreno
sonda = readOGR(paste0(dir.in,'/shapes'),'sonda') # Punto sonda
GAP_LAI_CI_110 = readOGR(paste0(dir.in,'/shapes'),'GAP_LAI')
FPAR_CI_110 = readOGR(paste0(dir.in,'/shapes'),'FPAR')

bio_febrero = mask(crop(stack('D:/MEMORIA_LRI/Materiales/IVs_y_modelos_regresion_IAF_fPAR/biophysical/S2B_MSIL2A_20190223T143749_N0211_R096_T19HBB_20190223T202555.tif'),ROI),ROI)
bio_marzo =  mask(crop(stack('D:/MEMORIA_LRI/Materiales/IVs_y_modelos_regresion_IAF_fPAR/biophysical/S2B_MSIL2A_20190325T143749_N0211_R096_T19HBB_20190325T185224.tif'),ROI),ROI)

LAI_S2_feb = bio_febrero$S2B_MSIL2A_20190223T143749_N0211_R096_T19HBB_20190223T202555.1
LAI_S2_mar = bio_marzo$S2B_MSIL2A_20190325T143749_N0211_R096_T19HBB_20190325T185224.1

FPAR_S2_feb = bio_febrero$S2B_MSIL2A_20190223T143749_N0211_R096_T19HBB_20190223T202555.3
FPAR_S2_mar = bio_marzo$S2B_MSIL2A_20190325T143749_N0211_R096_T19HBB_20190325T185224.3

coords_feb = coords[c(1:17),]
coords_abr = coords[c(18:34),]

# INDICES ESPECTRALES ----
{
  # NDVI (Rouse, 1974) ----
  ### (NIR-RED)/(NIR+RED)
  ### (B8A - B4)/(B8A + B4)
  
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
  }
  NDVI_stack = stack(NDVI_temp)
  
  # GNDVI (Gitelson, 1996) ----
  ### (NIR-GREEN)/(NIR+GREEN)
  ### (B8A - B3)/(B8A + B3)
  
  GNDVI_temp = list() # lista vacia para llenar con indices espectrales
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    green = mask(crop(raster(list.files()[3]),ROI),ROI) # banda 'green' Sentinel 2
    nir = mask(crop(raster(list.files()[10]),ROI),ROI) # banda '8A' Sentinel 2
    
    gndvi = ((nir/10000 - green/10000)/(nir/10000 + green/10000))
    
    GNDVI_temp[[i]] = gndvi
    names(GNDVI_temp)[i] = fechas[i]
  }
  GNDVI_stack = stack(GNDVI_temp)
  
  # WDVI (Clevers, 1989) ----
  ### (NIR - ((NIRsoil/REDsoil) * RED))
  ### (B8A - (B8Asoil/B4soil) * B4)
  
  WDVI_temp = list() # lista vacia para llenar con indices espectrales
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    red = crop(raster(list.files()[4]),ROI) # banda 'red' Sentinel 2
    nir = crop(raster(list.files()[10]),ROI) # banda '8A' Sentinel 2
    
    a = data.frame(red = red@data@values/10000, nir = nir@data@values/10000)
    soil_line = BSL(band3 = a$red, band4 = a$nir, method = "quantile", ulimit = 0.99, llimit = 0.005 )
    alpha = soil_line$BSL[2]
    wdvi = (nir/10000 - (alpha * red/10000))
    
    WDVI_temp[[i]] = wdvi
    names(WDVI_temp)[i] = fechas[i]
  }
  WDVI_stack = stack(WDVI_temp)
  WDVI_stack = stack(mask(crop(WDVI_stack,ROI),ROI))
  
  
  # SAVI (Huete, 1988) ----
  ### ((NIR - RED)/(NIR + RED + L))(1 + L)
  ### ((B8A - B4)/(B8A + B4 + L))(1 + L)
  
  SAVI_temp = list() # lista vacia para llenar con indices espectrales
  L = 0.5
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    red = mask(crop(raster(list.files()[4]),ROI),ROI) # banda 'green' Sentinel 2
    nir = mask(crop(raster(list.files()[10]),ROI),ROI) # banda '8A' Sentinel 2
    
    savi = ((nir/10000 - red/10000)/((nir/10000 + red/10000 + L))) * (1 + L)
    
    SAVI_temp[[i]] = savi
    names(SAVI_temp)[i] = fechas[i]
  }
  SAVI_stack = stack(SAVI_temp)
  
  # GSAVI (Li, 2010) ----
  ### ((NIR - GREEN)/(NIR + GREEN + L))(1 + L)
  ### ((B8A - B3)/(B8A + B3 + L))(1 + L)
  
  GSAVI_temp = list() # lista vacia para llenar con indices espectrales
  L = 0.5
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    green = mask(crop(raster(list.files()[3]),ROI),ROI) # banda 'green' Sentinel 2
    nir = mask(crop(raster(list.files()[10]),ROI),ROI) # banda '8A' Sentinel 2
    
    gsavi = ((nir/10000 - green/10000)/((nir/10000 + green/10000 + L))) * (1 + L)
    
    GSAVI_temp[[i]] = gsavi
    names(GSAVI_temp)[i] = fechas[i]
  }
  GSAVI_stack = stack(GSAVI_temp)
  
  # TSAVI (Baret, 1989) ----
  ### (a(NIR -aRED - b)/(RED + aNIR - ab))
  ### (a(B8A - aB4 - b)/(B4 + aB8A - ab))
  
  TSAVI_temp = list() # lista vacia para llenar con indices espectrales
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    red = crop(raster(list.files()[4]),ROI) # banda 'red' Sentinel 2
    nir = crop(raster(list.files()[10]),ROI) # banda '8A' Sentinel 2
    
    a = data.frame(red = red@data@values/10000, nir = nir@data@values/10000)
    soil_line = BSL(band3 = a$red, band4 = a$nir, method = "quantile", ulimit = 0.99, llimit = 0.005 )
    alpha = soil_line$BSL[2]
    beta = soil_line$BSL[1]
    tsavi = (alpha * (nir/10000 - alpha*red/10000 - beta)/(red/10000 + alpha*nir/10000 - alpha*beta))
    
    TSAVI_temp[[i]] = tsavi
    names(TSAVI_temp)[i] = fechas[i]
  }
  
  TSAVI_stack = stack(TSAVI_temp)
  TSAVI_stack = stack(mask(crop(TSAVI_stack,ROI),ROI))
  
  # PVI (Richardson and Wiegand, 1977) ----
  ### a * NIR - RED + b / (a^2 + 1)^(1/2)
  ### a * B8A - B4 + b / (a^2 + 1)^(1/2)
  
  PVI_temp = list() # lista vacia para llenar con indices espectrales
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    red = crop(raster(list.files()[4]),ROI) # banda 'green' Sentinel 2
    nir = crop(raster(list.files()[10]),ROI) # banda '8A' Sentinel 2
    
    a = data.frame(red = red@data@values/10000, nir = nir@data@values/10000)
    soil_line = BSL(band3 = a$red, band4 = a$nir, method = "quantile", ulimit = 0.99, llimit = 0.005 )
    alpha = soil_line$BSL[2]
    beta = soil_line$BSL[1]
    pvi = ((alpha * nir/10000 - red/10000 + beta)/((alpha^2 + 1)^(1/2)))
    
    PVI_temp[[i]] = pvi
    names(PVI_temp)[i] = fechas[i]
  }
  
  PVI_stack = stack(PVI_temp)
  PVI_stack = stack(mask(crop(PVI_stack,ROI),ROI))
  
  # NDWI (Gao, 1996) ----
  ### (NIR - SWIR)/(NIR + SWIR)
  ### (B8A - B11) / (B8A + B11)
  
  NDWI_temp = list() # lista vacia para llenar con indices espectrales
  for (i in 1:length(imgs_L2A)) {
    wd0 = (paste0(dir.in,'/imgs'))
    setwd(imgs_L2A[i])
    setwd(paste0(getwd(),'/GRANULE'))
    setwd(paste0(getwd(),'/',list.files()[1]))
    setwd(paste0(getwd(),'/','IMG_DATA/R20m'))
    
    swir = mask(crop(raster(list.files()[8]),ROI),ROI) # banda 'green' Sentinel 2
    nir = mask(crop(raster(list.files()[10]),ROI),ROI) # banda '8A' Sentinel 2
    
    ndwi = (nir/10000 - swir/10000) / (nir/10000 + swir/10000)
    
    NDWI_temp[[i]] = ndwi
    names(NDWI_temp)[i] = fechas[i]
  }
  NDWI_stack = stack(NDWI_temp)
}

# Extracts ----

NDVI_extract = c(extract(NDVI_stack[[1]],coords_feb),extract(NDVI_stack[[2]],coords_abr))
GNDVI_extract = c(extract(GNDVI_stack[[1]],coords_feb),extract(GNDVI_stack[[2]],coords_abr))
WDVI_extract = c(extract(WDVI_stack[[1]],coords_feb),extract(WDVI_stack[[2]],coords_abr))
SAVI_extract = c(extract(SAVI_stack[[1]],coords_feb),extract(SAVI_stack[[2]],coords_abr))
GSAVI_extract = c(extract(GSAVI_stack[[1]],coords_feb),extract(GSAVI_stack[[2]],coords_abr))
TSAVI_extract = c(extract(TSAVI_stack[[1]],coords_feb),extract(TSAVI_stack[[2]],coords_abr))
PVI_extract = c(extract(PVI_stack[[1]],coords_feb),extract(PVI_stack[[2]],coords_abr))
NDWI_extract = c(extract(NDWI_stack[[1]],coords_feb),extract(NDWI_stack[[2]],coords_abr))
LAI_S2_extract = c(extract(LAI_S2_feb,coords_feb),extract(LAI_S2_mar,coords_abr))
FPAR_S2_extract = c(extract(FPAR_S2_feb,coords_feb),extract(FPAR_S2_mar,coords_abr))
lai_CI_110 = GAP_LAI_CI_110$d__GAP_
fpar_CI_110 = FPAR_CI_110$d__FPAR

index_df = data.frame(lai_CI_110 = lai_CI_110,
                      fpar_CI_110 = fpar_CI_110,
                      LAI_S2 = LAI_S2_extract,
                      FPAR_S2 = FPAR_S2_extract,
                      NDVI = NDVI_extract,
                      GNDVI = GNDVI_extract,
                      WDVI = WDVI_extract,
                      SAVI = SAVI_extract,
                      GSAVI = GSAVI_extract,
                      TSAVI = TSAVI_extract,
                      PVI = PVI_extract,
                      NDWI = NDWI_extract)

#write.csv(index_df, paste0(dir.out, '/Index_df.csv'))

# ---- MODELOS REGRESION ----

setwd(paste0(dir.out,'/modelos_regresion'))

#pairs(index_df)
#cor(index_df)

samp_size <- floor(0.7 * nrow(index_df))
set.seed(123)
train_ind <- sample(seq_len(nrow(index_df)), size = samp_size)

train <- index_df[train_ind, ]
test <- index_df[-train_ind, ]

#setwd(paste0(dir.out, '_medios','/model_df'))
#write.csv(train, 'model_train_df.csv')
#write.csv(test, 'model_test_df.csv')


setwd('C:/SCRIPTS_MEMORIA/1/Resultados/modelos_regresion/LAI/lineal')
# GAP LAI linear
GAP_LAI_lineal_df = data.frame(r2 = 0, RMSE = 0, PBIAS = 0, MAPE = 0)
for (i in 3:12) {
  lm = lm(lai_CI_110~train[,i], data = train)
  r2 = summary(lm)$adj.r.squared
  obs = data.frame(test$lai_CI_110)
  pred = lm$coefficients[2]*test[i] + lm$coefficients[1]
  rmse = rmse(pred, obs)
  pbias = pbias(pred, obs)
  MAPE = MAPE(pred[[1]], obs[[1]])
  
  g1 = ggplot(train,aes(x=train[,i],y=lai_CI_110))+geom_point()+geom_smooth(method="lm", col = 'red', se = FALSE)+
    labs(title=paste0(colnames(train[i]), ' - GAP LAI '))+
    xlab(colnames(train[i]))+ylab("GAP LAI") +
    xlim(0,0.5) +
    annotate("text",x= 0.2, y=2.5,label=paste("GAP LAI =", round(lm$coefficients[2],2),"*", colnames(train[i]),round(lm$coefficients[1],2)))+
    annotate("text",x= 0.2, y=2.35,label=paste("r2 =", round(r2,3)))+
    annotate("text",x= 0.2, y=2.2,label=paste("RMSE =", round(rmse,3)))+
    annotate("text",x= 0.2, y=2.05,label=paste("PBIAS =",round(pbias,3),"%"))+
    annotate("text",x= 0.2, y=1.90,label=paste("MAPE =",round(MAPE,3)))+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
  
  
  df = data.frame(r2 = r2, RMSE = rmse, PBIAS = pbias, MAPE = MAPE)
  GAP_LAI_lineal_df[i-2,] = df 
  
  
  ggsave(paste0(colnames(train[i]),"_GAP_LAI.jpg"),plot=last_plot(),width=8,height=6,units="in",dpi=300)
  print(i)
}
write.csv(GAP_LAI_lineal_df, 'GAP_LAI_lineal_df.csv')


setwd('C:/SCRIPTS_MEMORIA/1/Resultados/modelos_regresion/FPAR/lineal')
#FPAR Lineal
rm(i)
FPAR_lineal_df = data.frame(r2 = 0, RMSE = 0, PBIAS = 0, MAPE = 0)
for (i in 3:12) {
  lm = lm(fpar_CI_110/100~train[,i], data = train)
  r2 = summary(lm)$adj.r.squared
  obs = data.frame(test$fpar_CI_110/100)
  pred = lm$coefficients[2]*test[i] + lm$coefficients[1]
  rmse = rmse(pred, obs)
  pbias = pbias(pred, obs)
  MAPE = MAPE(pred[[1]], obs[[1]])
  
  g1 = ggplot(train,aes(x=train[,i],y=fpar_CI_110/100))+geom_point()+geom_smooth(method="lm", col = 'red', se = FALSE)+
    labs(title=paste0(colnames(train[i]), ' - FPAR '))+
    xlab(colnames(train[i]))+ylab("FPAR") +
    xlim(0.2,1) + 
    annotate("text",x=0.485,y=1,label=paste("FPAR =", round(lm$coefficients[2],2),"*", colnames(train[i]),round(lm$coefficients[1],2)))+
    annotate("text",x=0.485,y=0.95,label=paste("r2 =", round(r2,3)))+
    annotate("text",x=0.485,y=0.9,label=paste("RMSE =", round(rmse,3)))+
    annotate("text",x=0.485,y=0.85,label=paste("PBIAS =",round(pbias,3),"%"))+
    annotate("text",x=0.485,y=0.8,label=paste("MAPE =",round(MAPE,3)))+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
  
  
  df = data.frame(r2 = r2, RMSE = rmse, PBIAS = pbias, MAPE = MAPE)
  FPAR_lineal_df[i-2,] = df
  
  ggsave(paste0(colnames(train[i]),"_FPAR.jpg"),plot=last_plot(),width=8,height=6,units="in",dpi=300)
  print(i)
}
write.csv(FPAR_lineal_df, 'FPAR_lineal_df.csv')

setwd('C:/SCRIPTS_MEMORIA/1/Resultados/modelos_regresion/LAI/log')
# GAP LAI log
GAP_LAI_log_df = data.frame(r2 = 0, RMSE = 0, PBIAS = 0, MAPE = 0)
for (k in 3:12) {
  logm = lm(lai_CI_110~log(train[,k]), data = train)
  r2 = summary(logm)$adj.r.squared
  obs = data.frame(test$lai_CI_110)
  pred = logm$coefficients[2]*log(test[k]) + logm$coefficients[1]
  rmse = rmse(pred, obs)
  pbias = pbias(pred, obs)
  MAPE = MAPE(pred[[1]], obs[[1]])
  
  g1 = ggplot(train,aes(x=train[,k],y=lai_CI_110))+geom_point()+geom_smooth(method="lm",formula = y ~ log(x) , col = 'red', se = FALSE)+
    labs(title=paste0(colnames(train[k]), '-GAP LAI'))+
    xlab(colnames(train[k]))+ylab("GAP LAI") +
    annotate("text",x=0.5,y=1.0,label=paste("GAP LAI =", round(logm$coefficients[2],2),"*",colnames(train[k]),'+',round(logm$coefficients[1],2)))+
    annotate("text",x=0.5,y=0.93,label=paste("r2 =", round(r2,2)))+
    annotate("text",x=0.5,y=0.86,label=paste("RMSE =", round(rmse,2)))+
    annotate("text",x=0.5,y=0.79,label=paste("PBIAS =",round(pbias,2),"%"))+
    annotate("text",x=0.5,y=0.72,label=paste("MAPE =",round(MAPE,2),"%"))+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
  
  df = data.frame(r2 = r2, RMSE = rmse, PBIAS = pbias, MAPE = MAPE)
  GAP_LAI_log_df[k-2,] = df
  
  ggsave(paste0(colnames(train[k]),"_GAP_LAI.jpg"),plot=last_plot(),width=8,height=6,units="in",dpi=300)
  print(k)
}
write.csv(GAP_LAI_log_df, 'GAP_LAI_log_df.csv')


setwd('C:/SCRIPTS_MEMORIA/1/Resultados/modelos_regresion/FPAR/log')
# FPAR exponencial log
FPAR_log_df = data.frame(r2 = 0, RMSE = 0, PBIAS = 0, MAPE = 0)
rm(k)
for (k in 3:12) {
  logm = lm(fpar_CI_110~log(train[,k]), data = train)
  r2 = summary(logm)$adj.r.squared
  obs = data.frame(test$fpar_CI_110)
  pred = logm$coefficients[2]*log(test[k]) + logm$coefficients[1]
  rmse = rmse(pred, obs)/100
  pbias = pbias(pred, obs)
  MAPE = MAPE(pred[[1]], obs[[1]])
  
  g1 = ggplot(train,aes(x=train[,k],y=fpar_CI_110))+geom_point()+geom_smooth(method="lm",formula = y ~ log(x) , col = 'red', se = FALSE)+
    labs(title=paste0(colnames(train[k]), '- FPAR'))+
    xlab(colnames(train[k]))+ylab("FPAR") +
    annotate("text",x=0.5,y=1.0,label=paste("FPAR =", round(logm$coefficients[2],2),"*",colnames(train[k]),'+',round(logm$coefficients[1],2)))+
    annotate("text",x=0.5,y=0.93,label=paste("r2 =", round(r2,2)))+
    annotate("text",x=0.5,y=0.86,label=paste("RMSE =", round(rmse,2)))+
    annotate("text",x=0.5,y=0.79,label=paste("PBIAS =",round(pbias,2),"%"))+
    annotate("text",x=0.5,y=0.72,label=paste("MAPE =",round(MAPE,2),"%"))+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
  
  df = data.frame(r2 = r2, RMSE = rmse, PBIAS = pbias, MAPE = MAPE)
  FPAR_log_df[k-2,] = df
  
  
  ggsave(paste0(colnames(train[k]),"_FPAR.jpg"),plot=last_plot(),width=8,height=6,units="in",dpi=300)
  print(k)
}
write.csv(FPAR_log_df, 'FPAR_log_df.csv')

#### FIN ####