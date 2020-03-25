####################################################################################################
## CODIGOS MEMORIA 
## 8 - ZONAS MANEJO - KMEANS
## Lucas Rivero Iribarne
####################################################################################################


#####################################################################################################
## Este codigo permite generar zonas de manejo a partir del algoritmo de clasificacion no supervisada
## "Kmeans" a partir de los datos de Kc
#####################################################################################################

rm(list = ls()) # limpia "Global Enviroment"
graphics.off() # limpia graficos

library(RStoolbox)
library(raster)

# ---- DEFINIR DIRECTORIOS ----
setwd('D:/MEMORIA_LRI') # directorio
dir.in = 'D:/MEMORIA_LRI/Resultados/Kx/CI_SAVI/Kc/tiff' # directorio entrada con imgs. Sentinel 2 .SAFE
dir.out = 'D:/MEMORIA_LRI/Resultados/zonas_manejo' # directorio salida


Kc_TS = stack(list.files(dir.in, pattern='\\.tif$', full.names = T))


input <- Kc_TS

## Run classification
set.seed(25)
unC <- unsuperClass(input, nSamples = 100, nClasses = 3, nStarts = 5)
unC
plot(unC)
## Plots
colors <- rainbow(5)
plot(unC$map, col = colors, legend = FALSE, axes = FALSE, box = FALSE)
legend(1,1, legend = paste0("C",1:5), fill = colors,
       title = "Classes", horiz = TRUE,  bty = "n")

par(olpar) # reset par


## Moving majority

a<-focal(unC$map, w=matrix(1,3,3), fun=modal)    # 3x3 moving window
plot(unC$map)
plot(a)

writeRaster(a, "Zonas_manejo_Kc.tiff")

###########################################################################################################################
#### SERIE TIEMPO LAI Y FPAR

TS_df = data.frame(read.csv('D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/SONDA/sonda_df.csv'))
TS_df = TS_df[,-1]
TS_df$fecha = as.Date(TS_df$fecha) 

TS_melt = reshape2::melt(TS_df, id.vars = c('fecha'), measure.vars = c('lai', 'fpar'))

g1 = ggplot(TS_melt, aes(x = fecha, y = value, color = variable)) +
  geom_line(size = 1.0) + 
  geom_point(size = 3.0, shape = 18, col = 'black') +
  geom_abline(intercept = 0, slope = 0) +
  ggtitle('Serie de tiempo IAF y FPAR') +
  scale_y_continuous(limits = c(0, 2.6), breaks = seq(0, 2.6, by = 0.2)) +
  xlab("Fecha") + ylab("IAF (m2/m2) y FPAR (%)") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3)



TS_df = data.frame(read.csv('D:/MEMORIA_LRI/Resultados/Kx/CI_SAVI/SONDA/sonda_Kx_CI_SAVI.csv'))
TS_df = TS_df[,-1]
TS_df$fecha = as.Date(TS_df$fecha) 

TS_melt = reshape2::melt(TS_df, id.vars = c('fecha'), measure.vars = c('Kcb', 'Ke', 'Kc', 'Kcmax'))

g2 = ggplot(TS_melt, aes(x = fecha, y = value, color = variable)) +
  geom_line(size = 1.0) + 
  geom_point(size = 3.0, shape = 18, col = 'black') +
  geom_abline(intercept = 0, slope = 0) +
  ggtitle('Serie de tiempo Kc, Kcb, Ke, Kcmax') +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, by = 0.1)) +
  xlab("Fecha") + ylab("Kc, Kcb, Ke, Kcmax (-)") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3)

######

ROI_ALTO = ROI[,1][c(1,3,4),]
ROI_MEDIO = ROI[,1][c(7,9,10),]
ROI_BAJO = ROI[,1][c(2,5,6,8),]

Kc_Alto = crop(mask(Kc_TS, ROI_ALTO), ROI_ALTO)
Kc_Medio = crop(mask(Kc_TS, ROI_MEDIO), ROI_MEDIO)
Kc_Bajo = crop(mask(Kc_TS, ROI_BAJO), ROI_BAJO)

Kc_ALTO = mean(cellStats(Kc_Alto, stat = 'mean', na.rm = T))
Kc_MEDIO = mean(cellStats(Kc_Medio, stat = 'mean', na.rm = T))
Kc_BAJO = mean(cellStats(Kc_Bajo, stat = 'mean', na.rm = T))


ROI_1 = ROI[,1][c(1),]
ROI_2 = ROI[,1][c(2),]
ROI_3 = ROI[,1][c(3),]
ROI_4 = ROI[,1][c(4),]
ROI_5 = ROI[,1][c(5),]
ROI_6 = ROI[,1][c(6),]
ROI_7 = ROI[,1][c(7),]
ROI_8 = ROI[,1][c(8),]
ROI_9 = ROI[,1][c(9),]
ROI_10 = ROI[,1][c(10),]

Kc_1 = crop(mask(Kc_TS,ROI_1 ), ROI_1)
Kc_2 = crop(mask(Kc_TS,ROI_2 ), ROI_2)
Kc_3 = crop(mask(Kc_TS,ROI_3 ), ROI_3)
Kc_4 = crop(mask(Kc_TS,ROI_4 ), ROI_4)
Kc_5 = crop(mask(Kc_TS,ROI_5 ), ROI_5)
Kc_6 = crop(mask(Kc_TS,ROI_6 ), ROI_6)
Kc_7 = crop(mask(Kc_TS,ROI_7 ), ROI_7)
Kc_8 = crop(mask(Kc_TS,ROI_8 ), ROI_8)
Kc_9 = crop(mask(Kc_TS,ROI_9 ), ROI_9)
Kc_10 = crop(mask(Kc_TS,ROI_10 ), ROI_10)

Kc_1 = (cellStats(Kc_1, stat = 'mean', na.rm = T)) # ALTO
Kc_2 = (cellStats(Kc_2, stat = 'mean', na.rm = T)) # BAJO
Kc_3 = (cellStats(Kc_3, stat = 'mean', na.rm = T)) # ALTO
Kc_4 = (cellStats(Kc_4, stat = 'mean', na.rm = T)) # ALTO
Kc_5 = (cellStats(Kc_5, stat = 'mean', na.rm = T)) # BAJO
Kc_6 = (cellStats(Kc_6, stat = 'mean', na.rm = T)) # BAJO
Kc_7 = (cellStats(Kc_7, stat = 'mean', na.rm = T)) # MEDIO
Kc_8 = (cellStats(Kc_8, stat = 'mean', na.rm = T)) # BAJO
Kc_9 = (cellStats(Kc_9, stat = 'mean', na.rm = T)) # MEDIO
Kc_10 = (cellStats(Kc_10, stat = 'mean', na.rm = T)) # MEDIO


###

Kc_df = data.frame(fecha = fecha,
                   Kc1 = Kc_1, # ALTO
                   Kc2 = Kc_2,# BAJO
                   Kc3 = Kc_3, # ALTO
                   Kc4 = Kc_4, # ALTO
                   Kc5 = Kc_5,# BAJO
                   Kc6 = Kc_6,# BAJO
                   Kc7 = Kc_7, # MEDIO
                   Kc8 = Kc_8,# BAJO
                   Kc9 = Kc_9, # MEDIO
                   Kc10 = Kc_10) # MEDIO


melt = reshape2::melt(Kc_df, id.vars = c('fecha'), measure.vars = c("Kc1","Kc2","Kc3","Kc4","Kc5","Kc6","Kc7","Kc8","Kc9","Kc10"))

g = ggplot(melt, aes(x = fecha, y = value, color = variable)) +
  geom_point(size = 2, shape = 16) +
  geom_line(size = 1.0) + 
  ggtitle("Kc") +
  xlab("Fecha") + ylab("Kc") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) 

### Level plot LAI & FPAR

# LAI

setwd('D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/tiff/LAI')
LAI_TS = stack(list.files(dir.in, pattern='\\.tif$', full.names = T))


cols <- colorRampPalette(brewer.pal(9,"RdYlGn"))
levelplot(LAI_TS[[c(1:20)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal IAF",
          names.attr=as.character(fecha[1:20]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(LAI_TS[[c(21:40)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal IAF",
          names.attr=as.character(fecha[21:40]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(LAI_TS[[c(41:52)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal IAF",
          names.attr=as.character(fecha[41:52]),
          scales=list(draw=FALSE )) # remove axes labels & ticks



# FPAR

setwd('D:/MEMORIA_LRI/Resultados/espacializacion_y_temporalizacion_IAF_y_fPAR/tiff/fPAR')
FPAR_TS = stack(list.files(getwd(), pattern='\\.tif$', full.names = T))

cols <- colorRampPalette(brewer.pal(9,"RdYlGn"))
levelplot(FPAR_TS[[c(1:20)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal FPAR",
          names.attr=as.character(fecha[1:20]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(FPAR_TS[[c(21:40)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal FPAR",
          names.attr=as.character(fecha[21:40]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(FPAR_TS[[c(41:52)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal FPAR",
          names.attr=as.character(fecha[41:52]),
          scales=list(draw=FALSE )) # remove axes labels & ticks



### Figuras Ke, Kcb y Kc

# Kcb

setwd('D:/MEMORIA_LRI/Resultados/Kx/CI_SAVI/Kcb/tiff')
kcb_TS = stack(list.files(getwd(), pattern='\\.tif$', full.names = T))

cols <- colorRampPalette(brewer.pal(9,"RdYlGn"))
levelplot(kcb_TS[[c(1:20)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal kcb",
          names.attr=as.character(fecha[1:20]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(kcb_TS[[c(21:40)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal kcb",
          names.attr=as.character(fecha[21:40]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(kcb_TS[[c(41:52)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal kcb",
          names.attr=as.character(fecha[41:52]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

writeRaster(kcb_TS[[30]], filename = "Kcb_2018_12_15.tiff", format = 'GTiff', overwrite = T)


# Ke

setwd('D:/MEMORIA_LRI/Resultados/Kx/CI_SAVI/Ke/tiff')
ke_TS = stack(list.files(getwd(), pattern='\\.tif$', full.names = T))

levelplot(ke_TS[[c(1:20)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal Ke",
          names.attr=as.character(fecha[1:20]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(ke_TS[[c(21:40)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal Ke",
          names.attr=as.character(fecha[21:40]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(ke_TS[[c(41:52)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal Ke",
          names.attr=as.character(fecha[41:52]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

writeRaster(ke_TS[[30]], filename = "Ke_2018_12_15.tiff", format = 'GTiff', overwrite = T)


# Kc

setwd('D:/MEMORIA_LRI/Resultados/Kx/CI_SAVI/Kc/tiff')
kc_TS = stack(list.files(getwd(), pattern='\\.tif$', full.names = T))

levelplot(kc_TS[[c(1:20)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal Kc",
          names.attr=as.character(fecha[1:20]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(kc_TS[[c(21:40)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal Kc",
          names.attr=as.character(fecha[21:40]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

levelplot(kc_TS[[c(41:52)]],
          layout=c(5, 4), # create a 13x4 layout for the data
          col.regions=cols, # add a color ramp
          main="Serie temporal Kc",
          names.attr=as.character(fecha[41:52]),
          scales=list(draw=FALSE )) # remove axes labels & ticks

writeRaster(kc_TS[[30]], filename = "Kc_2018_12_15.tiff", format = 'GTiff', overwrite = T)

### TS Kx sonda

sonda_TS = read.csv('D:/MEMORIA_LRI/Resultados/Kx/CI_SAVI/SONDA/sonda_Kx_CI_SAVI.csv')
sonda_TS_df = data.frame(fecha = as.Date(sonda_TS$fecha),
                         kcb = sonda_TS$Kcb,
                         ke = sonda_TS$Ke,
                         kc = sonda_TS$Kc)


TS_melt = reshape2::melt(sonda_TS_df, id.vars = c('fecha'), measure.vars = c('kcb', 'ke', 'kc'))

g1 = ggplot(TS_melt, aes(x = fecha, y = value, color = variable)) +
  geom_line(size = 1.0) + 
  geom_point(size = 3.0, shape = 18, col = 'black') +
  geom_abline(intercept = 0, slope = 0) +
  ggtitle('Serie de tiempo kcb, ke y kc') +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  xlab("Fecha") + ylab("kcb, ke y kc [-]") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3)

