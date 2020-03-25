###########################################################
## Comparacion Ke, Kcb y Kc HYDRUS 1D - Percepcion remota
###########################################################



################################################################################################
## ESTE CODIGO REQUIERE DE ENTRADAS:
#
# > OUTPUTs HYDRUS 1-D CON Ke, Kcb y Kc DERIVADOS DIRECTAMENTE DE HYDRUS 1-D
#
################################################################################################

######################################
## PLOTS SALIDAS HYDRUS 1D
######################################


################################################################################################
## Los datos de entrada tienen que cumplir:
#
# > Todos los archivos tienen que estar en la carpeta del directorio
#
# > Todos los archivos tienen que estar en formato .txt
#
# > Los nombres de los archivos deben de ser los siguienes:
#   - HYDRUS-1D_inputs.txt = inputs HYDRUS, de aqui salen los Kx satelitales
#   - Obs_Node.out.txt = Obs_Node salida HYDRUS, contenido de agua en cada punto de observacion
#   - sonda_&_pp_riego.txt = Contenido de agua sonda, precipitaciones y riegos en cm
#   - T_Level.out.txt = T_Level salida HYRUS, TR, TP, PERCO, EVAP
#
# > Esto ultimo es asi para que al llamar a los archivos entren en el orden correcto.
#   pueden tener otro nombre en la medida que su numero en el directorio sea el mismo
#
# > Obs_node hay que eliminarle las filas que no son info.
# > T_level hay que eliminarle las filas que no son info. y la fila en blanco bajo las unidades
# > Obs_node y T_level hay que "llevarlo a columnas" antes de guardarlo en .txt
################################################################################################

library(ggplot2)
library(ggpubr)
library(reshape)
library(hydroGOF)
library(gridExtra)
library(MLmetrics)
library(hydrusR)

t0 = Sys.time() # variable para medir tiempo de ejecucion codigo. NO RELEVANTE

d_cal = 365 # dias calibracion (sin hay gap de datos corresponde al ultimo dia con datos para calibrar)

setwd('D:/MEMORIA_LRI/HYDRUS_1D/salidas/HYDRUS_Kx')

files = list.files(pattern = '.txt')
HYDRUS_out = list.files(path = paste0(getwd(), '/HYDRUS_Kx'), pattern = '.out') 

Inputs_HYDRUS = read.delim(paste0(getwd(), '/', files[1]),header = TRUE, sep = '\t')
#Obs_Node = read.delim(paste0(getwd(), '/', files[2]),header = TRUE, sep = '\t')
sonda_pp_riego_ETo = read.delim(paste0(getwd(), '/', files[2]),header = TRUE, sep = '\t')
#T_level = read.delim(paste0(getwd(), '/', files[4]),header = TRUE, sep = '\t')
#T_level = T_level[- c(1),] # para quitar la fila de las unidades

Obs_Node = read.obs_node(project.path = 'D:/MEMORIA_LRI/HYDRUS_1D/salidas/HYDRUS_Kx/HYDRUS_Kx', out.file = "Obs_Node.out", obs.output = c('theta') ,obs.nodes = c(10,30,50))

header = read.table(paste0(getwd(), '/HYDRUS_Kx/', HYDRUS_out[9]), nrows = 1, header = FALSE, sep ='', stringsAsFactors = TRUE, skip = 5)
T_level = read.table(paste0(getwd(), '/HYDRUS_Kx/', HYDRUS_out[9]), header = FALSE, sep = '', skip = 8, nrows = 658)
colnames(T_level) = unlist(header)


day = substr(x = as.character(sonda_pp_riego_ETo$date),start = 1 ,stop = 2)
month = substr(x = as.character(sonda_pp_riego_ETo$date),start = 4 ,stop = 5) 
year = substr(x = as.character(sonda_pp_riego_ETo$date),start = 7 ,stop = 10)
fecha = paste0(year,'-',month,'-',day)
fecha = as.Date(fecha, format = '%Y-%m-%d')

hum = data.frame(fecha = fecha, # Data frame con contenidos de agua HYDRUS y sonda
                 hum10cm = Obs_Node$theta_10[1:length(fecha)],
                 hum30cm = Obs_Node$theta_30[1:length(fecha)],
                 hum50cm = Obs_Node$theta_50[1:length(fecha)],
                 sonda10cm = sonda_pp_riego_ETo$SondaHum_10cm...,
                 sonda30cm = sonda_pp_riego_ETo$SondaHum_30cm...,
                 sonda50cm = sonda_pp_riego_ETo$SondaHum_50cm...)

sum_EVAP = data.frame( sum_evap = as.numeric(paste(T_level$`sum(Evap)`[1:length(fecha)])))
EVAP = data.frame(EVAP = 1:length(sum_EVAP$sum_evap))

for (i in 1:length(sum_EVAP$sum_evap)) {# ciclo para calcular EVAPORACION DIARIA ya que
  if (i == 1) {                         # HYDRUS 1 D entrega la sumatoria de la EVAP
    EVAP$EVAP[i] = sum_EVAP$sum_evap[i]
  } else if (i != 1) {
    EVAP$EVAP[i] = sum_EVAP$sum_evap[i] - sum_EVAP$sum_evap[i - 1]
  }
  if (i == length(sum_EVAP$sum_evap)) {
    print(paste0('Calculando EVAP. --- Dia: ',i))
    print('EVAP. Listo!')
  } else {print(paste0('Calculando EVAP. --- Dia: ',i))}
}

TS = data.frame( fecha = fecha, # Data frame con variables HYDRUS 1D
                 TP_cm.dia = as.numeric(paste(T_level$rRoot))[1:length(fecha)], # tranpiracion real
                 TR_cm.dia = as.numeric(paste(T_level$vRoot))[1:length(fecha)], # transpiracion potencial
                 TR.TP.ratio = as.numeric(paste(T_level$vRoot))[1:length(fecha)]/as.numeric(paste(T_level$rRoot))[1:658], # TR/TP
                 PERCO = as.numeric(paste(T_level$vBot))[1:length(fecha)], # Percolacion
                 EVAP =  EVAP$EVAP, # Evaporacion
                 sum_EVAP = sum_EVAP$sum_evap,
                 Ke = Inputs_HYDRUS[,5][1:length(fecha)], # coeficiente de evaporacion desde el suelo
                 Kcb = Inputs_HYDRUS[,6][1:length(fecha)], # coeficiente basal de cultivo
                 Kc = Inputs_HYDRUS[,7][1:length(fecha)]) # Ke + Kcb

TS$ETo_cm = sonda_pp_riego_ETo$Eto.cm.# Evapotranspiracion de referencia
TS$Ke_HYDRUS = TS$EVAP/TS$ETo_cm # Ke HYDRUS
TS$Kcb_HYDRUS = TS$TR_cm.dia/TS$ETo_cm # Kcb HYDRUS
TS$Kcb_Pot_HYDRUS = TS$TP_cm.dia/TS$ETo_cm  # Kcb Potencial HYDRUS
TS$Kc_HYDRUS = TS$Ke_HYDRUS + TS$Kcb_HYDRUS# Kc HYDRUS
TS$Kc_Pot_HYDRUS = TS$Ke_HYDRUS + TS$Kcb_Pot_HYDRUS # Kc Potencial HYDRUS  

pp_df = data.frame(fecha = fecha,
                   pp = sonda_pp_riego_ETo$pp.cm.,
                   riego = sonda_pp_riego_ETo$riego.cm.,
                   pp_total = sonda_pp_riego_ETo$ppTotal.cm.)


ETc_df = data.frame(fecha = fecha,
                    ETo = sonda_pp_riego_ETo$Eto.cm.,
                    Kc_satelital = TS$Kc,
                    Kc_HYDRUS = TS$Kc_HYDRUS,
                    ETc_satelital = TS$Kc * sonda_pp_riego_ETo$Eto.cm.,
                    ETc_HYDRUS = TS$Kc_HYDRUS * sonda_pp_riego_ETo$Eto.cm.)

hum_sonda = data.frame(fecha = fecha, # Data frame con contenidos de agua sonda para water storage
                       sonda10cm = sonda_pp_riego_ETo$SondaHum_10cm...,
                       sonda30cm = sonda_pp_riego_ETo$SondaHum_30cm...,
                       sonda50cm = sonda_pp_riego_ETo$SondaHum_50cm...)

hum_sonda_transpuesta = t(hum_sonda)

WS_sonda = c()

for (i in 1:length(fecha)) {
  
  if (is.na(hum_sonda_transpuesta[2,i])) {
    
    h_h20 = NA
    WS_sonda = c(WS_sonda, h_h20)
    
  } else {
    
    humedad = data.frame(prof = c(0,10,30,50,70), hum = c(hum = as.numeric(hum_sonda_transpuesta[2,i]), as.numeric(hum_sonda_transpuesta[2:4,i]), as.numeric(hum_sonda_transpuesta[4,i])))
    ApproxHum = approxfun(x = humedad$prof, y = humedad$hum)
    Prof = seq(0, 70, by = 0.7)
    LinearFitHum = ApproxHum(Prof)
    h_h20 = sum(LinearFitHum * 0.715)
    WS_sonda = c(WS_sonda, h_h20)
    
  }
  if (i == length(fecha)) {
    print(paste0('Calculando WS_sonda --- Dia: ',i))
    print('WS_sonda Listo!')
  } else {print(paste0('Calculando WS_sonda --- Dia: ',i))}
}


WS_df = data.frame(fecha = fecha, # Water storage (cm)
                   WS_sonda = WS_sonda,
                   WS_HYDRUS = as.numeric(paste(T_level$Volume))[1:length(fecha)])

#### Data Frame fechas calibracion ####

hum_cal = hum[c(1:d_cal),] # humedad calibracion 
TS_cal = TS[c(1:d_cal),] # TS calibracion
ETc_df_cal = ETc_df[c(1:d_cal),] # ETc calibracion
WS_df_cal = WS_df[c(1:d_cal),] # Water Storage calibracion

#### Data Frame fechas sin calibracion ####

hum_val = hum[c((d_cal + 1):nrow(hum)),] # humedad Validacion 
TS_val = TS[c((d_cal + 1):nrow(TS)),] # TS validacion
ETc_df_val = ETc_df[c((d_cal + 1):nrow(ETc_df)),] # ETc validacion
WS_df_val = WS_df[c((d_cal + 1):nrow(WS_df)),] # Water Storage Validacion

#### PLOTS ####

#### Humedad sonda e HYDRUS en el tiempo ####

pp_df2 = data.frame(fecha = pp_df$fecha, pp = pp_df$pp/25, riego = pp_df$riego/25) # se multiplica para poder plotear junto a otro data frame con ggplot
pp_melt = reshape2::melt(pp_df2, id.vars = c('fecha'), measure.vars = c('pp', 'riego'))


hum_melt10 = reshape2::melt(hum, id.vars = c('fecha'), measure.vars = c("hum10cm","sonda10cm"))

g1 = ggplot(NULL, aes(x = fecha, y = value, color = variable)) +
  geom_bar(data = pp_melt[which(pp_melt$value>0),], stat = 'identity',position="stack") +
  geom_line(data = hum_melt10, size = 0.4) + 
  scale_colour_manual(values = c("red", "yellow", "green", 'blue')) +
  ggtitle("Humedad 10 cm") +
  xlab("Fecha") + ylab("Contenido agua (%)") +
  scale_y_continuous(sec.axis = sec_axis(~ .*25, name ="Precipitaciones y riegos (cm)")) + 
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.4) +
  geom_vline(xintercept = hum_melt10$fecha[d_cal], linetype="dashed", 
             color = "red", size=1)+
  scale_fill_manual(values = c("red", "yellow", "green", 'blue'))

hum_melt30 = reshape2::melt(hum, id.vars = c('fecha'), measure.vars = c("hum30cm","sonda30cm"))

g2 = ggplot(NULL, aes(x = fecha, y = value, color = variable)) +
  geom_bar(data = pp_melt[which(pp_melt$value>0),], stat = 'identity',position="stack") +
  geom_line(data = hum_melt30, size = 0.4) + 
  scale_colour_manual(values = c("red", "yellow", "green", 'blue')) +
  ggtitle("Humedad 30 cm") +
  xlab("Fecha") + ylab("Contenido agua (%)") +
  scale_y_continuous(sec.axis = sec_axis(~ .*25, name ="Precipitaciones y riegos (cm)")) + 
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.4) +
  geom_vline(xintercept = hum_melt30$fecha[d_cal], linetype="dashed", 
             color = "red", size=1)+
  scale_fill_manual(values = c("red", "yellow", "green", 'blue'))



hum_melt50 = reshape2::melt(hum, id.vars = c('fecha'), measure.vars = c("hum50cm","sonda50cm"))

g3 = ggplot(NULL, aes(x = fecha, y = value, color = variable)) +
  geom_bar(data = pp_melt[which(pp_melt$value>0),], stat = 'identity', position="stack") +
  geom_line(data = hum_melt50, size = 0.4) + 
  scale_colour_manual(values = c("red", "yellow", "green", 'blue')) +
  ggtitle("Humedad 50 cm") +
  xlab("Fecha") + ylab("Contenido agua (%)") +
  scale_y_continuous(sec.axis = sec_axis(~ .*25, name ="Precipitaciones y riegos (cm)")) + 
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.4) +
  geom_vline(xintercept = hum_melt50$fecha[d_cal], linetype="dashed", 
             color = "red", size=1)+
  scale_fill_manual(values = c("red", "yellow", "green", 'blue'))


# hum_melt80 = reshape2::melt(hum, id.vars = c('fecha'), measure.vars = c("hum80cm","sonda80cm"))
# 
# g4 = ggplot(hum_melt80, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   ggtitle("Humedad 80 cm") +
#   xlab("Fecha") + ylab("Contenido agua (%)") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = hum_melt80$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)
# 


#### REGRESION HYDRUS vs sonda validacion ####

lm = lm(hum$sonda10cm ~ hum$hum10cm, data = hum)
gof = gof(hum$hum10cm, hum$sonda10cm)
MAPE = MAPE(hum$hum10cm, hum$sonda10cm)

g5 = ggplot(hum, aes(x = hum10cm ,y = sonda10cm)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
  xlim(0.0,0.5) + 
  ylim(0.0,0.5) +
  ggtitle('HYDRUS vs sonda 10 cm') +
  xlab('humedad 10 cm HYDRUS') + ylab('humedad 10 cm sonda') +
  annotate("text",x= 0, y=0.5, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
  annotate("text",x= 0, y=0.45, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
  annotate("text",x= 0, y=0.40, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
  annotate("text",x= 0, y=0.35, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
  annotate("text",x= 0, y=0.30, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
  annotate("text",x= 0, y=0.25, hjust = 0, label=paste("KlingGupta =", round(gof[17],3))) +
  annotate("text",x= 0, y=0.20, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)

lm = lm(hum$sonda30cm ~ hum$hum30cm, data = hum)
gof = gof(hum$hum30cm, hum$sonda30cm)
MAPE = MAPE(hum$hum30cm, hum$sonda30cm)

g6 = ggplot(hum, aes(x = hum30cm ,y = sonda30cm)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
  xlim(0.25,0.45) + 
  ylim(0.25,0.45) +
  ggtitle('HYDRUS vs sonda 30 cm') +
  xlab('humedad 30 cm HYDRUS') + ylab('humedad 30 cm sonda') +
  annotate("text",x= 0.3, y=0.45, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
  annotate("text",x= 0.3, y=0.44, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
  annotate("text",x= 0.3, y=0.43, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
  annotate("text",x= 0.3, y=0.42, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
  annotate("text",x= 0.3, y=0.41, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
  annotate("text",x= 0.3, y=0.40, hjust = 0, label=paste("KlingGupta =", round(gof[17],3))) +
  annotate("text",x= 0.3, y=0.39, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


lm = lm(hum$sonda50cm~hum$hum50cm, data = hum)
gof = gof(hum$hum50cm, hum$sonda50cm)
MAPE = MAPE(hum$hum50cm, hum$sonda50cm)

g7 = ggplot(hum, aes(x = hum50cm ,y = sonda50cm)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
  xlim(0.25,0.45) + 
  ylim(0.25,0.45) +
  ggtitle('HYDRUS vs sonda 50 cm') +
  xlab('humedad 50 cm HYDRUS') + ylab('humedad 50 cm sonda') +
  annotate("text",x= 0.3, y=0.45, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
  annotate("text",x= 0.3, y=0.44, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
  annotate("text",x= 0.3, y=0.43, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
  annotate("text",x= 0.3, y=0.42, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
  annotate("text",x= 0.3, y=0.41, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
  annotate("text",x= 0.3, y=0.40, hjust = 0, label=paste("KlingGupta =", round(gof[17],3))) +
  annotate("text",x= 0.3, y=0.39, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)

# lm = lm(hum_val$sonda80cm ~ hum_val$hum80cm, data = hum)
# gof = gof(hum_val$hum80cm, hum_val$sonda80cm)
# MAPE = MAPE(hum_val$hum80cm, hum_val$sonda80cm)
# 
# g8 = ggplot(hum_val, aes(x = hum80cm ,y = sonda80cm)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.0,0.5) + 
#   ylim(0.0,0.5) +
#   ggtitle('HYDRUS vs sonda 80 cm') +
#   xlab('humedad 80 cm HYDRUS') + ylab('humedad 80 cm sonda') +
#   annotate("text",x= 0, y=0.5, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0, y=0.45, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0, y=0.40, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0, y=0.35, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0, y=0.30, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0, y=0.25, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0, y=0.20, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


#### SALIDAS HYDRUS Time Series ####

pp_df4 = data.frame(fecha = pp_df$fecha, pp = pp_df$pp/6, riego = pp_df$riego/6) # se multiplica para poder plotear junto a otro data frame con ggplot
pp_melt = reshape2::melt(pp_df4, id.vars = c('fecha'), measure.vars = c('pp', 'riego'))

TS_melt_1 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('TP_cm.dia', 'TR_cm.dia','PERCO', 'EVAP'))

g9 = ggplot(NULL, aes(x = fecha, y = value, color = variable)) +
  geom_bar(data = pp_melt[which(pp_melt$value>0),], stat = 'identity', position="stack") +
  geom_line(data = TS_melt_1, size = 1.0) + 
  geom_abline(intercept = 0, slope = 0) +
  ggtitle('Serie de tiempo') +
  xlab("Fecha") + ylab("TP, TR, EVAP, PERCO (cm/dia)") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  scale_y_continuous(sec.axis = sec_axis(~ .*6, name ="Precipitaciones y riegos (cm)")) + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
  geom_vline(xintercept = TS_melt_1$fecha[d_cal], linetype="dashed", 
             color = "red", size=1) 



# TS_melt_2 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Ke', 'Ke_HYDRUS'))
# 
# g10 = ggplot(TS_melt_2, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Ke') +
#   xlab("Fecha") + ylab("Ke, Ke_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_2$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Ke ~ TS_val$Ke_HYDRUS, data = TS_val)
# gof = gof(TS_val$Ke, TS_val$Ke_HYDRUS)
# MAPE = MAPE(TS_val$Ke, TS_val$Ke_HYDRUS)
# 
# g11 = ggplot(TS_val, aes(x = Ke_HYDRUS, y = Ke)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.003,0.005) + 
#   ylim(0.003,0.005) +
#   ggtitle('Ke HYDRUS vs Ke satelital') +
#   xlab('Ke HYDRUS') + ylab('Ke satelital') +
#   annotate("text",x= 0.003, y=0.005, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.003, y=0.0049, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.003, y=0.0048, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.003, y=0.0047, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.003, y=0.0046, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.003, y=0.0045, hjust = 0, label=paste("KlingGupta =", round(gof[17],3))) +
#   annotate("text",x= 0.003, y=0.0044, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


# TS_melt_3 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kcb', 'Kcb_HYDRUS'))
# 
# g12 = ggplot(TS_melt_3, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kcb') +
#   xlab("Fecha") + ylab("Kcb, Kcb_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_3$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Kcb ~ TS_val$Kcb_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kcb, TS_val$Kcb_HYDRUS)
# MAPE = MAPE(TS_val$Kcb, TS_val$Kcb_HYDRUS)
# 
# g13 = ggplot(TS_val, aes(x = Kcb_HYDRUS, y = Kcb)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.25,1.2) + 
#   ylim(0.25,1.2) +
#   ggtitle('Kcb HYDRUS vs Kcb satelital') +
#   xlab('Kcb HYDRUS') + ylab('Kcb satelital') +
#   annotate("text",x= 0.6, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.6, y=1.15, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.6, y=1.10, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.6, y=1.05, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.6, y=1.00, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.6, y=0.95, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.6, y=0.90, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)

# TS_melt_4 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kcb', 'Kcb_Pot_HYDRUS'))
# 
# g14 = ggplot(TS_melt_4, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kcb Potencial') +
#   xlab("Fecha") + ylab("Kcb, Kcb_Pot_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_4$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Kcb ~ TS_val$Kcb_Pot_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kcb, TS_val$Kcb_Pot_HYDRUS)
# MAPE = MAPE(TS_val$Kcb, TS_val$Kcb_Pot_HYDRUS)
# 
# g15 = ggplot(TS_val, aes(x = Kcb_Pot_HYDRUS, y = Kcb)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.7,1) + 
#   ylim(0.7,1) +
#   ggtitle('Kcb potencial  HYDRUS vs Kcb satelital') +
#   xlab('Kcb potencial HYDRUS') + ylab('Kcb satelital') +
#   annotate("text",x= 0.7, y=1.00, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.7, y=0.98, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.7, y=0.96, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.7, y=0.94, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.7, y=0.92, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.7, y=0.90, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.7, y=0.88, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)

# TS_melt_5 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kc', 'Kc_HYDRUS'))
# 
# g16 = ggplot(TS_melt_5, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kc') +
#   xlab("Fecha") + ylab("Kc, Kc_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_5$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)
# 
# 
# lm = lm(TS_val$Kc ~ TS_val$Kc_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kc, TS_val$Kc_HYDRUS)
# MAPE = MAPE(TS_val$Kc, TS_val$Kc_HYDRUS)
# 
# g17 = ggplot(TS_val, aes(x = Kc_HYDRUS, y = Kc)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.3,1.2) + 
#   ylim(0.3,1.2) +
#   ggtitle('Kc  HYDRUS vs Kc satelital') +
#   xlab('Kc HYDRUS') + ylab('Kc satelital') +
#   annotate("text",x= 0.7, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.7, y=1.15, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.7, y=1.10, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.7, y=1.05, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.7, y=1.00, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.7, y=0.95, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.7, y=0.90, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
# 
# TS_melt_6 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kc', 'Kc_Pot_HYDRUS'))
# 
# g18 = ggplot(TS_melt_6, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kc Potencial') +
#   xlab("Fecha") + ylab("Kc, Kc_Pot_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_6$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Kc ~ TS_val$Kc_Pot_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kc, TS_val$Kc_Pot_HYDRUS)
# MAPE = MAPE(TS_val$Kc, TS_val$Kc_Pot_HYDRUS)
# 
# g19 = ggplot(TS_val, aes(x = Kc_Pot_HYDRUS, y = Kc)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.6,1.2) + 
#   ylim(0.6,1.2) +
#   ggtitle('Kc potencial  HYDRUS vs Kc satelital') +
#   xlab('Kc potencial HYDRUS') + ylab('Kc satelital') +
#   annotate("text",x= 0.6, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.6, y=1.15, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.6, y=1.10, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.6, y=1.05, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.6, y=1.00, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.6, y=0.95, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.6, y=0.90, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)

#### PLOT pp y riego ####

pp_melt = reshape2::melt(pp_df, id.vars = c('fecha'), measure.vars = c('pp', 'riego'))

g20 = ggplot(pp_melt[which(pp_melt$value>0),], aes(x = fecha, y = value)) +
  geom_bar(stat="identity",position = 'stack', na.rm = TRUE, aes(fill = variable)) +
  xlab("Fecha") + ylab("pp (cm)") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 


#### PLOTs ETc

# ETc_melt = reshape2::melt(ETc_df, id.vars = c('fecha'), measure.vars = c('ETc_satelital', 'ETc_HYDRUS'))
# 
# g21 = ggplot(ETc_melt, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   ggtitle('Serie de tiempo ETc') +
#   xlab("Fecha") + ylab("ETc_Sat y ETc_HYDRUS (cm/dia)") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = ETc_melt$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1) #+
# #aes(group=rev(variable))


# lm = lm(ETc_df_val$ETc_satelital ~ ETc_df_val$ETc_HYDRUS, data = ETc_df_val)
# gof = gof(ETc_df_val$ETc_satelital, ETc_df_val$ETc_HYDRUS)
# MAPE = MAPE(ETc_df_val$ETc_satelital, ETc_df_val$ETc_HYDRUS)
# 
# g22 = ggplot(ETc_df_val, aes(x = ETc_HYDRUS, y = ETc_satelital)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0,1.2) + 
#   ylim(0,1.2) +
#   ggtitle('ETc HYDRUS vs ETc satelital') +
#   xlab('ETc HYDRUS') + ylab('ETc satelital') +
#   annotate("text",x= 0.0, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.0, y=1.13, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.0, y=1.06, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.0, y=0.99, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.0, y=0.92, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.0, y=0.85, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.0, y=0.78, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "lines"))

pp_df3 = data.frame(fecha = pp_df$fecha, pp = pp_df$pp*4.5, riego = pp_df$riego*4.5) # se multiplica para poder plotear junto a otro data frame con ggplot
pp_melt = reshape2::melt(pp_df3, id.vars = c('fecha'), measure.vars = c('pp', 'riego'))

WS_melt = reshape2::melt(WS_df,id.vars = c('fecha'), measure.vars = c('WS_sonda', 'WS_HYDRUS'))

g23 = ggplot(NULL, aes(x = fecha, y = value, color = variable)) +
  geom_bar(data = pp_melt[which(pp_melt$value>0),], stat = 'identity', position="stack") +
  geom_line(data = WS_melt, size = 0.4) + 
  scale_colour_manual(values = c("red", "yellow", "green", 'blue')) +
  ggtitle("Water storage") +
  xlab("Fecha") + ylab("water storage (cm)") +
  scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
  scale_y_continuous(sec.axis = sec_axis(~ ./4.5, name ="Precipitaciones y riegos (cm)")) + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
  geom_vline(xintercept = WS_melt$fecha[d_cal], linetype="dashed", 
             color = "red", size=1) + 
  scale_fill_manual(values = c("red", "yellow", "green", 'blue'))






lm = lm(WS_df_val$WS_sonda ~ WS_df_val$WS_HYDRUS, data = WS_df_val)
gof = gof(WS_df_val$WS_HYDRUS, WS_df_val$WS_sonda)
MAPE = MAPE(WS_df_val$WS_sonda, WS_df_val$WS_HYDRUS)

g24 = ggplot(WS_df_val, aes(x = WS_HYDRUS , y = WS_sonda )) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
  xlim(15,30) + 
  ylim(15,30) +
  ggtitle('WS HYDRUS vs WS sonda periodo post calibración') +
  xlab('WS HYDRUS (cm)') + ylab('WS sonda (cm)') +
  annotate("text",x= 15, y=30, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
  annotate("text",x= 15, y=29.6, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
  annotate("text",x= 15, y=29.2, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
  annotate("text",x= 15, y=28.8, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
  annotate("text",x= 15, y=28.4, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
  annotate("text",x= 15, y=28.0, hjust = 0, label=paste("KlingGupta =", round(gof[19],3))) +
  annotate("text",x= 15, y=27.6, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


lm = lm(WS_df_cal$WS_sonda ~ WS_df_cal$WS_HYDRUS, data = WS_df_cal)
gof = gof(WS_df_cal$WS_HYDRUS, WS_df_cal$WS_sonda)
MAPE = MAPE(WS_df_cal$WS_sonda, WS_df_cal$WS_HYDRUS)

g25 = ggplot(WS_df_cal, aes(x = WS_HYDRUS , y = WS_sonda )) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
  xlim(15,30) + 
  ylim(15,30) +
  ggtitle('WS HYDRUS vs WS sonda Periodo calibración') +
  xlab('WS HYDRUS (cm)') + ylab('WS sonda (cm)') +
  annotate("text",x= 15, y=30, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
  annotate("text",x= 15, y=29.6, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
  annotate("text",x= 15, y=29.2, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
  annotate("text",x= 15, y=28.8, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
  annotate("text",x= 15, y=28.4, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
  annotate("text",x= 15, y=28.0, hjust = 0, label=paste("KlingGupta =", round(gof[19],3))) +
  annotate("text",x= 15, y=27.6, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)

lm = lm(WS_df$WS_sonda ~ WS_df$WS_HYDRUS, data = WS_df)
gof = gof(WS_df$WS_HYDRUS, WS_df$WS_sonda)
MAPE = MAPE(WS_df$WS_sonda, WS_df$WS_HYDRUS)

g26 = ggplot(WS_df, aes(x = WS_HYDRUS , y = WS_sonda )) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
  xlim(15,30) + 
  ylim(15,30) +
  ggtitle('WS HYDRUS vs WS sonda periodo completo') +
  xlab('WS HYDRUS (cm)') + ylab('WS sonda (cm)') +
  annotate("text",x= 15, y=30, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
  annotate("text",x= 15, y=29.6, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
  annotate("text",x= 15, y=29.2, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
  annotate("text",x= 15, y=28.8, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
  annotate("text",x= 15, y=28.4, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
  annotate("text",x= 15, y=28.0, hjust = 0, label=paste("KlingGupta =", round(gof[19],3))) +
  annotate("text",x= 15, y=27.6, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


#### TODOS LOS PLOTS ####

# g1 # hum_10cm_TS
# g5 # hydrus_vs_sonda_10cm
# 
# g2 # hum_30cm_TS
# g6 # hydrus_vs_sonda_30cm
# 
# g3 # hum_50cm_TS
# g7 # hydrus_vs_sonda_50cm
# 
# g4 # hum_80cm_TS
# g8 # hydrus_vs_sonda_80cm
# 
# g9 # TP_TR_EVAP_PERCO_TS
# 
# g10 # Ke_TS
# g11 # Ke_hydrus_vs_Ke_sat
# 
# g12 # Kcb_TS
# g13 # Kcb_hydrus_vs_Kcb_sat
# 
# g14 # Kcb_TS_potencial
# g15 # Kcb_potencial_hydrus_vs_Kcb_sat
# 
# g16 # Kc_TS
# g17 # Kc_hydrus_vs_Kcb_sat
# 
# g18 # Kc_TS
# g19 # Kc_potencial_hydrus_vs_Kcb_sat
# 
# g20 # pp_riego
# 
# g21 # ETc TS
# g22 # ETc_hydrus_vs_ETc_sat
# 
# g23 # WS TS
# g24 # WS HYDRUS_vs_WS_sonda periodo post calibracion
# g25 # WS HYDRUS_vs_WS_sonda periodo calibracion
# g26 # WS HYDRUS_vs_WS_sonda todo el periodo

#### SAVE PLOTS ####

### individuales ####
ggsave(paste0(getwd(),'/plots/1_hum_10cm_TS.jpg') ,plot= g1 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/2_hum_30cm_TS.jpg') ,plot= g2 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/3_hum_50cm_TS.jpg') ,plot= g3 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/5_hydrus_vs_sonda_10cm.jpg') ,plot= g5 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/6_hydrus_vs_sonda_30cm.jpg') ,plot= g6 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/7_hydrus_vs_sonda_50cm.jpg') ,plot= g7 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/9_TP_TR_EVAP_PERCO_TS.jpg') ,plot= g9 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/10_Ke_TS.jpg') ,plot= g10 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/11_Ke_hydrus_vs_Ke_sat.jpg') ,plot= g11 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/12_Kcb_TS.jpg') ,plot= g12 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/13_Kcb_hydrus_vs_Kcb_sat.jpg') ,plot= g13 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/14_Kcb_potencial_TS.jpg') ,plot= g14 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/15_Kcb_potencial_hydrus_vs_Kcb_sat.jpg') ,plot= g15 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/16_Kc_TS.jpg') ,plot= g16 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/17_Kc_hydrus_vs_Kc_sat.jpg') ,plot= g17 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/18_Kc_potencial_TS.jpg') ,plot= g18 ,width = 8,height=6,units="in",dpi=300)
# ggsave(paste0(getwd(),'/plots/19_Kc_potencial_Kc_sat.jpg') ,plot= g19 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/20_pp_riego.jpg') ,plot= g20 ,width = 8,height=6,units="in",dpi=300)
#ggsave(paste0(getwd(),'/plots/21_ETc_TS.jpg') ,plot= g21 ,width = 8,height=6,units="in",dpi=300)
#ggsave(paste0(getwd(),'/plots/22_ETc_hydrus_vs_ETc_sat.jpg') ,plot= g22 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/23_WS_TS.jpg') ,plot= g23 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/24_WS_HYDRUSvs_WS_sonda_validacion.jpg') ,plot= g24 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/25_WS_HYDRUSvs_WS_sonda_calibracion.jpg') ,plot= g25 ,width = 8,height=6,units="in",dpi=300)
ggsave(paste0(getwd(),'/plots/26_WS_HYDRUSvs_WS_sonda_periodo_completo.jpg') ,plot= g26 ,width = 8,height=6,units="in",dpi=300)

tf = Sys.time()
print(tf - t0)
#### FIN ####




#### PLOTS

## Ke 

# TS_melt_2 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Ke', 'Ke_HYDRUS'))
# 
# g10 = ggplot(TS_melt_2, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Ke') +
#   xlab("Fecha") + ylab("Ke, Ke_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_2$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)

# lm = lm(TS_val$Ke ~ TS_val$Ke_HYDRUS, data = TS_val)
# gof = gof(TS_val$Ke, TS_val$Ke_HYDRUS)
# MAPE = MAPE(TS_val$Ke, TS_val$Ke_HYDRUS)
# 
# g11 = ggplot(TS_val, aes(x = Ke_HYDRUS, y = Ke)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.003,0.005) + 
#   ylim(0.003,0.005) +
#   ggtitle('Ke HYDRUS vs Ke satelital') +
#   xlab('Ke HYDRUS') + ylab('Ke satelital') +
#   annotate("text",x= 0.003, y=0.005, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.003, y=0.0049, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.003, y=0.0048, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.003, y=0.0047, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.003, y=0.0046, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.003, y=0.0045, hjust = 0, label=paste("KlingGupta =", round(gof[17],3))) +
#   annotate("text",x= 0.003, y=0.0044, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


## Kcb

# TS_melt_3 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kcb', 'Kcb_HYDRUS'))
# 
# g12 = ggplot(TS_melt_3, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kcb') +
#   xlab("Fecha") + ylab("Kcb, Kcb_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_3$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Kcb ~ TS_val$Kcb_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kcb, TS_val$Kcb_HYDRUS)
# MAPE = MAPE(TS_val$Kcb, TS_val$Kcb_HYDRUS)
# 
# g13 = ggplot(TS_val, aes(x = Kcb_HYDRUS, y = Kcb)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.25,1.2) + 
#   ylim(0.25,1.2) +
#   ggtitle('Kcb HYDRUS vs Kcb satelital') +
#   xlab('Kcb HYDRUS') + ylab('Kcb satelital') +
#   annotate("text",x= 0.6, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.6, y=1.15, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.6, y=1.10, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.6, y=1.05, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.6, y=1.00, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.6, y=0.95, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.6, y=0.90, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


# TS_melt_4 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kcb', 'Kcb_Pot_HYDRUS'))
# 
# g14 = ggplot(TS_melt_4, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kcb Potencial') +
#   xlab("Fecha") + ylab("Kcb, Kcb_Pot_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_4$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Kcb ~ TS_val$Kcb_Pot_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kcb, TS_val$Kcb_Pot_HYDRUS)
# MAPE = MAPE(TS_val$Kcb, TS_val$Kcb_Pot_HYDRUS)
# 
# g15 = ggplot(TS_val, aes(x = Kcb_Pot_HYDRUS, y = Kcb)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.7,1) + 
#   ylim(0.7,1) +
#   ggtitle('Kcb potencial  HYDRUS vs Kcb satelital') +
#   xlab('Kcb potencial HYDRUS') + ylab('Kcb satelital') +
#   annotate("text",x= 0.7, y=1.00, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.7, y=0.98, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.7, y=0.96, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.7, y=0.94, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.7, y=0.92, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.7, y=0.90, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.7, y=0.88, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


## Kc

# TS_melt_5 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kc', 'Kc_HYDRUS'))
# 
# g16 = ggplot(TS_melt_5, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kc') +
#   xlab("Fecha") + ylab("Kc, Kc_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_5$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)
# 
# 
# lm = lm(TS_val$Kc ~ TS_val$Kc_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kc, TS_val$Kc_HYDRUS)
# MAPE = MAPE(TS_val$Kc, TS_val$Kc_HYDRUS)
# 
# g17 = ggplot(TS_val, aes(x = Kc_HYDRUS, y = Kc)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.3,1.2) + 
#   ylim(0.3,1.2) +
#   ggtitle('Kc  HYDRUS vs Kc satelital') +
#   xlab('Kc HYDRUS') + ylab('Kc satelital') +
#   annotate("text",x= 0.7, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.7, y=1.15, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.7, y=1.10, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.7, y=1.05, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.7, y=1.00, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.7, y=0.95, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.7, y=0.90, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
# 
# TS_melt_6 = reshape2::melt(TS, id.vars = c('fecha'), measure.vars = c('Kc', 'Kc_Pot_HYDRUS'))
# 
# g18 = ggplot(TS_melt_6, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   geom_abline(intercept = 0, slope = 0) +
#   ggtitle('Serie de tiempo Kc Potencial') +
#   xlab("Fecha") + ylab("Kc, Kc_Pot_HYDRUS") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = TS_melt_6$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1)


# lm = lm(TS_val$Kc ~ TS_val$Kc_Pot_HYDRUS, data = TS_val)
# gof = gof(TS_val$Kc, TS_val$Kc_Pot_HYDRUS)
# MAPE = MAPE(TS_val$Kc, TS_val$Kc_Pot_HYDRUS)
# 
# g19 = ggplot(TS_val, aes(x = Kc_Pot_HYDRUS, y = Kc)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0.6,1.2) + 
#   ylim(0.6,1.2) +
#   ggtitle('Kc potencial  HYDRUS vs Kc satelital') +
#   xlab('Kc potencial HYDRUS') + ylab('Kc satelital') +
#   annotate("text",x= 0.6, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.6, y=1.15, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.6, y=1.10, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.6, y=1.05, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.6, y=1.00, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.6, y=0.95, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.6, y=0.90, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)


## ETc

# ETc_melt = reshape2::melt(ETc_df, id.vars = c('fecha'), measure.vars = c('ETc_satelital', 'ETc_HYDRUS'))
# 
# g21 = ggplot(ETc_melt, aes(x = fecha, y = value, color = variable)) +
#   geom_line(size = 1.0) + 
#   ggtitle('Serie de tiempo ETc') +
#   xlab("Fecha") + ylab("ETc_Sat y ETc_HYDRUS (cm/dia)") +
#   scale_x_date(breaks = '2 month', labels = fecha, date_labels = '%b %y') +
#   theme(plot.title = element_text(lineheight=.8, face="bold", size = 20)) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 0.3) +
#   geom_vline(xintercept = ETc_melt$fecha[d_cal], linetype="dashed", 
#              color = "red", size=1) #+
# #aes(group=rev(variable))


# lm = lm(ETc_df_val$ETc_satelital ~ ETc_df_val$ETc_HYDRUS, data = ETc_df_val)
# gof = gof(ETc_df_val$ETc_satelital, ETc_df_val$ETc_HYDRUS)
# MAPE = MAPE(ETc_df_val$ETc_satelital, ETc_df_val$ETc_HYDRUS)
# 
# g22 = ggplot(ETc_df_val, aes(x = ETc_HYDRUS, y = ETc_satelital)) +
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
#   geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2]) +
#   xlim(0,1.2) + 
#   ylim(0,1.2) +
#   ggtitle('ETc HYDRUS vs ETc satelital') +
#   xlab('ETc HYDRUS') + ylab('ETc satelital') +
#   annotate("text",x= 0.0, y=1.2, hjust = 0, label=paste("Slope =", round(lm$coefficients[2],2)))+
#   annotate("text",x= 0.0, y=1.13, hjust = 0, label=paste("R2 =", round(gof[17],3))) +
#   annotate("text",x= 0.0, y=1.06, hjust = 0, label=paste("PBIAS =", round(gof[6],3), '%')) +
#   annotate("text",x= 0.0, y=0.99, hjust = 0, label=paste("RMSE =", round(gof[4],3))) +
#   annotate("text",x= 0.0, y=0.92, hjust = 0, label=paste("NRMSE =", round(gof[5],3), '%')) +
#   annotate("text",x= 0.0, y=0.85, hjust = 0, label=paste("MAPE =", round(MAPE,3))) +
#   annotate("text",x= 0.0, y=0.78, hjust = 0, label=paste("NASH =", round(gof[9],3))) +
#   theme_bw() + theme( panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "lines"))




