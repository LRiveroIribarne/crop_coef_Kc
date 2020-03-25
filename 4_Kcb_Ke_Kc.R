#######################################################################################
## CODIGOS MEMORIA 
## 4 - Kcb, Ke, Kc
## Lucas Rivero Iribarne
#######################################################################################


############################################################################################################
## Este codigo calcula coeficiente basal de cultivo (Kcb), Coeficiende de evaporacion
##  desde el suelo (Ke) y coeficiente de cultivo (Kc) 
##  Requiere de:
##
##  1) Serie temporal de LAI y fPAR
##  2) Serie temporal de NDVI o SAVI
##  3) Serie temporal de u2 (velocidad del viento) diaria y HRmin (Humedad relativa minima) diaria
##  4) Altura promedio de las plantas 
##
## de la siguiente forma:
##
##  1) Se elige la combinacion de inputs (comb)
##  2) Se ingresan variables requeridas
##  3) Luego, en el apartado "No modificable":
##      1) Se calcula Kcb (Poças, 2015)
##      2) Se calcula Ke por metodo FAO 56 (Allen, 1998)
##      3) Se calcula Kc por metodo dual de FAO (Allen, 1998)
##      4) Se extraen Kcb, Ke y Kc en pixel sonda
##      5) Se interpolan valores entre fechas donde hay escenas
##
############################################################################################################

# ---- DEFINIR DIRECTORIOS ----
setwd('D:/MEMORIA_LRI') # directorio
dir.in = paste(getwd(),"/Materiales/Kx", sep = "") # directorio entrada con imgs. Sentinel 2 .SAFE
dir.out = paste(getwd(),"/Resultados/Kx", sep = "") # directorio salida

# ---- MODIFICABLE ----
# Combinaciones posibles:
#   (1) Kcb a partir de NDVI, LAI y fPAR a partir de CI 110
#   (2) Kcb a partir de NDVI, LAI y fPAR a partir de S2
#   (3) Kcb a partir de SAVI, LAI y fPAR a partir de CI 110
#   (4) Kcb a partir de SAVI, LAI y fPAR a partir de CI S2

comb = 3 # opciones (1), (2), (3), (4)

Kcmin = 0.15 # (Allen, 2009) modificado
h = 4 # altura maxima promedio de arboles en (m) para periodo estudiado
am = 3 # ancho mojamiento (m)
dist_sh = 2 # distancia sobre hilera (m)
dist_eh = 6 # distancia entre hilera (m)
fc = 2 # opciones (1) = fc (Allen, 1998), (2) = fc de Sentinel 2 SNAP
fw = 2 # opciones (1) = Riego por goteo segun (Allen, 1998), (2) = ecuacion (Lobos, 2017) 

# ---- NO MODIFICABLE ----
# ---- Kcb, Ke, Kc ---- 

t0 = Sys.time() # variable para medir tiempo de ejecucion codigo. NO RELEVANTE

NDVI = ndvi_stack_ROI # serie temporal de NDVI
SAVI = savi_stack_ROI # serie temporal de SAVI
fpar_CI = fpar_stack # serie temporal fPAR a partir de CI
lai_CI = lai_stack # serie temporal LAI a partir de CI
fpar_S2 = fpar_biophyscal_stack # serie temporal fPAR S2
lai_S2 = lai_biophysical_stack # serie temporal LAI S2
fcover_S2 = fcover_biophysical_stack # serie temporal FCOVER S2
nam2 = nam2 # fechas imgs "20171210" yyyymmdd

if (comb == 1) { # Kcb a partir de NDVI, LAI y fPAR a partir de CI 110
  
  IV = NDVI
  lai = lai_CI
  fpar = fpar_CI
  names = c('NDVI', 'LAI_CI', 'FPAR_CI', 'CI')
  
} else if (comb == 2) { # Kcb a partir de NDVI, LAI y fPAR a partir de S2
  
  IV = NDVI
  lai = lai_S2
  fpar = fpar_S2
  names = c('NDVI', 'LAI_S2', 'FPAR_S2', 'S2')
  
} else if (comb == 3) { # Kcb a partir de SAVI, LAI y fPAR a partir de CI 110
  
  IV = SAVI
  lai = lai_CI
  fpar = fpar_CI
  names = c('SAVI', 'LAI_CI', 'FPAR_CI', 'CI')
  
} else if (comb == 4) { # Kcb a partir de SAVI, LAI y fPAR a partir de S2
  
  IV = SAVI
  lai = lai_S2
  fpar = fpar_S2
  names = c('SAVI', 'LAI_S2', 'FPAR_S2', 'S2')
  
} else {
  cat('VALOR NO VALIDO! \nSeleccionar: 
      \n(1) Kcb a partir de NDVI, LAI y fPAR a partir de CI 110
      \n(2) Kcb a partir de NDVI, LAI y fPAR a partir de S2
      \n(3) Kcb a partir de SAVI, LAI y fPAR a partir de CI 110
      \n(4) Kcb a partir de SAVI, LAI y fPAR a partir de S2')
}

cat(paste('Se utilizara', names[1], 'para calcular Kcb,', names[2], 'y', names[3], 'para calcular Kd y fw respectivamente' ))
Sys.sleep(2)

# ---- Kcb ----
Kcmin = Kcmin 
Kd = 1 - exp(-0.7 * lai)
IVmax = cellStats(IV, max)
IVmin = cellStats(IV, min)
IVmax = max(IVmax)
IVmin = min(IVmin)

Kcb = Kcmin + Kd * ((IV - IVmin)/(IVmax - IVmin))

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Kcb/tiff'))
for (i in 1:nlayers(Kcb)) {
  
  writeRaster(Kcb[[i]], filename = paste0(nam2[i],'_Kcb_','.tif'),format = 'GTiff', overwrite = TRUE)
  print(i)
}

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Kcb/jpg'))
for (i in 1:nlayers(Kcb)) {
  
  jpeg(paste0(nam2[i],'_Kcb_','.jpg'))
  
  plot(Kcb[[i]])
  
  dev.off() 
  print(i)
}

print('Kcb LISTO!')
Sys.sleep(2)

# Ke 
h = h
meteo = read.csv('D:/MEMORIA_LRI/Materiales/meteo/meteo_imgs.csv')
Kr = 1 # coeficiente disminucion de evaporación desde el suelo
#u2 = c(1.0,0.6,0.7,0.6,0.7,1.0,0.8,0.7,0.8,0.8,0.5,0.4,0.7,0.5,0.8,0.3,0.2,0.3,0.5,0.3,0.3,0.6,0.6,0.4,0.5,0.5,0.6,0.5,0.8,0.6,0.6,0.9,0.6,0.8,0.6,0.8) # velocidad del viento a 2 m de altura en m/s
u2 = meteo$Velocidad.del.Viento..m.s.
u2 = Kd * 0 + u2 # cualquier raster por 0 + el valor que quiero asignar
#HRmin = c(24.1,21.2,20.1,20.1,36.2,46.9,29.4,26.6,37.8,34.2,16.6,32.1,17.6,25.4,12,39.7,34.6,41.2,22,47.1,40,28.6,37.8,17.9,20.9,39,17.7,13.2,25.9,19.2,21.2,34.6,18.2,34.8,22.4,34.9)
HRmin = meteo$HRmin....
HRmin = Kd * 0 + HRmin # Humedad relativa minima

Kcmax = list()
for (j in 1:length(imgs_L2A)) {
  
  Kcmax1 = (1.1 + ((0.04* (u2[[j]] - 2) - 0.004 * (HRmin[[j]] - 45)) * (h/3)^0.3))  
  Kcmax2 = (Kcb[[j]] + 0.05)
  
  KcmaxStack = stack(Kcmax1, Kcmax2)
  maxStack = max(KcmaxStack, na.rm = TRUE)
  
  Kcmax[[j]] = maxStack
  names(Kcmax)[j] = nam2[j]
  print(j)
}
Kcmax = stack(Kcmax)

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Kc_max/tiff'))
for (i in 1:nlayers(Kcmax)) {
  
  writeRaster(Kcmax[[i]], filename = paste0(nam2[i],'_Kcmax_','.tif'),format = 'GTiff', overwrite = TRUE)
  print(i)
}

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Kc_max/jpg'))
for (i in 1:nlayers(Kcmax)) {
  
  jpeg(paste0(nam2[i],'_Kcmax_','.jpg'))
  
  plot(Kcmax[[i]])
  
  dev.off() 
  print(i)
}

if (fc == 1) {
  fc = ((Kcb - Kcmin)/(Kcmax - Kcmin))^(1 + 0.5*h)
} else if (fc == 2) {
  fc = mask(crop(fcover_S2_stack, ROI), ROI)
} else { print('VALOR NO VALIDO de fc!')}

if (fw == 1) {
  fw = 0.35
} else if (fw == 2) {
  fw = (am * dist_sh * (1 - fpar))/(dist_eh * dist_eh) # (ancho.mojamiento * dist.sobre.hilera * (1 - FI))/ (marco.de.plantacion)
} else {print('VALOR NO VALIDO de fw!')}

few = list()
for (l in 1:length(imgs_L2A)) {
  
  few1 = 1 - fc
  few2 = fw
  
  fewminStack = stack(few1, few2)
  minStack = min(fewminStack, na.rm = TRUE)
  
  few[[l]] = minStack
  names(few)[l] = nam2[l]
  print(l)
}
few = stack(few)


Ke = list()
for (l in 1:length(imgs_L2A)) {
  
  ke1 = Kr*(Kcmax[[l]] - Kcb[[l]])
  ke2 = few[[l]] * Kcmax[[l]]
  
  keminStack = stack(ke1, ke2)
  minStack = min(keminStack, na.rm = TRUE)
  
  Ke[[l]] = minStack
  names(Ke)[l] = nam2[l]
  print(l)
}
Ke = stack(Ke)

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Ke/tiff'))
for (i in 1:nlayers(Ke)) {
  
  writeRaster(Ke[[i]], filename = paste0(nam2[i],'_Ke_','.tif'),format = 'GTiff', overwrite = TRUE)
  print(i)
}

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Ke/jpg'))
for (i in 1:nlayers(Ke)) {
  
  jpeg(paste0(nam2[i],'_Ke_','.jpg'))
  
  plot(Ke[[i]])
  
  dev.off() 
  print(i)
}

print('Ke LISTO!')
Sys.sleep(2)

# Kc
Kc = Kcb + Ke
names(Kc) = nam2

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Kc/tiff'))
for (i in 1:nlayers(Kc)) {
  
  writeRaster(Kc[[i]], filename = paste0(nam2[i],'_Kc_','.tif'),format = 'GTiff', overwrite = TRUE)
  print(i)
}

setwd(paste0(dir.out, '/', names[4],'_',names[1], '/Kc/jpg'))
for (i in 1:nlayers(Kc)) {
  
  jpeg(paste0(nam2[i],'_Kc_','.jpg'))
  
  plot(Kc[[i]])
  
  dev.off() 
  print(i)
}
print('Kc LISTO!')
Sys.sleep(2)

# ---- EXTRACT PIXEL SONDA ---- 

Kcb_sonda = extract(Kcb, sonda)
Ke_sonda = extract(Ke, sonda)
Kc_sonda = extract(Kc, sonda)
Kcmax_sonda = extract(Kcmax, sonda)


day = substr(x = nam2,start = 7 ,stop = 8)
month = substr(x = nam2,start = 5 ,stop = 6) 
year = substr(x = nam2,start = 1 ,stop = 4)
fecha = paste0(year,'-',month,'-',day)
fecha = as.Date(fecha, format = '%Y-%m-%d')
sonda_Kx_df = data.frame(fecha = fecha,
                         Kcb = c(Kcb_sonda),
                         Ke = c(Ke_sonda),
                         Kc = c(Kc_sonda),
                         Kcmax = c(Kcmax_sonda))

write.csv(sonda_Kx_df,paste0(dir.out, '/', names[4],'_',names[1], '/SONDA/sonda_Kx_',names[4],'_',names[1],'.csv'))

ApproxKcb_sonda = approxfun(x = sonda_Kx_df$fecha, y = sonda_Kx_df$Kcb)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitKcb_sonda = ApproxKcb_sonda(Dates)

ApproxKe_sonda = approxfun(x = sonda_Kx_df$fecha, y = sonda_Kx_df$Ke)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitKe_sonda = ApproxKe_sonda(Dates)

ApproxKc_sonda = approxfun(x = sonda_Kx_df$fecha, y = sonda_Kx_df$Kc)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitKc_sonda = ApproxKc_sonda(Dates)

ApproxKcmax_sonda = approxfun(x = sonda_Kx_df$fecha, y = sonda_Kx_df$Kcmax)
Dates = seq.Date(ymd("2017-12-10"), ymd("2019-11-20"), by = 1)
LinearFitKcmax_sonda = ApproxKcmax_sonda(Dates)

sonda_Kx_df_temporal = data.frame(fecha = Dates,
                                  Kcb = LinearFitKcb_sonda,
                                  Ke = LinearFitKe_sonda,
                                  Kc = LinearFitKc_sonda,
                                  Kcmax = LinearFitKcmax_sonda)

write.csv(sonda_Kx_df_temporal, paste0(dir.out, '/', names[4],'_',names[1], '/SONDA/sonda_Kx_',names[4],'_',names[1],'_temporal.csv'))

tf = Sys.time()
print(tf - t0)
#### FIN ####
