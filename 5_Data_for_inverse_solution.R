################################################
###### Data for inverse solution HYDRUS 1-D
################################################

###### MODIFICABLE ######

sonda = read.csv('D:/MEMORIA_LRI/HYDRUS_1D/salidas/MODELACION_DEFINITIVA_70_cm_tests/sonda_&_pp_riego.txt', sep = '\t')
dias = 658 # dias a modelar
obs.nodes =c(10, 30, 50, 80) # profundidades nodos de observacion (cm)
obs.nodes.position = seq(from = 1, to = length(obs.nodes), by = 1)
sonda = data.frame(fecha = as.vector(sonda$date), # data frame con fecha y humedad observada
                   obs.10 = sonda$SondaHum_10cm...,
                   obs.30 = sonda$SondaHum_30cm...,
                   obs.50 = sonda$SondaHum_50cm...,
                   obs.80 = sonda$SondaHum_80cm...)

type = 2 # water content
weight = 1

###### NO MODIFICABLE ######

obs.10 = sonda$obs.10 
obs.30 = sonda$obs.30
obs.50 = sonda$obs.50
obs.80 = sonda$obs.80


x = c()
for (i in 1:dias) {
  d = rep(i,length(obs.nodes))
  x = c(x, d)
  print(i)
}  

y = c() 
for (i in 1:dias) {
  hum.x = c(obs.10[i], obs.30[i], obs.50[i], obs.80[i])
  y = c(y, hum.x)
  print(i)
}


position = rep(obs.nodes.position, dias)

df = data.frame(x = x, y = y, type = type, position = position, weight = weight)

write.csv(df,paste0('D:/MEMORIA_LRI/HYDRUS_1D/salidas/MODELACION_DEFINITIVA_70_cm/Data_inverse_HYDRUS_',dias,'_days.csv'))
