library(runjags)
library(rjags)
library(coda)
library(ggplot2)
library(lme4)

require(lme4)

rm(list=ls())
setwd("/Volumes/Alehermosa/TESIS")

datafile = "COR_rl.txt" # PB

DATA_TOTAL = read.table(datafile,sep = '\t', header = T)

unique(DATA_TOTAL$SU)



DATA_TOTAL$nt=DATA_TOTAL$out_nt
DATA_TOTAL$buena =DATA_TOTAL$out_buena
DATA_TOTAL$Resp     = as.numeric(DATA_TOTAL$out_resp)
DATA_TOTAL$Bd1   = as.numeric(DATA_TOTAL$out_Bd1)
DATA_TOTAL$Bd2   = as.numeric(DATA_TOTAL$out_Bd2)
DATA_TOTAL$prob   = as.numeric(DATA_TOTAL$out_prob)
DATA_TOTAL$F   = as.numeric(DATA_TOTAL$out_F)

ID = as.numeric(DATA_TOTAL$SU)
ID = as.numeric(as.factor(DATA_TOTAL$SU))
sort(unique(ID))
nSU = max(ID) # número total de sujetos
nT = length(DATA_TOTAL$out_Bd1)


# detectar el primer trial y el shift de cada bloque 

firstT = as.numeric( (c(-1, diff(DATA_TOTAL$prob) )!=0) | (c(-1, diff(DATA_TOTAL$buena)!=0)) | DATA_TOTAL$nt==1) # detecta el primer trial 

library(insight)

if (sum(firstT) == 4*2*nSU) {print_color("Todo bien mi amo", "green")} else {warning("Esta mal la deteccion del los priemros trail")} # detecta todos los primeros trials (de cada bloque y los de shift)

indx = which(firstT == 1)
indxf = indx [ seq(1,length(indx),by=2) ]

DATA_TOTAL$firstT=0 # detecta el primer trial de cada bloque 
DATA_TOTAL$firstT[indxf]=1 

DATA_TOTAL$shift=0 # detecta el primer trial después del shift
DATA_TOTAL$shift[indx [seq(2,length(indx),by=2)]]=1 

firstT = DATA_TOTAL$firstT

# shift 
indx = which(DATA_TOTAL$shift == 1)

# end 
endT = which(firstT==1)-1
endT = c(endT[2:length(endT)], 
dim(DATA_TOTAL)[1])

nt_p = 7 # numero de ensayos para los episodios a analizar 

DATA_TOTAL$fase = 0

for (i in 1:nt_p) {
  DATA_TOTAL$fase[endT-i+1]  = 2  
  DATA_TOTAL$fase[indx+i-1] = 1
  DATA_TOTAL$fase[indx-i]   = -1
}

DATA_TOTAL$money = (1*(DATA_TOTAL$Resp==1)*(DATA_TOTAL$F==1))+
  (1*(DATA_TOTAL$Resp==2)*(DATA_TOTAL$F==1))


FX = aggregate(money ~ SU:fase,data=DATA_TOTAL,FUN = mean) # puntaje acumulado por fase 
# fase 0 pre shift
# fase 1 shift 
# fase 2 post shift 

FX <- FX %>% filter(fase %in% c(-1,1,2))

FXm = aggregate(money ~ fase,data=FX,FUN = mean)


