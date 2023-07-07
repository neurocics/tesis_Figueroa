library(runjags)
library(rjags)
library(coda)
library(ggplot2)
library(ggpubr)
library(lme4)
library(tidyverse)
require(lme4)
library(rstatix)

# Version 07.07.2023

rm(list=ls())
setwd("/Volumes/Alehermosa/TESIS")

datafile = "COR_rl.txt" #

DATA_TOTAL = read.table(datafile,sep = '\t', header = T)

unique(DATA_TOTAL$SU)
DATA_TOTAL$GR[DATA_TOTAL$GR==" controles"]="HC"
DATA_TOTAL$GR[DATA_TOTAL$GR==" pacientes"]="MS"

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

if (sum(firstT) == 4*2*nSU) {print_color("Todo bien mi amo", "green")} else 
  {warning("Esta mal la deteccion del los priemros trail")
    aggregate(firstT ~ SU , FUN=sum, data=DATA_TOTAL)
    
    } # detecta todos los primeros trials (de cada bloque y los de shift)




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

DATA_TOTAL$money = (DATA_TOTAL$Bd1*(DATA_TOTAL$Resp==1)*(DATA_TOTAL$F==1))+
  (DATA_TOTAL$Bd2*(DATA_TOTAL$Resp==2)*(DATA_TOTAL$F==1))


FX = aggregate(money ~ SU:fase:GR,data=DATA_TOTAL,FUN = mean) # puntaje acumulado por fase 
# fase 0 pre shift
# fase 1 shift 
# fase 2 post shift 

FX <- FX %>% filter(fase %in% c(-1,1,2))

FXm = aggregate(money ~ fase,data=FX,FUN = mean)

res.aov <- anova_test(
  data = FX, dv = money, wid = SU,
  between = GR, within = fase
)
get_anova_table(res.aov)

pwc <- FX %>% 
  group_by(GR) %>%
  emmeans_test(money ~ fase, p.adjust.method = "bonferroni") 
  #pairwise_wilcox_test(money ~ fase, p.adjust.method = "none")
pwc



pwc <- pwc %>% add_xy_position(x = "fase")
pwc <- pwc %>%
  mutate(y.position = c(1,1.1, 1.2,1.3, 1.4,1.5)*4,
         xmin = c(0.8,0.8,1.8,    1.2,1.2,2.2),
         xmax = c(1.8,2.8,2.8,    2.2,3.2,3.2),
  )
#pwc <- pwc %>%3
#  mutate(x.position = c(0.9,1.1,1.9,2.1,2.9,3.1))

#pwc.filtered <- pwc %>% filter(time != "t1")
bxp <- ggboxplot(
  FX, x = "fase", y = "money",
  color = "GR" ,#palette = "jco"
)

bxp + 
  ylab("Money")+
  scale_x_discrete("Stage", labels=c("-1"="preS","1"="Shift","2"="prostS" ))+
  stat_pvalue_manual(pwc,hide.ns = TRUE,tip.length = 0.02) + #tip.length = 0, hide.ns = TRUE,label = "p = {p.adj}"
  labs(
    title="Earning Modulation",
    subtitle = get_test_label(res.aov, detailed = TRUE,row=2),
    caption =  get_pwc_label(pwc)
  )


#
DATA_TOTAL$behavior_Shift <-  c(abs(diff( DATA_TOTAL$Resp)),0) 

# cambio adaptativo despues de una perdida
forIDX = DATA_TOTAL[((DATA_TOTAL$Resp==1 & DATA_TOTAL$prob<0.5) | 
                     (DATA_TOTAL$Resp==2 & DATA_TOTAL$prob>0.5)     )& (DATA_TOTAL$F==0) 
                 , ]
#forIDX$adap=0
forIDXm = aggregate(behavior_Shift ~ fase:SU:GR, data=forIDX,FUN=mean  )
forIDXm <- forIDXm %>% filter(fase %in% c(-1,1,2))

res.aov <- anova_test(
  data = forIDXm, dv = behavior_Shift, wid = SU,
  between = GR, within = fase
)
get_anova_table(res.aov)

pwc <- forIDXm %>% 
  group_by(GR) %>%
  emmeans_test(behavior_Shift ~ fase, p.adjust.method = "bonferroni") 
  #pairwise_wilcox_test(behavior_Shift ~ fase, p.adjust.method = "fdr")
pwc



pwc <- pwc %>% add_xy_position(x = "fase")
pwc <- pwc %>%
  mutate(y.position = c(1.1,1.1, 1.2,1.3, 1.4,1.5),
         xmin = c(0.8,0.8,1.8,    1.2,1.2,2.2),
         xmax = c(1.8,2.8,2.8,    2.2,3.2,3.2),
  )
#pwc <- pwc %>%3
#  mutate(x.position = c(0.9,1.1,1.9,2.1,2.9,3.1))

#pwc.filtered <- pwc %>% filter(time != "t1")
bxp <- ggboxplot(
  forIDXm, x = "fase", y = "behavior_Shift",
  color = "GR" ,#palette = "jco"
)

bxp + 
  ylab("Adaptative Shift")+
  scale_x_discrete("Stage", labels=c("-1"="preS","1"="Shift","2"="prostS" ))+
  stat_pvalue_manual(pwc,hide.ns = TRUE,tip.length = 0.02) + #tip.length = 0, hide.ns = TRUE,label = "p = {p.adj}"
  labs(
    title="Adaptative Shift Modulation",
    subtitle = get_test_label(res.aov, detailed = TRUE,row=2),
    caption =  get_pwc_label(pwc)
  )

