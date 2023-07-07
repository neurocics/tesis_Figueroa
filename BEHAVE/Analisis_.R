rm(list=ls())

# Instalar paquetes
#install.packages("runjags")
#install.packages("coda")
#install.packages("rjags")
library(runjags)
library(coda)
library(rjags)
library(ggplot2)
library(nlme)
require(lme4)
library(MASS)

library(rstatix)
library(tidyr)
#---------------------------------
datafile = "/Users/alejandra/Desktop/DOCTORADO/2021/TESIS\ /git_neurocovid/tesis_Figueroa/BEHAVE/COR_rl.txt" # AF
datafile = "~/Documents/GitHub/tesis_Figueroa/BEHAVE/COR_rl.txt" # PB

#---------------------------------
#AF
setwd("/Users/alejandra/Desktop/DOCTORADO/2021/TESIS\ /git_neurocovid/tesis_Figueroa/BEHAVE")
#PB
setwd("~/Documents/GitHub/tesis_Figueroa/BEHAVE")

DATA_RL <-  read.table(datafile, header = TRUE, sep="\t")
names(DATA_RL)
DATA_RL$Sub_ID=as.numeric(as.factor(DATA_RL$SU))
max(DATA_RL$Sub_ID)
length(unique(DATA_RL$Sub_ID[DATA_RL$GR==" controles"]))
length(unique(DATA_RL$Sub_ID[DATA_RL$GR==" pacientes"]))

# Prepare the data for JAGS

sub_ID = sort(unique(DATA_RL$Sub))






DATA_RL$nt=DATA_RL$out_nt
DATA_RL$buena =DATA_RL$out_buena
DATA_RL$resp     = as.numeric(DATA_RL$out_resp)
DATA_RL$Resp     = as.numeric(DATA_RL$out_resp)
DATA_RL$Bd1   = as.numeric(DATA_RL$out_Bd1)
DATA_RL$Bd2   = as.numeric(DATA_RL$out_Bd2)
DATA_RL$prob   = as.numeric(DATA_RL$out_prob)
DATA_RL$F   = as.numeric(DATA_RL$out_F)
ID = as.numeric(DATA_RL$Sub_ID)
nSUB = max(sub_ID)
nT = length(DATA_RL$out_Bd1)

firstT = as.numeric((c(-1, DATA_RL$prob[1:nT-1]-DATA_RL$prob[2:nT] )!=0) | (c(-1, DATA_RL$buena[1:nT-1]-DATA_RL$buena[2:nT] )!=0) | DATA_RL$nt==1)
indx = which(firstT == 1)
indxf = indx [ seq(1,length(indx),by=2) ]
DATA_RL$firstT=0 
DATA_RL$firstT[indxf]=1 
firstT = DATA_RL$firstT

shift = as.numeric((c(-1, DATA_RL$prob[1:nT-1]-DATA_RL$prob[2:nT] )!=0) | (c(-1, DATA_RL$buena[1:nT-1]-DATA_RL$buena[2:nT] )!=0) | DATA_RL$nt==1)
indx = which(shift == 1)
indxs = indx [ seq(2,length(indx),by=2) ]
DATA_RL$shift=0 
DATA_RL$shift[indxs]=1 

DATA_RL$adap =0 
DATA_RL$adap[indxs]=1
DATA_RL$adap[indxs+1]=1
DATA_RL$adap[indxs+2]=1
DATA_RL$adap[indxs+3]=1
DATA_RL$adap[indxs+4]=1
DATA_RL$adap[indxs+5]=1
DATA_RL$adap[indxs+6]=1
DATA_RL$adap[indxs+7]=1

#DATA_RL$adap =0 
DATA_RL$adap[indxs-1]=-1
DATA_RL$adap[indxs-2]=-1
DATA_RL$adap[indxs-3]=-1
DATA_RL$adap[indxs-4]=-1
DATA_RL$adap[indxs-5]=-1
DATA_RL$adap[indxs-6]=-1
DATA_RL$adap[indxs-7]=-1

endT = which(firstT==1)-1
endT = c(endT[2:length(endT)], length(DATA_RL$adap))

for (n in 1:length(endT)) {DATA_RL$adap[(indxs[n]+5):endT[n]]=2 }

DATA_RL$adap[endT]=3


DATA_RL$behavior_Shift <-  c(abs(diff( DATA_RL$Resp)),0) 

# cambio no asdaptativo despues de una perdida 
forIDX = DATA_RL[((DATA_RL$Resp==1 & DATA_RL$prob>0.5 )| (DATA_RL$Resp==2 & DATA_RL$prob<0.5 )) & (DATA_RL$F==0) 
                 & (DATA_RL$adap<0) , ]
forIDX$adap=0
forIDX0 = aggregate(behavior_Shift ~ adap:Sub_ID:GR, data=forIDX,FUN=mean  )


# cambio adaptativo despues de una perdida 
forIDX = DATA_RL[((DATA_RL$Resp==1 & DATA_RL$prob<0.5) | (DATA_RL$Resp==2 & DATA_RL$prob>0.5 ))& (DATA_RL$F==0) 
                 & (DATA_RL$adap<0) , ]
forIDX$adap=0
forIDX0m = aggregate(behavior_Shift ~ adap:Sub_ID:GR, data=forIDX,FUN=mean  )

# cambio adaptativo despues de una perdida en el shift 
forIDX = DATA_RL[((DATA_RL$Resp==1 & DATA_RL$prob<0.5 )|( DATA_RL$Resp==2 & DATA_RL$prob>0.5 ))& (DATA_RL$F==0) 
                 & (DATA_RL$adap==1), ]
forIDX1 = aggregate(behavior_Shift ~ adap:Sub_ID:GR, data=forIDX,FUN=mean  )

# cambio no adaptativo despues de una perdida en el shift 
forIDX = DATA_RL[((DATA_RL$Resp==1 & DATA_RL$prob>0.5 )|( DATA_RL$Resp==2 & DATA_RL$prob<0.5 ))& (DATA_RL$F==0) 
                 & (DATA_RL$adap==1), ]
forIDX1m = aggregate(behavior_Shift ~ adap:Sub_ID:GR, data=forIDX,FUN=mean  )



forIDX = DATA_RL[((DATA_RL$Resp==1 & DATA_RL$prob>0.5 )|( DATA_RL$Resp==2 & DATA_RL$prob<0.5 ))& (DATA_RL$F==0) 
                 & (DATA_RL$adap==2) , ]
forIDX2 = aggregate(behavior_Shift ~ adap:Sub_ID:GR, data=forIDX,FUN=mean  )
forIDX = DATA_RL[((DATA_RL$Resp==1 & DATA_RL$prob<0.5 )|( DATA_RL$Resp==2 & DATA_RL$prob>0.5 ))& (DATA_RL$F==0) 
                 & (DATA_RL$adap==2) , ]
forIDX2m = aggregate(behavior_Shift ~ adap:Sub_ID:GR, data=forIDX,FUN=mean  )


mean(forIDX0$behavior_Shift)  # good shift
mean(forIDX0m$behavior_Shift) # bad shift

mean(forIDX1$behavior_Shift)
mean(forIDX1m$behavior_Shift)


mean(forIDX2$behavior_Shift)
mean(forIDX2m$behavior_Shift)


forIDX = rbind(forIDX0,forIDX1,forIDX2)
friedman.test(behavior_Shift ~ adap | Sub_ID, data=forIDX)


mixed_anova <- aov_ez(
  id = "Sub_ID",
  dv = "behavior_Shift",
  data = forIDX,
  between = "GR",
  within = "adap"
)



DF = forIDX0$behavior_Shift-forIDX1$behavior_Shift

wilcox.test(forIDX0$behavior_Shift-forIDX1$behavior_Shift)
mean(forIDX1$behavior_Shift-forIDX0$behavior_Shift)

wilcox.test(forIDX0$behavior_Shift-forIDX0m$behavior_Shift)
wilcox.test(forIDX2$behavior_Shift-forIDX2m$behavior_Shift)

wilcox.test(forIDX0$behavior_Shift-forIDX2$behavior_Shift)
wilcox.test(forIDX0m$behavior_Shift-forIDX2m$behavior_Shift)


wilcox.test(forIDX1$behavior_Shift ~ forIDX1$GR)
wilcox.test(forIDX0$behavior_Shift ~ forIDX0$GR)


wilcox.test(DF ~ forIDX1$GR)

summary(glm(I(forIDX1$behavior_Shift-forIDX2$behavior_Shift) ~ forIDX1$GR))
summary(rlm(I(forIDX1$behavior_Shift-forIDX2$behavior_Shift) ~ forIDX1$GR))




wilcox.test((forIDX0$behavior_Shift-forIDX0m$behavior_Shift)-(forIDX1$behavior_Shift-forIDX1m$behavior_Shift))



forIDX = rbind(forIDX0,forIDX1,forIDX2)
friedman.test(behavior_Shift ~ adap|Sub_ID, data=forIDX)
aggregate(behavior_Shift ~ adap, data=forIDX, FUN=mean)

forIDX = rbind(forIDX0m,forIDX1m,forIDX2m)
friedman.test(behavior_Shift ~ adap|Sub_ID, data=forIDX)
aggregate(behavior_Shift ~ adap, data=forIDX, FUN=mean)

####





#cor.test(ACUg$accuracy[forMODEL$acu_g>0.6],forIDX0$behavior_Shift[forMODEL$acu_g>0.6], method = "spearman" )
#cor.test(ACU$accuracy[forMODEL$acu_g>0.6],forIDX0$behavior_Shift[forMODEL$acu_g>0.6] , method = "spearman" )

#cor.test(ACUg$accuracy[forMODEL$acu_g>0.6],forIDX1$behavior_Shift[forMODEL$acu_g>0.6], method = "spearman" )
#cor.test(ACU$accuracy[forMODEL$acu_g>0.6],forIDX1$behavior_Shift[forMODEL$acu_g>0.6] , method = "spearman" )

#cor.test(R.E.seqG[forMODEL$acu_g>0.6,3],forIDX0$behavior_Shift[forMODEL$acu_g>0.6] , method = "spearman" )


forMODEL = data.frame(shift_noA = (forIDX0$behavior_Shift),
                      shift_A = (forIDX1$behavior_Shift))
forMODEL$acu_ng= (ACU$accuracy)
forMODEL$acu_g= (ACUg$accuracy)
forMODEL$proactive= scale(as.numeric(R.E.seqG$Sub[,2]))
#forMODEL$proactive_c= scale(as.numeric(R.E.seqG[,2]))

modelo_lineal = lm( shift_noA ~ acu_ng + acu_g + proactive  , data=forMODEL[forMODEL$acu_g>0.2,])
summary(modelo_lineal)

modelo_lineal = lm( shift_A ~  acu_ng + acu_g + proactive, data=forMODEL[forMODEL$acu_g>0.2,])
summary(modelo_lineal)

scatter.smooth(forMODEL$proactive[forMODEL$acu_g>0.2],forIDX1$behavior_Shift[forMODEL$acu_g>0.2] ,span=1.5)   


RT = glmer(I(Resp==1) ~ Bd1 +Bd2  + prob  + (1 +  Bd1 + Bd2   + prob | Sub_ID), family = binomial, data = DATA_RL)
summary(RT)


y   = as.numeric(DATA_RL$Resp==1) # Left choice 

feedback_left = (DATA_RL$F*(DATA_RL$Resp==1)) + ((1-DATA_RL$F)*(1-(DATA_RL$Resp==1)))

dat <- dump.format(list(y=as.numeric(DATA_RL$Resp==1),
                        Bd1=DATA_RL$Bd1,
                        Bd2=DATA_RL$Bd2,
                        prob=DATA_RL$prob,
                        feedback=DATA_RL$F,
                        feedback_left = feedback_left,
                        ID=ID,
                        nT=nT,
                        s=ID,
                        nSUB=nSUB,
                        firstT=firstT
))

# Initialize chains
inits1 <- dump.format(list( mu.beta0=0.5,mu.beta1=-0.5,mu.beta2=0.5,mu.beta3=-0.5, .RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list( mu.beta0=5,mu.beta1=5,mu.beta2=-5,mu.beta3=-5, .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list( mu.beta0=-5,mu.beta1=-5,mu.beta2=5,mu.beta3=5,  .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))


# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.beta2", "mu.beta3","deviance")

# Run the function that fits the models using JAGS
results <- run.jags(model="modelo_logic_RL.R",
                    monitor=monitor, data=dat, n.chains=3, method="parallel",
                    inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results


# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columns represent the monitored variables
chains = rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
DIC = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC  # 3623.963

# Initialize chains
inits1 <- dump.format(list( mu.beta0=0.5,mu.beta1=-0.5,mu.beta2=0.5,mu.beta3=-0.5, alpha.mu = 0.5 ,.RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list( mu.beta0=5,mu.beta1=5,mu.beta2=-5,mu.beta3=-5, alpha.mu = 0.2 , .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list( mu.beta0=-5,mu.beta1=-5,mu.beta2=5,mu.beta3=5, alpha.mu = 0.7 ,  .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))


# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.beta2", "mu.beta3","deviance","alpha.mu")

# Run the function that fits the models using JAGS
results2 <- run.jags(model="modelo_logic_RL_learn.R",
                     monitor=monitor, data=dat, n.chains=3, method="parallel",
                     inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results2


# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columns represent the monitored variables
chains = rbind(results2$mcmc[[1]], results2$mcmc[[2]], results2$mcmc[[3]])
DIC2 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC2  # 2622.058 


# Initialize chains
inits1 <- dump.format(list( mu.beta0=-0.5,mu.beta1=0.1,mu.rew=0.8,mu.gamma=0.9, alpha.mu = 0.5 ,.RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list( mu.beta0=0,mu.beta1=-0.1,mu.rew=0.9,mu.gamma=1.1, alpha.mu = 0.2 , .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list( mu.beta0=-5,mu.beta1=0,mu.rew=1,mu.gamma=1, alpha.mu = 0.7 ,  .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))


# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.rew", "mu.gamma","deviance","alpha.mu")

# Run the function that fits the models using JAGS
results3 <- run.jags(model="modelo_prospect_RL_learn.R",
                     monitor=monitor, data=dat, n.chains=3, method="parallel",
                     inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results3
plot(results3)

# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columnsrepresent the monitored variables
chains = rbind(results3$mcmc[[1]], results3$mcmc[[2]], results3$mcmc[[3]])
DIC3 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC3  
# 2588.332
# 7811.224 --67s

# Initialize chains
inits1 <- dump.format(list( mu.beta0=-0.5,mu.beta1=0.1,mu.rew=0.8,mu.gamma=0.9, alpha.mu = 0.5 ,.RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list( mu.beta0=0,mu.beta1=-0.1,mu.rew=0.9,mu.gamma=1.1, alpha.mu = 0.2 , .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list( mu.beta0=-5,mu.beta1=0,mu.rew=1,mu.gamma=1, alpha.mu = 0.7 ,  .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))


# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.rew", "mu.gamma","deviance","alpha.mu")

# Run the function that fits the models using JAGS
results4 <- run.jags(model="modelo_prospect_RL_learn_2p.R",
                     monitor=monitor, data=dat, n.chains=3, method="parallel",
                     inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results4


# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columnsrepresent the monitored variables
chains = rbind(results4$mcmc[[1]], results4$mcmc[[2]], results4$mcmc[[3]])
DIC4 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC4  # 3932.88

# Initialize chains
inits1 <- dump.format(list( mu.beta0=-0.5,mu.beta1=0.1,mu.rew=0.8,mu.gamma=0.9, alpha.mu = 0.5 ,alpha_b.mu = 0.4 ,.RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list( mu.beta0=0,mu.beta1=-0.1,mu.rew=0.9,mu.gamma=1.1, alpha.mu = 0.2 ,alpha_b.mu = 0.3 , .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list( mu.beta0=-5,mu.beta1=0.2,mu.rew=0.7,mu.gamma=1, alpha.mu = 0.4 , alpha_b.mu = 0.2 , .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))


# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.rew", "mu.gamma","deviance","alpha.mu","alpha_b.mu")

# Run the function that fits the models using JAGS
results5 <- run.jags(model="modelo_prospect_RL_learn_belief.R",
                     monitor=monitor, data=dat, n.chains=3, method="parallel",
                     inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results5


# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columnsrepresent the monitored variables
chains = rbind(results5$mcmc[[1]], results5$mcmc[[2]], results5$mcmc[[3]])
DIC5 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC5  # 4013.697

# Initialize chains
inits1 <- dump.format(list( mu.beta0=-0.5,mu.beta1=0.1,mu.rew=0.8,mu.gamma=0.9, alpha_a.mu = 0.5 ,alpha_r.mu = 0.4 ,.RNG.name="base::Super-Duper", .RNG.seed=99999 ))
inits2 <- dump.format(list( mu.beta0=0,mu.beta1=-0.1,mu.rew=0.9,mu.gamma=1.1, alpha_a.mu = 0.2 ,alpha_r.mu = 0.3 , .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 ))
inits3 <- dump.format(list( mu.beta0=-5,mu.beta1=0.2,mu.rew=0.7,mu.gamma=1, alpha_a.mu = 0.4 , alpha_r.mu = 0.2 , .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 ))


# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.rew", "mu.gamma","deviance","alpha_a.mu","alpha_r.mu")

# Run the function that fits the models using JAGS
results6 <- run.jags(model="modelo_prospect_RL_learn_2a.R",
                     monitor=monitor, data=dat, n.chains=3, method="parallel",
                     inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results6
plot(results6)

# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columnsrepresent the monitored variables
chains = rbind(results6$mcmc[[1]], results6$mcmc[[2]], results6$mcmc[[3]])
DIC6 = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC6  # 2441.411 -- 67s





# Tell JAGS which latent variables to monitor
# In this case thet is a vector with Nsubjs positions (see JAGS code associated to this script)
monitor = c("mu.beta0", "mu.beta1", "mu.rew", "mu.gamma","deviance","alpha_a.mu","alpha_r.mu","alpha_A.mu","alpha_H.mu")

# Run the function that fits the models using JAGS
results6 <- run.jags(model="modelo_prospect_RL_learn_2a.R",
                     monitor=monitor, data=dat, n.chains=3, method="parallel",
                     inits=c(inits1, inits2, inits3), plots = FALSE, burnin=5000, sample=1000, thin=5)
results6
plot(results6)

# readout the 3 chains from the "results" structure and combine them into a single matrix
# each of the resulting matrix represent a single MCMC sample, the columnsrepresent the monitored variables
chains = rbind(results6ah$mcmc[[1]], results6ah$mcmc[[2]], results6ah$mcmc[[3]])
DIC6ah = mean(chains[,"deviance"]) + (sd(chains[,"deviance"])^2)/2
DIC6ah  # 



##### PARA PROBAR

dat_p <- dump.format(list(y=as.numeric(DATA_RL$Resp==1),
                        Bd1=DATA_RL$Bd1,
                        Bd2=DATA_RL$Bd2,
                        prob=DATA_RL$prob,
                        feedback=DATA_RL$F,
                        ID=ID,
                        nT=3000:7000,
                        s=ID,
                        nSUB=nSUB,
                        firstT=firstT,
                        Anosmia=DATA_RL$ANOSMIA,
                        Hospital=DATA_RL$HOSPITAL
))

