
rm(list= ls())

RS <- read.csv("~/R/Jitter/Pilot/data/final.csv")
RS$cond<- RS$cond.x
RS$cond.x<- NULL

library(simr)
library(tidyverse)
library(lmerTest)
library(ggeffects)

RS<-RS %>%
  mutate(corr_prob= ifelse(Rtn_sweep_type=='undersweep', 1, 0),
         cond= recode(cond, '1'= 'normal', '2'= 'left', '3'= 'right'),
         cond= as.factor(cond),
         cond= fct_relevel(cond, 'normal', 'left', 'right'))


RS%>%
  group_by(cond)%>%
  summarise(Corr_prob= mean(corr_prob),
            Land_pos= mean(char_line, na.rm= T),
            corr_sacc_land= mean(corr_sacc_land_char, na.rm=T))


#RS$cond<- as.factor(RS$cond)
levels(RS$cond)

M1<- glmer(corr_prob ~ cond +(cond|sub)+ (cond|item),
           data= RS, family = binomial)
summary(M1)

plot(ggeffect(M1, 'cond'))

## extract model coefficients:

b_CS <- coef(summary(M1))[,1] # fixed intercept and slopes
RE_CS <- VarCorr(M1) # random effects
#s_CS <- sigma(FFD) # residual sd


NSim= 10 # number of simulations per cell 
nsub = seq(6, 9, 3)  # number of subjects
nitems = 90 
data_loss<- 0.15 # percentage of data points to randomly remove



power<- NULL

for(i in 1:length(nsub)){
  
  LSQD<- c(rep(c("normal", "left", 'right'), nitems/3),
           rep(c("left", "right", 'normal'), nitems/3),
           rep(c("right", 'normal', "left"), nitems/3))
  
  sim_data <- data.frame(cond= rep(LSQD, nsub[i]/3))
  
  sim_data$sub <- rep(seq(1, nsub[i], 1), each= nitems)
  sim_data$item <- rep(1:nitems, times = nsub[i])
  
  #table(sim_data$item, sim_data$Cond)
  sim_data$rand<- rnorm(n= nrow(sim_data))
  sim_data$corr_prob<- ifelse(sim_data$rand>0, 1, 0)
  
  # simulate data loss:
  sim_data2<- sim_data[-sample(nrow(sim_data),
                                round(data_loss*nrow(sim_data))), ]
  
  sim_data2<-sim_data2 %>%mutate(
#           cond= recode(cond, '1'= 'normal', '2'= 'left', '3'= 'right'),
           cond= as.factor(cond),
           cond= fct_relevel(cond, 'normal', 'left', 'right'))
  
  contrasts(sim_data2$cond)
  
  model_CS <- makeGlmer(corr_prob~ cond+ (cond|sub)+ (cond|item),
                        fixef=b_CS, VarCorr= RE_CS, data=sim_data2,
                        family = binomial)
  
  summary(model_CS)
  
  ### Power:
  suppressMessages(
    p_CS<- powerSim(model_CS, nsim= NSim,
                     test = simr::fixed(xname = 'condright',
                     method = "z"), alpha =.05,
                    progress= T))
  
  sum_CS<- summary(p_CS)
  
  
}



sim_data$Subject <- rep(seq(1, nsub/2, 1),
                        each= nitems)
sim_data$item <- rep(1:nitems[j], times = nsub[i]/2)
sim_data$FFD<- rnorm(n= nrow(sim_data))
sim_data$SFD<- rnorm(n= nrow(sim_data))
sim_data$GD<- rnorm(n= nrow(sim_data))
sim_data$TVT<- rnorm(n= nrow(sim_data))












