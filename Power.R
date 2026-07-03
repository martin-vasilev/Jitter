
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


NSim      <- 200 # number of simulations per cell 
nsub      <- c(12, 24) # number of subjects
nitems    <- 90
data_loss <- 0.15
effects   <- c('condleft', 'condright')  

ES_reduction<- 0.25

if(ES_reduction>0){
  b_CS[2:3]<- b_CS[2:3]* (1- ES_reduction)
}

power <- NULL
for (i in seq_along(nsub)) {
  
  
  LSQD <- c(rep(c('normal','left','right'), nitems/3),
            rep(c('left','right','normal'), nitems/3),
            rep(c('right','normal','left'), nitems/3))
  
  sim_data <- data.frame(cond = rep(LSQD, nsub[i] / 3))
  sim_data$sub  <- rep(seq_len(nsub[i]), each = nitems)
  sim_data$item <- rep(1:nitems, times = nsub[i])
  sim_data$cond <- fct_relevel(as.factor(sim_data$cond),
                               'normal', 'left', 'right')
  
  
  sim_data2 <- sim_data[-sample(nrow(sim_data),
                                round(data_loss * nrow(sim_data))), ]
  
  
  model_CS <- makeGlmer(corr_prob ~ cond + (cond|sub) + (cond|item),
                        fixef = b_CS, VarCorr = RE_CS,
                        data = sim_data2, family = binomial)
  
  for (eff in effects) {
    p_CS <- suppressMessages(
      powerSim(model_CS, nsim = NSim, alpha = .05, progress = T,
               test = simr::fixed(xname = eff, method = "z")))
    s <- summary(p_CS)
    power <- rbind(power, data.frame(
      nsub   = nsub[i],
      effect = eff,
      power  = s$mean,
      lower  = s$lower,
      upper  = s$upper))
  }
  cat("done nsub =", nsub[i], "\n")
}

power 
print(power)


#landing position

#shift in characters
#data_choice <- "all" 

RS <- read.csv("Pilot/data/final.csv") %>% rename(cond = cond.x)
#RS <- RS %>% rename(cond = cond.x)
# RS <- switch(data_choice,
#              all = RS,
#              new = filter(RS, sample == "New"),
#              i90 = filter(RS, item <= 90),
#              stop("data_choice must be 'all', 'new', or 'i90'"))


RS <- RS %>%
  mutate(
    corr_prob = ifelse(Rtn_sweep_type == "undersweep", 1L, 0L),
    cond = recode(as.character(cond), "1"="normal", "2"="left", "3"="right"),
    cond = fct_relevel(factor(cond), "normal", "left", "right"),
    cs_land_abs = corr_sacc_land_char   
  )

RS_cs <- RS %>% filter(corr_prob == 1, !is.na(cs_land_abs))

RS %>% group_by(cond) %>%
  summarise(cs_land     = mean(cs_land_abs, na.rm = TRUE),
            n_corr_sacc = sum(corr_prob))


# model
M_land <- lmer(cs_land_abs ~ cond + (cond | sub) + (1 | item), data = RS_cs)
summary(M_land)

b_land  <- fixef(M_land)     # effect sizes (intercept, condleft, condright)
RE_land <- VarCorr(M_land)   # between-subject / between-item variance
s_land  <- sigma(M_land)     # residual SD

#power

NSim      <- 200                   
nsub_grid <- c(12, 24, 36, 48, 54)    
nitems    <- 90                   
data_loss <- 0.2    
ES_reduction<- 0.25

if(ES_reduction>0){
  b_land[2:3]<- b_land[2:3]* (1- ES_reduction)
}

us_rate <- RS %>% group_by(cond) %>% summarise(p = mean(corr_prob)) %>% deframe()

make_design <- function(nsub, nitems) {
  LSQD <- c(rep(c("normal","left","right"), nitems/3),
            rep(c("left","right","normal"), nitems/3),
            rep(c("right","normal","left"), nitems/3))
  d <- data.frame(cond = rep(LSQD, nsub/3))
  d$sub  <- rep(seq_len(nsub), each = nitems)
  d$item <- rep(seq_len(nitems), times = nsub)
  d$cond <- fct_relevel(factor(d$cond), "normal","left","right")
  d
}

results <- list()
for (n in nsub_grid) {
  
  cat(n);
  
  d <- make_design(n, nitems)
  d <- d[-sample(nrow(d), round(data_loss * nrow(d))), ]      # random trial loss
  
  # keep only trials with a corrective saccade (condition-specific rate)
  keep   <- runif(nrow(d)) < us_rate[as.character(d$cond)]
  d_land <- d[keep, ]
  d_land$cs_land_abs <- rnorm(n = nrow(d_land))                                     # placeholder; overwritten
  
  # the RE structure here MUST match M_land above -> (1|sub)+(1|item)
  m_land <- makeLmer(cs_land_abs ~ cond + (cond | sub) + (1 | item),
                     fixef = b_land, VarCorr = RE_land, sigma = s_land,
                     data = d_land)
  
  for (eff in c("condleft", "condright")) {
    ps <- powerSim(m_land, nsim = NSim, alpha = .05, progress = T,
                   test = simr::fixed(xname = eff, method = "z"))
    s  <- summary(ps)                                         # cols: mean, lower, upper
    results[[length(results) + 1]] <- data.frame(
      nsub = n, effect = eff,
      power = s$mean, lower = s$lower, upper = s$upper)
  }
}

power_table <- bind_rows(results)
print(power_table)









