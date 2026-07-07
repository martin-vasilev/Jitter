
rm(list= ls())


######## dependent variables analysis
#DV1: corrective saccade probability (Bernoulli GLMM)

library(tidyverse)

rs  <- read.csv("data/return-sweep_data.csv")
qst <- read.csv("data/quest.csv")              
fix <- read.csv("data/raw_fixations.csv")     


#calculate the pixels per character: 18.1 pixel/xharacter, 0.337 deg/char
l1  <- fix %>% filter(!is.na(SFIX), line == 1, outsideText == 0, !is.na(char_line))
PXC <- unname(coef(lm(xPos ~ char_line, data = l1))["char_line"])


dat <- rs %>%
  transmute(
    sub, item, cond,
    corr_sacc_prob = undersweep_prob,

    #   DV2 = landStart_RS + (x_corrective - x_RS) / PXC -> the corrective saccade exceeds the character margin
    corr_sacc_land_pos = if_else(undersweep_prob == 1,
                                 landStart + (nextX - xPos) / PXC,
                                 NA_real_),
    
    launch_site_raw = launchSite,
    # "return-sweep landing position" covariate = line-relative
    #   landStart (how far into line 2 the RS landed), NOT the word-relative
    #   land_pos column. Confirm。
    land_pos_raw    = landStart,
    DS_in_RS
  )


#relevel and dummy coding
dat$cond <- factor(dat$cond, levels = c(1, 2, 3),
                   labels = c("normal", "left", "right"))
dat$cond <- relevel(dat$cond, ref = "normal")   


#centering the covariates
dat$launch_site <- as.numeric(scale(dat$launch_site_raw, center = TRUE, scale = FALSE))
dat$land_pos    <- as.numeric(scale(dat$land_pos_raw,    center = TRUE, scale = FALSE))


##exclusion criteria
n_start <- nrow(dat)
#Display change must have occurred within the return-sweep window,
#not triggered outside the RS saccade, or completed >10 ms after the line-2 started
dat <- dat %>% filter(DS_in_RS == 1)

#Participant comprehension accuracy < 70% -> drop participant
acc      <- qst %>% group_by(sub) %>% summarise(acc = mean(accuracy, na.rm = TRUE))
keep_sub <- acc %>% filter(acc >= 0.70) %>% pull(sub)
dat      <- dat %>% filter(sub %in% keep_sub)

# Missing data: drop participants retaining < 50% of experimental trials.
retained <- dat %>% filter(item %in% 1:90) %>% count(sub, name = "n_kept")
drop_sub <- retained %>% filter(n_kept < 0.50 * 90) %>% pull(sub)
dat      <- dat %>% filter(!(sub %in% drop_sub))


dat$sub  <- factor(dat$sub)
dat$item <- factor(dat$item)




dat_prob <- dat                                            # DV1: all valid RS
dat_land <- dat %>% filter(!is.na(corr_sacc_land_pos))     # DV2: undersweeps only


cat(sprintf("px/char used            : %.2f\n", PXC))
cat(sprintf("RS before exclusions    : %d\n", n_start))
cat(sprintf("RS after exclusions     : %d  (DV1 / dat_prob)\n", nrow(dat_prob)))
cat(sprintf("Undersweeps (DV2 subset): %d  (dat_land)\n", nrow(dat_land)))
cat("\nP(corrective saccade) by condition (DV1):\n")
print(dat_prob %>% group_by(cond) %>%
        summarise(p = mean(corr_sacc_prob), n = n()))
cat("\nCorrective landing position by condition (DV2, chars):\n")
print(dat_land %>% group_by(cond) %>%
        summarise(mean = mean(corr_sacc_land_pos),
                  sd = sd(corr_sacc_land_pos), n = n()))


##### corrective saccade probability -- Buenoulli GLMM
library(tidyverse)
library(brms)

# 1 = control, 2 = left, 3 = right.

dat_prob$cond <- factor(dat_prob$cond,
                        levels = c("normal", "left", "right"))  # 'normal' = control
dat_prob$cond <- relevel(dat_prob$cond, ref = "normal")

f_prob <- bf(corr_sacc_prob ~ cond + launch_site + land_pos +
               (cond | sub) + (cond | item))

get_prior(f_prob, data = dat_prob, family = bernoulli())
#NAs： delete or filter ？？？
colSums(is.na(dat_prob[, c("corr_sacc_prob","launch_site","land_pos")]))
dat_prob[is.na(dat_prob$launch_site), 
         c("sub","item","cond","corr_sacc_prob","launch_site_raw")]


priors_prob <- c(
  set_prior("normal(0, 2)",    class = "b", coef = "condleft"),    
  set_prior("normal(0, 2)",    class = "b", coef = "condright"),   
  set_prior("normal(-0.3, 1)", class = "b", coef = "launch_site"), 
  set_prior("normal(0, 2)",    class = "b", coef = "land_pos"),
  set_prior("normal(0.6, 1)",  class = "Intercept")                
)

# ---- ---------------------------------------------------------------
m_prob <- brm(
  formula      = f_prob,
  data         = dat_prob,
  family       = bernoulli(),       
  prior        = priors_prob,
  sample_prior = TRUE,               
  warmup       = 1000,
  iter         = 6000,
  chains       = 4,
  cores        = 4,
  seed         = 1234,
  control      = list(adapt_delta = 0.95)  # separation -> reduce divergences
)


summary(m_prob)                
plot(m_prob)                   
pp_check(m_prob, ndraws = 100) 

# Inspect the condition contrasts (this is the DV1 result):
fixef(m_prob)[c("condleft", "condright"), ]


#####landing position -- Guassion GLMM

dat_land$cond <- factor(dat_land$cond, levels = c("normal", "left", "right"))
dat_land$cond <- relevel(dat_land$cond, ref = "normal")

f_land <- bf(corr_sacc_land_pos ~ cond + launch_site + land_pos +
               (cond | sub) + (cond | item))

get_prior(f_land, data = dat_land, family = gaussian())

priors_land <- c(
  set_prior("normal(0, 2)", class = "b", coef = "condleft"),    
  set_prior("normal(0, 2)", class = "b", coef = "condright"),  
  set_prior("normal(0, 1)", class = "b", coef = "launch_site"),
  set_prior("normal(0, 2)", class = "b", coef = "land_pos"),
  set_prior("normal(4, 2)", class = "Intercept") 
)


m_land <- brm(
  formula      = f_land,
  data         = dat_land,
  family       = gaussian(),
  prior        = priors_land,
  sample_prior = TRUE,              
  warmup       = 1000,
  iter         = 6000,
  chains       = 4,
  cores        = 4,
  seed         = 1234,
  control      = list(adapt_delta = 0.95)
)


summary(m_land)                 
plot(m_land)                   
pp_check(m_land, ndraws = 100)  

# The DV2 result — the condition contrasts on landing position:
fixef(m_land)[c("condleft", "condright"), ]

prior_summary(m_land)


### additional: latency










