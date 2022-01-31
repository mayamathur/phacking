
# PRELIMINARIES ---------------------------------------------------------------


# data-wrangling packages
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(here)
library(xtable)
library(testthat)

results.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2022-1-27 autoregressive processes and conservatism"


options(scipen=999)

# PART 2: LOOK AT LAST-DRAW TRY-TO-NMAX HACKING AND ITS STATIONARITY ---------------------------------------

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)")
source("helper_SAPH.R")



thetaiF = c()


# simulate from scratch
#**note I'm holding muin constant in order to condition on it
# since SE = 0 (almost), affirmative results are anything positive
dataset.name = "sim_data_10"
k = 5000
for ( i in 1:k ) {
  .d = sim_one_study_set(Nmax = 50,
                         Mu = -2,
                         T2 = 0, # keep all mui's the same for now
                         m = 50,
                         t2w = .5,
                         se = 0.0001,  # no within-study error for now
                         rho = 0.9,
                         hack = "no", #@@HANGED
                         return.only.published = FALSE)
  
  hack = .d$hack[1]
  
  # was this study successful?
  if ( hack == "no" ) {
    .d$study.success = any(.d$affirm == TRUE)
  } else {
    # when hack = "no", this should be equivalent to the above, I think
    .d$study.success = ( .d$affirm[.d$Di == 1] == TRUE )
  }
  
  .d$study = i
  .d$draw.index = 1:nrow(.d)
  
  # all studies
  if ( i == 1 ) d = .d else d = bind_rows(d, .d)
  
  thetaiF = c( thetaiF, .d$yi[.d$Di == 1] )
  
  # covEmp = c(covEmp, unique(d$covEmp))
  # rhoEmp = c(rhoEmp, unique(d$rhoEmp))
  # varEmp = c(varEmp, var(d$muin))
}
setwd(results.dir)


fwrite( d, paste(dataset.name, ".csv", sep = "") )


# read in existing dataset
setwd(results.dir)
d = fread( paste(dataset.name, ".csv", sep = "") )
dim(d)
table(d$Di)  # should match number of studies
table( d$study.success[ !duplicated(d$study) ] )


# make subsets
dn = d %>% filter(affirm == FALSE)
du = d %>% filter(study.success == FALSE)
length(unique(du$study))/length(unique(d$study))  # prop. of studies unsuccessful

# look at first few studies
View( d %>% select(study, draw.index, Di, study.success, affirm, mui, yi ) )

# sanity check for hack = "no" case
# look at some studies that had some draws affirmative and others not
temp = d %>% group_by(study) %>%
  mutate( study.prop.affirm = mean(affirm) ) %>%
  filter( study.prop.affirm != 0 & study.prop.affirm != 1 ) %>%
  select(study, draw.index, Di, study.success, study.prop.affirm, affirm, mui, yi )
View(temp)



# ~ Condition on nonaffirmative results ---------------------------------------

# ~~ Term 2 ---------------
# Cov(Fi, yi | mu_i^* = m, Ain^*=0)
cor.test(dn$Di, dn$yi) 
# yes, close to 0! 

# ~~ Term 2A ---------------
# P(all draws nonaffirm | yi=t, mu_i^* = m, Ain^*=0)
# should be decreasing in t for positive autocorrelation
mod = glm( (study.success==0) ~ yi,
             data = dn,
         family = binomial )
summary(mod)
exp(coef(mod))
# yes :)

# ~~ Term 2B ---------------
# P(this draw is favored | all draws nonaffirm, yi=t, mu_i^* = m, Ain^*=0)
# should also be decreasing in t
# **this one is actually very slightly increasing in t
mod = glm( Di ~ yi,
           data = dn[ dn$study.success == 0, ],
           family = binomial )
summary(mod)
exp(coef(mod))

# accounting for clustering w/in studies:
library(geepack)
geeMod = geeglm( Di ~ yi,
                 id = study,
                 family = binomial(link="logit"),
                 corstr = "exchangeable",
                 data = dn[ dn$study.success == 0, ] )
summary(geeMod)
exp(coef(geeMod))
# tiny p-value, but point estimate is 1?
# might be an issue with Wald inference

# singular
library(lme4)
glmmMod = glmer( Di ~ yi + (1|study),
                 family = binomial(link="logit"),
                 #corstr = "exchangeable",
                 data = dn[ dn$study.success == 0, ] )
summary(glmmMod)



# another way to look at size of favored draws among unsuccessful studies
du %>% group_by(Di) %>%
  summarise(mean(yi))




# ~ For no-hacking scenarios ---------------------------------------


# compare to all results among all studies
mod = glm( Di ~ yi,
           data = d,
           family = binomial )
summary(mod)
exp(coef(mod))

# compare to all nonaffirmative results among all studies
# **this removes the correlation!!!
# in hacked studies, this would be the same as the first regression above because unsuccessful studies ONLY have nonaffirmative results
# but when hack=0, this model includes nonaffirms from ultimately successful studies

# favored nonaffirm draws from studies that NEVER got an affirm
# + favored nonaffirm draws from studies that did get at least one, but not the favored one
mod = glm( Di ~ yi,
           data = dn,
           family = binomial )
summary(mod)
exp(coef(mod))

# counterpart to first model: nonaffirms among successful studies
mod = glm( Di ~ yi,
           data = dn[ dn$study.success == 1, ],
           family = binomial )
summary(mod)
exp(coef(mod))

# look at nonaffirms 

# nonaffirm draws from studies that NEVER got an affirm: 10th (nonaffirm) draw tends to be LARGER than those other nonaffirms
# nonaffirm draws from studies that DID get an affirm, but not the 10th one: 10th (nonaffirm) draw tends to be SMALLER than those other nonaffirms


as.data.frame( dn %>% group_by(study.success, draw.index==10) %>%
  summarise(mean(yi)) )

#*pattern holds even for the FIRST draw
as.data.frame( dn %>% group_by(study.success, draw.index==1) %>%
  summarise(mean(yi)) )

temp = d %>% group_by(study.success, affirm, draw.index) %>%
  summarise(mean(yi))

ggplot( data = temp,
        aes( x=draw.index,
            y=`mean(yi)`,
            color = affirm,
            lty = study.success) ) +
  geom_line() +
  ggtitle("All results from unhacked studies, Nmax=50, Mu = -2") +
  geom_hline(yintercept=0, color="gray") +
  theme_bw()



# ~ Condition on unsuccessful studies ---------------------------------------

dim(du)

# **the actual conservatism result
du %>% filter(Di == 1) %>%
  summarise(mean(yi))

# stationarity
# want mean(yi) to stay the same over draws - it does! :)
# trend in yi over successive draws
mod = lm( yi ~ draw.index,
          data = du )
summary(mod)

#**this doesn't make sense given the stationarity
mod = lm( yi ~ Di,
          data = du )
summary(mod)

# this should exactly agree with above by def'n of Di
# and does
mod = lm( yi ~ (draw.index == 10),
          data = du )
summary(mod)

table(du$Di, du$draw.index)


# ~ Compare to all studies ---------------------------------------
# obviously this is positive because this includes affirmative results
cov(d$Di, d$yi) 


# ~ Plot draws as a time series ---------------------------------------

# sample of studies
ind = sample( unique(du$study), 20)
dp = du %>% filter(study %in% ind)



ggplot( data = dp,
        aes( x = draw.index,
             y = yi,
            color = as.factor(Di) ) ) + 
  geom_point() + 
  geom_line() +
  # threshold for affirmative result, given that SEs are always 0
  geom_hline(yintercept = 0,
             color = "green") +
  theme_bw() +
  facet_wrap(~study) +
  ggtitle( paste( dataset.name,
                 ": draws in unsuccessful studies",
                 sep = "" ) )




# PART 1: PLOT TIME SERIES TO LOOK AT STATIONARITY ---------------------------------------

# Goal: Understand intuitively which kinds of time series are stationary (at least for their expectations)
# e.g., I don't understand why this holds for a moving-average process?

# Note about simts pkg: 
#  The error "Need to supply initial values within the ts.model object"
#  means you specified arguments that are out of possible range


# https://cran.r-project.org/web/packages/simts/vignettes/vignettes.html
library(simts)

n = 500


# White noise
x = gen_gts( n, WN(sigma2 = 1) )
mean(x)
plot(x)

# AR1, positive autocorrelation
# specification in package docs:
# X_t = φ X_{t-1} + \varepsilon_t
x = gen_gts( n, AR1(phi = -0.95, sigma2 = 1) )
mean(x); var(x)
plot(x)

# ARMA(1): special case with constant moving avg (should be just like AR(1))
# X_t = ∑_{j = 1}^p φ_j X_{t-j} + ∑_{j = 1}^q θ_j \varepsilon_{t-j} + \varepsilon_t
# "ar" argument is the phi's (provide only one to have an ARMA(1) process)
# "ma" is the theta's
x = gen_gts( n, ARMA(ar = .95, ma = 0, sigma2 = 1) )
mean(x)
plot(x)

# ARMA(1): special case with (almost) no autocorrelation
x = gen_gts( n, ARMA(ar = 0.01, ma = 0.9, sigma2 = 1) )
mean(x)
plot(x)
# @I don't understand why this doesn't increase its average over time? 


# ARMA(1): autocorrelated and moving average is positive
x = gen_gts( n, ARMA(ar = 0.95, ma = 0.9, sigma2 = 1) )
mean(x)
plot(x)















