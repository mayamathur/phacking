


# Goal: Try to get MAPs and MLEs for RTMA from TNE code itself


# ~ PRELIMINARIES ----------------------------------------------------

#rm(list=ls())


# data-wrangling packages
library(here)
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(xtable)
library(testthat)
library(Deriv)
library(mosaic)
library(hpa)
library(pracma)
library(truncnorm)
library(tmvtnorm)
library(Hmisc)
library(truncdist)
library(weightr)

options(scipen = 100)


prepped.data.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Prepped Hagger data (meta and replication)"
res.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Hagger comparison results"


setwd(prepped.data.dir)
dm = read.csv("prepped_hagger_meta_data.csv")

#code.dir = here("2021-10-7 check RTMA Jeffreys theory")

code.dir = here()
setwd(code.dir)
source("helper_SAPH.R")

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-11-23 repurpose TNE code")
#source("_helper.R")
source("fns_from_TNE.R")



# ~ SET PARAMETERS AND SIMULATE META-ANALYSIS ----------------------------------------------------

# Mu = 1
# T2t = 2
# m = 50
# se = 
# 
# Nmax = 1
# 
# d = sim_meta(Nmax = .p$Nmax,
#              Mu = .p$Mu,
#              T2 = .p$T2,
#              m = .p$m,
#              t2w = .p$t2w,
#              se = .p$se,
#              hack = .p$hack,
#              return.only.published = FALSE,
#              rho = .p$rho,
#              
#              k = .p$k,
#              k.hacked = .p$k.hacked )

# use Hagger replications, for which MLE was horrible

setwd(prepped.data.dir)
dm = read.csv("prepped_hagger_meta_data.csv")

# .obj = correct_meta_phack2(yi = dm$yi,
#                            vi = dm$vi)
# 
# # t-stat MLE is HUGE
# .obj$sanityChecks$tstatMeanMLE
# 
# # **plot shows that the MLE is so big because the nonaffirmative t-stats are actually left-skewed
# # very interesting
# plot_trunc_densities(.obj)


# ~ 2021-11-24: SIM META WITH SAME SE ----------------------------------------------------

# SUMMARY OF THIS SECTION:
# truth = (0.1, 0.5)

# at k=100:
# MLE = (0.21, 0.65) – too big 
# MAP = (0.31, 0.74) – too big

# unlike the sim meta below, now we allow different studies to have different SEs
# RESULTS DO NOT MAKE SENSE -- ESTIMATES ARE WAY TOO SMALL

Mu = .1
T2 = .25
t2w = 0.25
k = 1000
se = 0.5
m = 500

# from doParallel:
# scen.params = data.frame( scen = 1,
#                           Mu = 0.1,
#                           T2 = 0.25,
#                           t2w = 0.25,
#                           se = 0.5,
#                           m = 500,
#                           
#                           Nmax = 1,
#                           hack = "affirm", 
#                           rho = 0,
#                           
#                           k = 100,
#                           k.hacked = 0,
#                           
#                           get.CIs = TRUE)

# SAVE: this is how data were generated
d = sim_meta(Nmax = 1,
             Mu = Mu,
             T2 = T2,
             m = m,
             t2w = t2w,
             se = se,
             hack = "affirm",
             return.only.published = FALSE,
             rho = 0,

             k = k,
             k.hacked = 0 )


dn = d %>% filter(Di == 1 & affirm == FALSE)
dn$sei = sqrt(dn$vi)
# sanity check
all( dn$yi <= dn$tcrit * dn$sei )
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-10-13 numerical integration for RTMA Jeffreys prior")
fwrite(dn, "sim_meta_k1000.csv")

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-10-13 numerical integration for RTMA Jeffreys prior")
dn = fread("sim_meta_k1000.csv")



yi = dn$yi
sei = dn$sei
tcrit = dn$tcrit

Mu.start = 0
Tt.start = 1

### Version 1: MLE ###
res.MLE.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = FALSE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MLE.1$MuHat
res.MLE.1$TtHat


### Sanity check: optimize the nll directly

joint_nll2_simple = function(.Mu, .Tt) {
  joint_nll_2( .yi = yi, .sei = sei, .Mu = .Mu, .Tt = .Tt, .tcrit = tcrit )
}


res = mle( minuslogl = joint_nll2_simple,
           start = list( .Mu = Mu.start, .Tt = Tt.start),
           method = "Nelder-Mead" )
as.numeric(coef(res))
# yes, agrees with above :)

# vs. BFGS
res = mle( minuslogl = joint_nll2_simple,
           start = list( .Mu = Mu.start, .Tt = Tt.start),
           method = "BFGS" )
as.numeric(coef(res))
# yes, agrees with above :)



# does nlpost_jeffreys agree? YES
v1 = joint_nll_2( .yi = yi, .sei = sei, .Mu = Mu.start, .Tt = Tt.start, .tcrit = tcrit )
v2 = nlpost_jeffreys_RTMA( .pars = c(Mu.start, Tt.start),
                           .par2is = "Tt",
                           .yi = yi,
                           .sei = sei,
                           .tcrit = tcrit,
                           .usePrior = FALSE )
expect_equal(v1, v2)

### Version 2: MAP ###
res.MAP.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MAP.1$MuHat
res.MAP.1$TtHat


### MLE version 2: weightr ###

#  yes, agrees with MLE :)
( m1 = weightfunct( effect = yi,
                    v = sei^2,  
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate
res.MLE.1$MuHat; m1[[2]]$par[2]
# very close, but not exact because of t vs. z cutoff




# ~ 2021-11-23: WITH DIFFERENT SEs (HAGGER) ----------------------------------------------------

# SUMMARY OF THIS SECTION:
# MLE = 3.31 (seems very big, but agrees with weightr)
# MAP = 305!!! (horrible)

# scenario: different SEs as in Hagger meta, and we have to use the estimated SEs
#  rather than the true ones

zcrit = qnorm(.975)

dn = dm %>% filter( yi/sqrt(vi) < zcrit )

# set vars for all methods
yi = dn$yi
sei = sqrt(dn$vi)
kn = nrow(dn)
Mu.start = 0
Tt.start = 1

# sanity check 
all( dn$yi <= zcrit * sei )


### Version 1: MLE - agrees with weightr!! ###
res.MLE.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    tcrit = zcrit,
                                    
                                    usePrior = FALSE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MLE.1$MuHat
res.MLE.1$TtHat


### MLE version 2: weightr ###
# agrees with above :)
( m1 = weightfunct( effect = yi,
                    v = sei^2,  
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate
expect_equal( m1[[2]]$par[2], res.MLE.1$MuHat, tol = 0.001 ) 




### Sanity check: optimize the nll directly

joint_nll2_simple = function(.Mu, .Tt) {
  joint_nll_2( .yi = yi, .sei = sei, .Mu = .Mu, .Tt = .Tt, .tcrit = zcrit )
}


res = mle( minuslogl = joint_nll2_simple,
           start = list( .Mu = Mu.start, .Tt = Tt.start),
           method = "Nelder-Mead" )
as.numeric(coef(res))
# yes, agrees with above :)

# does nlpost_jeffreys agree? YES
v1 = joint_nll_2( .yi = yi, .sei = sei, .Mu = Mu.start, .Tt = Tt.start, .tcrit = zcrit )
v2 = nlpost_jeffreys_RTMA( .pars = c(Mu.start, Tt.start),
                      .par2is = "Tt",
                      .yi = yi,
                      .sei = sei,
                      .tcrit = zcrit,
                      .usePrior = FALSE )
expect_equal(v1, v2)

### Version 2: MAP (SD parametrization) ###
res.MAP.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    tcrit = zcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MAP.1$MuHat
res.MAP.1$TtHat
#bm: this one is huge!

### Version 3: MAP (var parametrization) ###
res.MAP.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "T2t",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    tcrit = zcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MAP.1$MuHat
res.MAP.1$TtHat





# ~ 2021-11-23: WITH ALL SES EQUAL ----------------------------------------------------

# SUMMARY OF THIS SECTION:
# MLE = 1.02 (excellent)
# MAP = 0.97 (also excellent)

# scenario: no hacking; all SEs equal
# this behaves reasonably with either MLE or MAP

Mu = 0.1
Tt = 0.5
k = 100
se = 1


# SAVE: this is how data were generated
d = sim_meta(Nmax = 1,
             Mu = Mu,
             T2 = Tt^2,
             m = 100,
             t2w = 0,
             se = se,
             hack = "affirm",
             return.only.published = FALSE,
             rho = 0,

             k = k,
             k.hacked = 0 )

# **get nonaffirms based on TRUE SEs
# temporary because otherwise they won't all have same se
tcrit = unique(dn$tcrit)
dn = d %>% filter( yi < tcrit * se )


# setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-10-13 numerical integration for RTMA Jeffreys prior")
# # fwrite(dn, "sim_meta_1.csv")
# dn = fread("sim_meta_1.csv")

# set vars for all methods
kn = nrow(dn)
Mu.start = 0
Tt.start = 1
sei.true = rep(se, kn)
tcrit = unique(dn$tcrit)
dn$sei = sqrt(dn$vi)

# sanity check 
max(dn$yi) <= tcrit * se

### Version 1: MLE from TNE ###
# because SEs are the same

# as in doParallel_TNE.R
p = data.frame(n = kn,
               a = -99,
               b = tcrit * se,
               stan.iter = 2000,  # default: 2000
               stan.maxtreedepth = 10, # default: 10
               stan.adapt_delta = 0.8,
               get.CIs = FALSE )

res.MLE.1 = estimate_mle(x = dn$yi,
                           p = p,
                           par2is = "var",  # NOTE: it prefers var parameterization
                           mu.start = Mu.start,
                           sigma.start = Tt.start,  # resulting estimate will INCLUDE se
                           get.CIs = p$get.CIs,
                           CI.method = "wald")

res.MLE.1$Mhat
res.MLE.1$Shat  # this includes se
( res.MLE.1$TtHat = sqrt( res.MLE.1$Shat^2 - se^2 ) )
# reasonable


### Version 1.1: Compare to MLE from estimate_jeffreys_RTMA

# should agree with above
# yes, very similar :)
res.MLE.2 = estimate_jeffreys_RTMA( yi = dn$yi,
                                      sei = se,
                                      par2is = "Tt",
                                      Mu.start = Mu.start,
                                      Tt.start = Tt.start,
                                      tcrit = tcrit,

                                      usePrior = FALSE,
                                      get.CIs = FALSE,
                                      CI.method = "wald" )

res.MLE.2$MuHat; res.MLE.1$Mhat
res.MLE.2$TtHat; res.MLE.1$TtHat



### Version 2: MAP from TNE ###
# makes sense
res.MAP.TNE = estimate_jeffreys( x = dn$yi,
                                 p = p,
                                 par2is = "sd",  # NOTE: it prefers var parameterization
                                 mu.start = Mu.start,
                                 sigma.start = Tt.start,  # resulting estimate will INCLUDE se
                                 get.CIs = p$get.CIs,
                                 CI.method = "wald" )

res.MAP.TNE$Mhat
res.MAP.TNE$Shat
( res.MAP.TNE$TtHat = sqrt( res.MAP.TNE$Shat^2 - se^2 ) )

# CI
res.MAP.TNE.CI = estimate_jeffreys( x = dn$yi,
                                    p = p,
                                    mu.start = Mu.start,
                                    sigma.start = Tt.start,
                                    get.CIs = TRUE)

res.MAP.TNE.CI$M.CI
res.MAP.TNE.CI$S.CI


### Version 3: MAP from SAPH code ###
res.MAP.2 = estimate_jeffreys_RTMA( yi = dn$yi,
                                    sei = se,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = FALSE,
                                    CI.method = "wald" )

res.MAP.2$MuHat; res.MAP.TNE$Mhat
res.MAP.2$TtHat; res.MAP.TNE$TtHat
# very close, but not exactly the same








