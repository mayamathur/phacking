


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


#code.dir = here("2021-10-7 check RTMA Jeffreys theory")

# code.dir = here()
# setwd(code.dir)
# source("helper_SAPH.R")

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-11-23 repurpose TNE code")
source("_helper.R")
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




# ~ 2021-11-23: WITH DIFFERENT SEs ----------------------------------------------------

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


### Version 1: MLE - makes sense! ###
res.MLE.1 = estimate_jeffreys_RTMA( yi = yi,
                                      sei = sei,
                                      par2is = "Tt",
                                      Mu.start = Mu.start,
                                      Tt.start = Tt.start,
                                      crit = zcrit,
                                      
                                      usePrior = FALSE,
                                      get.CIs = TRUE,
                                      CI.method = "wald" )

res.MLE.1$MuHat
res.MLE.1$TtHat


### Sanity check: optimize the nll directly

joint_nll2_simple = function(.Mu, .Tt) {
  joint_nll_2( .yi = yi, .sei = sei, .Mu = .Mu, .Tt = .Tt, .crit = zcrit )
}


res = mle( minuslogl = joint_nll2_simple,
           start = list( .Mu = Mu.start, .Tt = Tt.start),
           method = "Nelder-Mead" )
as.numeric(coef(res))
# yes, agrees with above :)

# does nlpost_jeffreys agree? YES
v1 = joint_nll_2( .yi = yi, .sei = sei, .Mu = Mu.start, .Tt = Tt.start, .crit = zcrit )
v2 = nlpost_jeffreys_RTMA( .pars = c(Mu.start, Tt.start),
                      .par2is = "Tt",
                      .yi = yi,
                      .sei = sei,
                      .crit = zcrit,
                      .usePrior = FALSE )
expect_equal(v1, v2)

### Version 2: MAP ###
res.MAP.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    Tt.start = Tt.start,
                                    crit = zcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MAP.1$MuHat
res.MAP.1$TtHat

# point estimates are similar (MAP vs. MLE) and they make sense! 


### MLE version 2: weightr ###

# this one does NOT agree
( m1 = weightfunct( effect = yi,
                    v = sei^2,  
                    steps = c(0.025, 1),
                    weights = c(0,1), # weight such that ONLY nonaffirmatives are published
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]


#bm


# ~ 2021-11-23: WITH ALL SES EQUAL ----------------------------------------------------

# scenario: no hacking; all SEs equal
# this behaves reasonably with either MLE or MAP

Mu = 1
Tt = 2
k = 1000
se = 1


# # SAVE: this is how data were generated
# d = sim_meta(Nmax = 1,
#              Mu = Mu,
#              T2 = Tt^2,
#              m = 100,
#              t2w = 0,
#              se = se,
#              hack = "affirm",
#              return.only.published = FALSE,
#              rho = 0,
#              
#              k = k,
#              k.hacked = 0 )
# 
# # **get nonaffirms based on TRUE SEs
# # temporary because otherwise they won't all have same se
# tcrit = unique(dn$tcrit)
# dn = d %>% filter( yi < tcrit * se )
# 
# 
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-10-13 numerical integration for RTMA Jeffreys prior")
# fwrite(dn, "sim_meta_1.csv")
dn = fread("sim_meta_1.csv")

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
res.MLE.1$Shat
# reasonable


### Version 1.1: Compare to estimate_jeffreys_RTMA

# should agree with above, I think, BUT DOESN'T
res.MLE.2 = estimate_jeffreys_RTMA( yi = dn$yi,
                                      sei = se,
                                      par2is = "Tt",
                                      Mu.start = Mu.start,
                                      Tt.start = Tt.start,
                                      crit = tcrit,

                                      usePrior = FALSE,
                                      get.CIs = FALSE,
                                      CI.method = "wald" )

res.MLE.2$Mhat
res.MLE.2$Shat


# do the nlposts (actually nll in this case) for a single observation agree?
# runs joint_nll_2 under the hood
nlpost_jeffreys_RTMA( .pars = c(Mu.start, Tt.start),
                      .par2is = "Tt",
                      .yi = yi[1],
                      .sei = sei[1],
                      .crit = tcrit,
                      .usePrior = FALSE )

nll(.pars = c(Mu.start, Tt.start),
    .x = yi,
    .a = -99,
    .b = sei[1] * tcrit,
    par2is = "sd")



### Version 2: MAP from TNE ###

res.MAP.TNE = estimate_jeffreys( x = dn$yi,
                                 p = p,
                                 par2is = "sd",  # NOTE: it prefers var parameterization
                                 mu.start = Mu.start,
                                 sigma.start = Tt.start,  # resulting estimate will INCLUDE se
                                 get.CIs = p$get.CIs,
                                 CI.method = "wald" )

res.MAP.TNE$Mhat
res.MAP.TNE$Shat

# CI
res.MAP.TNE.CI = estimate_jeffreys( x = dn$yi,
                                    p = p,
                                    mu.start = Mu.start,
                                    sigma.start = Tt.start,
                                    get.CIs = TRUE)

res.MAP.TNE.CI$M.CI
res.MAP.TNE.CI$S.CI



