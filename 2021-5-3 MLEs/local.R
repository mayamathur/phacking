
# This script is all about the NONaffirmatives.

# PRELIMINARIES ------------------------------

# borrow helper fns from previous sims
library(here)
setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
source("helper_SAPH.R")

# data-wrangling packages
library(dplyr)
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
# for this project
library(truncdist)
#library(ExtDist)
library(gmm)  # https://stackoverflow.com/questions/63511986/error-package-or-namespace-load-failed-for-gmm-in-dyn-loadfile-dllpath-dl
library(tmvtnorm)

results.dir = here("2021-5-3 MLEs/Numbered simulation experiments")
setwd(results.dir)




# LOOK AT GETTING MLES FROM NONAFFIRMATIVES -------------------------

# ~ EXPT 1 -------------------------

# all unhacked with Nmax = 10, so I expect the published nonaffirms
#  to exactly come from the trunc dist
p = data.frame( Mu = 1,
                T2 = 0.25,
                m = 500,
                t2w = 0.25,
                se = .5,

                Nmax = 10,
                hack = "affirm",

                k = 10000,
                k.hacked = 0,


                sim.name = "expt_1" )

# this is just to get the dataset, even though it also does some analyses
# x = quick_sim( .p = p,
#                .results.dir = results.dir )

# read back in
load("expt_1")
d = x$d

crit = unique(d$tcrit)


# published only (i.e., last draw of each set)
dp = d %>% filter(Di == 1)
expect_equal( 10000, length(unique(dp$study)) )



# ~~ Check that moments of published nonaffirms are as expected  -------------------------
tstats = d$tstat[ d$Di == 1 & d$affirm == 0 ]


# check that moments are what we expect
extrunc(spec = "norm",
       mean = p$Mu / p$se,
       #@doesn't yet use the delta-method thing
       sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
       b = crit )
mean(tstats)  # yes, matches :)

vartrunc(spec = "norm",
        mean = p$Mu / p$se,
        #@doesn't yet use the delta-method thing
        sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
        b = crit )
var(tstats)  # pretty close :)


# ~~ Get MLEs for Mu and T2 from the published nonaffirmatives  -------------------------

# has convergence issue, which might be why the MLEs are off
nonaffirm_mles( .tstats = tstats, 
                .t2w = t2w, 
                .se = se,
                .crit = crit )
unique(d$Mu); unique(d$T2)

# compare to LL evaluated at true parameters
-joint_ll( .tstats = tstats,
           .Mu = unique(d$Mu),
           .T2 = unique(d$T2),
           .t2w = t2w,
           .se = se,
           .crit = crit )

               
# ~~ Compare to using data drawn directly from trunc dist -------------------------

# **still has convergence issue
x2 = rtrunc(n = 5000,
           spec = "norm",
           mean = p$Mu / p$se,
           sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
           b = crit)

nonaffirm_mles( .tstats = x2, 
                .t2w = t2w, 
                .se = se,
                .crit = crit )
Mu; T2

# compare to LL evaluated at true parameters
-joint_ll( .tstats = tstats,
           .Mu = p$Mu,
           .T2 = p$T2,
           .t2w = t2w,
           .se = se,
           .crit = unique(d$tcrit) )

# ~~ Try a different package for truncated MLEs -------------------------

# note that the variance component will include t2w + se^2 here
# and mean will be Mu/se


### With data straight from trunc dist
# fn wants data given as a matrix
mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
summary(mle.fit)

mles = coef(mle.fit)
# pretty good
mles[1]; p$Mu/p$se
mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# # doesn't work
# also doesn't work in package example
# mle.profile1 <- profile(mle.fit, x2, method="BFGS", trace=TRUE)
# confint(mle.profile1)
# 
# par(mfrow=c(3,2))
# plot(mle.profile1)

### With real data
# **THIS SEEMS TO WORK PRETTY WELL
# try with my real data
mle.fit = mle.tmvnorm( as.matrix(tstats, ncol = 1), lower=-Inf, upper=crit)
summary(mle.fit)

mles = coef(mle.fit)
# pretty good :)
mles[1]; p$Mu/p$se
mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)



# ~ EXPT 2 -------------------------

# all unhacked with Nmax = 10 (as before), but now draws are correlated
p = data.frame( Mu = 1,
                T2 = 0.25,
                m = 500,
                t2w = 0.25,
                se = .5,
                
                Nmax = 10,
                hack = "affirm",
                
                
                k = 10000,
                k.hacked = 0,
                rho = 0.9,
                
                sim.name = "expt_2" )

# simulate from scratch
x = quick_sim( .p = p,
                .results.dir = results.dir )
d = x$d


# # read back in
# load("expt_2")
# d = x$d

#  look at autocorrelation of muin's
#  but note this will be less than rho because of 
#  small-sample bias of autocorrelation
d %>% filter( !duplicated(study) ) %>%
  summarise( sum(!is.na(rhoEmp)),
             sum(N > 1),
             mean(rhoEmp, na.rm = TRUE) )

crit = unique(d$tcrit)



# ~~ Check that moments of published nonaffirms are as expected  -------------------------
tstats = d$tstat[ d$Di == 1 & d$affirm == 0 ]

# check that moments are what we expect
extrunc(spec = "norm",
        mean = p$Mu / p$se,
        #@doesn't yet use the delta-method thing
        sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
        b = crit )
mean(tstats)  # yes, matches :)

vartrunc(spec = "norm",
         mean = p$Mu / p$se,
         #@doesn't yet use the delta-method thing
         sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
         b = crit )
var(tstats)  # pretty close :)


# ~~ MLEs  -------------------------

mle.fit = mle.tmvnorm( as.matrix(tstats, ncol = 1), lower=-Inf, upper=crit)
summary(mle.fit)

mles = coef(mle.fit)
# pretty good :)
mles[1]; p$Mu/p$se
mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)


           
               