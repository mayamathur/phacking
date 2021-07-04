
# Goal: I think that if I analyze only the *affirmatives*, the weightr package and LTN should coincide if I tell weightr that eta=infinity and all studies have same SE (because LTN doesn’t yet give us a way to handle different SEs). See if this is true. 

# Not concerned here with whether the MLEs are actually correct, but rather whether the selection model theory behaves
#  as I expect.

# Overall purpose is to see when selection models will work under p-hacking. 

# See project log for conclusions

# PRELIMINARIES ------------------------------

# borrow helper fns from previous sims
library(here)
setwd( here() )
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
library(weightr)
# other
library(xtable)
library(testthat)
# for this project
library(truncdist)
library(tmvtnorm)

results.dir = here( "2021-7-4 connection to selection models/Numbered simulation datasets" )

code.dir = here( "2021-7-4 connection to selection models" )
setwd(code.dir)
source("_helper.R")


# EXPT 1 ------------------------------

x = quick_sim( .p = data.frame( Mu = 0.1,
                                T2 = 0.25,
                                m = 500,
                                t2w = .25,
                                se = .5,

                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,

                                k = 200,
                                k.hacked = 200,


                                sim.name = "expt_1_res" ),
               .results.dir = results.dir )

setwd(data.dir)
load("expt_1_res")

# keep only affirmatives
d = x$d
d = d %>% filter( Di == 1 )
# should always be true
table(d$affirm)
# percent successful
sum(d$affirm)/200


# ~ Fit weightr with known weights -----------------------------

# ~~ Warmups -----

# warmup: uncorrected meta-analysis
m0 = rma.uni( yi = d$yi,
              vi = d$vi,
              method = "REML" )
m0$b

# warmup: don't specify weights
( m1 = weightfunct( effect = d$yi,
                    v = d$vi,
                    steps = c(0.025, 1),
                    #weights = c(1,0),
                    table = TRUE ) )

# "The order of parameters is as follows: variance component, mean or linear coefficients,
#  and weights"
m1[[2]]$par


# warmup: now specify that there's no pub bias
( m1 = weightfunct( effect = d$yi,
                    v = d$vi,
                    steps = c(0.025, 1),
                    weights = c(1,1),
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]

# sanity check: should match naive
expect_equal( as.numeric(m1[[2]]$par[2]), as.numeric(m0$b))



# ~~ For Real Now -----

### Method 1: weightr ###
# specify eta -> infinity
( m1 = weightfunct( effect = d$yi,
                    v = d$viTrue,  # use true variances as in LTN below
                    steps = c(0.025, 1),
                    weights = c(1,0),
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]


### Method 2: Compare to LTN on t-stats (rescaled) ###
crit = qnorm(.975)
mle.fit = mle.tmvnorm( X = as.matrix(d$tstat, ncol = 1),
                       lower = crit,
                       upper = Inf)
#summary(mle.fit)
mles = coef(mle.fit)


# APPROXIMATELY rescale MLEs to represent effect sizes rather than tstats
( Mhat = mles[1] * d$se[1] )

# this is NOT the same as weightr




### Method 3: Using my implementation of LTN incorporating SEs ###

# use the true weights
# MLEs should now be close to the true mu and V  
library(bbmle)


bbmle::mle2( minuslogl = function(.muhat, .t2) {
  minusLL_LTN(
    yi = d$yi,
    vi = d$vi,
    muhat = .muhat, 
    t2 = .t2 ) },
  start = list(.muhat = 0, .t2=0),
  method = "L-BFGS-B",
  lower = list(.muhat = -Inf, .t2=0) )
# MATCHES weightr :) 
#**So this confirms that my version of LTN that incorporates the SEs
#  coincides exactly with weightr if specifying eta=infty




# EXPT 2 ------------------------------

# now try an example where LTN should actually be unbiased: no heterogeneity at all

x = quick_sim( .p = data.frame( Mu = 0.1,
                                T2 = 0,
                                m = 500,
                                t2w = 0,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 5000,
                                k.hacked = 5000,
                                
                                
                                sim.name = "expt_2_res" ),
               .results.dir = results.dir )

setwd(data.dir)
load("expt_2_res")

# keep only affirmatives
d = x$d
d = d %>% filter( Di == 1 )
# should always be true
table(d$affirm)
# percent successful
sum(d$affirm)/5000


### Method 1: Weightr ###
# specify eta -> infinity
( m1 = weightfunct( effect = d$yi,
                    v = d$vi,  
                    steps = c(0.025, 1),
                    weights = c(1,0),
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]


### Method 3: Using my implementation of LTN incorporating SEs ###

# use the true weights
# MLEs should now be close to the true mu and V  
library(bbmle)


bbmle::mle2( minuslogl = function(.muhat, .t2) {
  minusLL_LTN(
    yi = d$yi,
    vi = d$vi,
    muhat = .muhat, 
    t2 = .t2 ) },
  start = list(.muhat = 0, .t2=0),
  method = "L-BFGS-B",
  lower = list(.muhat = -Inf, .t2=0) )


### Method 2: Compare to LTN on t-stats (rescaled) ###
crit = qnorm(.975)
mle.fit = mle.tmvnorm( X = as.matrix(d$tstat, ncol = 1),
                       lower = crit,
                       upper = Inf)
#summary(mle.fit)
mles = coef(mle.fit)


# APPROXIMATELY rescale MLEs to represent effect sizes rather than tstats
( Mhat = mles[1] * d$se[1] )


#bm: This example is exactly like row 2 in project log (unbiased using LTN), 
#  where I showed that the t-stats ARE LTN. But the MLEs aren't correct here.






# EXPT 3 ------------------------------

# now try an example where LTN should actually be unbiased: t2w>0, but rho=0

x = quick_sim( .p = data.frame( Mu = 0.1,
                                T2 = 0,
                                m = 500,
                                t2w = 0.25,
                                se = .5,
                                
                                Nmax = 10,
                                hack = "affirm",
                                rho = 0,
                                
                                k = 5000,
                                k.hacked = 5000,
                                
                                
                                sim.name = "expt_3_res" ),
               .results.dir = results.dir )

setwd(data.dir)
load("expt_3_res")

# keep only affirmatives
d = x$d
d = d %>% filter( Di == 1 )
# should always be true
table(d$affirm)
# percent successful
sum(d$affirm)/5000


### Method 1: Weightr ###
# specify eta -> infinity
( m1 = weightfunct( effect = d$yi,
                    v = d$vi,  
                    steps = c(0.025, 1),
                    weights = c(1,0),
                    table = TRUE ) )

# adjusted point estimate
m1[[2]]$par[2]


### Method 3: Using my implementation of LTN incorporating SEs ###

# use the true weights
# MLEs should now be close to the true mu and V  
library(bbmle)


bbmle::mle2( minuslogl = function(.muhat, .t2) {
  minusLL_LTN(
    yi = d$yi,
    vi = d$vi,
    muhat = .muhat, 
    t2 = .t2 ) },
  start = list(.muhat = 0, .t2=0),
  method = "L-BFGS-B",
  lower = list(.muhat = -Inf, .t2=0) )


### Method 2: Compare to LTN on t-stats (rescaled) ###
crit = qnorm(.975)
mle.fit = mle.tmvnorm( X = as.matrix(d$tstat, ncol = 1),
                       lower = crit,
                       upper = Inf)
#summary(mle.fit)
mles = coef(mle.fit)


# APPROXIMATELY rescale MLEs to represent effect sizes rather than tstats
( Mhat = mles[1] * d$se[1] )


#bm: was waiting for this to simulate :)

