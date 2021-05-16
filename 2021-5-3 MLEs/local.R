
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
                .t2w = p$t2w, 
                .se = p$se,
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
# MLE = 2.04 vs. theory = 2
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
# MLE = 2.12 vs. theory = 2 (seems biased upward?)
mle.fit = mle.tmvnorm( as.matrix(tstats, ncol = 1), lower=-Inf, upper=crit)
summary(mle.fit)

mles = coef(mle.fit)
# pretty good :)
mles[1]; p$Mu/p$se
mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

### Visually compare distribution of meta-analysis data vs. random generates straight from distribution

ggplot() + 
  geom_density( aes(x = x2),
                adjust = 0.3 ) + 
  
  geom_density( aes(x = tstats),
                adjust = 0.3,
                color = "red") + 
  theme_bw()

# looks damn similar to me...



# ~~ Bigger simulation with random generates straight from dist -----------------

sim.reps = 500

meanHat = c()
varHat = c()


for ( i in 1:sim.reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  x2 = rtrunc(n = 5000,
              spec = "norm",
              mean = p$Mu / p$se,
              sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
              b = crit)
  
  mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat = c(meanHat, mles[1])
  varHat = c(meanHat, mles[2])
}


# **THIS LOOKS JUST FINE
mean(meanHat); p$Mu/p$se
mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat)



# ~~ Very similar but using sim_meta -----------------

#bm

meanHat = c()
varHat = c()
# all tstats from all sims
tstats = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  dp = sim_meta(Nmax = p$Nmax,
                Mu = p$Mu,
                T2 = p$T2,
                
                
                m = p$m,
                t2w = p$t2w,
                se = p$se, 
                
                rho = 0,
                
                hack = p$hack,
                
                k = p$k,
                k.hacked = p$k.hacked,
                return.only.published = TRUE)
  
  dpn = dp[ dp$affirm == FALSE, ]
  
  mle.fit = mle.tmvnorm( as.matrix(dpn$tstats, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat = c(meanHat, mles[1])
  varHat = c(meanHat, mles[2])
  tstats = c(tstats, dpn$tstats)
}



mean(meanHat); p$Mu/p$se
mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat)



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


           
               


# EXPT 3.1: k=6000 ------------------------------------------

# Goal: Find out why the MLEs are biased by comparing to MLEs when data come
#  straight from distribution.
# Conclusion: Here the MLEs are PERFECT even when using data from sim_meta: 
#  could it be the very large k?

# set parameters
# exactly as in row 20 of "agg_rounded_annotated", where the tstat MLE was 0.24 instead of 0.2
p = data.frame( Mu = 0.1,
                T2 = 0.25,
                m = 500,
                t2w = 0.25,
                se = .5,
                
                Nmax = 1,
                hack = "affirm",
                
                k = 6000,  # to get about 5000 nonaffirmatives
                k.hacked = 0 )


# ~~ Random generates straight from dist -----------------

sim.reps = 500

meanHat = c()
varHat = c()
tstats = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  x2 = rtrunc(n = 5000,
              spec = "norm",
              mean = p$Mu / p$se,
              sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
              b = crit)
  
  mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat = c(meanHat, mles[1])
  varHat = c(varHat, mles[2])
  tstats = c(tstats, x2)
}


# **THIS LOOKS JUST FINE
mean(meanHat); p$Mu/p$se
mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat)

# trunc normal
hist(tstats)


# ~~ Very similar but using sim_meta -----------------

# fewer sim reps because takes longer
sim.reps = 250

meanHat2 = c()
varHat2 = c()
# all tstats from all sims
tstats2 = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 10 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  dp = sim_meta(Nmax = p$Nmax,
                Mu = p$Mu,
                T2 = p$T2,
                
                
                m = p$m,
                t2w = p$t2w,
                se = p$se, 
                
                rho = 0,
                
                hack = p$hack,
                
                k = p$k,
                k.hacked = p$k.hacked,
                return.only.published = TRUE)
  
  dpn = dp[ dp$affirm == FALSE, ]
  
  mle.fit = mle.tmvnorm( as.matrix(dpn$tstat, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat2 = c(meanHat2, mles[1])
  varHat2 = c(varHat2, mles[2])
  tstats2 = c(tstats2, dpn$tstat)
}

# THIS ONE IS ALSO FINE?????
# FASCINATING
# both are exactly right
mean(meanHat2); p$Mu/p$se
mean(varHat2); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat2)





# ~~ Visually compare tstats' dist -----------------

ggplot() +
  # straight from dist
  geom_density( aes(x = tstats),
                adjust = 0.3,
                lty = 2) +

  # from sim_meta
  geom_density( aes(x = tstats2),
                adjust = 0.3,
                color = "red",
                lty = 2) +
  theme_bw()


# EXPT 3.2: k=100 as on Sherlock ------------------------------------------

# Goal: Find out if lowering k explains discrepancy between Expt 3 and the Sherlock
#  results in agg_rounded_annotated. 
#  Here I'm reducing k to 100 to exactly match the Sherlock sims. 

# Now I get meanHat2 = 0.210, varHat2 = 2.96

# set parameters
# exactly as in row 20 of "agg_rounded_annotated", where the tstat MLE was 0.24 instead of 0.2
p = data.frame( Mu = 0.1,
                T2 = 0.25,
                m = 500,
                t2w = 0.25,
                se = .5,
                
                Nmax = 1,
                hack = "affirm",
                
                k = 100, 
                k.hacked = 0 )


# fewer sim reps because takes longer
sim.reps = 250

meanHat2 = c()
varHat2 = c()
# all tstats from all sims
tstats2 = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 10 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  dp = sim_meta(Nmax = p$Nmax,
                Mu = p$Mu,
                T2 = p$T2,
                
                
                m = p$m,
                t2w = p$t2w,
                se = p$se, 
                
                rho = 0,
                
                hack = p$hack,
                
                k = p$k,
                k.hacked = p$k.hacked,
                return.only.published = TRUE)
  
  dpn = dp[ dp$affirm == FALSE, ]
  
  mle.fit = mle.tmvnorm( as.matrix(dpn$tstat, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat2 = c(meanHat2, mles[1])
  varHat2 = c(varHat2, mles[2])
  tstats2 = c(tstats2, dpn$tstat)
}

# THIS ONE IS ALSO FINE?????
# FASCINATING
# both are exactly right
mean(meanHat2); p$Mu/p$se
mean(varHat2); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat2)




# EXPT 3.3: k = 20 ------------------------------------------

# Reduce k even more. Now the MLEs are horrible! 

# Now I get meanHat2 = 0.4829075, varHat2 = 3.433637

# set parameters
# exactly as in row 20 of "agg_rounded_annotated", where the tstat MLE was 0.24 instead of 0.2
p = data.frame( Mu = 0.1,
                T2 = 0.25,
                m = 500,
                t2w = 0.25,
                se = .5,
                
                Nmax = 1,
                hack = "affirm",
                
                k = 20, 
                k.hacked = 0 )


sim.reps = 1000


# ~~ From sim_meta -----------------

meanHat2 = c()
varHat2 = c()
# all tstats from all sims
tstats2 = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 10 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  dp = sim_meta(Nmax = p$Nmax,
                Mu = p$Mu,
                T2 = p$T2,
                
                
                m = p$m,
                t2w = p$t2w,
                se = p$se, 
                
                rho = 0,
                
                hack = p$hack,
                
                k = p$k,
                k.hacked = p$k.hacked,
                return.only.published = TRUE)
  
  dpn = dp[ dp$affirm == FALSE, ]
  
  mle.fit = mle.tmvnorm( as.matrix(dpn$tstat, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat2 = c(meanHat2, mles[1])
  varHat2 = c(varHat2, mles[2])
  tstats2 = c(tstats2, dpn$tstat)
}

# VERY, very wrong
# both are exactly right
mean(meanHat2); p$Mu/p$se
mean(varHat2); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat2)


ggplot() + 
  geom_density( aes(x = x2),
                adjust = 0.3 ) + 
  
  geom_density( aes(x = tstats),
                adjust = 0.3,
                color = "red") + 
  theme_bw()



# ~~ Random generates straight from dist -----------------

# this is also very biased!
# meanHat = 0.4442739
# varHat = 3.349863

sim.reps = 1000

meanHat = c()
varHat = c()
tstats = c()

for ( i in 1:sim.reps ) {
  
  if ( i %% 25 == 0 ) cat( paste("\n\nStarting rep", i ) )
  
  x2 = rtrunc(n = 20,
              spec = "norm",
              mean = p$Mu / p$se,
              sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
              b = crit)
  
  mle.fit = mle.tmvnorm( as.matrix(x2, ncol = 1), lower=-Inf, upper=crit)
  mles = coef(mle.fit)
  
  meanHat = c(meanHat, mles[1])
  varHat = c(varHat, mles[2])
  tstats = c(tstats, x2)
}


mean(meanHat); p$Mu/p$se
mean(varHat); (1/p$se^2) * (p$T2 + p$t2w + p$se^2)

# should be normal
hist(meanHat)

# trunc normal
hist(tstats)
