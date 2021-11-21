
# Conclusions:
# - Note that double.truncation:::NPMLE needs trunc points given as vectors
#   (i.e., it won't recycle)
#   o.w. will say "Error in if (Error < epsilon) { : missing value where TRUE/FALSE needed"

# PRELIMINARIES -----------------------------

library(here)
setwd(here())
source("helper_SAPH.R")


library(metafor)
library(weightr)
library(ggplot2)
library(dplyr)
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

# other
library(devtools)
library(double.truncation)
#devtools::install_github('jtleek/tidypvals')


# PLAIN TRUNCATED NORMAL -----------------------------

#bm: Try simulating data from plain normal dist, then truncate in different places
# see if you can recover the density
# if you end up extending Efron's thing to meta-analysis, it might look at lot like the calibrated estimates

# the more you push the b limits close to +Infinity, the better the MLE for mu becomes
# but it badly underestimates mu if the b limits aren't close to +Infinity



set.seed(1)
n = 1000
mu = 1
sd = 2

a = rep(-99, n)

# moderate truncation
# retains at least 50% of the distribution
#  => mean estimate is a dramatic underestimate
b = runif(n = n,
          min = mu,
          max = mu + 1*sd)  # causes bad behavior

# # minimal truncation
# # => mean estimate is still an underestimate, but is better
# b = runif(n = n,
#           min = mu + 2*sd,
#           max = mu + 3*sd)
# 
# # no truncation
# #  => mean estimate seems unbiased
# b = rep(99, n)

dt = data.frame(a = a,
                b = b)

dt = dt %>% rowwise() %>%
  mutate( x = rtrunc( n = 1,
                      spec = "norm",
                      mean = mu,
                      sd = sd,
                      a = a,
                      b = b ) )
hist(dt$x)
# sanity check
expect_equal( TRUE, all( dt$x < dt$b ) )

res = NPMLE( u.trunc = dt$a,
             y.trunc = dt$x,
             v.trunc = dt$b )


# NPMLE for mu
( Mhat = sum( sort(dt$x) * res$f ) )



# if you simulate from untruncated normal (a = -99, b = 99), 
#  f is same for all observations 
# but does NOT seem to work if data are truncated??
# plot( sort(dt$x), res$f )  
# plot( sort(dt$x), res$F )


# SIMULATED META-ANALYSIS -----------------------------

# no hacking
d = sim_meta(Nmax = 1,  
             Mu = 0,  
             T2 = 2,  
             
             # study parameters, assumed same for all studies:
             m = 50,  
             t2w = 0.25,  
             se = 0.5,  
             
             rho = 0,  
             
             hack = "affirm",  
             
             k = 1000,  
             k.hacked = 0 )

dp = d[ d$Di == TRUE, ]
table(dp$affirm)

dn = dp %>% filter(affirm == FALSE)



# ~ Try Efron nonparametric method  ---------------------

# lower trunc limit for each yi based on its own SE
dn$yi.cutoff = dn$tcrit * sqrt(dn$vi)

# small example for debugging
#dn = dn[1:10,]

library(double.truncation)
res = NPMLE(u.trunc = rep(-99, length(dn$yi)),
      y.trunc = dn$yi,
      v.trunc = dn$yi.cutoff )

plot( sort(dn$yi), res$f )

# package example from Efron paper
y.trunc=c(0.75, 1.25, 1.50, 1.05, 2.40, 2.50, 2.25)
u.trunc=c(0.4, 0.8, 0.0, 0.3, 1.1, 2.3, 1.3)
v.trunc=c(2.0, 1.8, 2.3, 1.4, 3.0, 3.4, 2.6)
res = NPMLE(u.trunc,y.trunc,v.trunc)
plot( sort(y.trunc), res$f )

NPMLE( u.trunc, y.trunc, rep(99, length(y.trunc) ) )



