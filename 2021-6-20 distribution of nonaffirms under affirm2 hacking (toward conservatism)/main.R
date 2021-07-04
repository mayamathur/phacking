
# PRELIMINARIES -----------------------------

library(here)

# general helpers
setwd(here())
source("helper_SAPH.R")

# helpers for the theory herein
setwd(here("2021-5-18 analytic nonaffirm MLE bias corrections"))
source("2021-5-18 helper.R")

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
library(Hmisc)
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
library(TruncatedNormal)
library(mle.tools)
library(fitdistrplus)
library(Deriv)

options(scipen=99)


# SIMULATE -----------------------------

# want to look at nonaffirmatives from hacked vs. unhacked studies


# same parameters as in Sherlock sims
Mu = 1
crit = qnorm(.975)
k = 10^4
k.hacked = k/2
d1 = sim_meta( Nmax = 10,
              Mu = Mu,
              T2 = 0.25,
              m = 500,
              t2w = 0.25,
              se = 0.5,
              hack = "affirm2",
              rho = 0,
              
              k = k,
              k.hacked = k.hacked,
              return.only.published = FALSE )

# keep only published nonaffirmatives
d = d1 %>% filter( Di == 1 &
                     affirm == FALSE )

# how many of each type?
table(d$hack)

# proportion of hacked studies that were unsuccessful 
# when this proportion is high, hacked nonaffirms should closely resemble unhacked ones 
#  because there's little selection on mui
( unsuccessful =  round( 100 * sum(d$hack == "affirm2") / k.hacked ) )

# distribution of t-stats in published nonaffirmatives
string = paste("Mu = ", Mu, "; ", unsuccessful, "% of hacked studies were unsuccessful", sep = "")
ggplot( data = d, 
        aes(x = tstat,
            color = hack) ) + 
  geom_density() + 
  theme_bw() + 
  ggtitle(string)


# distribution of mui in published nonaffirmatives
ggplot( data = d, 
        aes(x = mui,
            color = hack) ) + 
  geom_density() + 
  
  theme_bw() + 
  geom_vline(xintercept = Mu) +
  ggtitle(string)


# try fitting RTMA via MLE (unhacked nonaffirms only)
.obj = correct_meta_phack2(yi = d$yi[ d$hack == "no" ],
                           vi = d$vi[ d$hack == "no" ],
                           crit = qnorm(.975))

( MLE = round(.obj$sanityChecks$tstatMeanMLE, 2 ) )  # now matches the mean nicely
p = plot_trunc_densities(.obj,
                     showAffirms = FALSE)

string = paste("Mu = ", Mu, "; Muhat = ", MLE, "; unhacked nonaffirms", sep ="")
p + ggtitle(string)


# QQ plot vs. fitted truncated normal (should match)
qqtrunc( d$tstat[ d$hack == "no" ],
         spec = "norm",
         a = -Inf,
         mean = .obj$sanityChecks$tstatMeanMLE,
         sd = sqrt(.obj$sanityChecks$tstatVarMLE),
         b = crit )




# try fitting RTMA via MLE (HACKED nonaffirms only)
.obj = correct_meta_phack2(yi = d$yi[ d$hack == "affirm2" ],
                           vi = d$vi[ d$hack == "affirm2" ],
                           crit = qnorm(.975))


( MLE = round(.obj$sanityChecks$tstatMeanMLE, 2 ) )  # now matches the mean nicely
p = plot_trunc_densities(.obj,
                         showAffirms = FALSE)

string = paste("Mu = ", Mu, "; Muhat = ", MLE, "; hacked nonaffirms", sep ="")
p + ggtitle(string)



# QQ plot vs. fitted truncated normal (probably won't match because misspecified?)
qqtrunc( d$tstat[ d$hack == "affirm2" ],
         spec = "norm",
         a = -Inf,
         mean = .obj$sanityChecks$tstatMeanMLE,
         sd = sqrt(.obj$sanityChecks$tstatVarMLE),
         b = crit )

#bm: wanted to try this with many more hacked studies (10K)


