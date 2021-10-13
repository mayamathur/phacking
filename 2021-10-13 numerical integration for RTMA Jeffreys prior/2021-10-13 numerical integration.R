
# Goal: Try to get Jeffreys prior numerically because the expectations will be hard


# ~ PRELIMINARIES ----------------------------------------------------

rm(list=ls())


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


prepped.data.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Prepped Hagger data (meta and replication)"
res.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Hagger comparison results"


#code.dir = here("2021-10-7 check RTMA Jeffreys theory")

code.dir = here()
setwd(code.dir)
source("helper_SAPH.R")



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

yi = dm$yi
sei = sqrt(dm$vi)


# ~ GET EXPECTED FISHER NUMERICALLY ----------------------------------------------------


.Mu = 0.5
.Tt = 0.2
.yi = yi
.sei = sei

# expectation of -Dij wrt yi
integrand_Dij = function( i, j, .yi, .sei, .Mu, .Tt, .crit = qnorm(0.975) ) {
  
  # function of yi whose expectation we want
  if ( i == 1 & j == 1 ) {
    term1 = -get_D11(.yi = .yi,
                     .sei = .sei,
                     .Mu = .Mu,
                     .Tt = .Tt)
  } else if ( ( i == 2 & j == 1 ) | ( i == 1 & j == 2 ) ) {
    term1 = -get_D21(.yi = .yi,
                     .sei = .sei,
                     .Mu = .Mu,
                     .Tt = .Tt)
  } else if ( i == 2 & j == 2 ) {
    term1 = -get_D22(.yi = .yi,
                     .sei = .sei,
                     .Mu = .Mu,
                     .Tt = .Tt)
  }
  
  # density of yi
  term2 = dtruncnorm( x = .yi,
                      b = .crit * .sei )
  
  term1 * term2
}


# for a single choice of .Mu and .Tt and a vector of .yi, .sei
E_fisher_RTMA = function( .sei, .Mu, .Tt, .crit = qnorm(0.975) ) {
  
  # dataframe to store the 4 integrals for each observation
  .d = data.frame( sei = .sei )
  
  .d = .d %>% rowwise() %>%
    mutate( F11 = integrate( function(x) integrand_Dij(i = 1, 
                                                       j = 1, 
                                                       .yi = x,
                                                       .sei = sei,
                                                       .Mu = .Mu,
                                                       .Tt = .Tt),
                             lower = 0,
                             upper = Inf )$value,
            
            F21 = integrate( function(x) integrand_Dij(i = 2, 
                                                       j = 1,
                                                       .yi = x,
                                                       .sei = sei,
                                                       .Mu = .Mu,
                                                       .Tt = .Tt),
                             lower = 0,
                             upper = Inf )$value,

            F22 = integrate( function(x) integrand_Dij(i = 2, 
                                                       j = 2,
                                                       .yi = x,
                                                       .sei = sei,
                                                       .Mu = .Mu,
                                                       .Tt = .Tt),
                             lower = 0,
                             upper = Inf )$value )
  
  return( matrix( c( sum(.d$F11), sum(.d$F21), sum(.d$F21), sum(.d$F22) ),
                  byrow = TRUE,
                  nrow = 2 ) )
  
  
}

# for ONE observation
integrate( function(x) integrand_D11(.yi = x, .sei = .sei, .Mu = .Mu, .Tt = .Tt),
           lower = 0,
           upper = Inf )











