
# Goal: Try to get Jeffreys prior numerically because the expectations will be hard


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


prepped.data.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Prepped Hagger data (meta and replication)"
res.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Hagger comparison results"


#code.dir = here("2021-10-7 check RTMA Jeffreys theory")

code.dir = here()
setwd(code.dir)
source("helper_SAPH.R")

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

yi = dm$yi
sei = sqrt(dm$vi)



# ~ SANITY CHECK: SUM OF FISHER VS. FISHER OF SUM ----------------------------------------------------

# PASSES :)
# If all sei are the same, then the yi's are iid, so summing expected Fishers
#  for each observation should be the same as getting the joint Fisher directly

Mu = 1
Tt = 2
k = 1000
se = 1.5

# by setting se super small, we're effectively just truncating at yi > 0 
#  since all studies with yi > 0 are significant 
dp = sim_meta(Nmax = 1,
              Mu = Mu,
              T2 = Tt^2,
              m = 100,
              t2w = 0,
              se = se,
              hack = "affirm",
              return.only.published = TRUE,
              rho = 0,
              
              k = k,
              k.hacked = k )

hist(dp$yi)
kn = nrow(dp)
tcrit = unique(dp$tcrit)

### Version 1: Joint EFish as in TNE
( EFish1 = E_fisher_TNE( .mu = Mu,
                         .sigma = sqrt(Tt^2 + se^2), 
                         .n = kn,
                         .a = -99,
                         #*note that .b here needs to be on yi scale, not Zi scale
                         .b = tcrit*se ) )

### Version 2: Sum individual Efish for each observation
( EFish2 = E_fisher_RTMA( .sei = sei.true,
                          .Mu = Mu,
                          .Tt = Tt,
                          .crit = tcrit ) )

expect_equal( Efish1, Efish2 )


# ~ SANITY CHECK: AGREEMENT WITH MLE ----------------------------------------------------

# check that not using the prior agrees with mle

# you already have the usePrior argument :)




# ~ SANITY CHECK: PLOT PRIOR ----------------------------------------------------


if ( redo.prior.plots == TRUE ) {
  
  
  dp = expand.grid( .kn = c(50), # changing this doesn't matter much
                    .Mu = seq(-4, 4, .1),
                    .Tt = seq(0.5, 4, .1),
                    .sei = 1.5 )
  
  
  # make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
  dp = dp %>%
    rowwise() %>%
    mutate( lprior = lprior( .sei = rep(.sei, .kn),
                             .Mu = .Mu,
                             .Tt = .Tt,
                             # just use z-approx for now
                             .crit = 1.96 )  )
  
  # check again for NA values occurring when EFisher is NA and remove them
  table(is.na(dp$lprior))
  dp = dp %>% filter( !is.na(lprior) )
  
  # set up colors for contours
  get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
  myColors = get_colors(n=15)  # chose number based on errors from ggplot if it was fewer
  
  title.string = paste( "Log-prior for 50 studies all with SE=", unique(dp$.sei))
  
  # make plot
  p = ggplot( data = dp,
              aes(x = .Mu,
                  y = .Tt,
                  z = lprior) ) +
    
    geom_contour_filled() +
    
    # close, but not enough colors
    scale_fill_manual(values = myColors) +
    
    geom_contour(color = "white") +
    
    xlab( bquote(mu) ) +
    ylab( bquote(tau[t]) ) +
    
    geom_vline( xintercept = 0, lty = 2 ) +
    
    ggtitle(title.string) +
    
    facet_wrap( .kn ~.,
                scales = "fixed" ) +
    
    # scale_y_continuous(breaks = seq( min(dp$.sigma), max(dp$.sigma), 0.5),
    #                    limits = c( min(dp$.sigma), max(dp$.sigma) ) ) +
    
    theme_bw(base_size = 16) +
    theme(text = element_text(face = "bold"),
          axis.title = element_text(size=20),
          legend.position = "none")
  
  
} # end "if redo.prior.plots == TRUE"


# ~ PLOT PRIOR FOR HAGGER ----------------------------------------------------

if ( redo.prior.plots == TRUE ) {
  yi = dm$yi
  sei = sqrt(dm$vi)
  # just as an experiment, what happens if all have same sei?
  sei = rep( mean(sqrt(dm$vi)), length(dm$vi) )
  
  
  dp = expand.grid( .Mu = seq(-2, 2, .1),
                    .Tt = seq(0.1, 2, .1) )
  
  
  # make plotting dataframe by calculating log-prior for a grid of values (mu, sigma)
  dp = dp %>%
    rowwise() %>%
    mutate( lprior = lprior( .sei = sei,
                             .Mu = .Mu,
                             .Tt = .Tt,
                             # just use z-approx for now
                             .crit = 1.96 )  )
  
  # check again for NA values occurring when EFisher is NA and remove them
  table(is.na(dp$lprior))
  dp = dp %>% filter( !is.na(lprior) )
  
  # set up colors for contours
  get_colors = colorRampPalette( c("lemonchiffon1", "chocolate4") )
  myColors = get_colors(n=15)  # chose number based on errors from ggplot if it was fewer
  
  title.string = paste( "Log-prior for Hagger data, setting all studies' sei to their mean")
  
  # make plot
  p = ggplot( data = dp,
              aes(x = .Mu,
                  y = .Tt,
                  z = lprior) ) +
    
    geom_contour_filled() +
    
    # close, but not enough colors
    scale_fill_manual(values = myColors) +
    
    geom_contour(color = "white") +
    
    xlab( bquote(mu) ) +
    ylab( bquote(tau[t]) ) +
    
    geom_vline( xintercept = 0, lty = 2 ) +
    
    ggtitle(title.string) +
    
    
    # scale_y_continuous(breaks = seq( min(dp$.sigma), max(dp$.sigma), 0.5),
    #                    limits = c( min(dp$.sigma), max(dp$.sigma) ) ) +
    
    theme_bw(base_size = 16) +
    theme(text = element_text(face = "bold"),
          axis.title = element_text(size=20),
          legend.position = "none")
  
} # end "if redo.prior.plots == TRUE"



# ~ GET NLPOSTERIOR FOR SIM META ----------------------------------------------------

Mu = 1
Tt = 2
k = 1000
se = 1
tcrit = unique(dn$tcrit)

# SAVE: this is how data were generated
# by setting se super small, we're effectively just truncating at yi > 0
#  since all studies with yi > 0 are significant
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

# get nonaffirms based on TRUE SEs
# temporary because otherwise they won't all have same se
dn = d %>% filter( yi < tcrit * se )


setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/2021-10-13 numerical integration for RTMA Jeffreys prior")
fwrite(dn, "sim_meta_1.csv")
dn = fread("sim_meta_1.csv")

# set vars for all methods
kn = nrow(dn)
Mu.start = 0
Tt.start = 1
sei.true = rep(se, kn)


### Version 1: MAP ###
res.MAP = estimate_jeffreys_RTMA(yi = dn$yi,
                                 
                                 sei = sei.true,  # true SEs
                                 # sei = sqrt(dn$vi),  # estimated SEs
                                 
                                 par2is = "Tt",
                                 Mu.start,
                                 Tt.start = Tt.start,
                                 crit = tcrit,
                                 
                                 usePrior = TRUE,
                                 get.CIs = FALSE,
                                 CI.method = "wald")
res.MAP$MuHat

### Version 1: MLE ###
# without using prior
# doesn't work at all ("Lapack routine dgesv: system is exactly singular: U[1,1] = 0")
res.MLE = estimate_jeffreys_RTMA(yi = dn$yi,
                                 
                                 sei = sei.true,  # true SEs
                                 # sei = sqrt(dn$vi),  # estimated SEs
                                 
                                 par2is = "Tt",
                                 Mu.start,
                                 Tt.start = Tt.start,
                                 crit = tcrit,
                                 
                                 usePrior = FALSE,  # MLE
                                 get.CIs = FALSE,
                                 CI.method = "wald")

# this does give a value
( nll.start.RTMA = nlpost_jeffreys_RTMA( .pars = c(0,1),
                                            .par2is = "Tt",  # "Tt" or "Tt2"
                                            .yi = dn$yi,
                                            .sei = sei.true,
                                            .crit = qnorm(.975),
                                            
                                            # if .usePrior = FALSE, will just be the MLE
                                            .usePrior = FALSE ) )


### Version 3: Directly use TNE ###
# because SEs are the same

# as in doParallel_TNE.R
p = data.frame(n = kn,
               a = -99,
               b = tcrit * se,
               stan.iter = 2000,  # default: 2000
               stan.maxtreedepth = 10, # default: 10
               stan.adapt_delta = 0.8,
               get.CIs = FALSE )

res.MLE.TNE = estimate_mle(x = dn$yi,
                           p = p,
                           par2is = "var",  # NOTE: it prefers var parameterization
                           mu.start = Mu.start,
                           sigma.start = Tt.start,  # resulting estimate will INCLUDE se
                           get.CIs = p$get.CIs,
                           CI.method = "wald")

res.MLE.TNE$Mhat



( nll.start.TNE = nll(.pars = c(Mu.start, Tt.start),
                    par2is = "sd",
                    .x = dn$yi,
                    .a = -99,
                    .b = tcrit*se) )


### Figure out why NLLs disagree ###

### from inside nll()
( term1 = dnorm(x = dn$yi,
                mean = Mu.start,
                sd = Tt.start,  
                log = TRUE) )

( term2 = length(dn$yi) * log( pmvnorm(lower = -99,
                                    upper = tcrit*se,
                                    mean = Mu.start,
                                    # note use of sigma^2 here because of pmvnorm's different parameterization:
                                    sigma = Tt.start^2 ) ) )


-( sum(term1) - term2 )

# sum(term1) = -2601.394
# term2 = -5.377309

### from inside nll2() (SAPH)

joint_nll_2(.yi = dn$yi,
            .sei = sei.true,
            .Mu = Mu.start,
            .Tt = Tt.start,
            .crit = tcrit)
# inside this fn:
# term1 = -1482.057
# term2 = 1.078192


#bm: MAP one is giving something huge and negative
# MLE one 
# **would be good to look at case with all SEs the same, because that should agree exactly with TNE estimate_jeffreys
# YOU GOT THIS





