
# Goal: Do some simple sanity checks with the weighting-based idea on Overleaf.

# PRELIMINARIES ------------------------------

library(here)
setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
source("helper.R")

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

# SIMULATE DATA ------------------------------

# bm: for now maybe don't try to force unhackeds to be nonaffirmative; just assume we know which studies are unhacked
#  then try applying my weighting estimators and see how they compare to estimate from meta-analyzing all studies, incl
#  underlying ones? :)

# parameters also needed later
Mu = 1
T2 = 0.1
m = 50
t2w = .1
se = 1


# # simulate a huge dataset, including unpublished ones
# # don't need to redo this, so commented out
# d = sim_meta(Nmax = 20,
#              Mu = Mu,
#              T2 = T2,
#              m = m,
#              t2w = t2w,
#              se = se,
#              hack = "affirm",
#              return.only.published = FALSE,
#              
#              k = 200,
#              k.hacked = 100 )


# dataset of only published results
dp = d %>% filter(Di == 1 )
dim(dp)

# save for later
# setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
# fwrite(d, "sim_meta_all_studies.csv")
# fwrite(dp, "sim_meta_published_studies.csv")


# read back in
setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
d = fread("sim_meta_all_studies.csv")
dp = fwrite("sim_meta_published_studies.csv")



# DATA SIMULATION SANITY CHECKS ------------------------------

table(d$hack)

nrow(d)
nrow(dp)

length(unique(d$study))


# all results, sorted by hacking status and publication status
t = d %>%
  group_by(hack, Di) %>%
  summarise( n(),
             k = length(unique(study)),
             mean(affirm),
             mean(mui),
             var(mui),
             mean(yi) )
t

# unhacked, published results only
# so only one per study set
# same as second row of above table
duh = d %>% filter(hack == "no" & Di == 1)


# ANALYZE ------------------------------

# ~ Simple meta-analyses ------------------------------

# unbiased meta-analysis of all studies, even unpublished ones
# account for clustering of draws within studies
# *the tau^2 estimate will be close to T2
( modAll = rma.mv( yi = yi,
                   V = vi,
                   data = d,
                   method = "REML",
                   random = ~1 | study ) )
# looks good :)

# biased meta-analysis of only published studies
( modPub = rma( yi = dp$yi,
                vi = dp$vi,
                method = "REML",
                knha = TRUE ) )

# ~ Weighted estimator ------------------------------

# ~~ 1. Meta-analyze observed, unhacked studies to unbiasedly estimate mu, t2 ------------------------------
# this one can be really imprecise because t2w is still in there
# is it because tau^2 estimate includes both T2 and t2w?
( modUH = rma( yi = yi,
               vi = vi,
               data = duh,
               method = "REML",
               knha = TRUE ) )

Mhat.UH = modUH$b
# *important: since t2w is a sensitivity parameter, we can just subtract it off
T2.UH = modUH$tau2 - t2w

# # debug: try without filtering on pub status
# fake = d %>% filter(hack == "no")
# rma.mv( yi = yi,
#         V = vi,
#         data = fake,
#         method = "REML",
#         random = ~1 | study )
# 
# mean(fake$mui)


# ~~ 2. Bias-correct each hacked result using above estimates in the truncated distribution ------------------------------

# *estimate* the bias of the hacked affirmatives using the expectation of truncated normal

# estimate the noncentrality parameter
# uses estimated mean, tau^2, and t2w treated as known
ncp = Mhat.UH / sqrt( T2.UH + t2w + se^2 )

affirmExp = extrunc( spec = "t",
                     ncp = ncp,
                     df = m-1,
                     # since I simulated iid studies here, the truncation cutoff is always the same
                     a = d$tcrit[1] )

estBias = affirmExp - Mhat.UH

# # sanity check: compare this estimated expectation to real one:
# extrunc( spec = "t",
#          ncp = Mu / sqrt( T2 + t2w + se^2 ),
#          df = m-1,
#          # since I simulated iid studies here, the truncation cutoff is always the same
#          a = d$tcrit[1] )
# #...and to empirical one
# t$`mean(yi)`[ t$hack == "affirm" & t$Di == 1 ]
# # all quite close, even with T2 estimate pretty off in this sample! 

# bias-correct the affirmatives in published data
dp$yiCorr = dp$yi
dp$yiCorr[ dp$hack == "affirm" ] = dp$yi[ dp$hack == "affirm" ] - c(estBias)

dp %>% group_by(hack) %>%
  summarise( n(),
             k = length(unique(study)),
             mean(affirm),
             mean(mui),
             var(mui),
             mean(yi),
             mean(yiCorr)
             )
# seems to have worked pretty well! 


# ~~ 3. Adjust all studies' variances ------------------------------

# **not actually adjusting because issue doesn't seem that bad
# for now, don't adjust nonaffirmatives' variances because uncorrelated p-hacking

# adjust affirmatives' variances

# how bad is the variance estimation?
d %>% group_by(hack, Di) %>%
  summarise( mean(sqrt(vi)) )

# doesn't seem very severe

# # not quite sure how to do the bias-correction, actually
# ncp2 = Mhat.UH / se
# affirmVarExp = vartrunc( spec = "t",
#                      ncp = ncp2,
#                      df = m-1,
#                      # since I simulated iid studies here, the truncation cutoff is always the same
#                      a = d$tcrit[1] )
# 
# varBias = affirmVarExp - Mhat.UH


# ~~ 4. Meta-analyze these new things ------------------------------

( modCorr = rma( yi = dp$yiCorr,
                vi = dp$vi,
                method = "REML",
                knha = TRUE ) )


# hmm! pretty good!
# c.f. uncorrected model:
modPub







correct_studies = function() {
  
  # this one can be really imprecise because t2w is still in there
  # is it because tau^2 estimate includes both T2 and t2w?
  ( modUH = rma( yi = yi,
                 vi = vi,
                 data = duh,
                 method = "REML",
                 knha = TRUE ) )
  
  Mhat.UH = modUH$b
  # *important: since t2w is a sensitivity parameter, we can just subtract it off
  T2.UH = modUH$tau2 - t2w
  
  # # debug: try without filtering on pub status
  # fake = d %>% filter(hack == "no")
  # rma.mv( yi = yi,
  #         V = vi,
  #         data = fake,
  #         method = "REML",
  #         random = ~1 | study )
  # 
  # mean(fake$mui)
  
  
  # ~~ 2. Bias-correct each hacked result using above estimates in the truncated distribution ------------------------------
  
  # *estimate* the bias of the hacked affirmatives using the expectation of truncated normal
  
  # estimate the noncentrality parameter
  # uses estimated mean, tau^2, and t2w treated as known
  ncp = Mhat.UH / sqrt( T2.UH + t2w + se^2 )
  
  affirmExp = extrunc( spec = "t",
                       ncp = ncp,
                       df = m-1,
                       # since I simulated iid studies here, the truncation cutoff is always the same
                       a = d$tcrit[1] )
  
  estBias = affirmExp - Mhat.UH
  
  # # sanity check: compare this estimated expectation to real one:
  # extrunc( spec = "t",
  #          ncp = Mu / sqrt( T2 + t2w + se^2 ),
  #          df = m-1,
  #          # since I simulated iid studies here, the truncation cutoff is always the same
  #          a = d$tcrit[1] )
  # #...and to empirical one
  # t$`mean(yi)`[ t$hack == "affirm" & t$Di == 1 ]
  # # all quite close, even with T2 estimate pretty off in this sample! 
  
  # bias-correct the affirmatives in published data
  dp$yiCorr = dp$yi
  dp$yiCorr[ dp$hack == "affirm" ] = dp$yi[ dp$hack == "affirm" ] - c(estBias)
  
  dp %>% group_by(hack) %>%
    summarise( n(),
               k = length(unique(study)),
               mean(affirm),
               mean(mui),
               var(mui),
               mean(yi),
               mean(yiCorr)
    )
  # seems to have worked pretty well! 
  
}








