
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


# simulate a huge dataset, including unpublished ones
# don't need to redo this, so commented out
d = sim_meta(Nmax = 20,
             Mu = 0.1,
             T2 = 0.1,
             m = 50,
             t2w = .1,
             se = 1,
             hack = "affirm",
             return.only.published = FALSE,

             k = 500,
             k.hacked = 0 )


# dataset of only published results
dp = d %>% filter(Di == 1 )
dim(dp)

# setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
# fwrite(d, "sim_meta_all_studies.csv")
# fwrite(dp, "sim_meta_published_studies.csv")
# 
# 
# # read back in
# setwd( here("2021-4-20 simple sims with partial try-to-Nmax hacking") )
# d = fread("sim_meta_all_studies.csv")
# dp = fwrite("sim_meta_published_studies.csv")



# DATA SIMULATION SANITY CHECKS ------------------------------

table(d$hack)

nrow(d)
nrow(dp)

length(unique(d$study))


# all results (regardless of publication)
d %>%
  group_by(hack) %>%
  summarise( n(),
             kUnique = length(unique(study)),
             mean(affirm),
             mean(mui),
             var(mui),
             mean(yi) )

# look at the published results only
dp %>%
  group_by(hack) %>%
  summarise( n(),
             kUnique = length(unique(study)),
             mean(affirm),
             mean(mui),
             var(mui),
             mean(yi) )

# why is yi seemingly biased downward here???

# DEBUGGING
temp = d[ d$study ==1,]
# END DEBUGGING

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

# 1. Meta-analyze observed, unhacked studies to unbiasedly estimate mu, t2
# this one can be really imprecise because t2w is still in there
# is it because tau^2 estimate includes both T2 and t2w?
( modUH = rma( yi = yi,
               vi = vi,
               data = duh,
               method = "REML",
               knha = TRUE ) )

Mhat.UH = modUH$b; 
# *important: since t2w is a sensitivity parameter, we can just subtract it off
T2.UH = modUH$tau2 - 0.1

# # debug: try without filtering on pub status
# fake = d %>% filter(hack == "no")
# rma.mv( yi = yi,
#         V = vi,
#         data = fake,
#         method = "REML",
#         random = ~1 | study )
# 
# mean(fake$mui)


# 2. Bias-correct each hacked result using above estimates in the truncated distribution

#bm: stopped here:
# need to use Mhat.UH and T2.UH :)

# 3. Adjust all studies' variances

# 4. Meta-analyze these new things

















