
# PRELIMINARIES ----------------------------------------------------

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
library(RColorBrewer)


# ~~ Set global vars needed throughout -------------------------
# no sci notation
options(scipen=999)

# control which results should be redone and/or overwritten
redo.plots = TRUE


# variables that define the scenarios (though some aren't actually manipulated, like stan.iter)
param.vars = c("unique.scen",  
               "method",
               "boot.reps",
               "stan.iter",
               "stan.adapt_delta",
               "stan.maxtreedepth",
               "trunc.type",
               "prop.retained",
               "mu",
               "V",
               "n")

# used later to create plots and tables, but needed to check var types 
#  upon reading in data
estNames = c("Mhat", "Shat")

outcomeNames = c("Bias", "RMSE", "EmpSE", "EstFail",
                 "Cover", "Width", "CIFail", "Rhat")


# ~~ Set directories -------------------------
code.dir = here("Sherlock code")


data.dir = str_replace( string = here(),
                        pattern = "Code \\(git\\)",
                        replacement = "Sherlock simulation results/Pilot simulations/2022-2-27 one scen" )

results.dir = str_replace( string = here(),
                           pattern = "Code \\(git\\)",
                           replacement = "Sherlock simulation results/Pilot simulations/2022-2-27 one scen" )

# overleaf.dir.general = "/Users/mmathur/Dropbox/Apps/Overleaf/TNE (truncated normal estimation)/R_objects"
# overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/TNE (truncated normal estimation)/R_objects/figures"

setwd(code.dir)
source("helper_SAPH.R")


# ~~ Get aggregated data -------------------------

# setwd(data.dir)
# agg = fread("agg_dataset_clean.csv")
# # check when the dataset was last modified to make sure we're working with correct version
# file.info("agg_dataset_clean.csv")$mtime
# 
# # sanity checks
# expect_equal( nuni(agg$scen.name), 437 )
# 
# 
# ##@check for var type problems for all variables in estNames, outcomeNames
# 
# # temporary fixes for columns with values that are outside machine precision:
# #  https://stackoverflow.com/questions/21752121/fread-from-data-table-package-cant-read-small-numbers
# agg$ShatRMSE = as.numeric(agg$ShatRMSE)
# agg$MhatRMSE = as.numeric(agg$MhatRMSE)
# # if combining separate Sherlock runs
# agg = bind_rows( fread("*2021-7-28 agg_dataset.csv"),  # mle and boot-mle
#                  fread("*2021-8-1 agg_dataset.csv"), # trunc-SP
#                  fread("*2021-8-8 agg_dataset.csv"), )  # jeffreys, boot-jeffreys, mle again



# ~~ Get iterate-level data -------------------------

setwd(data.dir)
# check when the dataset was last modified to make sure we're working with correct version
s = fread( "stitched.csv")

#@TEMP
#s = fread("long_results_job_1_.csv")

file.info("stitched.csv")$mtime

dim(s)


s %>% filter(method == "jeffreys-mcmc-max-lp-iterate") %>% 
  select(Mhat, optimx.Mhat.winner, overall.error, MhatRhat)


#**Important: getting a lot of weird errors for MCMC
s %>% group_by(method) %>%
  summarise( meanNA(is.na(Mhat)))

#@@a lot of those cryptic cluster errors
s$overall.error[ s$method == "jeffreys-mcmc-pmean"]

# 2022-2-28: many of the MCMC had this error:
# ReadItem: unknown type 0, perhaps written by later version of R


# ~~ Aggregate locally -------------------------

agg = make_agg_data(s)


# ~~ Main results -------------------------

# look at just certain cols
t = agg %>% select(method, 
                   sim.reps.actual,
                   all_of(names_with(agg, "Mhat")) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))

data.frame(t)
View(t)


# ~~ Look at sanity checks -------------------------

# scenario diagnostics for scenario
agg.checks = agg %>% select( 
  all_of(names_with(agg, "sancheck.")) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))


t(agg.checks)


# EXPLORE BAD ITERATES -------------------------

# add info that will later be collected by helper
#@but note this DOESN'T include kkt criteria
keepers = names_with(s, "convcode")

temp = s %>% select(keepers)
s$optimx.Nconvcode0 = rowSums(temp == 0)


# sanity check
keepers = names_with(s, "convcode")

temp = s %>% filter(Mhat > 9 & method == "jeffreys-sd") %>%
  select(keepers)
data.frame(temp)




keepers = c("Mhat", "Shat", names_with(s, "optimx."))

temp = s %>% filter(Mhat > 9 & method == "jeffreys-sd") %>%
  select(keepers)

data.frame(temp)



# look at results among only iterates that 
N.optimizers = length( names_with(s, "nll") )
summary(s$optimx.Nconvcode0)
summary(s$optimx.Pagree.of.convergers.Mhat.winner)

meanNA( s$Mhat[ s$method == "jeffreys-sd"] > 4 )




# *very informative plots:
ggplot( s %>% filter( method == "jeffreys-sd" &
                       !is.na(Mhat) ),
        aes(x = as.factor(optimx.Nconvcode0),
            y = Mhat) ) +
  geom_violin() +
  theme_classic()


ggplot( s %>% filter( method == "jeffreys-sd" &
                        !is.na(Mhat) ),
        aes(x = as.factor(optimx.Pagree.of.convergers.Mhat.winner == 1 & optimx.Nconvcode0 > 3),
            y = Mhat) ) +
  geom_violin() +
  theme_classic()


# ~~ Restrict to iterates with better convergence properties ---------------------

# with this threshold, we retain 92% of iterates; MLE-sd is unbiased with coverage 92%
#  both jeffreys have Mhat = 0.22 and coverage 96%
thresh = 6
meanNA(s$optimx.Nconvcode0 > thresh & s$optimx.Pagree.of.convergers.Mhat.winner == 1)
agg2 = make_agg_data(s %>%
                       filter( optimx.Nconvcode0 > thresh &
                                 optimx.Pagree.of.convergers.Mhat.winner == 1) )

# look at just certain cols
t2 = agg2 %>% select(method, 
                   sim.reps.actual,
                   all_of(names_with(agg2, "Mhat")) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))

data.frame(t2)
View(t2)


