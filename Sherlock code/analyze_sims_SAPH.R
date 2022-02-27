
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
redo.contour.plot = FALSE

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

# choose methods to show in main text (all others will only be in Supplement) 
mainTextEstMethods = c("Jeffreys mode", "MLE")  
mainTextInfMethods = c("Jeffreys posterior quantiles", "MLE Wald", "MLE profile")  


# ~~ Set directories -------------------------
code.dir = here("Sherlock code")

# SAVE: general data and results directories
# data.dir and results.dir need to be saved outside what here() considers the top-level dir
#  because too big for git
# data.dir = paste( str_remove( string = here(),
#                               pattern = "Git/Code" ),
#                   "Sherlock simulation results/Stitched data",
#                   sep = "" )
# results.dir = paste( str_remove( string = here(),
#                                  pattern = "Git/Code" ),
#                      "Sherlock simulation results/Analysis output",
#                      sep = "" )

# directories specific to this simulation run
#@ PREVIOUS FULL SIMS:
# data.dir = paste( str_remove( string = here(),
#                               pattern = "Git/Code" ),
#                   "Sherlock simulation results/All simulations in manuscript/*2021-10-10 full sims (5K) with all scens/Data/Prepped data",
#                   sep = "" )
# results.dir = paste( str_remove( string = here(),
#                                  pattern = "Git/Code" ),
#                      "Sherlock simulation results/All simulations in manuscript/*2021-10-10 full sims (5K) with all scens/Results",
#                      sep = "" )

# NEW FULL SIMS:
data.dir = paste( str_remove( string = here(),
                              pattern = "Git/Code" ),
                  "Sherlock simulation results/All simulations in manuscript/2021-11-11 full sims with all scens/Data/Prepped data",
                  sep = "" )
results.dir = paste( str_remove( string = here(),
                                 pattern = "Git/Code" ),
                     "Sherlock simulation results/All simulations in manuscript/2021-11-11 full sims with all scens/Results",
                     sep = "" )

#@CHANGED TO TEMP DIR TO PRESERVE THE EXISTING STATS_FOR_PAPER.CSV
overleaf.dir.general = "/Users/mmathur/Dropbox/Apps/Overleaf/TNE (truncated normal estimation)/R_objects"
overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/TNE (truncated normal estimation)/R_objects/figures"

setwd(code.dir)
source("helper_TNE.R")


# ~~ Get aggregated data -------------------------
setwd(data.dir)
agg = fread("agg_dataset_clean.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg_dataset_clean.csv")$mtime

# sanity checks
expect_equal( nuni(agg$scen.name), 437 )


##@check for var type problems for all variables in estNames, outcomeNames

# temporary fixes for columns with values that are outside machine precision:
#  https://stackoverflow.com/questions/21752121/fread-from-data-table-package-cant-read-small-numbers
agg$ShatRMSE = as.numeric(agg$ShatRMSE)
agg$MhatRMSE = as.numeric(agg$MhatRMSE)

# # if combining separate Sherlock runs
# agg = bind_rows( fread("*2021-7-28 agg_dataset.csv"),  # mle and boot-mle
#                  fread("*2021-8-1 agg_dataset.csv"), # trunc-SP
#                  fread("*2021-8-8 agg_dataset.csv"), )  # jeffreys, boot-jeffreys, mle again



# ~~ Get iterate-level data -------------------------

# use the small iterate-level data (fewer unique n and prop.retained) to 
#  make local manipulations computationally feasible

# this retains 14% of all iterates as follows (from stitch_on_sherlock.R):
# s.small = s %>% filter( prop.retained %in% c(0.1, 0.3, 0.5, 0.9) &
#                             n %in% c(10, 50, 80, 500, 1000) )


setwd(data.dir)
# check when the dataset was last modified to make sure we're working with correct version
s.small = fread( "stitched_small_subset.csv")
file.info("stitched_small_subset.csv")$mtime

# based on its dimensions on 2021-11-11
expect_equal(nrow(s.small), 2194986)