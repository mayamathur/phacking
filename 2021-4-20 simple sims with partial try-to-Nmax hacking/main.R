
# Goal: Do some simple sanity checks with the weighting-based idea on Overleaf.

# PRELIMINARIES ---------------------------------------------------------------

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

# SIMULATE DATA ---------------------------------------------------------------

# data simulation:
# - for each affirmative, simulate infinite draws for affirmatives (have option of correlated vs. uncorrelated to check if still trunc normal under correlated draws)

# - for each nonaffirmative, simulate Nmax total underlying nonaffirmatives and take the last one to be the observed one
#  - I think we need to retain the whole study set, though

