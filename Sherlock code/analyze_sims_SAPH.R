
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

# # save
# data.dir = str_replace( string = here(),
#                         pattern = "Code \\(git\\)",
#                         replacement = "Sherlock simulation results/Pilot simulations/2022-2-27 one scen" )

data.dir = str_replace( string = here(),
                        pattern = "Code \\(git\\)",
                        replacement = "Sherlock simulation results/Pilot simulations/2022-3-2 same scen; rm prior from mcmc" )


results.dir = data.dir

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


#**Check for MCMC errors and similar
s %>% group_by(method) %>%
  summarise( meanNA(is.na(Mhat)))

#@@a lot of those cryptic cluster errors
s$overall.error[ s$method == "jeffreys-mcmc-pmean"]

table(s$method, s$overall.error)

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

setwd(results.dir)
fwrite(t, "results_all_iterates.csv")


# ~~ Look at sanity checks -------------------------

# scenario diagnostics for scenario
agg.checks = agg %>% select( 
  all_of(names_with(agg, "sancheck.")) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))

# transpose it
t = tibble( t(agg.checks[1,]) )
t = t %>% add_column( variable = names(agg.checks), .before = 1 )

setwd(results.dir)
fwrite(t, "sanity_check_stats.csv")


# EXPLORE BAD ITERATES -------------------------


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
        aes(x = as.factor(optimx.Nconvergers),
            y = Mhat) ) +
  geom_hline(yintercept = scen.params$Mu) +
  geom_violin(draw_quantiles = TRUE) +
  theme_classic() +
  scale_y_continuous( limits = c(-2, 5))


ggplot( s %>% filter( method == "jeffreys-sd" &
                        !is.na(Mhat) ),
        aes(x = as.factor(optimx.Nagree.of.convergers.Mhat.winner),
            y = Mhat) ) +
  geom_hline(yintercept = scen.params$Mu) +
  geom_violin(draw_quantiles = TRUE) +
  theme_classic() +
  scale_y_continuous( limits = c(-2, 5))


ggplot( s %>% filter( method == "jeffreys-sd" &
                        !is.na(Mhat) ),
        aes(x = as.factor(optimx.Pagree.of.convergers.Mhat.winner == 1 & optimx.Nconvergers > 4),
            y = Mhat) ) +
  geom_violin(draw_quantiles = TRUE) +
  theme_classic() +
  scale_y_continuous( limits = c(-2, 5))


# ~~ **Restrict to iterates with better convergence properties ---------------------

### Explore possible optimx agreement thresholds ###

# note that these look at each method separately, rather than keeping entire sim iterates
# prop.table( table(s$optimx.Nagree.of.convergers.Mhat.winner) )
# meanNA(s$optimx.Nagree.of.convergers.Mhat.winner > 5)
# 
# thresh = 6
# meanNA(s$optimx.Nagree.of.convergers.Mhat.winner > thresh &
#          s$optimx.Pagree.of.convergers.Mhat.winner == 1)



### Mark iterates with good optimx properties for ALL methods ###
non.optimx.methods = c("naive", "gold-std", "maon", "2psm", "jeffreys-mcmc-pmean", 
            "jeffreys-mcmc-pmed", "jeffreys-mcmc-max-lp-iterate")


s$job.rep.name = paste( s$job.name, s$rep.name, sep = "_" )

# for the optimx methods, mark whether the iterate should be retained for ALL methods
#  "min" below is to consider performance across all methods that used optimx
s = s %>% group_by(job.rep.name) %>%
  mutate( keep.iterate = min(optimx.Nagree.of.convergers.Mhat.winner, na.rm = TRUE) > thresh &
                            min(optimx.Pagree.of.convergers.Mhat.winner, na.rm = TRUE) == 1 )

table(s$keep.iterate)

good.iterates = unique( s$job.rep.name[ s$keep.iterate == TRUE ] )
length(good.iterates) / nuni(s$job.rep.name)  # percent of iterates that were "good"

# sanity check
View( s %>% select(job.rep.name, method, keep.iterate, optimx.Nagree.of.convergers.Mhat.winner,
             optimx.Pagree.of.convergers.Mhat.winner) )


# should be the same for all methods
table(s$keep.iterate, s$method)

# sanity check
# within each job.rep.name, the decision to keep the iterate or not should be the same
#  for all methods and should never be NA
t = s %>% group_by(job.rep.name) %>%
  summarise( SD = sd(keep.iterate) )
expect_equal( TRUE, all(t$SD == 0) )



### Performance among iterates with good optimx properties for ALL methods ###

s2 = s %>% filter(keep.iterate == TRUE)
dim(s2)

# look at individual iterates
View( s2 %>% select(job.rep.name, method, Mhat, keep.iterate, optimx.Nagree.of.convergers.Mhat.winner,
                   optimx.Pagree.of.convergers.Mhat.winner) )

agg2 = make_agg_data(s2)

# look at just certain cols
t2 = agg2 %>% select(method, 
                   sim.reps.actual,
                   all_of(names_with(agg2, "Mhat")) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))

data.frame(t2)
View(t2)

setwd(results.dir)
fwrite(t2, "results_optimx_thresh_iterates.csv")


# ~~ Sanity checks ---------------------

# maxLPiterate should be close to MAP
# it's not...hmmm
# the discrepancy is NOT due to just a few crazy iterates because this is considering
#  only the thresholded dataset AND we can see that even the median discrepancy is 0.14
summary( abs( s2$Mhat[s2$method == "jeffreys-mcmc-max-lp-iterate"] - 
                s2$Mhat[s2$method == "jeffreys-sd"] ) )

# for the case where prior is commented out, instead it should agree with MLE
summary( abs( s2$Mhat[s2$method == "jeffreys-mcmc-max-lp-iterate"] - 
                s2$Mhat[s2$method == "mle-sd"] ) )
# it does NOT agree very closely

