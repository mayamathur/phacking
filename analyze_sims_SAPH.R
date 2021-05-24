

# PRELIMINARIES --------------------------------------------

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

setwd(here())
source("analyze_sims_helper_SAPH.R")

setwd(here("*2021-5-17 nonaffirm MLEs with large k and affirm2 hacking"))
s = fread("stitched.csv")

length(unique(s$scenName))

# check reps run vs. expected
expect_equal( unique(table(s$scenName)), 500 )

# TEMP ONLY
# remove old files
s = s[ s$hack == "affirm2", ]


paramVars = names(s)[ 1 : ( which( names(s) == "MhatAll" ) - 1 ) ]
( paramVars = paramVars[ !paramVars %in% c("V1", "scenName", "repName") ] )

outcomeVars = names(s)[ which( names(s) == "MhatAll") : length(names(s)) ]
( outcomeVars = outcomeVars[ !outcomeVars %in% c("repSeconds") ] )
  
# non-parameter vars to keep only first entry per scenario because static
firstOnly = "scenName"

# sanity check:
# make sure we listed all the param vars
t = s %>% group_by_at(paramVars) %>%
  summarise( scenName = scenName[1],
             reps = n(),
             .groups = "keep" )
# 500 such that each scenario is uniquely defined by the param vars
expect_equal( unique(t$reps), 500 )

# if some scenarios timed out, might have fewer than 500 reps
table(t$reps)




# AGGREGATE  --------------------------------------------


# aggregate by scenario
agg = s %>% 
  
  # take just first entry of non-parameter variables that are static within scenarios
  group_by_at(paramVars) %>%
  mutate_at( firstOnly, 
             function(x) x[1] ) %>%
  
  #group_by(calib.method.pretty) %>%
  group_by_at(paramVars) %>%
  summarise_at( .vars = outcomeVars,
                .funs = meanNA )


# I expect k.hacked scenarios to be unbiased
# order from smallest to largest corrected estimate
View( agg[ order(agg$MhatCorr), ] )

agg %>% group_by( k.hacked, T2, t2w ) %>%
  summarise( MhatAll = meanNA(MhatAll),
             MhatCorr = meanNA(MhatCorr),
             MhatCoverCorr = meanNA(MhatCoverCorr) )

fwrite(agg, "agg.csv")

# rounded version
aggRounded = agg %>% mutate_at( .vars = outcomeVars,
                                .fun = function(x) round(x,2) )
fwrite(aggRounded, "agg_rounded.csv")


# reminder from genSbatch:
#
# main scenarios of interest:
# 1. Nmax = 1, k.hacked = 0, rho = 0 (basically a sanity check)
# 2. Nmax > 1, k.hacked = 0, rho = 0.9 
# 3. Nmax > 1, k.hacked = 50, rho = 0 or 0.9 (conservative?)
# 
# scen.params = expand_grid( Mu = 0.1,
#                            T2 = c(0, 0.25),
#                            m = 500,
#                            t2w = c(0, 0.25),
#                            se = 0.5,
#                            
#                            Nmax = c(1, 10),
#                            hack = "affirm",
#                            rho = c(0, 0.9),
#                            
#                            k = 100,
#                            k.hacked = c(0, 50) )
# 
# # remove nonsense combinations
# # rho > 0 is pointless if there's only 1 draw
# scen.params = scen.params %>% dplyr::filter( !(rho > 0 & Nmax == 1) )
# 
# scen.params = scen.params %>% add_column( scen = 1:nrow(scen.params),
#                                           .before = 1 )

