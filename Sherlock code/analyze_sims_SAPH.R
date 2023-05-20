
# PRELIMINARIES ----------------------------------------------------

#rm(list=ls())

# data-wrangling packages
library(here)
library(plotly)  # must be BEFORE dplyr or else plotly::select will take over
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(tidyverse)
library(fastDummies)
library(xlsx)
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
library(sjmisc)

# prevent masking
select = dplyr::select


# ~~ User-specified global vars -------------------------
# no sci notation
options(scipen=999)

# control which results should be redone and/or overwritten
#@ not all fns respect this setting
overwrite.res = FALSE


# ~~ Set directories -------------------------
code.dir = here("Sherlock code")

data.dir = str_replace( string = here(),
                        pattern = "Code \\(git\\)",
                        replacement = "Simulation results/2023-05-08 full Stefan sims (as in RSM_1)/Data" )



results.dir = str_replace( string = here(),
                           pattern = "Code \\(git\\)",
                           replacement = "Simulation results/2023-05-08 full Stefan sims (as in RSM_1)/Results" )


overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH)/figures_SAPH/sims"

# # alternative for running new simulations
# data.dir = str_replace( string = here(),
#                         pattern = "Code \\(git\\)",
#                         replacement = "Simulation results" )
# 
# results.dir = data.dir


setwd(code.dir)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")




# TEMP ----------------

setwd(data.dir)
s = fread("stitched.csv")

aggo = make_agg_data(s)
setwd(data.dir)
fwrite(aggo, "agg.csv")


# ~~ Get agg data -------------------------

# if only analyzing a single set of sims (no merging):
setwd(data.dir)
aggo = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

sim.env = unique(aggo$sim.env)

dim(aggo)  # will exceed number of scens because of multiple methods
s

if (sim.env == "mathur") expect_equal( 84, nuni(aggo$scen.name) ) 
if (sim.env == "stefan") expect_equal( 70, nuni(aggo$scen.name) ) 


# prettify variable names
agg = wrangle_agg_local(aggo)

# drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
agg = agg %>% filter( !is.na(method.pretty) )

# look at number of actual sim reps
table(agg$sim.reps.actual)



# ~~ List variable names -------------------------

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()







# 2022-4-4: EFFECT OF SCEN PARAMS ON DATASETS -------------------------

#@not updated

param.vars.manip2 = drop_vec_elements(param.vars.manip, "method")

t = agg %>% group_by_at( param.vars.manip2 ) %>%
  # keep only scens in main text
  #filter(Mu == 0.5 & true.sei.expr == "0.02 + rexp(n = 1, rate = 3)") %>%
  select( all_of( contains("sancheck") ) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

# setwd(results.dir)
# setwd("Sanity checks")
# write.xlsx( as.data.frame(t), "table_sanchecks.xlsx")


# ONE-OFF STATS FOR PAPER  -------------------------

#@pull in from MBMA


# ******** WINNER TABLES -------------------------

# ~ Winner tables shown in paper ------------------------------

# can toggle output of fn below by changing the default arg of 
#  make_winner_table between display = "dataframe" (easy viewing)
#  and display = "xtable" (Overleaf)

# 1: all scenarios
make_both_winner_tables(.agg = agg)
# here, reason SM-step appears unbiased is that it's positively biased under evil.selection=0
#   but negatively biased under evil.selection=1


# 2: by k
make_both_winner_tables(.agg = agg %>% filter(k.pub.nonaffirm == 10) )
make_both_winner_tables(.agg = agg %>% filter(k.pub.nonaffirm == 100) )


# OLD

# # ******** PLOTS (BIG AND NOT PRETTIFIED) -------------------------
# 
# Ynames = rev(MhatYNames)
# 
# # alternatively, run just a subset:
# # Ynames = c("MhatWidth", "MhatCover", "MhatBias",
# #            "MhatEstFail",
# #            # last 2 are useful for looking at MAN
# #            "Mhat", "MhatTestReject")
# 
# # to help decide which vars to include in plot:
# param.vars.manip2
# 
# 
# # in case you want to filter scens:
# # full set for reference:
# # c("naive", "gold-std", "maon", "2psm", "pcurve", "jeffreys-mcmc-pmean", 
# #   "jeffreys-mcmc-pmed", "jeffreys-mcmc-max-lp-iterate", "jeffreys-sd", 
# #   "jeffreys-var", "mle-sd", "csm-mle-sd", "mle-var", "2psm-csm-dataset", 
# #   "prereg-naive", "ltn-mle-sd")
# ( all.methods = unique(agg$method) )
# #toDrop = c("jeffreys-mcmc-pmean", "jeffreys-mcmc-max-lp-iterate")
# toDrop = NULL
# method.keepers = all.methods[ !all.methods %in% toDrop ]
# 
# 
# # for each hacking method and Mu, make facetted plotly
# 
# for ( .hack in unique(agg$hack) ) {
#   
#   
#   for ( .Mu in unique(agg$Mu) ) {
#     
#     cat( paste("\n\n -------- STARTING Mu=", .Mu, ", hack=", .hack, sep = "") )
#     
#     aggp = agg %>% filter(method %in% method.keepers &
#                             Mu == .Mu &
#                             hack == .hack)
#     # to label the plots
#     prefix = paste( "2022-5-4 sims; ",
#                     "Mu=", .Mu,
#                     "; hack=", .hack, 
#                     sep = "")
#     
#     
#     # temporarily set wd
#     # results.dir.temp = paste(results.dir,
#     #                          "/Big unprettified plots/",
#     #                          .Mu,
#     #                          "/hack=",
#     #                          .hack,
#     #                          sep = "")
#     
#     results.dir.temp = paste(results.dir,
#                              "/Big unprettified plots",
#                              sep = "")
#     
#     
#     # set one of the two facetting variables for plots
#     aggp$tempFacetVar2 = paste( "t2a=", aggp$t2a, "; t2w=", aggp$t2w, sep = "")
#     table(aggp$tempFacetVar2)
#     
#     
#     for ( Yname in Ynames) {
#       
#       # to run "manually"
#       #Yname = "MhatBias"
#       #Yname = "MhatCover"
#       
#       y.breaks = NULL
#       if ( Yname == "MhatBias") y.breaks = seq(-0.5, 0.5, 0.1)
#       if ( Yname == "MhatWidth") y.breaks = seq(0, 10, 0.5)
#       
#       p  = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
#                                .Yname = Yname,
#                                .colorVarName = "method",
#                                .facetVar1Name = "true.sei.expr.pretty",
#                                .facetVar2Name = "tempFacetVar2",
#                                .dat = aggp,
#                                .ggtitle = prefix,
#                                .y.breaks = y.breaks,
#                                .writePlot = FALSE,
#                                .results.dir = results.dir.temp)
#       
#       
#       pl = ggplotly(p)
#       
#       # in filename, mark the most important plots with asterisk
#       if ( Yname %in% c("MhatBias", "MhatCover", "MhatWidth") ){
#         new.prefix = paste("*", prefix, sep = "")
#       } else {
#         new.prefix = prefix
#       }
#       
#       # how to save a plotly as html
#       # https://www.biostars.org/p/458325/
#       setwd(results.dir.temp)
#       string = paste(new.prefix, Yname, "plotly.html", sep="_")
#       htmlwidgets::saveWidget(pl, string)
#       
#     }
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
# # ******** PLOTS (SIMPLE AND PRETTY FOR MAIN TEXT) -------------------------
# 
# 
# 
# # for each hack type, arrange plots so each facet row is an outcome
# ( all.methods = unique(agg$method.pretty) )
# ( method.keepers = all.methods[ !is.na(all.methods) &
#                                   all.methods != "Gold standard"] )
# 
# 
# # outcomes to show in main text figures
# YnamesMain = c("MhatBias", "MhatCover", "MhatWidth")
# 
# # outcomes to show in supplement figures
# YnamesSupp = c("MhatBias", "MhatCover", "MhatWidth",
#                "MhatTestReject")
# 
# # this dataset will be one full-page figure in main text or Supp depending on hack type
# # by default, these write only to Overleaf dir
# pl1 = sim_plot_multiple_outcomes(.hack = "favor-best-affirm-wch",
#                                  .ggtitle = bquote( "SWS favors best affirmative; stringent SAS;" ~ mu ~ "= 0.5" ),
#                                  .local.results.dir = results.dir )
# 
# 
# pl2 = sim_plot_multiple_outcomes(.hack = "affirm",
#                                  .ggtitle = bquote( "SWS favors first affirmative; stringent SAS; " ~ mu ~ "= 0.5" ),
#                                  .local.results.dir = results.dir)
# 
# 
# 
# pl3 = sim_plot_multiple_outcomes(.hack = "affirm2",
#                                  .ggtitle = bquote( "SWS favors first affirmative; no SAS; " ~ mu ~ "= 0.5" ),
#                                  .local.results.dir = results.dir)






