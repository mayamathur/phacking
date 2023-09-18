

# NOTES ----------------------------------------------------

# Only results from sim.env = stefan are in paper per reviewers' comments.
#  Analyses with sim.env = mathur retained here for completeness.


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


stitch.from.scratch = TRUE

# control which results should be redone and/or overwritten
#@ not all fns respect this setting
overwrite.res = TRUE


# ~~ Set directories -------------------------
code.dir = here("Sherlock code")

# dataset prepped by prep_sims_SAPH.R
data.dir = str_replace( string = here(),
                           pattern = "Code \\(git\\)",
                           replacement = "Simulation results/2023-06-21 aggregated simulations (as in RSM_1)" )


overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/P-hacking (SAPH)/figures_SAPH/sims"


setwd(code.dir)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")




# ~~ Get agg data -------------------------

# if only analyzing a single sim environment (no merging):
setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg)


# drop any "NA" methods (i.e., ones that didn't get labeled in wrangle_agg_local)
agg = agg %>% filter( method.pretty != "" )

# look at number of actual sim reps
table(agg$sim.reps.actual)


#@important: stefan and mathur share scen.names, so relabel them
agg$scen.name2 = paste(agg$sim.env, agg$scen.name)
# stefan: 40 scens, mathur: 48
expect_equal( nuni(agg$scen.name2), 40 + 48)
 


# ~~ List variable names -------------------------

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()


# ~~ Make data subsets -------------------------

# stefan only
aggs = agg %>% filter(sim.env == "stefan")
# mathur only
aggm = agg %>% filter(sim.env == "mathur")

# ~~ Convergence stats by method -------------------------

summary(aggs$MhatEstConverge)
summary(aggs$MhatCIFail)

aggs %>% group_by(method)


# convergence rates
t = aggs %>% group_by(method) %>%
  summarise( mean(1-MhatEstFail), 
             min(1-MhatEstFail),
             
             mean(1-MhatCIFail),
             min(1-MhatCIFail) )

View(t)



# ******** WINNER TABLES -------------------------

# ~ Winner tables shown in paper ------------------------------

# can toggle output of fn below by changing the default arg of 
#  make_winner_table between display = "dataframe" (easy viewing)
#  and display = "xtable" (Overleaf)



# ~~~ Stefan ------------------------------
# 1: all scenarios
make_both_winner_tables(.agg = aggs)

# by k
make_both_winner_tables(.agg = aggs %>% filter(k.pub.nonaffirm == 10) )


# by two- vs. one-tailed selection
make_both_winner_tables(.agg = aggs %>% filter(alternative.stefan == "two.sided") )
make_both_winner_tables(.agg = aggs %>% filter(alternative.stefan == "greater") )

# by strategy (favor-first vs. favor-best)
make_both_winner_tables(.agg = aggs %>% filter(strategy.stefan == "firstsig") )
make_both_winner_tables(.agg = aggs %>% filter(strategy.stefan == "smallest") )



# ~~~ Mathur  ------------------------------

# not in paper per reviewers

# 1: all scenarios
make_both_winner_tables(.agg = aggm)
# here, reason SM-step appears unbiased is that it's positively biased under evil.selection=0
#   but negatively biased under evil.selection=1
#bm: wrong number of scens

# 2: by k
make_both_winner_tables(.agg = aggm %>% filter(k.pub.nonaffirm == 10) )
make_both_winner_tables(.agg = aggm %>% filter(k.pub.nonaffirm == 100) )



# 3: evil.selection (misspecified)
make_both_winner_tables(.agg = aggm %>% filter(rtma.misspec == TRUE))

make_both_winner_tables(.agg = aggm %>% filter(rtma.misspec == FALSE))



# # ~~~ Mathur, heterogeneity subsets  ------------------------------
# 
# make_both_winner_tables(.agg = aggm %>% filter(t2a == 0) )
# 





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






