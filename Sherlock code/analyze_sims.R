
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
library(sjmisc)
library(plotly)

# ~~ Set global vars needed throughout -------------------------
# no sci notation
options(scipen=999)

# control which results should be redone and/or overwritten
redo.plots = TRUE


# variables that define the scenarios (though some aren't actually manipulated, like stan.iter)
#@THIS IS ONLY THE ONES THAT I CHOSE TO VARY
param.vars.manipulated = c("k.pub.nonaffirm",
                           "rho",
                           "true.sei.expr")

# used later to create plots and tables, but needed to check var types 
#  upon reading in data
estNames = c("Mhat", "Shat")

outcomeNames = c("Bias", "RMSE", "EmpSE", "EstFail",
                 "Cover", "Width", "CIFail", "Rhat")


# ~~ Set directories -------------------------
code.dir = here("Sherlock code")


data.dir = str_replace( string = here(),
                        pattern = "Code \\(git\\)",
                        replacement = "Sherlock simulation results/Pilot simulations" )


results.dir = data.dir

# overleaf.dir.general = "/Users/mmathur/Dropbox/Apps/Overleaf/TNE (truncated normal estimation)/R_objects"
# overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/TNE (truncated normal estimation)/R_objects/figures"

setwd(code.dir)
source("helper_SAPH.R")
source("analyze_sims_helper_SAPH.R")



# ~~ Get iterate-level data -------------------------

setwd(data.dir)
# check when the dataset was last modified to make sure we're working with correct version
s = fread( "stitched.csv")

file.info("stitched.csv")$mtime

dim(s)


#**Check for MCMC errors and similar
#*# frequent errors for both jeffreys-var and mle-var, but not the corresponding SD param'zations
s %>% group_by(method) %>%
  summarise( meanNA(is.na(Mhat)))

#@@a lot of those cryptic cluster errors
unique( s$overall.error[ s$method == "jeffreys-var"] )

table(s$method, s$overall.error)


# ~~ Aggregate locally or read in aggregated data -------------------------

#agg = make_agg_data(s)


setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg)



# MAYBE SAVE? THESE WORK WHEN CONSIDERING 1 SCEN.
# # ~~ Main results -------------------------
# 
# t = agg %>% group_by(k.pub.nonaffirm, true.sei.expr, rho, method) %>%
#   select(method, 
#          sim.reps.actual,
#          all_of(names_with(agg, "Mhat")) ) %>%
#   mutate_if(is.numeric, function(x) round(x,2))
# 
# data.frame(t)
# View(t)
# 
# # look at just certain cols
# t = agg %>% select(method, 
#                    sim.reps.actual,
#                    all_of(names_with(agg, "Mhat")) ) %>%
#   mutate_if(is.numeric, function(x) round(x,2))
# 
# data.frame(t)
# View(t)
# 
# setwd(results.dir)
# fwrite(t, "results_all_iterates.csv")
# 
# 
# # ~~ Look at sanity checks -------------------------
# 
# # scenario diagnostics for scenario
# agg.checks = agg %>% select( 
#   all_of(names_with(agg, "sancheck.")) ) %>%
#   mutate_if(is.numeric, function(x) round(x,2))
# 
# # transpose it
# t = tibble( t(agg.checks[1,]) )
# t = t %>% add_column( variable = names(agg.checks), .before = 1 )
# 
# setwd(results.dir)
# fwrite(t, "sanity_check_stats.csv")



# 2022-3-5: SIMPLE PLOTS -------------------------

# **this is great
Ynames = rev( c("MhatBias", "MhatRMSE", "MhatCover", "MhatWidth",
                "MhatEstFail", "MhatCIFail",
                "MhatRhatGt1.01", "OptimxPropAgreeConvergersMhatWinner") )

for ( Yname in Ynames) {
  
  # to run "manually"
  #Yname = "MhatBias"
  
  p = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
                          .Yname = Yname,
                          .colorVarName = "method",
                          .facetVar1Name = "rho.pretty",
                          .facetVar2Name = "true.sei.expr.pretty",
                          .dat = agg,
                          .ggtitle = "",
                          .writePlot = FALSE,
                          .results.dir = NULL)
  
  # this is a great way to view plots!!
  pl = ggplotly(p)
  
  # how to save a plotly as html
  # https://www.biostars.org/p/458325/
  setwd(results.dir)
  string = paste(Yname, "_plotly.html", sep="")
  htmlwidgets::saveWidget(pl, string)
  
}




# 2022-3-5: INVESTIGATE 2PSM RESULTS -------------------------

# look at discrepancy in yi published affirms from hacked vs. unhacked studies
#  to help understand 2PSM results

# Conclusion:
# - rho has little effect on discrepancy (unexpected)
# - but true.sei.expr matters: when seis are tightly clustered around 0.55, the discrepancy is in the unexpected direction (i.e., hacked affirms are smaller)!

# hacked vs. unhacked yi for published affirmatives
t = s %>% group_by( k.pub.nonaffirm, true.sei.expr, rho ) %>%
  summarise(              Power = meanNA(sancheck.prob.unhacked.udraws.affirm),
                          Mean.yi.hacked = meanNA(sancheck.mean.yi.hacked.pub.affirm),
                          Mean.yi.unhacked = meanNA(sancheck.mean.yi.unhacked.pub.affirm),
                          Discrep = Mean.yi.hacked - Mean.yi.unhacked,
                          Discrep.ratio = (Mean.yi.hacked - Mean.yi.unhacked) / abs(Mean.yi.unhacked) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

View(t)
setwd(results.dir)
fwrite(t, "table_hacked_vs_unhacked_pub_affirms.csv")

# effect of true.sei.expr on power
t = s %>% group_by( true.sei.expr, rho ) %>%
  summarise_at( all_of( c("sancheck.dp.meanN.hacked",
                          "sancheck.prob.hacked.udraws.affirm",
                          "sancheck.prob.unhacked.udraws.affirm",
                          "sancheck.prob.hacked.ustudies.published") ),
                function(x) meanNA(x) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

View(t)
setwd(results.dir)
fwrite(t, "table_underlying_draw_power.csv")
#**save this after re-running

# questions:
# - why is sancheck.prob.hacked.udraws.affirm != sancheck.prob.unhacked.udraws.affirm even when rho = 0? It's because hacked studies make more draws when they have smaller mui, so small-mui hacked studies are overrepresented.
# - with these settings, first draws have 13-17% power 
# - rho has little effect on anything because it hardly changes how many draws the hacked studies make

# 2022-3-6: EXPLORE NMAX -------------------------

# how hard is it to get affirmative result?
t = s %>% filter()
group_by( true.sei.expr, rho ) %>%
  summarise_at( all_of( c("sancheck.dp.meanN.hacked",
                          "sancheck.prob.hacked.udraws.affirm",
                          "sancheck.prob.unhacked.udraws.affirm",
                          "sancheck.prob.hacked.ustudies.published") ),
                function(x) meanNA(x) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

View(t)



# 2022-3-2: EXPLORE BAD ITERATES -------------------------


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
thresh = 6
meanNA(s$optimx.Nagree.of.convergers.Mhat.winner > thresh &
         s$optimx.Pagree.of.convergers.Mhat.winner == 1)



### Mark iterates with good optimx properties for ALL methods ###
non.optimx.methods = c("naive", "gold-std", "maon", "2psm", "jeffreys-mcmc-pmean", 
                       "jeffreys-mcmc-pmed", "jeffreys-mcmc-max-lp-iterate")


s$job.rep.name = paste( s$job.name, s$rep.name, sep = "_" )



# for the optimx methods, mark whether the iterate should be retained for ALL methods
#  "min" below is to consider performance across all methods that used optimx
s = s %>% group_by(job.rep.name) %>%
  mutate( keep.iterate = min(optimx.Nagree.of.convergers.Mhat.winner, na.rm = TRUE) > thresh &
            min(optimx.Pagree.of.convergers.Mhat.winner, na.rm = TRUE) == 1 )

mean(s$keep.iterate)

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