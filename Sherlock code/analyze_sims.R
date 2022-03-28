
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
redo.plots = TRUE


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
# MAY TAKE A LONG TIME!
#s = fread( "stitched.csv")

file.info("stitched.csv")$mtime

dim(s)
nuni(s$scen.name)

#**Check for MCMC errors and similar
#*# frequent errors for both jeffreys-var and mle-var, but not the corresponding SD param'zations
s %>% group_by(method) %>%
  summarise( meanNA(is.na(Mhat)))

# # look at errors
# unique( s$overall.error[ s$method == "jeffreys-var"] )
# table(s$method, s$overall.error)


# ~~ Get agg data -------------------------

#agg = make_agg_data(s)


setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime

dim(agg)
nuni(agg$scen.name)

agg = wrangle_agg_local(agg)



# ~~ List variable names -------------------------

### Names of statistical metrics ###
# used later to create plots and tables, but needed to check var types 
#  upon reading in data
estNames = c("Mhat", "Shat")

# blank entry is to get Mhat itself, which is useful for 
#  looking at whether MAON>0
mainYNames = c("Bias", "", "RMSE", "Cover", "Width", "EmpSE")

otherYNames = c("EstFail", "CIFail", "RhatGt1.01", "RhatGt1.05")

# these ones don't fit in nicely because the "Mhat" is in the middle of string
#"OptimxPropAgreeConvergersMhatWinner", "OptimxNAgreeOfConvergersMhatWinner"
MhatMainYNames = paste( "Mhat", c(mainYNames), sep = "" )
MhatYNames = c( paste( "Mhat", c(mainYNames, otherYNames), sep = "" ),
                "OptimxPropAgreeConvergersMhatWinner", "OptimxNAgreeOfConvergersMhatWinner" )


### Names of parameter variables ###
# figure out which scen params were actually manipulated
#@this assumes that "Nmax" is always the first param var and "method" is always the last
( param.vars = names(agg)[ which( names(agg) == "Nmax" ) : which( names(agg) == "method" ) ] )

# how many levels does each param var have in dataset?
( n.levels = agg %>% dplyr::select(param.vars) %>%
    summarise_all( function(x) nuni(x) ) )

( param.vars.manip = names(n.levels)[ n.levels > 1 ] )


# eliminate redundant ones
if ( "t2a" %in% param.vars.manip ) param.vars.manip = drop_vec_elements( param.vars.manip, c("S", "V") )


( param.vars.manip2 = drop_vec_elements(param.vars.manip, "method") )


# # sanity check: SDs of all analysis variables should be 0 within unique scenarios
# t = data.frame( s3 %>% group_by(unique.scen) %>%
#                   summarise_at( analysis.vars, sd ) )





# SORT ROWS BY PERFORMANCE -------------------------


# ~ Coverage -------------------------
t.sort = sort_agg(MhatCover, desc = FALSE)

# sort by method, then performance within method
t.sort = t.sort %>% arrange(method, MhatCover)

View(t.sort)

setwd(results.dir)
write.xlsx( as.data.frame(t.sort), "agg_sorted_by_method_and_coverage.xlsx")

# ~ Bias (Signed) -------------------------

t.sort = sort_agg(MhatBias, desc = TRUE)

# sort by method, then performance within method
t.sort = t.sort %>% arrange(method, desc(MhatBias))

View(t.sort)

setwd(results.dir)
write.xlsx( as.data.frame(t.sort), "agg_sorted_by_method_and_coverage.xlsx")

# REGRESS PERFORMANCE ON MANIPULATED SCEN PARS  -------------------------

.method = "jeffreys-mcmc-pmed"

# ~ MhatCover --------------------
RHS = paste( param.vars.manip2, collapse = " + " )
formula = paste( "MhatCover", RHS, sep = " ~ ")
mod = lm( eval( parse( text = formula) ),
          data = agg %>% filter(method == .method) )

summary(mod)


# ~ MhatBias --------------------
RHS = paste( param.vars.manip2, collapse = " + " )
formula = paste( "abs(MhatBias)", RHS, sep = " ~ ")
mod = lm( eval( parse( text = formula) ),
          data = agg %>% filter(method == .method) )


summary(mod)



# Conclusions
# Things that HURT 2PSM performance
# -	Wider range of SEs (e.g., U[1, 3] vs. U[0.1,3]; exponential is especially bad)
# -	Larger t2a
# -	Higher prob.hacked
# -	k.pub.nonaffirm doesn’t matter


# ******** PLOTS -------------------------


Ynames = rev(MhatYNames)

# alternatively, run just a subset:
Ynames = c("MhatWidth", "MhatCover", "MhatBias",
           # last 2 are useful for looking at MAON
           "Mhat", "MhatTestReject")

#@temp if not running optimx:
#Ynames = Ynames[3:11]

# to help decide which vars to include in plot:
param.vars.manip2


# in case you want to filter scens:
# full set for reference:
# c("naive", "gold-std", "maon", "2psm", "pcurve", "jeffreys-mcmc-pmean", 
#   "jeffreys-mcmc-pmed", "jeffreys-mcmc-max-lp-iterate", "jeffreys-sd", 
#   "jeffreys-var", "mle-sd", "csm-mle-sd", "mle-var", "2psm-csm-dataset", 
#   "prereg-naive", "ltn-mle-sd")
( all.methods = unique(agg$method) )
toDrop = c("jeffreys-mcmc-pmean", "jeffreys-mcmc-max-lp-iterate")
method.keepers = all.methods[ !all.methods %in% toDrop ]



aggp = agg %>% filter(method %in% method.keepers &
                        Mu == 0.5 &
                        prob.hacked == 0.8 &
                        hack == "favor-best-affirm-wch")
# to label the plots
prefix = "hack=favor-best; pr.hack=0.8"


# for 2022-3-25 sims
aggp$tempFacetVar2 = paste( "t2a=", aggp$t2a, "; t2w=", aggp$t2w, sep = "")
table(aggp$tempFacetVar2)

# aggp$tempFacetVar1 = paste( "pr.hack=", aggp$prob.hacked, sep = "")
# table(aggp$tempFacetVar1)


for ( Yname in Ynames) {
  
  # to run "manually"
  #Yname = "MhatBias"
  #Yname = "MhatCover"
  
  y.breaks = NULL
  if ( Yname == "MhatBias") y.breaks = seq(-0.5, 0.5, 0.1)
  if ( Yname == "MhatWidth") y.breaks = seq(0, 10, 0.5)
  
  p = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
                          .Yname = Yname,
                          .colorVarName = "method",
                          #.facetVar1Name = "tempFacetVar1",
                          .facetVar1Name = "true.sei.expr.pretty",
                          .facetVar2Name = "tempFacetVar2",
                          .dat = aggp,
                          .ggtitle = prefix,
                          .y.breaks = y.breaks,
                          .writePlot = FALSE,
                          .results.dir = NULL)
  
  # # for 2022-3-16
  # p = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
  #                         .Yname = Yname,
  #                         .colorVarName = "method",
  #                         .facetVar1Name = "Nmax",
  #                         .facetVar2Name = "hack",
  #                         .dat = aggp,
  #                         .ggtitle = "",
  #                         .writePlot = FALSE,
  #                         .results.dir = NULL)
  
  
  # for 2022-3-8 sims
  # p = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
  #                         .Yname = Yname,
  #                         .colorVarName = "method",
  #                         .facetVar1Name = "tempFacetVar",
  #                         .facetVar2Name = "true.sei.expr.pretty",
  #                         .dat = aggp,
  #                         .ggtitle = "",
  #                         .writePlot = FALSE,
  #                         .results.dir = NULL)
  
  # # SAVE: this was for the 2022-3-7 and earlier sims (based on what they manipulated)
  # p = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
  #                         .Yname = Yname,
  #                         .colorVarName = "method",
  #                         .facetVar1Name = "rho.pretty",
  #                         .facetVar2Name = "true.sei.expr.pretty",
  #                         .dat = agg,
  #                         .ggtitle = "",
  #                         .writePlot = FALSE,
  #                         .results.dir = NULL)
  
  # this is a great way to view plots!!
  pl = ggplotly(p)
  pl
  
  # how to save a plotly as html
  # https://www.biostars.org/p/458325/
  setwd(results.dir)
  string = paste(prefix, Yname, "plotly.html", sep="_")
  htmlwidgets::saveWidget(pl, string)
  
}




# EFFECT OF SCEN PARAMS ON DATASETS -------------------------

# look at discrepancy in yi published affirms from hacked vs. unhacked studies
#  to help understand 2PSM results

# Conclusion:
# - rho has little effect on discrepancy (unexpected)
# - but true.sei.expr matters: when seis are tightly clustered around 0.55, the discrepancy is in the unexpected direction (i.e., hacked affirms are smaller)!

### Published affirmatives: hacked vs. unhacked yi ###
param.vars.manip2 = drop_vec_elements(param.vars.manip, "method")

t.affirm = s %>% group_by_at( param.vars.manip2 ) %>%
  
  summarise( scen.name = as.numeric( paste( unique(scen.name) ) ),  # should only be 1 scen
             Power = meanNA(sancheck.prob.unhacked.udraws.affirm),
             Mean.yi.hacked = meanNA(sancheck.mean.yi.hacked.pub.affirm),
             Mean.yi.unhacked = meanNA(sancheck.mean.yi.unhacked.pub.affirm),
             Discrep = Mean.yi.hacked - Mean.yi.unhacked,
             Discrep.ratio = (Mean.yi.hacked - Mean.yi.unhacked) / abs(Mean.yi.unhacked) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

View(t.affirm)
setwd(results.dir)
write.xlsx( as.data.frame(t.affirm), "table_hacked_vs_unhacked_pub_affirms.xlsx")


### Published NONaffirmatives: hacked vs. unhacked yi ###
t.nonaffirm = s %>% group_by_at( param.vars.manip2 ) %>%
  
  summarise( scen.name = as.numeric( paste( unique(scen.name) ) ),  # should only be 1 scen
             Power = meanNA(sancheck.prob.unhacked.udraws.affirm),
             Mean.yi.hacked = meanNA(sancheck.mean.yi.hacked.pub.nonaffirm),
             Mean.yi.unhacked = meanNA(sancheck.mean.yi.unhacked.pub.nonaffirm),
             Discrep = Mean.yi.hacked - Mean.yi.unhacked,
             Discrep.ratio = (Mean.yi.hacked - Mean.yi.unhacked) / abs(Mean.yi.unhacked) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

View(t.nonaffirm)
setwd(results.dir)
write.xlsx( as.data.frame(t.nonaffirm), "table_hacked_vs_unhacked_pub_nonaffirms.xlsx")


### 2022-3-8: look at scen where 2PSM was especially an underestimate ###
# find scenario names 
agg %>% filter( true.sei.expr == "0.1 + rexp(n = 1, rate = 1.5)" &
                  prob.hacked == 0.8 &
                  Mu == 0.1 &
                  t2a == 1 ) %>%
  select(scen.name, method, MhatBias, MhatCover)

data.frame( t.affirm %>% filter(scen.name == 32 ) )

data.frame( t.nonaffirm %>% filter(scen.name == 32 ) )
#bm
#***VERY INTERESTING THAT IN SCEN 32, THE 2PSM IS BIASED DOWNWARD, YET
# mean.yi.hacked affirm is LARGER in hacked studies than unhacked (2 vs. 1.51)


### **Effect of true.sei.expr on power ###
t = s %>% group_by( true.sei.expr, rho ) %>%
  summarise_at( all_of( c("sancheck.dp.meanN.hacked",
                          #"sancheck.dp.q90N.hacked",  #@add after next round of sims
                          "sancheck.prob.hacked.udraws.affirm",
                          "sancheck.prob.unhacked.udraws.affirm",
                          "sancheck.prob.hacked.ustudies.published"
                          #"sancheck.prob.published.affirm.is.hacked"
  ) ),
  function(x) meanNA(x) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

View(t)
setwd(results.dir)
write.xlsx( as.data.frame(t), "table_underlying_draw_power.xlsx")


# questions:
# - why is sancheck.prob.hacked.udraws.affirm != sancheck.prob.unhacked.udraws.affirm even when rho = 0? It's because hacked studies make more draws when they have smaller mui, so small-mui hacked studies are overrepresented.
# - with these settings, first draws have 13-17% power 
# - rho has little effect on anything because it hardly changes how many draws the hacked studies make


# 2022-3-24: AGAIN TRY TO REPLICATION 2022-3-8 -------------------------

# direct replication of 2022-3-8
agg %>% filter( Mu == 0.5 & 
                  t2a == 1.5 & 
                  k.pub.nonaffirm == 50 &
                  prob.hacked == 0.5 &
                  method == "jeffreys-mcmc-pmed" &
                  true.sei.expr == "0.1 + rexp(n = 1, rate = 1.5)") %>%
  
  select(scen.name,
         hack,
         Mhat,
         MhatBias,
         MhatCover,
         MhatWidth) %>%
  mutate_if(is.numeric, function(x) round(x,2))


# look at other true.sei.expr as well
agg %>% filter( Mu == 0.5 & 
                  t2a == 1.5 & 
                  k.pub.nonaffirm == 50 &
                  prob.hacked == 0.5 &
                  method == "jeffreys-mcmc-pmed") %>%
  
  select(scen.name,
         true.sei.expr,
         hack,
         Mhat,
         MhatBias,
         MhatCover,
         MhatWidth) %>%
  mutate_if(is.numeric, function(x) round(x,2))



# 2022-3-22: WHY JEFFREYS-MCMC-PMED SUDDENLY BAD FOR HACK TYPE AFFIRM (WCH)? -------------------------

agg %>% filter( Mu == 0.5 & 
                  t2a == 1.5 & 
                  k.pub.nonaffirm == 50 &
                  prob.hacked == 0.5 &
                  method == "jeffreys-mcmc-pmed") %>%
  select(scen.name,
         hack,
         Mhat,
         MhatBias,
         MhatCover,
         MhatWidth) %>%
  mutate_if(is.numeric, function(x) round(x,2))


# compare to 3-8 sims
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Sherlock simulation results/Pilot simulations/2022-3-8")
aggo = fread("agg.csv")

aggo %>% filter( Mu == 0.5 & 
                   t2a == 1.5 & 
                   k.pub.nonaffirm == 50 &
                   prob.hacked == 0.5 &
                   method == "jeffreys-mcmc-pmed") %>%
  select(scen.name,
         true.sei.expr,
         hack,
         Mhat,
         MhatBias,
         MhatCover,
         MhatWidth) %>%
  mutate_if(is.numeric, function(x) round(x,2))


# One issue is that you were looking at plots for a diffent hack type.
# Also, now the CIs are super wide. 
# Otherwise, I think the difference must be in the sei's.
# Current sei distribution is much smaller and also less variable.
# Example:

x1 = 0.1 + rexp(n = 1000, rate = 1.5)
x2 = rbeta(n = 1000, 2, 5) 

summary(x1); sd(x1)
summary(x2); sd(x2)

hist(x1)
hist(x2)

#BM: IDEA TO HELP DIAGNOSE:
# Add back the previous sei.expr 

# PROBABLY UNUSED (MOVE TO NEW FILE?) -------------------------


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



# ISOLATE SCENS THAT WERE BAD FOR 2PSM -------------------------

# I picked scens using the sorted table:
View(t.sort)


# ~ Scen with large Mu, large t2a, and expo SEs -------------------------
scen = 144
t %>% filter(scen.name == scen)
# **this is a great scenario because Jeffreys isn't even THAT wide (width 2.54) and has coverage 94%, but 2PSM has coverage 10%!

# look at power characteristics of this scen
data.frame( t.affirm %>% filter(scen.name == scen) )



# ~ Scen with large Mu, large t2a, and less heterogeneity than above -------------------------
scen = 125
t %>% filter(scen.name == scen)
# **this is a great scenario because Jeffreys isn't even THAT wide (width 2.54) and has coverage 94%, but 2PSM has coverage 10%!

# look at power characteristics of this scen
data.frame( t.affirm %>% filter(scen.name == scen) )


# ~ Scen with large Mu, large t2a, and unif SEs -------------------------
scen = 150
t %>% filter(scen.name == scen)

# **also a great scen!

# look at power characteristics of this scen
data.frame( t.affirm %>% filter(scen.name == scen) )

#**Q: Are there cases where 2PSM is badly biased downward that even MAON is less conservative than 2PSM?



# 2022-3-16: Specific to these sims -------------------------

t = s %>% group_by(scen.name, method) %>%
  summarise(Nmax = Nmax[1],
            hack = hack[1],
            MhatMn = meanNA(Mhat),
            MhatMed = median(Mhat),
            Mhat2.5 = quantile(Mhat, 0.025),
            Mhat97.5 = quantile(Mhat, 0.975) ) %>%
  mutate_if(is.numeric, function(x) round(x,2))


View(t)

setwd(results.dir)
write.xlsx( as.data.frame(t), "summary_results.xlsx")

Mhats = s$Mhat[s$scen.name == 4 & s$method == "csm-mle-sd"]
mean(Mhats)

mean(Mhats < 0.4)

# the issue is NOT that there are too few affirmatives
s %>% filter(scen.name == 4 &
               method == "csm-mle-sd" &
               Mhat < -10 ) %>%
  summarise(sancheck.dp.k.affirm)



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