

# IMPORTANT NOTES -----------------------------

# This script goes with 2022-3-11_applied_doParallel_SAPH.R.
#  That script is in the Sherlock folder b/c it's designed to be run via sbatch file. 

# SAPB-E documentation for these datasets:
# codebook_for_b2_data_prepped_step2.csv

# PRELIMINARIES ----------------------------------------

toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "boot",
           "metafor",
           "robumeta",
           "data.table",
           "purrr",
           "metRology",
           "fansi",
           "MetaUtility",
           "ICC",
           "cfdecomp",
           "tidyr",
           "truncdist",
           "tibble",
           "tmvtnorm",
           "testthat",
           "truncreg",
           "truncnorm",
           "rstan",
           "optimx",
           "weightr",
           "DataEditR",
           "plotly")
lapply( toLoad,
        require,
        character.only = TRUE)


# helper fns
code.dir = here("Sherlock code")
setwd(code.dir)
source("helper_SAPH.R")

results.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 SAPBE results from Sherlock/sapbe"

# ~ PREP DATA -----------------------------

# these datasets are directly copied from this dir:
# "~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Empirical benchmarks/Data collection/Prepped data for analysis"
setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 prep SAPBE dataset for SAPH/Datasets from SAPBE")
b2 = fread("**b2_data_prepped_step2.csv")
f2 = fread("**f2_data_aggregated_step2.csv")

# exclude duplicated Metalab sensitivity estimates
unique( b2$meta.name[ grepl(pattern = "[.]sens", b2$meta.name) == TRUE ] )
b2 = b2 %>% filter( grepl(pattern = "[.]sens", meta.name) == FALSE ) 

expect_equal(63, nuni(b2$meta.name))

# how much clustering?
# average number of studies per cluster
t = b2 %>% group_by(meta.name) %>%
  summarise( n.randomly.chosen = sum(randomly.chosen),
             n.total = n(),
             k.per.cluster = round(n.total/n.randomly.chosen, 3),
             # to see if dataset is public vs. scraped
             Data.source = Data.source[1] )
View(t)

# only analyze ones with relatively little clustering
meta.keepers = t$meta.name[ t$k.per.cluster < 1.2 ]

b2 = b2 %>% filter(meta.name %in% meta.keepers)
f2 = f2 %>% filter(meta.name %in% meta.keepers) 

nuni(b2$meta.name)

#**quick look at these metas' characteristics
t = f2 %>% select(meta.name, group, discipline,
                  k.all, k.rc, k.nonaffirm.rc,
                  Mhat, Mhat.Lo, Mhat.Hi,
                  Mhat.Worst, Mhat.Worst.Lo,
                  LogEta) %>%
  mutate_if(is.numeric, function(x) round(x,2))

View(t)

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 prep SAPBE dataset for SAPH/Datasets from SAPBE")

fwrite(t, "quick_summary_sapbe_meta_subset.csv")


# WRITE PREPPED DATA AND ALSO SEND TO CLUSTER -----------------------------

setwd("~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-11 prep SAPBE dataset for SAPH/Datasets prepped for SAPH")

fwrite(b2, "b2_long_prepped.csv")
fwrite(f2, "f2_short_prepped.csv")

# send to cluster
scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/Applied\ examples/2022-3-11\ prep\ SAPBE\ dataset\ for\ SAPH/Datasets\ prepped\ for\ SAPH/* mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/applied_examples/data/sapbe



# ANALYZE RESULTS FROM CLUSTER -----------------------------

setwd(results.dir)

r = fread("results_sapbe_all.csv")


# ~ View results in DataEditR ---------------------

temp = r %>% select(meta.name, SAPBE.group, method, Mhat, MLo, MHi,
                    SAPBE.Eta,
                    MhatRhat, ShatRhat, optim.converged,
                    optimx.Nagree.of.convergers.Mhat.winner) %>%
  mutate_if(is.numeric, function(x) round(x, 2))

data_edit(temp, viewer = "browser", code = TRUE)


# ~ Add Fit Diagnostics for RTMA  ---------------------

# need to work with estimate-level data to do ks test
temp = b2 %>% group_by(meta.name) %>%
  mutate( ks.pval = my_ks_test_RTMA( yi = EstF,
                                     sei = SE,
                                     Mhat = r$Mhat[ r$meta.name == meta.name &
                                                      r$method == "jeffreys-mcmc-pmed"],
                                     Shat = r$Shat[ r$meta.name == meta.name &
                                                      r$method == "jeffreys-mcmc-pmed" ] ) )


# # SAVE :)
# # sanity check for my_ks_test_RTMA: reproduce for a few metas
# .meta.name = "24988220"
# Mhat = r$Mhat[ r$meta.name == .meta.name &
#                  r$method == "jeffreys-mcmc-pmed" ]
# Shat = r$Shat[ r$meta.name == .meta.name &
#                  r$method == "jeffreys-mcmc-pmed" ]
# 
# dpn = b2 %>% filter( meta.name == .meta.name )
# mine = my_ks_test_RTMA( yi = dpn$EstF,
#                  sei = dpn$SE,
#                  Mhat = Mhat,
#                  Shat = Shat )
# 
# expect_equal( mine, unique( temp$ks.pval[ temp$meta.name == .meta.name ] ) )

# add the ks p-value to "r" dataset
r = r %>% rowwise() %>%
  mutate( ks.pval = unique( temp$ks.pval[ temp$meta.name == meta.name ] ) )

# new meta.name variable containing the p-val
r$meta.name.pretty = paste( r$meta.name, " (fit p=", round(r$ks.pval, 3), ")", sep = "" )

# Metas with strange results:
# pb_5: p < 0.0001 (bad fit)

# Metas with reasonable results:
# 24988220: p = 0.90


# ~ Plotly ---------------------

r = r %>% group_by(meta.name) %>%
  mutate( psm.gt.naive = ifelse( Mhat[ method == "2psm" ] > Mhat[ method == "naive" ],
                                 "2psm > naive",
                                 "2psm <= naive" ),
          pmed.gt.naive = ifelse( Mhat[ method == "jeffreys-mcmc-pmed" ] > Mhat[ method == "naive" ],
                                  "pmed > naive",
                                  "pmed <= naive") )

# **in half the metas, 2PSM exceeds naive
table( r$psm.gt.naive[ !duplicated(r$meta.name) ],
       r$pmed.gt.naive[ !duplicated(r$meta.name) ] )



# to choose axis limits
summary(r$Mhat)
xmin = -1
xmax = 4

my.shapes = c(16, 2)

p = ggplot( data = r,
            aes(x = Mhat,
                y = meta.name.pretty, 
                color = method,
                shape = (method == "naive" ) ) ) +
  
  geom_vline(xintercept = 0,
             lty = 2,
             color = "gray") +
  
  geom_point(size = 2,
             position = position_dodge(width = 0.5) ) +
  
  geom_errorbarh( aes(xmax = MHi, xmin = MLo),
                  width = 0,
                  position=position_dodge(width=0.5) ) +
  
  coord_cartesian( xlim = c(xmin, xmax) )+
  #scale_x_continuous(breaks=c(-3,0,3)) +
  
  scale_shape_manual(values = my.shapes) +
  
  theme_bw() +
  facet_wrap(pmed.gt.naive ~ psm.gt.naive,
             scales = "free",
             drop = TRUE)


pl = ggplotly(p)
pl


# how to save a plotly as html
# https://www.biostars.org/p/458325/
setwd(results.dir)
string = paste("sapbe_all_metas_plotly.html", sep="_")
htmlwidgets::saveWidget(pl, string)



# ~ Explore best-fitting meta ---------------------

.meta.name = "24988220"
Mhat = r$Mhat[ r$meta.name == .meta.name &
                 r$method == "jeffreys-mcmc-pmed" ]
Shat = r$Shat[ r$meta.name == .meta.name &
                 r$method == "jeffreys-mcmc-pmed" ]

dp = b2 %>% filter( meta.name == .meta.name )
r2 = r %>% filter( meta.name == .meta.name )



# ~~ Density plot *including* affirms -----
dp$yi = dp$EstF
dp$vi = dp$SE^2
plot_trunc_densities_RTMA(d = dp,
                          Mhat = Mhat,
                          Shat = Shat,
                          showAffirms = TRUE)

# ~~ One-tailed p-value histogram from PublicationBias ----

# not helpful compared to PublicationBias version
library(PublicationBias)
pval_plot(yi = dp$EstF, vi = dp$SE^2)

# ~~ Simple 2-tailed p-value density plot -----

xmin = 0
xmax = 1

p = ggplot(data = data.frame(x = c(xmin, 3)),
           aes(x)) +
  
  geom_vline(xintercept = 0.05,
             lwd = 2,
             color = "red") +
  
  geom_histogram(data = dp,
                 aes(x = Pval.Two))

# estimated density of estimates
geom_density( data = dp,
              aes(x = Pval.Two),
              adjust = .3 ) +
  
  # # estimated density from meta-analysis
  # stat_function( fun = dtrunc,
  #                n = 101,
  #                args = list( spec = "norm",
  #                             mean = 0,
  #                             sd = 1,
  #                             b = qnorm(.975) ),
  #                #aes(y = .25 * ..count..),  # doesn't work
  #                lwd = 1.2,
  #                color = "red") +


ylab("") +
  #scale_x_continuous( breaks = seq(xmin, 3, 0.5)) +
  xlab("Across-study Z-score using (Mhat, Shat)") +
  theme_minimal() +
  scale_y_continuous(breaks = NULL) +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16))


# ~~ Sort the p-values ----

sort(dp$Pval.Two)



# ~~ Info about dataset ----

r2 %>% select(SAPBE.Analysis.Scale,
              SAPBE.group,
              SAPBE.Data.source)

# topic: polymorphism vs. MI risk
# author conclusion: significant pub bias based on funnel plot & Egger

