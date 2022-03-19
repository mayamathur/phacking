

# IMPORTANT NOTES -----------------------------

# This script goes with 2022-3-11_applied_doParallel_SAPH.R.
#  That script is in the Sherlock folder b/c it's designed to be run via sbatch file. 

# SAPB-E documentation for these datasets:
# codebook_for_b2_data_prepped_step2.csv

# PRELIMINARIES ----------------------------------------


# ~ User-Specified Parameters --------------------

# "sapbe", "kvarven"
dataset.name = "sapbe"

# ~ Packages, Etc.  --------------------
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
           "plotly",
           "readxl")
lapply( toLoad,
        require,
        character.only = TRUE)


# helper fns
code.dir = here("Sherlock code")
setwd(code.dir)
source("helper_SAPH.R")


if ( dataset.name == "sapbe" ) {
  raw.data.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/SAPBE/Datasets from SAPBE"
  
  prepped.data.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/SAPBE/Prepped datasets for SAPH"
  
  results.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/SAPBE/Results from Sherlock"
}

if ( dataset.name == "kvarven" ) {
  data.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-13 prep Kvarven dataset for SAPH"
  results.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-12 Kvarven results from Sherlock"
  
}

if ( dataset.name == "olsson" ) {
  data.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/Olsson-Collentine"
  results.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/Olsson-Collentine/Results from Sherlock"
  
}



# ~ PREP DATA -----------------------------


# ~ Prep SAPB-E Metas -----------------------------


if ( dataset.name == "sapbe" ) {
  # these datasets are directly copied from this dir:
  # "~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Empirical benchmarks/Data collection/Prepped data for analysis"
  setwd(raw.data.dir)
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
  
  # make columns with standardized names to match doParallel from sim study
  b2$yi = b2$EstF  # **uses the direction-flipped estimates
  b2$vi = b2$SE^2
  b2$sei = b2$SE
  b2$Zi = b2$yi / sqrt(b2$vi)
  #**note we're using all estimates, not just randomly-chosen ones, because we're focusing on metas with little clustering
  
  #**quick look at these metas' characteristics
  t = f2 %>% select(meta.name, group, discipline,
                    k.all, k.rc, k.nonaffirm.rc,
                    Mhat, Mhat.Lo, Mhat.Hi,
                    Mhat.Worst, Mhat.Worst.Lo,
                    LogEta) %>%
    mutate_if(is.numeric, function(x) round(x,2))
  
  View(t)
  
  setwd(prepped.data.dir)
  fwrite(t, "quick_summary_sapbe_meta_subset.csv")
  
  ### Write Prepped Data ###
  setwd(prepped.data.dir)
  fwrite(b2, "b2_long_prepped.csv")
  fwrite(f2, "f2_short_prepped.csv")
}




# ~ Prep Kvarven Metas -----------------------------

# this code is mostly taken from our RSOS paper code: analysis_code/analysis_2/analyze.R
if ( dataset.name == "kvarven" ) {
  
  setwd(data.dir)
  setwd("Datasets from Kvarven")
  
  # Kvarven's summary results containing replications
  agg = read_xls("Dataset.xls")
  # remove blank rows
  agg = agg[ !is.na(agg$metaanalysis), ]
  # for later merging joy: identify the vars that came from Kvarven's own dataset
  names(agg) = paste( "kvarven.", names(agg), sep = "" )
  # make merger variable with just first author's last name
  agg$meta = unlist( lapply( agg$kvarven.metaanalysis, FUN = function(x) strsplit( x, " " )[[1]][1] ) )
  
  
  # list all meta-analysis datasets
  setwd(data.dir)
  setwd("Datasets from Kvarven/Meta")
  files = list.files()[ grepl( pattern = ".csv", x = list.files() ) ]
  
  # separation character different for different meta-analyses
  sep.vec = rep( ";", length(files) )
  sep.vec[1] = ","
  sep.vec[10] = ","
  
  
  ##### Read In Each Meta and Add to Combined Dataset #####
  
  for ( i in 1:length(files) ) {
    
    cat( paste("\n ****Starting", files[i] ) )
    
    setwd(data.dir)
    setwd("Datasets from Kvarven/Meta")
    
    # d and var are always the first 2 columns, even though they are named
    #  differently
    # EXCEPT Schimmack, which uses Cohen's q and has cols in different order
    if ( files[i] != "Schimmack.csv" ) {
      d = read.csv( files[i], sep = sep.vec[i] )[,1:2]
      print( names(d) )  # as a sanity check
      # standardize names
      names(d) = c("d", "var")
    } 
    
    if ( files[i] == "Schimmack.csv" ) {
      d = read.csv( files[i], sep = ";" )
      d = d %>% select( q, Var.q. )
      # not actually Cohen's d, but irrelevant for the analyses we're going to do
      names(d) = c("d", "var")
    }
    
    # extract name of meta-analysis
    d = d %>% add_column( meta = str_remove_all( files[i], ".csv" ),
                          .before = 1 )
    
    # merge in other variables from Kvarven's own dataset
    d = merge( d, agg, all.x = TRUE, by = "meta" ) 
    
    # add it to combined dataset, named for consistency with SAPB-E
    if ( i == 1 ) b2 = d else b2 = rbind(b2, d) 
    
    
    cat( paste("\n ****Done", files[i] ) )
  }
  
  table(b2$meta)
  
  
  ### Sanity Check ###
  
  # have all pooled estimates been coded as positive?
  expect_equal( any( b2$kvarven.meta_s < 0 ), FALSE )
  
  ### Write Prepped Data ###
  setwd(data.dir)
  setwd("Datasets prepped for SAPH")
  
  fwrite(b2, "b2_long_prepped_kvarven.csv")
  
} # end "if(dataset.name == "kvarven")



# send to cluster (specific to SAPB-E)
# scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/Applied\ examples/2022-3-11\ prep\ SAPBE\ dataset\ for\ SAPH/Datasets\ prepped\ for\ SAPH/* mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/applied_examples/data/sapbe




# (NOW RUN ANALYSIS CODE ON CLUSTER) -----------------------------






# ANALYZE RESULTS FROM CLUSTER -----------------------------



if ( dataset.name == "sapbe" ) {
  # dataset with all methods prior to csm and ltn
  setwd(results.dir)
  setwd("2022-3-11")
  r1 = fread("2022-3-11_results_sapbe_all.csv")
  
  # dataset with just those 2 methods
  setwd(results.dir)
  setwd("2022-3-16 just CSM and LTMA")
  r2 = fread("2022-3-16_results_sapbe_all.csv")
  
  r = bind_rows(r1, r2) %>% arrange(meta.name)
  
  
  ### Get Prepped Meta-Level Data ###
  setwd(prepped.data.dir)
  b2 = fread("b2_long_prepped.csv")
}
  
if ( dataset.name == "kvarven" ) {
  setwd(results.dir)
  r = fread("results_kvarven_all.csv")
}


# ~ View results in DataEditR ---------------------

if ( dataset.name == "sapbe" ) {
  temp = r %>% select(meta.name, SAPBE.group, method, Mhat, MLo, MHi,
                      SAPBE.Eta,
                      MhatRhat, ShatRhat, optim.converged,
                      optimx.Nagree.of.convergers.Mhat.winner) %>%
    mutate_if(is.numeric, function(x) round(x, 2))
}

if ( dataset.name == "kvarven" ) {
  temp = r %>% select(meta.name, Mhat, MLo, MHi,
                      MhatRhat, ShatRhat, optim.converged,
                      optimx.Nagree.of.convergers.Mhat.winner) %>%
    mutate_if(is.numeric, function(x) round(x, 2))
}

data_edit(temp, viewer = "browser", code = TRUE)


# ~ Fit Diagnostics: QQ Plot  ---------------------

plot.method = "jeffreys-mcmc-pmed"

r2 = r %>% filter( method == plot.method & 
                     !is.na(Mhat) )

b3 = b2 %>% filter( )

for ( .m in unique(r2$meta.name) ) {
  
  yi = b2$yi[ b2$meta.name == .m & b2$affirm == FALSE ]
  sei = sqrt(b2$vi[ b2$meta.name == .m & b2$affirm == FALSE ])
  Mhat = r$Mhat[ r$meta.name == .m & r$method == plot.method ]
  Shat = r$Shat[ r$meta.name == .m & r$method == plot.method ]
  
  p1 = yi_qqplot(yi = yi,
            sei = sei, 
            Mhat = Mhat,
            Shat = Shat)
  
  
  
  ### SEEMS CRAZY!
  p2 = plot_trunc_densities_RTMA(d = b2 %>% filter(meta.name == .m),
                                Mhat = Mhat,
                                Shat = Shat,
                                showAffirms = FALSE)
  
  setwd(sherlock.results.dir)
  ggsave( filename = paste( "density_plot", i, ".pdf", sep="_" ),
          width = 10, 
          height = 10)
  
  
  
}




# ~ Add Fit Diagnostics for RTMA  ---------------------

#bm: maybe add ks test for 2PSM here instead?
# need to work with estimate-level data to do ks test

if ( dataset.name == "sapbe" ) {
  #@2022-3-14: THIS IS FROM BEFORE I UPDATED MY_KS_TEST; NOW WILL NEED TO PASS ONLY NONAFFIRMS
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
  
}


if ( dataset.name == "kvarven" ) {
  
  temp = suppressWarnings(  b2 %>% group_by(meta.name) %>%
                              filter(affirm == FALSE) %>%
                              mutate( ks.pval = my_ks_test_RTMA( yi = yi,
                                                                 sei = sqrt(vi),
                                                                 Mhat = r$Mhat[ r$meta.name == meta.name &
                                                                                  r$method == "jeffreys-mcmc-pmed"],
                                                                 Shat = r$Shat[ r$meta.name == meta.name &
                                                                                  r$method == "jeffreys-mcmc-pmed" ] ) ) )
  
  #bm: need to finish checking this and then need to edit above to get ks.pval for each method (jeffreys and 2psm);
  #  jeffreys will only have p-value for nonaffirms; 2psm will have p-value for affirms and nonaffirms
  # SAVE :)
  # sanity check for my_ks_test_RTMA: reproduce for a few metas
  .meta.name = "Meissner"
  Mhat = r$Mhat[ r$meta.name == .meta.name &
                   r$method == "jeffreys-mcmc-pmed" ]
  Shat = r$Shat[ r$meta.name == .meta.name &
                   r$method == "jeffreys-mcmc-pmed" ]
  
  b3 = b2 %>% filter( meta.name == .meta.name ) %>% filter(affirm == FALSE)
  nrow(b3)
  mine = my_ks_test_RTMA( yi = b3$yi,
                          sei = sqrt(b3$vi),
                          Mhat = Mhat,
                          Shat = Shat )
  
  expect_equal( mine, unique(temp$ks.pval[ temp$meta.name == .meta.name ]) )
}


# add the ks p-value to "r" dataset
r = r %>% rowwise() %>%
  mutate( ks.pval = temp$ks.pval[ temp$meta.name == meta.name ][1] )

# new meta.name variable containing the p-val
r$meta.name.pretty = paste( r$meta.name, " (fit p=", round(r$ks.pval, 3), ")", sep = "" )




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

if ( dataset.name == "sapbe" ) {
  xmin = -1
  xmax = 8
}

if ( dataset.name == "kvarven" ) {
  xmin = -0.5
  xmax = 4
}

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
string = paste("plotly_all_metas.html", sep="_")
htmlwidgets::saveWidget(pl, string)



# SPECIFIC TO KVARVEN ----------------------

# ~ Explore a horrible meta -------------------

b3 = b2 %>% filter(meta.name == "Hagger")


View( r %>% filter(meta.name == "Hagger") %>% select(method, Mhat, MhatRhat, optimx.Mhat.winner, optimx.Shat.winner, optimx.Nagree.of.convergers.Mhat.winner) )
 
# look at 2PSM's fitted density (reasonable) instead of RTMA's
Mhat = r$Mhat[ r$meta.name == "Hagger" & r$method == "2psm" ]
Shat = r$Shat[ r$meta.name == "Hagger" & r$method == "2psm" ]

plot_trunc_densities_RTMA(d = b3,
                          Mhat = Mhat,
                          Shat = Shat,
                          showAffirms = TRUE)

# KS test for 2PSM fit
my_ks_test_RTMA( yi = b3$yi,
                 sei = sqrt(b3$vi),
                 Mhat = Mhat,
                 Shat = Shat )

# ~~ One-tailed p-value histogram from PublicationBias ----

# not helpful compared to PublicationBias version
library(PublicationBias)
pval_plot(yi = b3$yi, vi = b3$vi)


# SPECIFIC TO SAPBE -----------------------

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

