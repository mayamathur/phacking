
# IMPORTANT NOTES -----------------------------

# This script is designed to be run in interactive Sherlock session,
#  NOT via an sbatch file. It analyzes ONE meta-analysis.
#  To analyze multiple metas, use the script 2022-3-11 applied doParallel SAPH.R

# To quickly run this script in high-mem interactive session:
# setwd("/home/groups/manishad/SAPH/applied_examples/code"); source("2_analyze_lodder_sherlock.R")

# To update this code on Sherlock, need to upload it manually since it's not in the dir
#  that pushes with my existing command.


# This is very similar to doParallel.R.
# The only real additions/changes are:
#  - plot_trunc_densities_RTMA


# "show rep res"
# quickly look at results when running locally
srr = function() {
  
  if( "optimx.Mhat.winner" %in% names(rep.res) ) {
    cat("\n")
    print( rep.res %>% select(method, Mhat, MLo, MHi,
                              Shat
                              # optim.converged,
                              # optimx.Mhat.winner,
                              # optimx.Nconvergers,
                              # optimx.Pagree.of.convergers.Mhat.winner
                              ) %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  } else {
    cat("\n")
    print( rep.res %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  }
}



# ~~ User-Specified Parameters -----------------------------------------------

run.local = FALSE

sherlock.data.dir = "/home/groups/manishad/SAPH/applied_examples/data"
sherlock.code.dir = "/home/groups/manishad/SAPH"
sherlock.results.dir = "/home/groups/manishad/SAPH/applied_examples/results/lodder"


local.data.dir = "~/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/Lodder/Prepped data"
local.results.dir = NULL

# specify which methods to run, as in doParallel
# but obviously can't run gold-std on a non-simulated meta-analysis
all.methods = "naive ; maon ; pcurve ; 2psm ; jeffreys-mcmc ; 2psm-csm-dataset ; csm-mcmc ; prereg-naive"
#all.methods = "naive ; maon ; 2psm ; mle-sd; csm-mle-sd"
# parse methods string
all.methods = unlist( strsplit( x = all.methods,
                                split = " ; " ) )

run.optimx = TRUE
stan.adapt_delta = 0.98
stan.maxtreedepth = 20
# hacky because estimate_jeffreys_mcmc_RTMA looks for p as global var
p = data.frame( stan.adapt_delta = 0.98,
                stan.maxtreedepth = 20 )


# for labelling the results file
meta.name = "lodder"


# ~~ Load Data and Packages -----------------------------------------------


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
           "here")



# load dataset and code
if ( run.local == FALSE ) {
  # set up dataset (currently specific to Hagger)
  setwd(sherlock.code.dir)
  source("helper_SAPH.R")
  
  setwd(sherlock.data.dir)
  # "dp" because this is published studies only
  dp = read.csv("lodder_prepped.csv")

  
  
  # load packages with informative messages if one can't be installed
  # **Common reason to get the "could not library" error: You did ml load R/XXX using an old version
  any.failed = FALSE
  for (pkg in toLoad) {
    
    cat( paste("\nAbout to try loading package", pkg) )
    
    tryCatch({
      # eval below needed because library() will otherwise be confused
      # https://www.mitchelloharawild.com/blog/loading-r-packages-in-a-loop/
      eval( bquote( library( .(pkg) ) ) )
    }, error = function(err) {
      cat( paste("\n*** COULD NOT LIBRARYIZE PACKAGE:", pkg) )
      any.failed <<- TRUE
    })
    
  }
  if ( any.failed == TRUE ) stop("Some packages couldn't be installed. See outfile for details of which ones.")
  
  # helper code
  setwd(sherlock.code.dir)
  source("helper_SAPH.R")
}


head(dp)
# for the full dataset that includes all interaction rows
expect_equal( nrow(dp), 287 )


# MAKE DATA SUBSETS ------------------------------

# published nonaffirmatives only
dpn = dp[ dp$affirm == FALSE, ]

# special dataset for CSM: 
# throws away affirmatives from hacked studies (i.e., all non-prereg studies)
dp.csm = dp %>% filter( Preregistered == TRUE | affirm == FALSE )

dp.csm %>% group_by(Preregistered, affirm) %>%
  summarise(n())

# unhacked only
dp.prereg = dp %>% filter(Preregistered == TRUE)
dpn.prereg = dpn %>% filter(Preregistered == TRUE)


# COMPILE STAN MODEL ONCE AT BEGINNING ------------------------------

if ( "jeffreys-mcmc" %in% all.methods ) {
  setwd(sherlock.code.dir)
  source("init_stan_model_applied_SAPH.R")
}



# RUN ANALYSIS ------------------------------

# initialize rep.res st run_method_safe and other standalone estimation fns
#  will correctly recognize it as having 0 rows
rep.res = data.frame()

# ~~ Start Values ------------------------------
#@UNLIKE in doParallel, here we start at (0,1)
#  by default, but if running method-naive, those will be the start values instead
Mhat.start = 0
Shat.start = 1

# in case we're not doing jeffreys-mcmc or it fails
Mhat.MaxLP = NA
Shat.MaxLP = NA

Mhat.MAP = NA
Shat.MAP = NA

# ~~ Fit Naive Meta-Analysis (All PUBLISHED Draws) ------------------------------

if ( "naive" %in% all.methods ) {
  rep.res = run_method_safe(method.label = c("naive"),
                            method.fn = function() {
                              mod = rma( yi = dp$yi,
                                         vi = dp$vi,
                                         method = "REML",
                                         knha = TRUE )
                              
                              report_meta(mod, .mod.type = "rma")
                            },
                            .rep.res = rep.res )
  
  Mhat.naive = rep.res$Mhat[ rep.res$method == "naive" ]
  Shat.naive = rep.res$Shat[ rep.res$method == "naive" ]
}

srr()

# ~~ Change Starting Values -----
if ( !is.na(Mhat.naive) ) Mhat.start = Mhat.naive 
if ( !is.na(Shat.naive) ) Shat.start = Shat.naive 



# ~~ Fit MAON (Nonaffirmative Published Draws) ------------------------------

if ( "maon" %in% all.methods ) {
  
  rep.res = run_method_safe(method.label = c("maon"),
                            method.fn = function() {
                              mod = robu( yi ~ 1, 
                                          data = dpn, 
                                          studynum = 1:nrow(dpn),
                                          var.eff.size = vi,
                                          small = TRUE )
                              
                              report_meta(mod, .mod.type = "robu")
                            },
                            .rep.res = rep.res )
  
}


srr()

# ~~ Fit 2PSM (All Published Draws) ------------------------------

if ( "2psm" %in% all.methods ) {
  
  rep.res = run_method_safe(method.label = c("2psm"),
                            method.fn = function() {
                              mod = weightfunct( effect = dp$yi,
                                                 v = dp$vi,
                                                 steps = c(0.025, 1),
                                                 table = TRUE ) 
                              
                              H = mod[[2]]$hessian
                              ses = sqrt( diag( solve(H) ) )
                              
                              # follow the same return structure as report_meta
                              list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                        MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                        MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                        
                                                        # could definitely get these from weightr
                                                        # I didn't even try
                                                        Shat = NA,
                                                        SLo = NA,
                                                        SHi = NA ) )
                            },
                            .rep.res = rep.res )
  
}

srr()


# ~~ Fit P-Curve (Published Affirmatives) ------------------------------

if ( "pcurve" %in% all.methods ) {
  # since pcurve.opt internally retains only 
  rep.res = run_method_safe(method.label = c("pcurve"),
                            method.fn = function() {
                              #@later, revisit the decision to use df_obs = 1000
                              #   to effectively treat yi/sei z-scores
                              
                              # pass all significant studies (that's what pcurve.opt does internally):
                              Mhat = pcurve.opt( t_obs = dp$Zi,
                                                 df_obs = rep(1000, length(dp$Zi)),
                                                 dmin = -5, #@HARD-CODED and arbitrary
                                                 dmax = 5)
                              
                              
                              # # using only published affirmatives:
                              # Mhat = pcurve.opt( t_obs = dpa$yi/dpa$sei,
                              #                    df_obs = rep(1000, length(dpa$yi)),
                              #                    dmin = -5, #@HARD-CODED and arbitrary
                              #                    dmax = 5)
                              
                              
                              return( list( stats = data.frame( Mhat = Mhat) ) )
                            },
                            .rep.res = rep.res )
  
}

# ~~ MCMC ------------------------------

if ( "jeffreys-mcmc" %in% all.methods ) {
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
  #  returns 2 lines of output, one for posterior mean and one for posterior median
  # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
  rep.res = run_method_safe(method.label = c("jeffreys-mcmc-pmean",
                                             "jeffreys-mcmc-pmed",
                                             "jeffreys-mcmc-max-lp-iterate"),
                            method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi,
                                                                               .sei = sqrt(dpn$vi),
                                                                               .tcrit = qnorm(0.975),
                                                                               .Mu.start = Mhat.start,
                                                                               .Tt.start = Shat.start,
                                                                               .stan.adapt_delta = stan.adapt_delta,
                                                                               .stan.maxtreedepth = stan.maxtreedepth),
                            .rep.res = rep.res )
  
  
  Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
  Shat.MaxLP = rep.res$Shat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
  
  cat("\n doParallel flag: Done jeffreys-mcmc if applicable")
}

cat("\n")
rep.res
cat("\n")

# ~~ Change Starting Values -----
if ( !is.na(Mhat.MaxLP) ) Mhat.start = Mhat.MaxLP 
if ( !is.na(Shat.MaxLP) ) Shat.start = Shat.MaxLP 

# ~~ MAP (SD param) ------------------------------

if ( "jeffreys-sd" %in% all.methods ) {
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  rep.res = run_method_safe(method.label = c("jeffreys-sd"),
                            method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                          sei = sqrt(dpn$vi),
                                                                          par2is = "Tt",
                                                                          tcrit = qnorm(0.975), 
                                                                          Mu.start = Mhat.start,
                                                                          par2.start = Shat.start,
                                                                          usePrior = TRUE,
                                                                          get.CIs = TRUE,
                                                                          CI.method = "wald",
                                                                          
                                                                          run.optimx = run.optimx),
                            .rep.res = rep.res )
  
  Mhat.MAP = rep.res$Mhat[ rep.res$method == "jeffreys-sd" ]
  Shat.MAP = rep.res$Shat[ rep.res$method == "jeffreys-sd" ]
}

cat("\n")
rep.res
cat("\n")


# ~~ MAP (Var param) ------------------------------

if ( "jeffreys-var" %in% all.methods ) {
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  rep.res = run_method_safe(method.label = c("jeffreys-var"),
                            method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                          sei = sqrt(dpn$vi),
                                                                          par2is = "T2t",
                                                                          tcrit = qnorm(0.975), 
                                                                          Mu.start = Mhat.start,
                                                                          par2.start = Shat.start^2,
                                                                          usePrior = TRUE,
                                                                          get.CIs = TRUE,
                                                                          CI.method = "wald",
                                                                          
                                                                          run.optimx = run.optimx),
                            .rep.res = rep.res )
  
  Mhat.MAP = rep.res$Mhat[ rep.res$method == "jeffreys-var" ]
  Shat.MAP = rep.res$Shat[ rep.res$method == "jeffreys-var" ]
}

cat("\n")
rep.res
cat("\n")

# ~~ Change Starting Values -----

if ( !is.na(Mhat.MAP) ) Mhat.start = Mhat.MAP 
if ( !is.na(Shat.MAP) ) Shat.start = Shat.MAP 

# ~~ MLE (SD param) ------------------------------

if ( "mle-sd" %in% all.methods ) {
  
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  rep.res = run_method_safe(method.label = c("mle-sd"),
                            method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                          sei = sqrt(dpn$vi),
                                                                          par2is = "Tt",
                                                                          tcrit = qnorm(0.975), 
                                                                          Mu.start = Mhat.start,
                                                                          par2.start = Shat.start,
                                                                          usePrior = FALSE,
                                                                          get.CIs = TRUE,
                                                                          CI.method = "wald",
                                                                          run.optimx = run.optimx),
                            .rep.res = rep.res )
  
  
  
  cat("\n doParallel flag: Done mle-sd if applicable")
}

srr()


# ~~ Conditional Selection Model: MLE *with* prereg affirms (SD param) ------------------------------

if ( "csm-mle-sd" %in% all.methods ) {
  
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  rep.res = run_method_safe(method.label = c("csm-mle-sd"),
                            method.fn = function() estimate_jeffreys_RTMA(yi = dp.csm$yi,
                                                                          sei = sqrt(dp.csm$vi),
                                                                          par2is = "Tt",
                                                                          tcrit = qnorm(0.975), 
                                                                          Mu.start = Mhat.start,
                                                                          par2.start = Shat.start,
                                                                          usePrior = FALSE,
                                                                          get.CIs = TRUE,
                                                                          CI.method = "wald",
                                                                          run.optimx = run.optimx),
                            .rep.res = rep.res )
  
  
  
  cat("\n doParallel flag: Done csm-mle-sd if applicable")
}

srr()

### LEFT-truncated normal
# include only affirmatives
rep.res = run_method_safe(method.label = c("ltn-mle-sd"),
                          method.fn = function() estimate_jeffreys_RTMA(yi = dp$yi[ dp$affirm == TRUE ],
                                                                        sei = sqrt(dp$vi[ dp$affirm == TRUE ]),
                                                                        par2is = "Tt",
                                                                        tcrit = qnorm(0.975), 
                                                                        Mu.start = Mhat.start,
                                                                        par2.start = Shat.start,
                                                                        usePrior = FALSE,
                                                                        get.CIs = TRUE,
                                                                        CI.method = "wald",
                                                                        run.optimx = run.optimx),
                          .rep.res = rep.res )



# ~~ CSM MCMC ------------------------------

if ( "csm-mcmc" %in% all.methods ) {
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
  #  returns 2 lines of output, one for posterior mean and one for posterior median
  # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
  rep.res = run_method_safe(method.label = c("csm-mcmc-pmean",
                                             "csm-mcmc-pmed",
                                             "csm-mcmc-max-lp-iterate"),
                            method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dp.csm$yi,
                                                                               .sei = sqrt(dp.csm$vi),
                                                                               .tcrit = dp.csm$tcrit,
                                                                               
                                                                               .Mu.start = Mhat.start,
                                                                               .Tt.start = max(0.01, Shat.start),
                                                                               .stan.adapt_delta = stan.adapt_delta,
                                                                               .stan.maxtreedepth = stan.maxtreedepth),
                            .rep.res = rep.res )
  
  cat("\n doParallel flag: Done csm-mcmc if applicable")
}


# ~~ 2PSM With CSM Dataset ------------------------------

if ( "2psm-csm-dataset" %in% all.methods &
     nrow(dp.csm) > 0 ) {
  
  rep.res = run_method_safe(method.label = c("2psm-csm-dataset"),
                            method.fn = function() {
                              mod = weightfunct( effect = dp.csm$yi,
                                                 v = dp.csm$vi,
                                                 steps = c(0.025, 1),
                                                 table = TRUE ) 
                              
                              H = mod[[2]]$hessian
                              ses = sqrt( diag( solve(H) ) )
                              
                              # follow the same return structure as report_meta
                              list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                        MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                        MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                        
                                                        Shat = sqrt( mod[[2]]$par[1] ),
                                                        # truncate lower limit at 0
                                                        SLo = sqrt( max( 0, mod[[2]]$par[1] - qnorm(.975) * ses[1] ) ),
                                                        SHi = sqrt( mod[[2]]$par[1] + qnorm(.975) * ses[1] ) ) ) 
                            },
                            .rep.res = rep.res )
  
  cat("\n doParallel flag: Done 2psm-csm-dataset if applicable")
  
}



# ~~ MLE (Var param) ------------------------------

if ( "mle-var" %in% all.methods ) {
  rep.res = run_method_safe(method.label = c("mle-var"),
                            method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                          sei = sqrt(dpn$vi),
                                                                          par2is = "T2t",
                                                                          tcrit = qnorm(0.975), 
                                                                          Mu.start = Mhat.start,
                                                                          par2.start = Shat.start^2,
                                                                          usePrior = FALSE,
                                                                          get.CIs = TRUE,
                                                                          CI.method = "wald",
                                                                          run.optimx = run.optimx),
                            .rep.res = rep.res )
  
  
  
  cat("\n doParallel flag: Done mle-var if applicable")
  
}

srr()



# ANALYSES OF PREREG STUDIES (SPECIFIC TO LODDER) ------------------

# ~~ Fit Naive (Prereg Only) ------------------------------

if ( "prereg-naive" %in% all.methods ) {
  rep.res = run_method_safe(method.label = c("prereg-naive"),
                            method.fn = function() {
                              mod = rma( yi = dp.prereg$yi,
                                         vi = dp.prereg$vi,
                                         method = "REML",
                                         knha = TRUE )
                              
                              report_meta(mod, .mod.type = "rma")
                            },
                            .rep.res = rep.res )
  
  Mhat.naive = rep.res$Mhat[ rep.res$method == "prereg-naive" ]
  Shat.naive = rep.res$Shat[ rep.res$method == "prereg-naive" ]
}

srr()

Mhat.start = Mhat.naive
Shat.start = Shat.naive

# ~~ MCMC (Prereg Only) ------------------------------

if ( "prereg-jeffreys-mcmc" %in% all.methods ) {
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
  #  returns 2 lines of output, one for posterior mean and one for posterior median
  # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
  rep.res = run_method_safe(method.label = c("prereg-jeffreys-mcmc-pmean",
                                             "prereg-jeffreys-mcmc-pmed",
                                             "prereg-jeffreys-mcmc-max-lp-iterate"),
                            method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dpn.prereg$yi,
                                                                               .sei = sqrt(dpn.prereg$vi),
                                                                               .tcrit = qnorm(0.975),
                                                                               .Mu.start = Mhat.start,
                                                                               .Tt.start = Shat.start,
                                                                               .stan.adapt_delta = stan.adapt_delta,
                                                                               .stan.maxtreedepth = stan.maxtreedepth),
                            .rep.res = rep.res )
  
  
  Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
  Shat.MaxLP = rep.res$Shat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
  
  cat("\n doParallel flag: Done jeffreys-mcmc if applicable")
}


# ~~ Fit 2PSM (Prereg Only) ------------------------------

if ( "prereg-2psm" %in% all.methods ) {
  
  rep.res = run_method_safe(method.label = c("prereg-2psm"),
                            method.fn = function() {
                              mod = weightfunct( effect = dp.prereg$yi,
                                                 v = dp.prereg$vi,
                                                 steps = c(0.025, 1),
                                                 table = TRUE ) 
                              
                              H = mod[[2]]$hessian
                              ses = sqrt( diag( solve(H) ) )
                              
                              # follow the same return structure as report_meta
                              list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                        MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                        MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                        
                                                        # could definitely get these from weightr
                                                        # I didn't even try
                                                        Shat = NA,
                                                        SLo = NA,
                                                        SHi = NA ) )
                            },
                            .rep.res = rep.res )
  
}

srr()



# ~~ Conditional Selection Model: all prereg studies (SD param) ------------------------------

if ( "prereg-csm-mle-sd" %in% all.methods ) {
  
  # # temp for refreshing code
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # source("helper_SAPH.R")
  
  rep.res = run_method_safe(method.label = c("prereg-csm-mle-sd"),
                            method.fn = function() estimate_jeffreys_RTMA(yi = dp.prereg$yi,
                                                                          sei = sqrt(dp.prereg$vi),
                                                                          par2is = "Tt",
                                                                          tcrit = qnorm(0.975), 
                                                                          Mu.start = Mhat.start,
                                                                          par2.start = Shat.start,
                                                                          usePrior = FALSE,
                                                                          get.CIs = TRUE,
                                                                          CI.method = "wald",
                                                                          run.optimx = run.optimx),
                            .rep.res = rep.res )
  
  
  
  cat("\n doParallel flag: Done csm-mle-sd if applicable")
}

srr()

# ### LEFT-truncated normal
# # include only affirmatives
# rep.res = run_method_safe(method.label = c("ltn-mle-sd"),
#                           method.fn = function() estimate_jeffreys_RTMA(yi = dp$yi[ dp$affirm == TRUE ],
#                                                                         sei = sqrt(dp$vi[ dp$affirm == TRUE ]),
#                                                                         par2is = "Tt",
#                                                                         tcrit = qnorm(0.975), 
#                                                                         Mu.start = Mhat.start,
#                                                                         par2.start = Shat.start,
#                                                                         usePrior = FALSE,
#                                                                         get.CIs = TRUE,
#                                                                         CI.method = "wald",
#                                                                         run.optimx = run.optimx),
#                           .rep.res = rep.res )



# OTHER NEW ANALYSES --------

# ~~ Fit 2PSM (CSM dataset) -----------------------------

# because this is estimating eta, I expect it to *not* agree with CSM
# because we've discarded some affirmatives, it should underestimate eta and overestimate Mu

if ( "csm-dataset-2psm" %in% all.methods ) {
  
  rep.res = run_method_safe(method.label = c("csm-dataset-2psm"),
                            method.fn = function() {
                              mod = weightfunct( effect = dp.csm$yi,
                                                 v = dp.csm$vi,
                                                 steps = c(0.025, 1),
                                                 table = TRUE ) 
                              
                              H = mod[[2]]$hessian
                              ses = sqrt( diag( solve(H) ) )
                              
                              # follow the same return structure as report_meta
                              list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                        MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                        MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                        
                                                        # could definitely get these from weightr
                                                        # I didn't even try
                                                        Shat = NA,
                                                        SLo = NA,
                                                        SHi = NA ) )
                            },
                            .rep.res = rep.res )
  
}

# sanity check with metafor


# DATASET INFO AND SANITY CHECKS ------------------------------

# Info About Dataset ------------

rep.res$k.pub = nrow(dp)
rep.res$k.pub.affirm = sum(dp$affirm == TRUE)
rep.res$k.pub.nonaffirm = sum(dp$affirm == FALSE)
rep.res$prob.pub.study.affirm = rep.res$k.pub.affirm / rep.res$k.pub

# Show Results --------------

cat("\n\n***REP.RES at the end:")
rep.res %>% select(method, Mhat, MLo, MHi, MhatRhat, optimx.Nagree.of.convergers.Mhat.winner) %>%
  mutate_if(is.numeric, function(x) round(x,2) )


# # QQ Plot: RTMA ------------------
# 
# p = yi_qqplot(yi = dpn$yi,
#               sei = dpn$sei,
#               Mhat = Mhat,
#               Shat = Shat)
# 
# 
# # QQ Plot: CSM -------------------
# 
# # this will be a bad fit because of hacking
# p = yi_qqplot(yi = dp.csm$yi,
#               sei = dp.csm$sei,
#               Mhat = Mhat,
#               Shat = Shat)
# 
# # QQ Plot: All Results -------------------
# 
# # this will be a bad fit because of hacking
# p = yi_qqplot(yi = dp$yi,
#               sei = dp$sei,
#               Mhat = Mhat,
#               Shat = Shat)
# 
# 
# ### Make RTMA Density Plot 
# 
# plot.method = "jeffreys-mcmc-pmed"
# 
# # catch possibility that we didn't run this method
# if ( length( rep.res$Mhat[ rep.res$method == plot.method ] ) > 0 ) {
#   p = plot_trunc_densities_RTMA(d = dp,
#                                 Mhat = rep.res$Mhat[ rep.res$method == plot.method ],
#                                 Shat = rep.res$Shat[ rep.res$method == plot.method ],
#                                 showAffirms  = FALSE)
#   
#   setwd(sherlock.results.dir)
#   ggsave( filename = paste( "density_plot", meta.name, ".pdf", sep="_" ),
#           width = 10, 
#           height = 10)
#   
# }


# ~ WRITE LONG RESULTS ------------------------------

setwd(sherlock.results.dir)
fwrite( rep.res, paste( "results_", meta.name, ".csv", sep="" ) )

