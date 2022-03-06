


# IMPORTANT NOTES -----------------------------

# This is very similar to doParallel.R from 2022-3-6.
# The only real additions/changes are:
#  - plot_trunc_densities_RTMA

# for interactive Sherlock:
# ml load v8
# ml load R/4.1.2
# srun --mem=32G --time=2:00:00 --pty bash
# R

# ~~ User-Specified Parameters -----------------------------------------------

run.local = FALSE

sherlock.data.dir = "/home/groups/manishad/SAPH/applied_examples/data"
sherlock.code.dir = "/home/groups/manishad/SAPH"
sherlock.results.dir = "/home/groups/manishad/SAPH/applied_examples/results"

#@SPECIFIC TO HAGGER
local.data.dir = "~/Dropbox/Personal computer/Miscellaneous/Journal article library/Reproducibility:replicability/Kvarven comparing meta-analysis and multisite replications/Reanalyze Kvarven/Prepped Hagger data (meta and replication)"
local.results.dir = NULL

# specify which methods to run, as in doParallel
# but obviously can't run gold-std
# all.methods = "naive ; maon ; 2psm ; jeffreys-mcmc ; jeffreys-sd ; mle-sd ; mle-var"
all.methods = "naive ; maon ; 2psm ; jeffreys-sd ; mle-sd ; mle-var"
run.optimx = TRUE
stan.adapt_delta = 0.98
stan.maxtreedepth = 20


# for labelling the results file
meta.name = "hagger"


# ~~ Load Data and Packages -----------------------------------------------

# parse methods string
all.methods = unlist( strsplit( x = all.methods,
                                split = " ; " ) )
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
           "weightr")



#@SPECIFIC TO HAGGER
if ( run.local == FALSE ) {
  # set up dataset (currently specific to Hagger)
  setwd(sherlock.code.dir)
  source("helper_SAPH.R")
  
  setwd(sherlock.data.dir)
  # "dp" because this is published studies only
  dp = read.csv("prepped_hagger_meta_data.csv")
  dp$yi = dp$d
  dp$vi = dp$var
  
  
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

if ( run.local == TRUE ) {
  lapply( toLoad,
          require,
          character.only = TRUE)
  
  
  # helper fns
  code.dir = here("Sherlock code")
  setwd(code.dir)
  source("helper_SAPH.R")
  
  # "dp" because this is published studies only
  setwd(local.data.dir)
  dp = read.csv("prepped_hagger_meta_data.csv")
}


head(dp)


# PREP DATASET ------------------------------

# exclude missing data because estimate_jeffreys_RTMA, etc., aren't designed to handle that
dp = dp %>% filter( !is.na(yi) & !is.na(vi) )


# published nonaffirmatives only
dpn = dp[ dp$affirm == FALSE, ]



# COMPILE STAN MODEL ONCE AT BEGINNING ------------------------------

if ( "jeffreys-mcmc" %in% all.methods ) {
  setwd(sherlock.code.dir)
  source("init_stan_model_applied_SAPH.R")
}

#@BM: this step breaks


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

cat("\n")
rep.res
cat("\n")

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


cat("\n")
rep.res
cat("\n")

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


cat("\n")
rep.res
cat("\n")

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

cat("\n")
rep.res
cat("\n")

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

cat("\n")
rep.res
cat("\n")

# ~~ Info About Dataset and Sanity Checks ------------------------------

### Info About Dataset

rep.res$k.pub = nrow(dp)
rep.res$k.pub.affirm = sum(dp$affirm == TRUE)
rep.res$k.pub.nonaffirm = sum(dp$affirm == FALSE)
rep.res$prob.pub.study.affirm = rep.res$k.pub.affirm / rep.res$k.pub



### Make Plot 

plot.method = "jeffreys-mcmc-pmed"

# catch possibility that we didn't run this method
if ( length( rep.res$Mhat[ rep.res$method == plot.method ] ) > 0 ) {
  p = plot_trunc_densities_RTMA(d = dp,
                                Mhat = rep.res$Mhat[ rep.res$method == plot.method ],
                                Shat = rep.res$Shat[ rep.res$method == plot.method ],
                                showAffirms  = FALSE)
  
  setwd(results.dir)
  ggsave( filename = paste( "density_plot", meta.name, ".pdf", sep="_" ),
          width = 10, 
          height = 10)
  
}





# ~ WRITE LONG RESULTS ------------------------------

setwd(results.dir)
fwrite( rep.res, paste( "results", meta.name, ".csv", sep="_" ) )











