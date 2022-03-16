


# IMPORTANT NOTES -----------------------------

# This script is designed to be run interactively in Sherlock. 
#  However, copy-pasting the foreach loop doesn't work because the Sherlock
#   console gets confused about the indentation, etc. 

#  Instead, run this script as follows in high-mem interactive session:
# setwd("/home/groups/manishad/SAPH"); source("2022-3-11_applied_doParallel_SAPH.R")

# To run via sbatch (not tested):
# sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/2022-3-11_applied.sbatch

# To bring back all results locally:

#SAPB
# scp -r mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/applied_examples/results/sapbe /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/Applied\ examples/2022-3-11\ SAPBE\ results\ from\ Sherlock

#Kvarven
# scp -r mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/applied_examples/results/kvarven /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/Applied\ examples

# This script is very similar to doParallel.R from 2022-3-11.
# The only real additions/changes are:
#  - plot_trunc_densities_RTMA

# # FOR LOCAL TESTING
# setwd("/Users/mmathur/Dropbox/Personal computer/Independent studies/2021/Sensitivity analysis for p-hacking (SAPH)/Code (git)/Applied examples/2022-3-13 prep Kvarven dataset for SAPH/Datasets prepped for SAPH")
# b2 = read.csv("b2_long_prepped_kvarven.csv")




# ~~ User-Specified Parameters -----------------------------------------------

# to run via sbatch file, must be FALSE 
#  affects which version of init_stan_model is used
run.interactive = TRUE

# should we recompile stan model?
# set to FALSE if you've already compiled in a given interactive session
compile.stan.model = TRUE

# "sapbe" or "kvarven"
dataset.name = "kvarven"

if ( dataset.name == "sapbe" ) {
  sherlock.data.dir = "/home/groups/manishad/SAPH/applied_examples/data/sapbe"
  sherlock.code.dir = "/home/groups/manishad/SAPH"
  sherlock.results.dir = "/home/groups/manishad/SAPH/applied_examples/results/sapbe"
}

if ( dataset.name == "kvarven" ) {
  sherlock.data.dir = "/home/groups/manishad/SAPH/applied_examples/data/kvarven"
  sherlock.code.dir = "/home/groups/manishad/SAPH"
  sherlock.results.dir = "/home/groups/manishad/SAPH/applied_examples/results/kvarven"
}

# specify which methods to run, as in doParallel
# but obviously can't run gold-std on a non-simulated meta-analysis
all.methods = "naive ; maon ; 2psm ; pcurve ; jeffreys-mcmc ; jeffreys-sd ; mle-sd"
#all.methods = "naive ; maon ; 2psm"
#@temp:
run.optimx = TRUE
stan.adapt_delta = 0.98
stan.maxtreedepth = 20
# hacky because estimate_jeffreys_mcmc_RTMA looks for p as global var
p = data.frame( stan.adapt_delta = 0.98,
                stan.maxtreedepth = 20 )


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
           "weightr",
           "data.table")


setwd(sherlock.code.dir)
source("helper_SAPH.R")



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


# read in dataset
setwd(sherlock.data.dir)

if ( dataset.name == "sapbe" ){
  b2 = fread("b2_long_prepped.csv")
  f2 = fread("f2_short_prepped.csv")
  
  # make columns with standardized names to match doParallel from sim study
  b2$yi = b2$EstF  # **uses the direction-flipped estimates
  b2$vi = b2$SE^2
  b2$Zi = b2$yi / sqrt(b2$vi)
  #**note we're using all estimates, not just randomly-chosen ones, because we're focusing on metas with little clustering
}

if ( dataset.name == "kvarven" ){
  b2 = fread("b2_long_prepped_kvarven.csv")
  
  # make columns with standardized names to match doParallel from sim study
  # per sanity check in prep code, all metas' pooled ests have already 
  #   been coded as >0
  b2$yi = b2$d  
  b2$vi = b2$var
  b2$Zi = b2$yi / sqrt(b2$vi)
  
  # affirm status
  b2$pval.two = 2 * ( 1 - pnorm( abs(b2$Zi) ) ) 
  b2$affirm = (b2$pval.two <= 0.05) & (b2$yi > 0)
  expect_equal( b2$affirm, b2$Zi > qnorm(0.975) ) 
  
  # match SAPB-E naming convention
  b2 = b2 %>% rename( meta.name = meta)

}





# COMPILE STAN MODEL ONCE AT BEGINNING ------------------------------

if ( "jeffreys-mcmc" %in% all.methods &
     compile.stan.model == TRUE ) {
  setwd(sherlock.code.dir)

  # this version of init_stan doesn't use parallelization
  if ( run.interactive == TRUE ) source("init_stan_model_applied_SAPH.R")

  # this version of init_stan DOES use parallelization
  if ( run.interactive == FALSE ) {
    registerDoParallel(cores=16)
    source("init_stan_model_SAPH.R")
  }
}




# RUN ANALYSIS FOR MULTIPLE METAS ------------------------------



# ~~ Beginning of ForEach Loop -----------------------------

# to analyze all of them
( meta.names.to.analyze = unique(b2$meta.name) )

#@analyze only a subset of them
#meta.names.to.analyze = unique(b2$meta.name)[1:2]

if ( exists("rs") ) rm(rs)

# system.time is in seconds
doParallel.seconds = system.time({
  rs = foreach( i = meta.names.to.analyze, .combine = bind_rows ) %dopar% {
    # for debugging (so that outfile will contain all printed things):
    #for ( i in meta.names.to.analyze ) {
    
    cat("\n\n~~~~~~~~~~~~~~~~ BEGIN ", i, "~~~~~~~~~~~~~~~~")
    
    # results for just this simulation rep
    if ( exists("rep.res") ) suppressWarnings( rm(rep.res) )
    
   
    # ~~ Get Subset Data for This Meta -----------------
    dp = b2 %>% filter(meta.name == i)
    if ( dataset.name == "sapbe" ) f2.subset = f2 %>% filter(meta.name == i)

    # published nonaffirmatives only
    dpn = dp[ dp$affirm == FALSE, ]
    
    # published affirmatives only
    dpa = dp[ dp$affirm == TRUE, ]
    
   cat("\n\n ****** PARTIAL HEAD OF DP:\n")
   print( head(dp %>% select(meta.name, yi)) )

    
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
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
    
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
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
                                                            Shat = sqrt( mod[[2]]$par[1] ),
                                                            # truncate lower limit at 0
                                                            SLo = sqrt( max( 0, mod[[2]]$par[1] - qnorm(.975) * ses[1] ) ),
                                                            SHi = sqrt( mod[[2]]$par[1] + qnorm(.975) * ses[1] ) ) )
                                },
                                .rep.res = rep.res )
      
      
      # #bm: KS test
      # # I think you should incorporate this into the 2psm fn itself
      # Mhat.2PSM = rep.res$Mhat[ rep.res$method == "2psm"]
      # Shat.2PSM = rep.res$Shat[ rep.res$method == "2psm"]
      # 
      # if ( !is.na(Mhat.2PSM) & !is.na(Shat.2PSM) ) {
      #   my_ks_test_RTMA( yi = dp$yi,
      #                    sei = sqrt(dp$vi),
      #                    Mhat = rep.res$Mhat[ rep.res$method == "2psm"],
      #                    Shat = rep.res$Shat[ rep.res$method == "2psm" ] )
      #   
      #   # nonaffirms only
      #   my_ks_test_RTMA( yi = dp$yi[ dp$affirm == FALSE ],
      #                    sei = sqrt(dp$vi[ dp$affirm == FALSE ]),
      #                    Mhat = rep.res$Mhat[ rep.res$method == "2psm"],
      #                    Shat = rep.res$Shat[ rep.res$method == "2psm" ] )
      #   
      #   # affirms only
      #   my_ks_test_RTMA( yi = dp$yi[ dp$affirm == TRUE ],
      #                    sei = sqrt(dp$vi[ dp$affirm == TRUE ]),
      #                    Mhat = rep.res$Mhat[ rep.res$method == "2psm"],
      #                    Shat = rep.res$Shat[ rep.res$method == "2psm" ] )
      # }
      
    }
    
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
    # ~~ Fit P-Curve (Published Affirmatives) ------------------------------
    
    if ( "pcurve" %in% all.methods ) {
      # since pcurve.opt internally retains only 
      rep.res = run_method_safe(method.label = c("pcurve"),
                                method.fn = function() {
                                  #@later, revisit the decision to use df_obs = 1000
                                  #   to effectively treat yi/sei z-scores
                                  
                                  # using all significant studies (that's what pcurve.opt does internally):
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
    
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
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
    
    # cat("\n")
    # print(rep.res)
    # cat("\n")
    
    # ~~ Info About Dataset and Sanity Checks ------------------------------
    
    ### Info About Dataset
    rep.res = rep.res %>% add_column(.before = 1, meta.name = i)
    rep.res$k.pub = nrow(dp)
    rep.res$k.pub.affirm = sum(dp$affirm == TRUE)
    rep.res$k.pub.nonaffirm = sum(dp$affirm == FALSE)
    rep.res$prob.pub.study.affirm = rep.res$k.pub.affirm / rep.res$k.pub
    
    ### Info Specific to Dataset
    
    if ( dataset.name == "sapbe" ) {
      f2.keepers = c("author.year", "meta.study.type", "LogEta", 
                     "VarLogEta", "shapiro.pval", "Mhat", "Mhat.Lo", "Mhat.Hi", "Mhat.Pval", 
                     "Analysis.Scale", "Sval.Est", "Sval.CI", "Sval.Est.Num", "Sval.CI.Num", 
                     "Mhat.Worst", "Mhat.Worst.Lo", "Mhat.Worst.Hi", "Mhat.Worst.Pval", 
                     "group", "group.pretty", "discipline", "Data.source", "k.all", 
                     "k.rc", "k.affirm.rc", "k.nonaffirm.rc",
                     "meta.med.year", "k.year.coded", "Pdisaffirm.ratio", "PZeroEst", 
                     "excl.Pdisaffirm.ratio", "excl.PZeroEst", "excl.shapiro", "excl.other")
      f2.subset = f2.subset %>% select( all_of(f2.keepers) )
      f2.subset = f2.subset %>% add_column(.after = "LogEta",
                                           Eta = exp(f2.subset$LogEta) )
      
      # for clarity, name all columns from SAPBE
      # since, for example, those use the randomly chosen estimates
      names(f2.subset) = paste("SAPBE.", names(f2.subset), sep = "")
      
      # add f2 info to rep.res
      rep.res = cbind( rep.res, f2.subset )
    }
    

    
    ### Show Results
    
    cat( "\n\n***** META", i, "PARTIAL REP.RES at the end:" )
    print( rep.res %>% select(method, Mhat, MLo, MHi, Shat) %>%
      mutate_if(is.numeric, function(x) round(x,2) ) )
    
    ### Make Plot 
    plot.method = "jeffreys-mcmc-pmed"
    
    # catch possibility that we didn't run this method
    if ( length( rep.res$Mhat[ rep.res$method == plot.method ] ) > 0 ) {
      
      cat("\n About to make plot")
      
      p = plot_trunc_densities_RTMA(d = dp,
                                    Mhat = rep.res$Mhat[ rep.res$method == plot.method ],
                                    Shat = rep.res$Shat[ rep.res$method == plot.method ],
                                    showAffirms = FALSE)
      
      setwd(sherlock.results.dir)
      ggsave( filename = paste( "density_plot", i, ".pdf", sep="_" ),
              width = 10, 
              height = 10)
      
      cat("\n Done writing plot")
      
      
    }
    
    rep.res
    
  }  ### end foreach loop
  
} )[3]  # end system.time


# quick look
rs %>% select(meta.name, method, Mhat, MLo, MHi) %>%
  mutate_if(is.numeric, function(x) round(x, 2))




# ~ WRITE LONG RESULTS ------------------------------

setwd(sherlock.results.dir)
if ( dataset.name == "sapbe" ) fwrite( rs, paste( "results_sapbe_all.csv", sep="" ) )
if ( dataset.name == "kvarven" ) fwrite( rs, paste( "results_sapbe_all.csv", sep="" ) )


