
# IMPORTANT NOTES -----------------------------

# MUST USE ml load R/4.0.2!!!!

# The returned Vhat is an estimate of T2 + t2w, *not* T2 itself
#
# correct_meta_phack1 will NOT work well for small m. I already tried this.
#   This is because, for computational convenience, I'm using the truncated normal dist
#   for the tstats, and it won't hold well when m is small. 
# 
# Key hypotheses to check in simulation results:
#  - T2 = 0, t2w > 0, Nmax > 1: Nonaffirms should fit the trunc distribution because
#    there's no selection on mui. So this scenario should look just like T2 = 0, t2w > 0, Nmax = 1.

#@To do before running:
# - Make sure sbatches uses R/4.0.2


# 2021-5-8: scens that are unbiased, correct coverage, etc., locally:
# scen  Mu   T2   m  t2w  se Nmax   hack rho   k k.hacked
# 1    1 0.1 0.25 500 0.25 0.5   10 affirm 0.9 100        0

# because Sherlock 2.0 restores previous workspace
rm( list = ls() )


# are we running locally?
run.local = FALSE


allPackages = c("here",
                "magrittr",
                "dplyr",
                "data.table",
                "tidyverse",
                "tidyr",
                "metafor",
                "robumeta",
                "testthat",
                "truncdist",
                "gmm",
                "tmvtnorm",
                "doParallel",
                "foreach")




# FOR CLUSTER USE ------------------------------




if (run.local == FALSE) {
  
  # load command line arguments
  args = commandArgs(trailingOnly = TRUE); print(args)
  jobname = args[1]
  scen = args[2]  # this will be a number
  
  # # install any missing packages
  # # find the ones that aren't already installed
  # libDir = "/home/users/mmathur/Rpackages/"
  # ( packagesNeeded = allPackages[ !( allPackages %in% installed.packages(lib.loc = libDir)[,"Package"] ) ] )
  # if( length(packagesNeeded) > 0 ) install.packages(packagesNeeded, lib = libDir)
  # 
  # # load all packages
  # lapply( allPackages,
  #         require,
  #         character.only = TRUE,
  #         lib.loc = libDir)
  
  # version without libDir
 #  cd $HOME
 #  rm -rf R/
 #    rm -rf .Rprofile
 #  rm -rf R_libs
 #  ml R/4.0.2
 # R
  
  ( packagesNeeded = allPackages[ !( allPackages %in% installed.packages()[,"Package"] ) ] )
  if( length(packagesNeeded) > 0 ) install.packages(packagesNeeded)
  
  # load all packages
  lapply( allPackages,
          require,
          character.only = TRUE)
  
  #**you need to see all "TRUE" printed by this in order for the package to actually be loaded
  
  # get scen parameters (made by genSbatch.R)
  path = "/home/groups/manishad/SAPH"
  setwd(path)
  scen.params = read.csv( "scen_params.csv" )
  p <<- scen.params[ scen.params$scen == scen, ]
  print(p)
  
  # # alternatively, generate a simple scen.params in order to run doParallel manually in
  # #  Sherlock as a test
  # scen.params = data.frame( scen = 1,
  #                           Mu = 0.1,
  #                           T2 = 0.25,
  #                           m = 500,
  #                           t2w = .25,
  #                           se = 0.5,
  # 
  #                           Nmax = 10,
  #                           hack = "affirm2", # **
  #                           rho = 0.9,
  # 
  #                           k = 100,
  #                           k.hacked = 100 )  # all published nonaffirms are from hacked studies
  # scen = 1

  # helper code
  setwd(path)
  source("helper_SAPH.R")
  
  # locally, with total k = 100, Nmax = 10, and sim.reps = 250, took 93 min total
  
  # simulation reps to run within this job
  # **this need to match n.reps.in.doParallel in the genSbatch script
  sim.reps = 5  #@update this 
  

  # set the number of cores
  registerDoParallel(cores=16)
  
}



# FOR LOCAL USE  ------------------------------
if ( run.local == TRUE ) {
  #rm(list=ls())
  
  
  lapply( allPackages,
          require,
          character.only = TRUE)
  
  # helper fns
  code.dir = here()
  setwd(code.dir)
  source("helper_SAPH.R")
  
  scen.params = data.frame( scen = 1,
                            Mu = 0.1,
                            T2 = 0.25,
                            m = 500,
                            t2w = .25,
                            se = 0.5,
                            
                            Nmax = 10,
                            hack = "affirm2", # **
                            rho = 0.9,
                            
                            k = 10,
                            k.hacked = 10 )  # all published nonaffirms are from hacked studies
  
  
  
  
  sim.reps = 2  # reps to run in this iterate
  
  # set the number of local cores
  registerDoParallel(cores=8)
  
  scen = 1
  # data.frame(scen.params %>% filter(scen.name == scen))
}



# RUN SIMULATION ------------------------------


#for ( scen in scen.params$scen.name ) {  # can't use this part on the cluster

# system.time is in seconds
doParallelTime = system.time({
  #@change this back before running for real
  #rs = foreach( i = 1:sim.reps, .combine=rbind ) %dopar% {
    # for debugging (out file will contain all printed things):
    for ( i in 1:sim.reps ) {
      
    cat("\n\n~~~~~~~~~~~~~~~~ BEGIN SIM REP", i, "~~~~~~~~~~~~~~~~")
    
    # results for just this simulation rep
    if ( exists("repRes") ) rm(repRes)
    
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    p = scen.params[ scen.params$scen == scen, names(scen.params) != "scen"]
    
    cat("\n\nHEAD OF SCEN.PARAMS:")
    print(p)
    
    #bm
    # ~~ Simulate Dataset ------------------------------
    # includes unpublished studies
    d = sim_meta( Nmax = p$Nmax,
                 Mu = p$Mu,
                 T2 = p$T2,
                 m = p$m,
                 t2w = p$t2w,
                 se = p$se,
                 hack = p$hack,
                 rho = p$rho,
                 
                 k = p$k,
                 k.hacked = p$k.hacked,
                 return.only.published = FALSE )
    
    
    # dataset of only published results
    dp = d %>% filter(Di == 1)
    
    cat("\n\nHEAD OF DP:")
    print(head(dp))
    
    # # published hacked ones only
    # dph = d %>% filter(hack == p$hack & Di == 1)
    # 
    # # published unhacked ones only
    # # so only one per study set
    # # same as second row of above table
    # dpu = d %>% filter(hack == "no" & Di == 1)
    
    
    
    # ~~ Fit Gold-Standard Meta-Analysis to All Results  ------------------------------
    # unbiased meta-analysis of all studies, even unpublished ones
    # account for clustering of draws within studies
    # *the tau^2 estimate will be close to T2
    
    # this takes forever and gets hung up if k =10^3
    # even if I omit the random intercept
    
    # only for non-huge k to prevent hangups
    if ( p$k < 500 ) {
      modAll = rma.mv( yi = yi,
                       V = vi,
                       data = d,
                       method = "REML",
                       random = ~1 | study ) 
    }
    
    
    # ~~ Fit Naive Meta-Analysis on Published Studies  ------------------------------
    # biased meta-analysis of only published studies
    
    # only for non-huge k to prevent hangups
    if ( nrow(dp) < 500 ) {
      modPub = rma( yi = dp$yi,
                    vi = dp$vi,
                    method = "REML",
                    knha = TRUE ) 
    }
    
    
    cat("\n\nSURVIVED MODPUB STEP")
    
    # ~~ Bias-Corrected Estimator #1: Nonaffirms Only ------------------------------
    
    #bm: next try rm(list=ls()) and see if below fn still works :)
    # also think about any diagnostics to add
    #  like the tstat expectation and variance among affirms from hacked vs. 
    #  unhacked studies
    #  then should be able to try running doParallel locally!
    
    modCorr = correct_meta_phack1( .p = p,
                                   .dp = dp )
    
    cat("\n\nSURVIVED MODCORR STEP")
    #print(head(dp))
    
    # add to results
    repRes = add_method_result_row(repRes = NA,
                                   corrObject = modCorr,
                                   methName = "AllNonaffirms")
    
    # # SAVE 
    # # methods from earlier simulations where I was bias-correcting the affirmatives
    # # # ~~ Bias-Corrected Estimator #1 ------------------------------
    # # omniscient version: we know which studies are unhacked
    # #  (also includes some affirmatives)
    # 
    # corr1 = correct_dataset_phack(.dp = dp,
    #                               .p = p,
    #                               hackAssumption = "omniscient")
    # 
    # expect_equal( corr1$sanityChecks$kAssumedHacked,
    #               sum(dp$hack == p$hack) )
    # 
    # expect_equal( corr1$modUH$k,
    #               sum(dp$hack == "no") )
    # 
    # # add to results
    # repRes = add_method_result_row(repRes = NA,
    #                                corrObject = corr1,
    #                                methName = "omniscient")
    # 
    # # ~~ Bias-Corrected Estimator #2 ------------------------------
    # # omniscient version: we know which studies are unhacked
    # #  (also includes some affirmatives)
    # 
    # corr2 = correct_dataset_phack(.dp = dp,
    #                               .p = p,
    #                               hackAssumption = "allAffirms")
    # 
    # # this check assumes p$hack = "affirm" instead of "signif"
    # expect_equal( corr2$sanityChecks$kAssumedHacked,
    #               sum(dp$affirm) )
    # 
    # expect_equal( corr2$modUH$k,
    #               sum(dp$affirm == FALSE) )
    # 
    # 
    # # add to results
    # repRes = add_method_result_row(repRes = repRes,
    #                                corrObject = corr2,
    #                                methName = "allAffirms")
    
    
    
    # ~ Write Results ------------------------------
    
    # add stats for the simulated dataset that don't change based on correction method
    # these are prefaced by dataset name for clarity
    # **note: all columns before methName don't depend on the correction method used
    repRes = repRes %>% add_column(repName = i,
                                   
                                   # dp.k = nrow(dp),
                                   # dp.kAffirm = sum(dp$affirm == TRUE),
                                   # dp.kNonaffirm = sum(dp$affirm == FALSE),
                                   # dp.Nrealized = mean(dp$N),
                                   
                                   report_rma(modAll,
                                              .Mu = p$Mu,
                                              .suffix = "All"),
                                   report_rma(modPub,
                                              .Mu = p$Mu,
                                              .suffix = "Naive"),
                                   
                                   .before = 1 )
    
    
    # add in scenario parameters
    repRes = repRes %>% add_column( scenName = scen,
                                    p,
                                    .after = 1 )
    
    
    repRes
    
    cat("\nSURVIVED MAKING REPRES:")
    print(repRes)
    
  }  ### end foreach loop
  
} )[3]  # end system.time



#@commented out temporarily

# ### Add meta-variables to dataset ###
# 
# # estimated time for 1 simulation rep
# # use NAs for additional methods so that the SUM of the rep times will be the
# #  total computational time
# nMethods = length(unique(rs$methName))
# rs$repSeconds = rep( c( doParallelTime / sim.reps,
#                         rep( NA, nMethods - 1 ) ), sim.reps )
# 
# expect_equal( as.numeric( sum(rs$repSeconds, na.rm = TRUE) ),
#               as.numeric(doParallelTime) )
# 
# 
# 
# ### LOCAL ONLY: Quick look at results ###
# 
# if ( run.local == TRUE ) {
#   takeMean = names(rs)[ !names(rs) %in% c( "repName",
#                                            "scenName",
#                                            "methName",
#                                            names(scen.params) ) ]
#   
#   
#   resTable = rs %>% group_by(methName) %>%
#     #mutate(simReps = n()) %>%
#     summarise_at( takeMean,
#                   function(x) round( mean(x, na.rm = TRUE), 2 ) )
#   View(resTable)
#   
#   # should be similar:
#   # *does NOT match with method = "affirm2"
#   resTable$TheoryExpTstat; resTable$MeanTstatUnhacked
#   # the other ones:
#   resTable$MeanTstat; resTable$MeanTstatHacked; resTable$MeanTstatUnhacked
#   
#   resTable$TheoryVarTstat; resTable$EstVarTstatUnhacked
#   
#   
#   # bias:
#   scen.params$Mu
#   resTable$MhatAll
#   resTable$MhatNaive
#   resTable$MhatCorr
#   
#   
#   
#   # CI coverage:
#   resTable$MhatCoverAll  # should definitely be good
#   resTable$MhatCoverNaive
#   resTable$MhatCoverCorr
#   
#   # setwd("Results")
#   # fwrite(resTable, "all_hacked_affirm2_Nmax10.csv")
#   # fwrite(resTable, "resTable.csv")
# }
# 
# 
# 
# # ~ WRITE LONG RESULTS ------------------------------
# if ( run.local == FALSE ) {
#   setwd("/home/groups/manishad/SAPH/long_results")
#   fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
# }