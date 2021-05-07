

# **note that the returned Vhat is an estimate of T2 + t2w, not T2 itself



# because Sherlock 2.0 restores previous workspace
rm( list = ls() )

# are we running locally?
run.local = FALSE



# # FOR CLUSTER USE ------------------------------
# if (run.local == FALSE) {
#   
#   # load command line arguments
#   args = commandArgs(trailingOnly = TRUE)
#   jobname = args[1]
#   scen = args[2]  # this will be a letter
#   
#   # get scen parameters
#   setwd("/home/groups/manishad/MRM")
#   scen.params = read.csv( "scen_params.csv" )
#   p = scen.params[ scen.params$scen.name == scen, ]
#   
#   print(p)
#   
#   
#   # simulation reps to run within this job
#   # this need to match n.reps.in.doParallel in the genSbatch script
#   sim.reps = 500  # for main sims
#   # used boot.reps=5,000 in NPPhat but have reduced to 1,000
#   # MR bt mn both correct still times out with 5:00:00 at 1000 boot.reps
#   # JUST TO LOOK AT TIMEOUT ISSUE:
#   # for largest scenario (k=150) with method MR and boot.reps=50 and sim.reps=100, one sbatch took 15 min
#   #  so boot.reps=1,000 should be about 5 hrs
#   # reducing sim.reps to 50 should be about 2.5 hrs
#   boot.reps = 1000
#   
#   
#   # EDITED FOR C++ ISSUE WITH PACKAGE INSTALLATION
#   library(crayon, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(dplyr, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(foreach, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(doParallel, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(boot, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(metafor, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(robumeta, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(data.table, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(purrr, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(metRology, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(fansi, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(MetaUtility, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(ICC, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(cfdecomp, lib.loc = "/home/groups/manishad/Rpackages/")
#   library(tidyr, lib.loc = "/home/groups/manishad/Rpackages/")
#   
#   # for use in ml load R
#   # install.packages( c("ICC", "cfdecomp", "tidyr"), lib = "/home/groups/manishad/Rpackages/" )
#   
#   path = "/home/groups/manishad/MRM"
#   setwd(path)
#   source("bootfuns.R")
#   source("helper_MRM.R")
#   
#   # set the number of cores
#   registerDoParallel(cores=16)
# 
# }



# FOR LOCAL USE  ------------------------------
if ( run.local == TRUE ) {
  #rm(list=ls())
  
  library(here)
  # data-wrangling packages
  library(dplyr)
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
  # for this project
  library(truncdist)
  #library(ExtDist)
  library(gmm)  # https://stackoverflow.com/questions/63511986/error-package-or-namespace-load-failed-for-gmm-in-dyn-loadfile-dllpath-dl
  library(tmvtnorm)
  
  # helper fns
  code.dir = here("2021-4-20 simple sims with partial try-to-Nmax hacking")
  setwd(code.dir)
  source("helper_SAPH.R")
  
  # # doesn't match expectations:
  # scen.params = data.frame( scen = 1,
  #                           Mu = 1,
  #                           T2 = 0.25,
  #                           m = 50,
  #                           t2w = .25,
  #                           se = 0.5,
  #                           
  #                           Nmax = 5,
  #                           hack = "affirm",
  #                           rho = 0.9,
  #                           
  #                           k = 100,
  #                           k.hacked = 50 )
  
  # # try to recreate successful Expt 2:
  # # this works pretty well for 1 sim rep
  # #  e.g., TheoryExpTstat = 0.5954914 vs. MeanTstatUnhacked = 0.5765191
  # #  TheoryVarTstat = 1.076894 vs. EstVarTstatUnhacked = 1.149218
  # scen.params = data.frame( scen = 1,
  #                           Mu = 1,
  #                           T2 = 0.25,
  #                           m = 500,
  #                           t2w = .25,
  #                           se = 0.5,
  #                           
  #                           Nmax = 10,
  #                           hack = "affirm",
  #                           rho = 0.9,
  #                           
  #                           k = 10^3,
  #                           k.hacked = 0 )
  
  # same as above, but m=50 instead of 500
  # theoretical moments differ from above because tcrit changes based on m:
  #  tcrit = qt(0.975, df = m-1)
  # TheoryExpTstat = 0.6241134 vs. MeanTstatUnhacked = 0.5702767
  # TheoryVarTstat = 1.093762 vs. EstVarTstatUnhacked = 1.259577
  scen.params = data.frame( scen = 1,
                            Mu = 1,
                            T2 = 0.25,
                            m = 50,
                            t2w = .25,
                            se = 0.5,

                            Nmax = 10,
                            hack = "affirm",
                            rho = 0.9,

                            k = 10^3,
                            k.hacked = 0 )
  
  # without heterogeneity:
  # TheoryExpTstat = 1.208201 vs. MeanTstatUnhacked = 1.233518
  # TheoryVarTstat = 0.365473 vs. EstVarTstatUnhacked = 0.3996075
  scen.params = data.frame( scen = 1,
                            Mu = 1,
                            T2 = 0,
                            m = 50,
                            t2w = 0,
                            se = 0.5,
                            
                            Nmax = 10,
                            hack = "affirm",
                            rho = 0.9,
                            
                            k = 10^3,
                            k.hacked = 0 )
  
  # Nmax = 1 with heterogeneity
  # TheoryExpTstat = 0.6241134 vs. MeanTstatUnhacked = 0.5515947
  # TheoryVarTstat = 1.093762 vs. EstVarTstatUnhacked = 1.198789
  #bm: trying to understand this one
  scen.params = data.frame( scen = 1,
                            Mu = 1,
                            T2 = 0.25,
                            m = 50,
                            t2w = .25,
                            se = 0.5,
                            
                            Nmax = 1,
                            hack = "affirm",
                            rho = 0.9,
                            
                            k = 10^3,
                            k.hacked = 0 )
  
  # reduce m even more to see if empirical moments change
  # TheoryExpTstat = 1.071456 vs. MeanTstatUnhacked = 0.9725942
  # TheoryVarTstat = 1.416841 vs. EstVarTstatUnhacked = 1.601244
  scen.params = data.frame( scen = 1,
                            Mu = 1,
                            T2 = 0.25,
                            m = 5,
                            t2w = .25,
                            se = 0.5,
                            
                            Nmax = 10,
                            hack = "affirm",
                            rho = 0.9,
                            
                            k = 10^3,
                            k.hacked = 0 )
  
  #bm: seems like consistently, the mean is lower than trunc mean and var is higher
  #  except when T2 = t2w = 0
  
  
  sim.reps = 1  # reps to run in this iterate
  
  
  library(foreach)
  library(doParallel)
  
  # set the number of local cores
  registerDoParallel(cores=8)
  
  scen = 1
  # data.frame(scen.params %>% filter(scen.name == scen))
}



# ~ RUN SIMULATION ------------------------------


#for ( scen in scen.params$scen.name ) {  # can't use this part on the cluster

# system.time is in seconds
doParallelTime = system.time({
  rs = foreach( i = 1:sim.reps, .combine=rbind ) %dopar% {
    # for debugging:
    #for ( i in 1:sim.reps ) {
    
    # results for just this simulation rep
    if ( exists("repRes") ) rm(repRes)
    
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    p = scen.params[ scen.params$scen == scen, names(scen.params) != "scen"]
    
    
    # ~~ Simulate Dataset ------------------------------
    # includes unpublished studies
    d = sim_meta(Nmax = p$Nmax,
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
    # ( modAll = rma.mv( yi = yi,
    #                    V = vi,
    #                    data = d,
    #                    method = "REML",
    #                    random = ~1 | study ) )
    
    
    
    # ~~ Fit Naive Meta-Analysis on Published Studies  ------------------------------
    # biased meta-analysis of only published studies
    # ( modPub = rma( yi = dp$yi,
    #                 vi = dp$vi,
    #                 method = "REML",
    #                 knha = TRUE ) )
    
 
    # ~~ Bias-Corrected Estimator #1: Nonaffirms Only ------------------------------
    
    #bm: next try rm(list=ls()) and see if below fn still works :)
    # also think about any diagnostics to add
    #  like the tstat expectation and variance among affirms from hacked vs. 
    #  unhacked studies
    #  then should be able to try running doParallel locally!
    
    modCorr = correct_meta_phack1( .p = p,
                                   .dp = dp )
    
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
                                   
                                   dp.k = nrow(dp),
                                   dp.kAffirm = sum(dp$affirm == TRUE),
                                   dp.kNonaffirm = sum(dp$affirm == FALSE),
                                   
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
    
  }  ### end foreach loop
  
} )[3]  # end system.time




### Put in dataset: Estimated time for 1 sim rep ###
# use NAs for additional methods so that the SUM of the rep times will be the
#  total computational time
nMethods = length(unique(rs$methName))
rs$repSeconds = rep( c( doParallelTime / sim.reps,
                        rep( NA, nMethods - 1 ) ), sim.reps )

expect_equal( as.numeric( sum(rs$repSeconds, na.rm = TRUE) ),
              as.numeric(doParallelTime) )


### Local only: Quick look at results ###
takeMean = names(rs)[ !names(rs) %in% c( "repName",
                                         "scenName",
                                         "methName",
                                         names(scen.params) ) ]


resTable = rs %>% group_by(methName) %>%
  #mutate(simReps = n()) %>%
  summarise_at( takeMean,
                function(x) round( mean(x, na.rm = TRUE), 2 ) )
View(resTable)

# should be similar:
# *does NOT match with method = "affirm2"
resTable$TheoryExpTstat; resTable$MeanTstatUnhacked
# the other ones:
resTable$MeanTstat; resTable$MeanTstatHacked; resTable$MeanTstatUnhacked

resTable$TheoryVarTstat; resTable$EstVarTstatUnhacked

#bm
# bias:
p$Mu
resTable$MhatAll
resTable$MhatNaive
resTable$MhatCorr


# CI coverage:
resTable$MhatCoverAll  # should definitely be good
resTable$MhatCoverNaive
resTable$MhatCoverCorr

# setwd("Results")
# fwrite(rs, "rs_one_scenario.csv")
# fwrite(resTable, "resTable.csv")


# WRITE LONG RESULTS ------------------------------
if ( run.local == FALSE ) {
  setwd("/home/groups/manishad/MRM/sim_results/long")
  write.csv( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}