
# IMPORTANT NOTES -----------------------------

# MUST USE ml load R/4.0.2!!!!

# Before running, search "#@"

# Important things to remember: 
#
# - The returned Vhat is an estimate of T2 + t2w, *not* T2 itself
#
# - correct_meta_phack1 will NOT work well for small m. I already tried this.
#   This is because, for computational convenience, I'm using the truncated normal dist
#   for the tstats, and it won't hold well when m is small. 
# 
# - Key hypotheses to check in simulation results: T2 = 0, t2w > 0, Nmax > 1: Nonaffirms should fit the trunc distribution because
#    there's no selection on mui. So this scenario should look just like T2 = 0, t2w > 0, Nmax = 1.

# Debugging help:
# 
# - The jobs may fail before fitting modAll with no apparent errors if 
#   k is too large for rma.mv. In that case, try setting p$k < 500 for modAll
#  and modPub small to prevent those models from being fit. 


# because Sherlock 2.0 restores previous workspace
rm( list = ls() )


# are we running locally?
run.local = TRUE


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
           "optimx")

if ( run.local == TRUE ) toLoad = c(toLoad, "here")


# allPackages = c("here",
#                 "magrittr",
#                 "dplyr",
#                 "data.table",
#                 "tidyverse",
#                 "tidyr",
#                 "metafor",
#                 "robumeta",
#                 "testthat",
#                 "truncdist",
#                 "gmm",
#                 "tmvtnorm",
#                 "doParallel",
#                 "foreach")




# FOR CLUSTER USE ------------------------------




if (run.local == FALSE) {
  
  # load command line arguments
  args = commandArgs(trailingOnly = TRUE); print(args)
  jobname = args[1]
  scen = args[2]  # this will be a number
  
  # install packages with informative messages if one can't be installed
  # **Common reason to get the "could not library" error: You did ml load R/XXX using an old version
  any.failed = FALSE
  for (pkg in toLoad) {
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
  
  
  # FOR AUTOMATIC CLUSTER RUN
  # # get scen parameters (made by genSbatch.R)
  # path = "/home/groups/manishad/SAPH"
  # setwd(path)
  # scen.params = read.csv( "scen_params.csv" )
  # p <<- scen.params[ scen.params$scen == scen, ]
  # print(p)
  # 
  # # helper code
  # setwd(path)
  # source("helper_SAPH.R")
  
  # FOR INTERACTIVE CLUSTER RUN
  # # alternatively, generate a simple scen.params in order to run doParallel manually in
  # #  Sherlock as a test
  path = "/home/groups/manishad/SAPH"
  setwd(path)
  source("helper_SAPH.R")
  scen.params = data.frame(scen = 1,
                           
                           # args from sim_meta_2
                           Nmax = 1,
                           Mu = 0.1,
                           t2a = 0.25,
                           t2w = 0.25,
                           m = 500,
                           true.sei.expr = "runif(n = 1, min = 0.1, max = 1)",
                           hack = "affirm",
                           rho = 0,
                           k.pub.nonaffirm = 500,
                           prob.hacked = 0,
                           
                           # Stan control args
                           stan.maxtreedepth = 20,
                           stan.adapt_delta = 0.98,
                           
                           get.CIs = TRUE)
  scen = 1
  
  # locally, with total k = 100, Nmax = 10, and sim.reps = 250, took 93 min total
  # for that I did sim.reps = 100 per doParallel
  
  # simulation reps to run within this job
  # **this need to match n.reps.in.doParallel in the genSbatch script
  sim.reps = 5  #@update this 
  

  # set the number of cores
  registerDoParallel(cores=16)
  
}



# FOR LOCAL USE  ------------------------------
if ( run.local == TRUE ) {
  #rm(list=ls())
  
  
  lapply( toLoad,
          require,
          character.only = TRUE)

  
  # helper fns
  code.dir = here()
  setwd(code.dir)
  source("helper_SAPH.R")
  
  
  scen.params = data.frame(scen = 1,
                           
                           # args from sim_meta_2
                           Nmax = 10,
                           Mu = 0.1,
                           t2a = 0.25,
                           t2w = 0.25,
                           m = 500,
                           true.sei.expr = "runif(n = 1, min = 0.1, max = 1)",
                           hack = "affirm",
                           rho = 0,
                           k.pub.nonaffirm = 50,
                           prob.hacked = 0.5,
                           
                           # Stan control args
                           stan.maxtreedepth = 20,
                           stan.adapt_delta = 0.98,
                           
                           get.CIs = TRUE)
  
  
  
  
  sim.reps = 1  # reps to run in this iterate
  
  # set the number of local cores
  registerDoParallel(cores=8)
  
  scen = 1
  # data.frame(scen.params %>% filter(scen.name == scen))
}



# RUN SIMULATION ------------------------------

if ( exists("rs") ) rm(rs)

#for ( scen in scen.params$scen.name ) {  # can't use this part on the cluster

# system.time is in seconds
doParallelTime = system.time({
  #@change this back before running for real
  rs = foreach( i = 1:sim.reps, .combine=rbind ) %dopar% {
    # for debugging (out file will contain all printed things):
    #for ( i in 1:sim.reps ) {
      
    # only print info for first sim rep for visual clarity
    if ( i == 1 ) cat("\n\n~~~~~~~~~~~~~~~~ BEGIN SIM REP", i, "~~~~~~~~~~~~~~~~")
    
    # results for just this simulation rep
    if ( exists("rep.res") ) rm(rep.res)
    
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    p = scen.params[ scen.params$scen == scen, names(scen.params) != "scen"]
    
    # calculate TOTAL heterogeneity
    p$V = p$t2a + p$t2w
    p$S = sqrt(p$V)
    
    if ( i == 1 ) cat("\n\nHEAD OF SCEN.PARAMS:\n")
    if ( i == 1 ) print(p)
    
    # ~~ Simulate Dataset ------------------------------
    # includes unpublished studies
    d = sim_meta_2( Nmax = p$Nmax,
                 Mu = p$Mu,
                 t2a = p$t2a,
                 m = p$m,
                 t2w = p$t2w,
                 true.sei.expr = p$true.sei.expr,
                 hack = p$hack,
                 rho = p$rho,
                 
                 k.pub.nonaffirm = p$k.pub.nonaffirm,
                 prob.hacked = p$prob.hacked,
                 return.only.published = FALSE)
    
    
    # dataset of only published results
    dp = d %>% filter(Di == 1)
    
    if ( i == 1 ) cat("\n\nHEAD OF DP:\n")
    if ( i == 1 ) print(head(dp))
    
    # # published hacked ones only
    # dph = d %>% filter(hack == p$hack & Di == 1)
    # 
    # # published unhacked ones only
    # # so only one per study set
    # # same as second row of above table
    # dpu = d %>% filter(hack == "no" & Di == 1)
    
    # initialize rep.res st run_method_safe and other standalone estimation fns
    #  will correctly recognize it as having 0 rows
    rep.res = data.frame()
    
    
    # ~~ Start Values ------------------------------
    Mhat.start = p$Mu
    Shat.start = p$S
    
    # in case we're not doing jeffreys-mcmc or it fails
    Mhat.MaxLP = NA
    Shat.MaxLP = NA
    
    Mhat.MAP = NA
    Shat.MAP = NA
    
    
    # ~~ Fit Naive Meta-Analysis (All PUBLISHED Draws) ------------------------------

    rep.res = run_method_safe(method.label = c("naive"),
                              method.fn = function() {
                                mod = rma( yi = dp$yi,
                                               vi = dp$vi,
                                               method = "REML",
                                               knha = TRUE )
                                
                                report_rma(mod)
                              },
                              .rep.res = rep.res )
    
    
    # ~~ Fit Gold-Standard Meta-Analysis (ALL FIRST Draws) ------------------------------
    
    # keep first draws only
    #@NOT YET TEST ON DATASET WITH MULTIPLE DRAWS PER STUDY
    d.first = d[ !duplicated(d$study), ]
    
    rep.res = run_method_safe(method.label = c("gold-std"),
                              method.fn = function() {
                                mod.all = rma( yi = d.first$yi,
                                               vi = d.first$vi,
                                               method = "REML",
                                               knha = TRUE )
                                
                                report_rma(mod.all)
                              },
                              .rep.res = rep.res )
 
    
    # ~~ MCMC ------------------------------
    
    # published nonaffirmatives only
    dpn = dp[ dp$affirm == FALSE, ]
    
    # # temp for refreshing code
    # path = "/home/groups/manishad/SAPH"
    # setwd(path)
    # source("helper_SAPH.R")
    

    # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
    #  returns 2 lines of output, one for posterior mean and one for posterior median
    # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
    # rep.res = run_method_safe(method.label = c("jeffreys-mcmc-pmean",
    #                                            "jeffreys-mcmc-pmed",
    #                                            "jeffreys-mcmc-max-lp-iterate"),
    #                           method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi,
    #                                                                              .sei = sqrt(dpn$vi),
    #                                                                              .tcrit = dpn$tcrit, 
    #                                                                              .Mu.start = Mhat.start,
    #                                                                              .Tt.start = Shat.start,
    #                                                                              .stan.adapt_delta = p$stan.adapt_delta,
    #                                                                              .stan.maxtreedepth = p$stan.maxtreedepth),
    #                           .rep.res = rep.res )
    
    
    #@TEMP TO AVOID SLOW MCMC PROCESS
    # example of rep.res at this point (to avoid having to re-run the method; it's slow):
    new.rows = structure(list(method = c("jeffreys-mcmc-pmean", "jeffreys-mcmc-pmed","jeffreys-mcmc-max-lp-iterate"), Mhat = c(0.196768652710757,0.19216059632037, 0.178221722070505), Shat = c(0.791416690654329, 0.790397445863146, 0.77856755485064), MhatSE = c(0.00266867786965314, 0.00266867786965314, 0.00266867786965314), ShatSE = c(0.00190121931174173, 0.00190121931174173, 0.00190121931174173), MLo = c(0.0368586243904465, 0.0368586243904465, 0.0368586243904465), MHi = c(0.384655875223785, 0.384655875223785, 0.384655875223785), SLo = c(0.672480635939751, 0.672480635939751, 0.672480635939751), SHi = c(0.921641744315143, 0.921641744315143, 0.921641744315143), stan.warned = c(0, 0, 0), stan.warning = c(NA, NA, NA), MhatRhat = c(1.00428369239268, 1.00428369239268, 1.00428369239268), ShatRhat = c(1.00170913954403, 1.00170913954403, 1.00170913954403), overall.error = c(NA,NA,NA)), row.names = 2:4, class = "data.frame")
    rep.res = bind_rows(rep.res, new.rows)
    # END TEMP
    
    Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
    Shat.MaxLP = rep.res$Shat[ rep.res$method == "jeffreys-mcmc-max-lp-iterate" ]
    
    cat("\n doParallel flag: Done jeffreys-mcmc if applicable")
    

    
    
    # ~~ Change Starting Values -----
    if ( !is.na(Mhat.MaxLP) ) Mhat.start = Mhat.MaxLP 
    if ( !is.na(Shat.MaxLP) ) Shat.start = Shat.MaxLP 
    
    # ~~ MAP (SD param) ------------------------------
    
    # # temp for refreshing code
    # path = "/home/groups/manishad/SAPH"
    # setwd(path)
    # source("helper_SAPH.R")
    
    rep.res = run_method_safe(method.label = c("jeffreys-sd-param"),
                              method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                            sei = sqrt(dpn$vi),
                                                                            par2is = "Tt",
                                                                            tcrit = dpn$tcrit, 
                                                                            Mu.start = Mhat.start,
                                                                            par2.start = Shat.start,
                                                                            usePrior = TRUE,
                                                                            get.CIs = p$get.CIs,
                                                                            CI.method = "wald"),
                              .rep.res = rep.res )
    
    Mhat.MAP = rep.res$Mhat[ rep.res$method == "jeffreys-sd-param" ]
    Shat.MAP = rep.res$Shat[ rep.res$method == "jeffreys-sd-param" ]
    
    
    # ~~ Change Starting Values -----
    
    if ( !is.na(Mhat.MAP) ) Mhat.start = Mhat.MAP 
    if ( !is.na(Shat.MAP) ) Shat.start = Shat.MAP 

    # ~~ MLE (SD param) ------------------------------
    
    # # temp for refreshing code
    # path = "/home/groups/manishad/SAPH"
    # setwd(path)
    # source("helper_SAPH.R")
    
    rep.res = run_method_safe(method.label = c("mle-sd"),
                              method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                            sei = sqrt(dpn$vi),
                                                                            par2is = "Tt",
                                                                            tcrit = dpn$tcrit, 
                                                                            Mu.start = Mhat.start,
                                                                            par2.start = Shat.start,
                                                                            usePrior = FALSE,
                                                                            get.CIs = p$get.CIs,
                                                                            CI.method = "wald"),
                              .rep.res = rep.res )
    
  
    
    cat("\n doParallel flag: Done mle-sd if applicable")
    
   
    # ~~ MLE (Var param) ------------------------------
    
    rep.res = run_method_safe(method.label = c("mle-var"),
                              method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi,
                                                                            sei = sqrt(dpn$vi),
                                                                            par2is = "T2t",
                                                                            tcrit = dpn$tcrit, 
                                                                            Mu.start = Mhat.start,
                                                                            par2.start = Shat.start^2,
                                                                            usePrior = FALSE,
                                                                            get.CIs = p$get.CIs,
                                                                            CI.method = "wald"),
                              .rep.res = rep.res )
    
    
    
    cat("\n doParallel flag: Done mle-var if applicable")
    
   
    
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
    
    
    # # SAVE: OLD RETRUN STRUCTURE
    # # ~ Write Results ------------------------------
    # 
    # # add stats for the simulated dataset that don't change based on correction method
    # # these are prefaced by dataset name for clarity
    # # **note: all columns before methName don't depend on the correction method used
    # repRes = repRes %>% add_column(repName = i,
    #                                
    #                                # dp.k = nrow(dp),
    #                                # dp.kAffirm = sum(dp$affirm == TRUE),
    #                                # dp.kNonaffirm = sum(dp$affirm == FALSE),
    #                                # dp.Nrealized = mean(dp$N),
    #                                
    #                                report_rma(modAll,
    #                                           .Mu = p$Mu,
    #                                           .suffix = "All"),
    #                                report_rma(modPub,
    #                                           .Mu = p$Mu,
    #                                           .suffix = "Naive"),
    #                                
    #                                .before = 1 )
    # 
    # # add in scenario parameters
    # repRes = repRes %>% add_column( scenName = scen,
    #                                 p,
    #                                 .after = 1 )
    # 
    # 
    # repRes
    # 
    # cat("\nSURVIVED MAKING REPRES:")
    # print(repRes)
    
    # ~ Write Results ------------------------------

    # add in scenario parameters
    # do NOT use rbind here; bind_cols accommodates possibility that some methods' rep.res
    #  have more columns than others
    rep.res = p %>% bind_cols( rep.res )
    
    # add more info
    rep.res = rep.res %>% add_column( rep.name = i, .before = 1 )
    rep.res = rep.res %>% add_column( scen.name = scen, .before = 1 )
    rep.res = rep.res %>% add_column( job.name = jobname, .before = 1 )
    
    # add info about simulated datasets
    
    perc.studies.published = 100 * mean( d.first$study %in% unique(dp$study) )
    expect_equal( perc.studies.published, nrow(dp)/nrow(d.first) )
    perc.hacked.studies.published = 100 * mean( d.first$study[ d.first$hack != "no" ] %in% unique( dp$study[ dp$hack != "no" ] ) )
    perc.unhacked.studies.published = 100 * mean( d.first$study[ d.first$hack == "no" ] %in% unique( dp$study[ dp$hack == "no" ] ) )
    #bm
    
    
    
    perc.hacked.draws.published = 100 * mean( d.first$study[ d.first$hack != "no" ] %in% unique( dp$study[ dp$hack != "no" ] ) )
    perc.unhacked.draws.published = 100 * mean( d.first$study[ d.first$hack == "no" ] %in% unique( dp$study[ dp$hack == "no" ] ) )
    
    # a benchmark for average power
    perc.unhacked.draws.affirm = 100 * mean( d$affirm[ d$hack == "no" ] )
    
    
    rep.res = rep.res %>% add_column(   dp.k = nrow(dp),
                                        dp.kAffirm = sum(dp$affirm == TRUE),
                                        dp.kNonaffirm = sum(dp$affirm == FALSE),
                                        
                                        # means draws per HACKED, published study
                                        dp.Nrealized = mean( dp$N[dp$hack == "affirm"] ),
                                        
                                        perc.studies.published = ,
                                        perc.hacked.published = 100 * mean( d.first$study %in% unique(dp$study) ) 
                                          
                                          )
    
    rep.res
    
  }  ### end foreach loop
  
} )[3]  # end system.time


# quick look locally
rs %>% mutate_if(is.numeric, function(x) round(x,2) )

#@if using the for-loop for debugging instead of foreach, need to comment
#  out these parts that look for "rs"

### Add meta-variables to dataset ###

# estimated time for 1 simulation rep
# use NAs for additional methods so that the SUM of the rep times will be the
#  total computational time
nMethods = length(unique(rs$methName))
rs$doParallelTime = doParallelTime
rs$repSeconds = rep( c( doParallelTime / sim.reps,
                        rep( NA, nMethods - 1 ) ), sim.reps )

expect_equal( as.numeric( sum(rs$repSeconds, na.rm = TRUE) ),
              as.numeric(doParallelTime) )



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



# ~ WRITE LONG RESULTS ------------------------------
if ( run.local == FALSE ) {
  setwd("/home/groups/manishad/SAPH/long_results")
  fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}