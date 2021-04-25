

# simulate meta-analysis with hacking

# bias-correction #1 (omniscient):
#  we know which studies are hacked

# bias-correction #2 (assume all affirmatives hacked):
#  should be conservative




# audited 2020-6-17

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
  rm(list=ls())
  
  library(here)
  
  # helper fns
  code.dir = here("2021-4-20 simple sims with partial try-to-Nmax hacking")
  setwd(code.dir)
  
  source("helper_SAPH.R")
  
  scen.params = data.frame( scen = 1,
                            Mu = 1,
                            T2 = 0.1,
                            m = 50,
                            t2w = .1,
                            se = 1,
                            
                            Nmax = 20,
                            hack = "affirm",
                            
                            k = 50,
                            k.hacked = 20 )
  
  
  sim.reps = 5  # reps to run in this iterate; leave this alone!
  
  
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
rep.time = system.time({
  rs = foreach( i = 1:sim.reps, .combine=rbind ) %dopar% {
    # for debugging:
    #for ( i in 1:sim.reps ) {
    
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
                 return.only.published = FALSE,
                 
                 k = p$k,
                 k.hacked = p$k.hacked )
    
    
    # dataset of only published results
    dp = d %>% filter(Di == 1 )
    
    # published hacked ones only
    dph = d %>% filter(hack == p$hack & Di == 1)
    
    # published unhacked ones only
    # so only one per study set
    # same as second row of above table
    dpu = d %>% filter(hack == "no" & Di == 1)
    
    # for sanity checks
    # all results, sorted by hacking status and publication status
    t = d %>%
      group_by(hack, Di) %>%
      summarise( n(),
                 k = length(unique(study)),
                 mean(affirm),
                 mean(mui),
                 var(mui),
                 mean(yi),
                 mean(vi) )
    
    
    # ~~ Fit Gold-Standard Meta-Analysis to All Results  ------------------------------
    # unbiased meta-analysis of all studies, even unpublished ones
    # account for clustering of draws within studies
    # *the tau^2 estimate will be close to T2
    ( modAll = rma.mv( yi = yi,
                       V = vi,
                       data = d,
                       method = "REML",
                       random = ~1 | study ) )
    
    
    # ~~ Fit Naive Meta-Analysis on Published Studies  ------------------------------
    # biased meta-analysis of only published studies
    ( modPub = rma( yi = dp$yi,
                    vi = dp$vi,
                    method = "REML",
                    knha = TRUE ) )
    
    # ~~ Bias-Corrected Estimator #1 ------------------------------
    # omniscient version: we know which studies are unhacked
    #  (also includes some affirmatives)
    # meta-analyze only the observed, unhacked studies
    ( modUH = rma( yi = yi,
                   vi = vi,
                   data = dpu,
                   method = "REML",
                   knha = TRUE ) )
    Mhat.UH = modUH$b
    # *important: since t2w is a sensitivity parameter, we can just subtract it off
    T2.UH = modUH$tau2 - p$t2w
    
    ### 1.1: *Estimate* bias of each hacked result ###
    # i.e., truncated t expectation
    # estimate the noncentrality parameter
    # uses estimated mean, estimated tau^2, and known t2w 
    dph$ncp = c(Mhat.UH) / sqrt( c(T2.UH) + p$t2w + dph$vi )
    
    # estimated expectation of each hacked result
    # extrunc always seems to give warnings about precision not
    #  being achieved
    dph = dph %>% filter(hack == "affirm") %>%
      rowwise() %>%
      mutate( affirmExp =  suppressWarnings( extrunc( spec = "t",
                                                      ncp = ncp,
                                                      df = p$m-1,
                                                      a = tcrit ) ) )
    
    # estimated bias of hacked results
    dph$estBias = dph$affirmExp - c(Mhat.UH)
    
    # sanity check:
    # also calculate the REAL truncated expectations
    #  (i.e., using the real T2, Mu, and se rather than sample estimates)
    #@to do in helper code: include SE and m (sample size) in returned dataset
    #  so we can easily have them vary across study sets
    dph = dph %>% filter(hack == "affirm") %>%
      rowwise() %>%
      mutate( affirmExpTrue =  suppressWarnings( extrunc( spec = "t",
                                                          ncp = p$Mu / sqrt( p$T2 + p$t2w + p$se^2 ),
                                                          df = p$m-1,
                                                          a = tcrit ) ) )
    
    
    
    # # another sanity check:and to empirical one
    # t$`mean(yi)`[ t$hack == "affirm" & t$Di == 1 ]
    # # all quite close, even with T2 estimate pretty off in this sample!
    
    ### 1.2: Bias-correct the published, hacked results ###
    dph$yiCorr = dph$yi
    dph$yiCorr = dph$yi - dph$estBias
    
    # put in big dataset as well
    dp$yiCorr = dp$yi
    dp$yiCorr[ dp$hack == p$hack ] = dph$yiCorr
    
    # save for sanity checks
    tab2 = dp %>% group_by(hack) %>%
      summarise( n(),
                 k = length(unique(study)),
                 mean(affirm),
                 mean(mui),
                 var(mui),
                 mean(yi),
                 mean(yiCorr) )
    
    ### 1.3: Meta-analyze corrected estimates ### 
    ( modCorr = rma( yi = dp$yiCorr,
                     vi = dp$vi,
                     method = "REML",
                     knha = TRUE ) )
    
    
    
    
    
 
    # how bad is the variance estimation?
    d %>% group_by(hack, Di) %>%
      summarise( mean(sqrt(vi)) )
    
    

    # ~~ Bias-Corrected Estimator #2 ------------------------------
    # omniscient version: conservatively assume only nonaffirmatives are unhacked
    
    
    
    # ~ Write Results ------------------------------
    
    
    # add in scenario parameters
    rows$scen.name = scen
    rows = as.data.frame( merge(rows, scen.params,
                                by = "scen.name") )
    rows
    
  }  ### end foreach loop
  
} )[3]  # end timer


rs$repTime = rep.time



# WRITE LONG RESULTS ------------------------------
if ( run.local == FALSE ) {
  setwd("/home/groups/manishad/MRM/sim_results/long")
  write.csv( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}