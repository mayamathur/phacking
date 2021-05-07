
# NOTES ---------------------------------------------------------------

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study

# If Nmax is small, rhoEmp (empirical autocorrelation of muin's) will be smaller
#  than rho. That's okay because it reflects small-sample bias in autocorrelation
# estimate itself, not a problem with the simulation code
# For more about the small-sample bias: # https://www.jstor.org/stable/2332719?seq=1#metadata_info_tab_contents


# 2021-5-3 ANALYSIS FNS ---------------------------------------------------------------

# log-likelihood
joint_ll = function(.tstats, .Mu, .T2, .t2w, .se, .crit) {
  sum( log( dtrunc(x = .tstats,
                   spec = "norm",
                   mean = .Mu / .se,
                   #@doesn't yet use the delta-method thing
                   sd = sqrt( (1/.se^2) * (.T2 + .t2w + .se^2) ),
                   b = .crit ) ) )
}

# joint_ll( .tstats = c(1.3, 1.2, 1, 0.9),
#     .Mu = 1,
#     .T2 = 0.25,
#     .t2w = 0.25,
#     .se = 0.5)

# get a corrected MLE from nonaffirms only
# .t2w and .se are taken to be fixed and known
nonaffirm_mles = function(.tstats, .t2w, .se, .crit) {
  
  # par given as (Mu, T2)
  optim( par = c( mean(.tstats * .se), 0),  # starting point for optimization
         f = function(.x) # the negative ll
           -joint_ll( .tstats = .tstats,
                      .Mu = .x[1],
                      .T2 = .x[2],
                      .t2w = .t2w,
                      .se = .se,
                      .crit = .crit),
         lower = c(-Inf, 0),
         method = "L-BFGS-B" )
  
}



# ANALYSIS FNS ---------------------------------------------------------------

# **this is the major workhorse fn
# use the studies assumed to be unhacked to corrected the ones assumed to be hacked
# hackAssumed: "omniscient" (we know which are hacked) or "allAffirms"
correct_dataset_phack = function( .dp,  # published studies
                                  .p,  # parameters as dataframe
                                  hackAssumption ) {
  
  # indicator for whether a study is ASSUMED to be hacked for analysis purposes
  if (hackAssumption == "omniscient" ) .dp$assumedHacked = (.dp$hack == .p$hack)
  if (hackAssumption == "allAffirms" ) .dp$assumedHacked = .dp$affirm
  
  ### Make filtered datasets ###
  # published ASSUMED-hacked ones only
  dph = .dp %>% filter(assumedHacked == TRUE )
  
  # published unhacked ones only
  # so only one per study set
  # same as second row of above table
  dpu = .dp %>% filter(assumedHacked == FALSE )
  
  # meta-analyze only the assumed-unhacked studies
  ( modUH = rma( yi = yi,
                 vi = vi,
                 data = dpu,
                 method = "REML",
                 knha = TRUE ) )
  Mhat.UH = modUH$b
  # *important: since t2w is a sensitivity parameter, we can just subtract it off
  T2.UH = max(0, modUH$tau2 - .p$t2w )
  
  ### *Estimate* bias of each assumed-hacked result ###
  # i.e., truncated t expectation
  # estimate the noncentrality parameter
  # uses estimated mean, estimated tau^2, and known t2w 
  dph$ncp = c(Mhat.UH) / sqrt( c(T2.UH) + .p$t2w + dph$vi )
  
  # estimated expectation of each hacked result
  # extrunc always seems to give warnings about precision not
  #  being achieved
  dph = dph %>%
    rowwise() %>%
    mutate( hackedExp =  extrunc( spec = "t",
                                  ncp = ncp,
                                  df = m-1,
                                  a = tcrit ) ) 
  
  # estimated bias of hacked results
  dph$estBias = dph$hackedExp - c(Mhat.UH)
  
  # sanity check:
  # also calculate the REAL truncated expectations
  #  (i.e., using the real T2, Mu, and se rather than sample estimates)
  dph = dph %>%
    rowwise() %>%
    mutate( hackedExpTrue =  extrunc( spec = "t",
                                      ncp = .p$Mu / sqrt( .p$T2 + .p$t2w + viTrue ),
                                      df = m-1,
                                      a = tcrit ) ) 
  
  
  
  # # another sanity check:and to empirical one
  # t$`mean(yi)`[ t$hack == "affirm" & t$Di == 1 ]
  # # all quite close, even with T2 estimate pretty off in this sample!
  
  ### Bias-correct the published, hacked results ###
  dph$yiCorr = dph$yi - dph$estBias
  
  # put in big dataset (with assumed-unhacked ones) as well
  .dp$yiCorr = .dp$yi
  .dp$yiCorr[ .dp$assumedHacked == TRUE ] = dph$yiCorr
  
  ### Corrected meta-analysis ###
  modCorr = rma( yi = .dp$yiCorr,
                 vi = .dp$vi,
                 method = "REML",
                 knha = TRUE )
  
  ### Return all the things ###
  return( list(data = .dp,  # corrected dataset
               metaCorr = report_rma(modCorr, 
                                     .Mu = .p$Mu,
                                     .suffix = "Corr"),
               sanityChecks = data.frame( Mhat.UH = Mhat.UH,
                                          T2.UH = T2.UH,
                                          
                                          kAssumedHacked = sum(.dp$assumedHacked),
                                          
                                          # these 3 should be similar
                                          meanHackedExp = mean(dph$hackedExp),
                                          meanHackedExpTrue = mean(dph$hackedExpTrue),
                                          
                                          # for "omniscient" mode, this should be similar to the 2 above
                                          yiMeanAssumedHacked = mean( .dp$yi[ .dp$assumedHacked == TRUE ] ),
                                          
                                          
                                          yiCorrMeanAssumedHacked = mean( .dp$yiCorr[ .dp$assumedHacked == TRUE ] ),
                                          
                                          
                                          yiMeanAssumedUnhacked = mean( .dp$yiCorr[ .dp$assumedHacked == FALSE ] ),
                                          
                                          # this is for looking at how biased the SEs become in the hacked studies
                                          # (when using omniscient mode)
                                          viMeanAssumedHacked = mean( .dp$vi[ .dp$assumedHacked == TRUE ] ) ),
               
               
               
               modUH = modUH,
               modCorr = modCorr ) )
  
}




# use only the nonaffirms to get trunc MLE
#  and throws away the affirms
# **note that the returned Vhat is an estimate of T2 + t2w, not T2 itself
correct_meta_phack1 = function( .dp,  # published studies
                                  .p  # parameters as dataframe
                                   ) {
  
  # published affirmatives only
  dpa = .dp[ .dp$affirm == FALSE, ]
  
  
  ### MLEs from trunc normal ###
  # these are the MLEs of the *t-stats*
  mle.fit = mle.tmvnorm( as.matrix(dpa$tstat, ncol = 1), lower=-Inf, upper=crit)
  #summary(mle.fit)
  mles = coef(mle.fit)
  
  #@can't get CIs yet because of package issue
  #  emailed Wilhelm about this on 2021-5-7
  # CIs.tstats = confint( profile(mles, dpa$tstat, method="BFGS", trace=TRUE) )
  
  # get Wald CI a different way
  tstat.mu.SE = attr( summary(mle.fit), "coef" )[ "mu_1", "Std. Error" ]
  tstat.mu.CI = c( mles[1] - tstat.mu.SE * qnorm(0.975),
                   mles[1] + tstat.mu.SE * qnorm(0.975) )

  # # pretty good :)
  # mles[1]; p$Mu/p$se
  # mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)
  
  # rescale MLEs to represent effect sizes rather than tstats
  Mhat = mles[1] * p$se
  # **use Vhat to represent MARGINAL heterogeneity (i.e., T2 + t2w)
  Vhat = ( mles[2] * p$se^2 ) - p$se^2
  
  # and rescale CI limits
  MhatCI = tstat.mu.CI * p$se
  
  ### Sanity checks: Moments of published nonaffirms vs. theory ###
  # check that moments are what we expect

  theoryExpTstat = extrunc(spec = "norm",
                          mean = p$Mu / p$se,
                          #@doesn't yet use the delta-method thing
                          sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
                          b = crit )
  
  theoryVarTstat = vartrunc(spec = "norm",
           mean = p$Mu / p$se,
           #@doesn't yet use the delta-method thing
           sd = sqrt( (1/p$se^2) * (p$T2 + p$t2w + p$se^2) ),
           b = crit )

 
  
  ### Return all the things ###
  return( list( metaCorr = data.frame( Mhat = Mhat,
                                       MhatLo = MhatCI[1],
                                       MhatHi = MhatCI[2],
                                       MhatCover = ( MhatCI[1] <= p$Mu ) & ( MhatCI[2] >= p$Mu ),
                                       Vhat = Vhat),

               
               # note that all of these stats pertain to only published nonaffirmatives
               sanityChecks = data.frame( kNonaffirmPub = sum( dpa$affirm == FALSE ),
                                          
                                          # check that moments of published nonaffirms are what we expect
                                          # these should be similar
                                          theoryExpTstat = theoryExpTstat,
                                          empExpTstat = mean(dpa$tstat),
                                          
                                          theoryVarTstat = theoryVarTstat,
                                          empVarTstat = var(dpa$tstat),
                                        
                                          
                                          yiMeanNonaffirmPub = mean(dpa$yi),
                                          viMeanNonaffirmPub = mean(dpa$vi) ) ) )
  
}


# nicely report a metafor object with optional suffix to denote which model it is
# includes coverage
report_rma = function(.mod,
                      .Mu,  # true mean (to get coverage)
                      .suffix = "") {
  
  .res = data.frame( .mod$b,
                     .mod$ci.lb,
                     .mod$ci.ub,
                     (.mod$ci.lb <= .Mu) & (.mod$ci.ub >= .Mu), 
                     .mod$pval )
  
  names(.res) = paste( c("Mhat", "MhatLo", "MhatHi", "MhatCover", "MhatPval"), .suffix, sep = "" )
  return(.res)
}



# DATA SIMULATION ---------------------------------------------------------------


# runs a simple simulation to compare empirical moments to theoretical ones 
#  that I expected to match the empirical ones
# writes the dataset to .results.dir if it isn't NA (in which case needs to have a variable, .p$sim.name)
quick_sim = function(.p,
                     .results.dir = NA,
                     printRes = FALSE ) {
  
  # IMPORTANT: IF YOU ADD ARGS TO SIM_ONE_STUDY OR MAKE_ONE_DRAW, MUST ADD THEM 
  #  HERE OR ELSE QUICK_SIM WILL SILENTLY NOT PASS THEM
  # simulate a huge dataset, all finitely hacked
  d = sim_meta(Nmax = .p$Nmax,
               Mu = .p$Mu,
               T2 = .p$T2,
               m = .p$m,
               t2w = .p$t2w,
               se = .p$se,
               hack = .p$hack,
               return.only.published = FALSE,
               rho = .p$rho,
               
               k = .p$k,
               k.hacked = .p$k.hacked )
  
  
  # add in the parameters that aren't already in dataset
  shortParams = .p[ , !names(.p) %in% names(d) ]
  d = cbind( d, shortParams )
  
  # dataset of only published results
  dph = d[ d$Di == TRUE, ]
  
  Mu = unique(d$Mu)
  T2 = unique(d$T2)
  t2w = unique(d$t2w)
  m = unique(d$m)
  se = unique(d$se)
  
  library(msm)
  correctedSE = deltamethod( g = ~ x1/sqrt(x2),
                             mean = c(Mu, se^2),
                             cov = matrix( c( T2 + t2w + se^2, 0, 0, var(d$vi) ),
                                           nrow = 2 ) )
  
  crit = unique(d$tcrit)
  
  
  res = data.frame( matrix( NA, nrow = 2, ncol = 2) )
  names(res) = c("affirms", "nonaffirms")
  row.names(res) = c("theoryExp", "empiricalExp")
  
  
  ## Expectation of affirmatives
  res[ "theoryExp", "affirms" ] = extrunc( spec = "norm",
                                           mean = Mu / se,
                                           sd = correctedSE,
                                           a = crit )
  
  res[ "empiricalExp", "affirms" ] = mean( dph$tstat )
  
  
  res[ "theoryVar", "affirms" ] = vartrunc( spec = "norm",
                                            mean = Mu / se,
                                            sd = correctedSE,
                                            a = crit )
  
  res[ "empiricalVar", "affirms" ] = var( dph$tstat )
  
  
  # would be hard to look at affirms because of duplication within studies (messes up var)
  
  if ( printRes == TRUE ) print(res)
  
  library(Hmisc)
  returnList = llist(d, res, correctedSE)
  
  # save dataset
  if ( !is.na(.results.dir) & !is.na(.p$sim.name) ) {
    setwd(.results.dir)
    save( returnList, file=.p$sim.name )
  }
  
  
  return( returnList )
  
}


# *note that the number of reported, hacked studies might be less than k.hacked
#  if all Nmax draws are unsuccessful

# also note that for the unhacked study sets, the single published result could be 
# affirmative OR nonaffirmative

# simulate meta-analysis in which the hacking follows a mixture distribution:
# some studies are unhacked, in which case we always make Nmax draws and then report the last one (which is equivalent to only making 1 draw)
# and some studies are hacked, in which case we make UP TO Nmax draws and stop
# either when we get the first affirmative OR when we reach Nmax draws
sim_meta = function(Nmax,  # max draws to try
                    Mu,  # overall mean for meta-analysis
                    T2,  # across-study heterogeneity
                    
                    # study parameters, assumed same for all studies:
                    m,  # sample size for this study
                    t2w,  # within-study heterogeneity
                    se,  # TRUE SE
                    
                    rho = 0,  # autocorrelation of muin's
                
                    hack,  # mechanism of hacking for studies that DO hack (so not "no")
                    
                    # args not passed to sim_one_study_set:
                    k,  # number of studies
                    k.hacked,  # number of hacked studies
                    
                    return.only.published = FALSE
) {
  
  # # test only
  # Nmax = 20
  # Mu = 0.1
  # T2 = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "affirm"
  # return.only.published = FALSE
  # rho=0
  # k = 30
  # k.hacked = 20
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k", "k.hacked")]
  
  
  if ( hack == "no" ) stop("hack should only be 'affirm' or 'signif' for this fn")
  
  k.unhacked = k - k.hacked
  
  
  ### Simulate the unhacked studies ###
  if ( k.unhacked > 0 ) {
    for ( i in 1:(k - k.hacked) ) {
      
      if ( i %% 50 == 0 ) cat("\nSimulating study #", i)
      
      # to generate unhacked studies, need to change argument "hack"
      .argsUnhacked = .args
      .argsUnhacked[ names(.args) == "hack" ] = "no"
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsUnhacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
  }
  
  ### Simulate hacked studies ###
  if ( k.hacked > 0 ) {
    if ( exists(".dat") ) startInd = max(.dat$study) + 1 else startInd = 1
    
    for ( i in startInd:(startInd + k.hacked - 1) ) {
      
      
      if ( i %% 50 == 0 ) cat("\nSimulating study #", i)
      
      # for unhacked studies, no need to change argument "hack"
      .argsHacked = .args
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsHacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
  }
  
  return(.dat)
  
}


# ### Example 1
# d = sim_meta(Nmax = 20,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1,
#              hack = "affirm",
#              return.only.published = FALSE,
# 
#              k = 30,
#              k.hacked = 10
# 
# )
# 
# 
# 
# 
# nrow(d)
# 
# # look at the published results only
# d %>% filter(Di == 1 ) %>%
#   group_by(hack) %>%
#   summarise( n(),
#              mean(affirm),
#              mean(mui),
#              mean(yi) )
# 
# ### Example 2: Affirm2 hacking
# #bm
# d = sim_meta(Nmax = 5,
#              Mu = 0.1,
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = 1,
#              hack = "affirm2",
#              return.only.published = FALSE,
# 
#              k = 100,
#              k.hacked = 100
# 
# )
# 
# table(d$N)
# 
# # there should be some published nonaffirms
# d %>% filter(Di == 1) %>%
#   summarise( mean(affirm) )
# 
# # all are hacked, so every published nonaffirm should come from a study set
# #  that reached Nmax
# expect_equal( unique( d$N[ d$Di == 1 & d$affirm == FALSE] ),
#               5 )
# 
# 
# 
# nrow(d)
#
# ### Correlated draws
# # this is slow
# d = sim_meta(Nmax = 100,
#              Mu = -0.5,  # make it hard to get affirms
#              T2 = 0.1,
#              m = 50,
#              t2w = .5,
#              se = .5,
#              rho = 0.9,
#              hack = "affirm",
#              return.only.published = FALSE,
# 
#              k = 1000,
#              k.hacked = 500 )
# 
# # all results, even unpublished ones
# # primarily to check autocorrelation of muin's
# d %>% filter( !duplicated(study) ) %>%
#   group_by(hack) %>%
#   summarise( sum(!is.na(rhoEmp)),
#              sum(N > 1),
#              mean(rhoEmp, na.rm = TRUE) )


# ~ Simulate a single study ----------------- 
# simulated study from potentially heterogeneous meta-analysis distribution
#  + its own heterogeneous distribution

# for hack = "no":
# draws exactly Nmax results and treats the last one as the reported one
#  so the final result could be affirmative or nonaffirmative

# for hack = "affirm" or "signif":
# draws until affirm or signif result is obtained or Nmax is reached
# then reports the last result

# hack options:
#  - "no": no hacking (report all Nmax results)
#  - "affirm": hack until you get an affirmative result or you reach Nmax,
#    but if you reach Nmax, do NOT report any result at all
#  - "signif": similar to "affirm", but hack to significance
#  - "affirm2": similar to "affirm", but you always report the last draw, even
#    if it was nonaffirm (no file drawer)

# NOTE: If you add args here, need to update quick_sim as well
sim_one_study_set = function(Nmax,  # max draws to try
                             Mu,  # overall mean for meta-analysis
                             T2,  # across-study heterogeneity
                             m,  # sample size for this study
                             t2w,  # within-study heterogeneity
                             se,  # TRUE SE for this study
                             return.only.published = FALSE,
                             hack, # should this study set be hacked? ("no", "affirm","affirm2", "signif")
                             
                             # for correlated draws; see make_one_draw
                             rho = 0
                              ) {  
  
  
  # # test only
  # Nmax = 20
  # Mu = 0.1
  # T2 = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "affirm"
  # rho=0
  
  # mean for this study set
  # doesn't have t2w because that applies to results within this study set
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(T2),
                   n = 1)

  # TRUE SD (not estimated)
  sd.y = se * sqrt(m)

  # collect all args from outer fn, including default ones
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  .args$mui = mui
  .args$sd.y = sd.y


  stop = FALSE  # indicator for whether to stop drawing results
  N = 0  # counts draws actually made


  # we use this loop whether there's hacking or not
  while( stop == FALSE & N < Nmax ) {

    if ( rho == 0 ) {
      # make uncorrelated draw
      newRow = do.call( make_one_draw, .args )
    } else {
      # make correlated draw
      if ( N == 0 ) .args$last.muin = NA  # on first draw, so there's no previous one
      if ( N > 0 ) .args$last.muin = d$muin[ nrow(d) ]
      newRow = do.call( make_one_draw, .args )
    }


    # number of draws made so far
    N = N + 1

    # add new draw to dataset
    if ( N == 1 ) d = newRow
    if ( N > 1 ) d = rbind( d, newRow )

    # check if it's time to stop drawing results
    if (hack == "signif") {
      stop = (newRow$pval < 0.05)
    } else if ( hack %in% c("affirm", "affirm2") ) {
      stop = (newRow$pval < 0.05 & newRow$yi > 0)
    } else if ( hack == "no" ) {
      # if this study set is unhacked, then success is just whether
      #  we've reached Nmax draws
      if ( hack == "no") stop = (N == Nmax)
    } else {
      stop("No stopping criterion implemented for your chosen hack mechanism")
    }
    
  }  # end while-loop

  # record info in dataset
  d$N = N
  d$hack = hack

  # empirical correlation of muin's
  #  but note this will be biased for rho in small samples (i.e., Nmax small)
  if ( nrow(d) > 1 ) {
    # get lag-1 autocorrelation
    d$rhoEmp = cor( d$muin[ 2:length(d$muin) ],
                    d$muin[ 1: ( length(d$muin) - 1 ) ] )

    # mostly for debugging; could remove later
    d$covEmp = cov( d$muin[ 2:length(d$muin) ],
                    d$muin[ 1: ( length(d$muin) - 1 ) ] )

  } else {
    d$rhoEmp = NA
    d$covEmp = NA
  }


  
  # convenience indicators for significance and affirmative status
  d$signif = d$pval < 0.05
  d$affirm = (d$pval < 0.05 & d$yi > 0)
  
  # publication indicator
  # in the first 2 cases, Di=1 for only the last draw IF we got an affirm result
  #  but if we didn't, then it will always be 0
  if ( hack == "signif" ) d$Di = (d$signif == TRUE)
  if (hack == "affirm") d$Di = (d$affirm == TRUE)
  # if no hacking or affirmative hacking without file drawer, assume only last draw is published
  if ( hack %in% c("no", "affirm2") ) {
    d$Di = 0
    d$Di[ length(d$Di) ] = 1
  }
  
  
  if ( return.only.published == TRUE ) d = d[ d$Di == 1, ]
  
  return(d)
  
}


# ### example
# d = sim_one_study_set(Nmax = 5,
#                       Mu = 0.1,
#                       T2 = 0.1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "affirm2",
#                       return.only.published = FALSE)
# nrow(d)
# d


# ### sanity check by simulation
# for ( i in 1:2000 ) {
#   newRows = sim_one_study_set(Nmax = 20,
#                          Mu = 0.1,
#                          T2 = 0.1,
#                          m = 50,
#                          t2w = .1,
#                          se = 1,
#                          hack = "no",
#                          return.only.published = FALSE )
# 
#   if ( i == 1 ) .d = newRows else .d = rbind(.d, newRows)
# 
# }
# 
# # all studies
# # note that conditional on Di == 0, variance might be off because of repeated rows
# # but means should be correct
# .d %>% group_by(Di == 1) %>%
#   summarise(n(),
#             mean(mui),
#             mean(muin),
#             var(muin),
#             mean(yi) )
# # seems fine
# 
# ### check correlated draws: large Nmax
# # no hacking so that we can get a lot of draws
# d = sim_one_study_set(Nmax = 1000,
#                       Mu = 0,
#                       T2 = 0.1,
#                       m = 50,
#                       t2w = .5,
#                       se = 0.5,
#                       rho = 0.9,
#                       hack = "no",
#                       return.only.published = FALSE)
# nrow(d)
# # should match each other and should be close to rho
# table(d$rhoEmp)
# acf(d$muin, lag = 1)$acf
#
# ### check correlated draws: small Nmax
# #  wrong because of small-sample bias in autocorrelation (not my fn's fault)
# #  about the small-sample bias: # https://www.jstor.org/stable/2332719?seq=1#metadata_info_tab_contents
# # issue is that the sample variance and autocorrelation estimate are non-independent
# #  in small samples
# rhoEmp = c()
# covEmp = c()
# varEmp = c()
# for ( i in 1:250 ) {
#   d = sim_one_study_set(Nmax = 10,
#                         Mu = 0,
#                         T2 = 0.1,
#                         m = 50,
#                         t2w = .5,
#                         se = 0.5,
#                         rho = 0.9,
#                         hack = "no",
#                         return.only.published = FALSE)
#   
#   covEmp = c(covEmp, unique(d$covEmp))
#   rhoEmp = c(rhoEmp, unique(d$rhoEmp))
#   varEmp = c(varEmp, var(d$muin))
# }
# 
# # TOO LOW when Nmax is small! 
# mean(rhoEmp)  # when Nmax = 10: 0.46; when Nmax = 200: 0.88
# mean(covEmp)  # when Nmax = 10: 0.09; when Nmax = 200: 0.40
# mean(varEmp)  # when Nmax = 10: 0.15; when Nmax = 200: 0.47 (should equal t2w = 0.5)
# # sample variance is biased downward in small samples







# # ~ Sanity check  ---------------------------------------------------------------
# #  if Nmax -> Inf and we hack until affirmative, 
# #   published results should follow truncated t distribution
# # also need to set heterogeneity to 0?
# d = data.frame( matrix(nrow = 500, ncol = 1))
# Mu = 1
# T2 = 0.5
# t2w = 0.3
# se = 1
# 
# d = d %>% rowwise() %>%
#   mutate( sim_one_study_set( Nmax = Inf,
#                              Mu = Mu,
#                              T2 = T2,
#                              m = 50,
#                              t2w = t2w,
#                              se = se,
#                              hack = "affirm",
#                              return.only.published = TRUE ) )
# 
# summary(d$N)
# 
# # calculate noncentrality parameter
# ncp = Mu / sqrt( T2 + t2w + se^2 )
# 
# # compare to truncated t
# qqtrunc(x = d$tstat,
#         spec = "t",
#         ncp = ncp,
#         df = 50-1,
#         # since I simulated iid studies here, the truncation cutoff is always the same
#         a = d$tcrit[1] )
# 
# # vs. actually drawing directly from truncated t
# x = rtrunc( n = 500,
#             spec = "t",
#             ncp = ncp,
#             df = 50-1,
#             a = d$tcrit[1])
# 
# plot( density(x) ) +
#   lines( density(d$tstat), col = "red")
# # looks great! :)



# ~ Draw one unbiased result within one study ------------------
# muin should be NA if either we want uncorrelated draws OR it's the first of a series of correlated draws
make_one_draw = function(mui,  # mean for this study set
                         t2w,
                         sd.y,  # TRUE SD
                         m,
                         
                         # for making correlated draws
                         rho = 0,  # autocorrelation of muin's (not yi's)
                         last.muin = NA,  
                         ...) {
  
  
  # true mean for draw n (based on within-study heterogeneity)
  # either we want uncorrelated draws OR it's the first of a series of correlated draws
  if ( rho == 0 | is.na(last.muin) ) {
    muin = rnorm(mean = mui,
                 sd = sqrt(t2w),
                 n = 1)
  }
  
  # make correlated draw
  if ( rho != 0 & !is.na(last.muin) ) {
    # draw from BVN conditional, given muin from last draw
    # conditional moments given here:
    #  https://www.amherst.edu/system/files/media/1150/bivarnorm.PDF
    muin = rnorm(mean = mui + rho * (last.muin - mui),
                 sd = sqrt( t2w * (1 - rho^2) ),
                 n = 1)
  }
  
  
  # draw subject-level data from this study's population effect
  y = rnorm( mean = muin,
             sd = sd.y,
             n = m)
  
  # run a one-sample t-test
  test = t.test(y,
                alternative = "two.sided")
  
  pval = test$p.value
  tstat = test$statistic
  vi = test$stderr^2  # ESTIMATED variance
  
  
  # if (hack == "signif") success = (pval < 0.05)
  # if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  
  return( data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = m-1),
                     mui = mui,
                     muin = muin,
                     yi = mean(y),
                     vi = vi,
                     viTrue = sd.y^2 / m,  # true variance; will equal p$se^2
                     m = m ) )
  #success = success,
  #N = Nmax ) )
}

# make_one_draw(mui = 0.1,
#               t2w = 0,
#               sd.y = 0.3,
#               m = 30 )
# 
# # rho = 1, so should get exactly the same muin again
# make_one_draw(mui = 0.1,
#               t2w = 0.1,
#               sd.y = 0.3,
#               m = 30,
#               rho = 1,
#               last.muin = 0.8)
# 
# # sanity check by simulation: uncorrelated draws
# for ( i in 1:5000 ) {
#   
#   # get last draw
#   last.muin = NA
#   if ( i > 1 ) last.muin = .d$muin[ nrow(.d) ]
#   
#   newRow = make_one_draw(mui = 0.1,
#                   t2w = 0.1,
#                   sd.y = 0.3,
#                   m = 30,
#                   rho = 0.9,
#                   last.muin = last.muin)
# 
#   if ( i == 1 ) .d = newRow else .d = rbind(.d, newRow)
# 
# }
# 
# .d %>% summarise( mean(mui),
#                   mean(muin),
#                   var(muin),
#                   mean(yi) )
# 
# 
# # sanity check by simulation: correlated draws
# for ( i in 1:5000 ) {
#   newRow = make_one_draw(mui = 0.1,
#                          t2w = 0.1,
#                          sd.y = 0.3,
#                          m = 30,
#                          rho = 1,
#                          last.muin = 0.8)
#   
#   if ( i == 1 ) .d = newRow else .d = rbind(.d, newRow)
#   
# }
# 
# .d %>% summarise(
#   mean(mui),
#   mean(muin),
#   var(muin),
#   mean(yi) )
# 
# # look at autocorrelation of muin's (should match rho above)
# acf(.d$muin, lag = 1)$acf
# 
# # look at autocorrelation of yi's (<rho because of SE>0)
# acf(.d$yi, lag = 1)$acf




# DATA WRANGLING ---------------------------------------------------------------

# corrObject: something returned by correct_dataset_phack
# looks for (or makes) global object, "res"
add_method_result_row = function(repRes = NA,
                                 corrObject,
                                 methName) {
  
  newRow = cbind( corrObject$metaCorr,
                  corrObject$sanityChecks )
  
  newRow = newRow %>% add_column(.before = 1,
                                 methName = methName )
  
  
  # "if" condition is hacky way to deal with repRes = NA case
  if ( is.null( nrow(repRes) ) ) repRes = newRow else repRes = rbind(repRes, newRow)
  return(repRes)
}


