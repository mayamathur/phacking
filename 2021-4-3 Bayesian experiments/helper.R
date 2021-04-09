
run_sampler = function() {
  
  
  # at each step, compute log of joint posterior, jumping dist, and ratio 
  # jumping dist must be symmetric
  for (i in 1:B) {
    
    # monitor progress
    if ( i %% 50 == 0 ) cat( "\nStarting iterate ", i )
    
    ### Draw "Primary" Parameters
    # i.e., from marginals
    # last iteration's draws for scalars
    # if ( i == 1 ) lastScalarPars = initPars  # init values
    # if ( i > 1 ) lastScalarPars = triedPars$scalars[i-1,]
    lastScalarPars = triedPars$scalars[i,]
    
    # generate new candidates
    # jumping dist is independent normals centered around last iteration's draw
    
    candZ = rnorm( n = 1,
                   mean = lastScalarPars$Z,
                   sd = jump.SD )
    
    # draw variance components from truncated normal
    # asymmetric proposal dist is fine if you use Metropolis-Hastings rather than Metropolis
    # see Gelman pg 279
    # the difference is the way the acceptance ratio is calculated later
    candtau = rtrunc( n = 1,
                      spec = "norm",
                      mean = lastScalarPars$tau,
                      sd = jump.SD,
                      a = 0)
    
    candtauw = rtrunc( n = 1,
                       spec = "norm",
                       mean = lastScalarPars$tauw,
                       sd = jump.SD,
                       a = 0 )
    
    # for the acceptance ratio
    # log-prob of making the jump we actually made
    termA.J = log( dnorm( candZ,
                          mean = lastScalarPars$Z,
                          sd = jump.SD ) ) +
      
      log( dtrunc( candtau,
                   spec = "norm",
                   mean = lastScalarPars$tau,
                   sd = jump.SD,
                   a = 0) ) + 
      
      log( dtrunc( candtauw,
                   spec = "norm",
                   mean = lastScalarPars$tauw,
                   sd = jump.SD,
                   a = 0) )
    
    #...vs log-prob of making the opposite jump
    termB.J = log( dnorm( lastScalarPars$Z,
                          mean = candZ,
                          sd = jump.SD ) ) +
      
      log( dtrunc( lastScalarPars$tau,
                   spec = "norm",
                   mean = candtau,
                   sd = jump.SD,
                   a = 0) ) + 
      
      log( dtrunc( lastScalarPars$tauw,
                   spec = "norm",
                   mean = candtauw,
                   sd = jump.SD,
                   a = 0) )
    
    
    ### Draw "Intermediate" Parameters (Zi, Ni) Conditional on Primaries
    # random intercept for each study
    ( candZi = rnorm( n = k,
                      mean = candZ,
                      sd = candtau ) )
    
    # Ni
    SDi = sqrt(candtauw^2 + 1)  # marginal SD of the Z-statistics
    crit = qnorm(.975)
    Psig.pos = 1 - pnorm( (crit - candZi) / SDi )
    Psig.neg = pnorm( (-crit - candZi) / SDi )
    
    if ( hack == "signif" ) power_i = Psig.pos + Psig.neg
    if ( hack == "affirm" ) power_i = Psig.pos
    
    # is the final draw affirmative or not?
    if ( hack == "signif" ) favored = abs(d$Zhat) > qnorm(0.975)
    if ( hack == "affirm" ) favored = d$Zhat > qnorm(0.975)
    
    
    # truncated geometric (max = maxN)
    # dtrunc isn't vectorized
    candNi = vapply( X = seq( 1:k ),
                     FUN = function(i) rtrunc( n = 1,
                                               spec = "geom",
                                               prob = power_i[i],
                                               b = N.max ),
                     FUN.VALUE = -99)
    #@ to avoid impossible choices of Ni, anytime the final draw was not favored, 
    # Ni should be equal to N.max
    candNi[ favored == FALSE ] = N.max
    
    # sanity check
    #table(candNi[d$Zhat < 1.96])
    
    
    # calculate new posterior using the candidates
    # use log densities to avoid craziness
    newPost = log_joint_post( .Z = candZ,
                              .tau = candtau,
                              .tauw = candtauw,
                              .Zi = candZi,
                              .Ni = candNi,
                              
                              # **here's where the observed data come in
                              .Zhat = d$Zhati,
                              
                              
                              
                              hack = "affirm",
                              N.max = 20 )
    
    # posterior of previous values
    oldPost = log_joint_post( .Z = lastScalarPars$Z,
                              .tau = lastScalarPars$tau,
                              .tauw = lastScalarPars$tauw,
                              .Zi = as.numeric( triedPars$Zi[i,] ),
                              .Ni = as.numeric( triedPars$Ni[i,] ),
                              
                              # **here's where the observed data come in
                              .Zhat = d$Zhati,
                              
                              
                              
                              hack = "affirm",
                              N.max = 20 )
    
    # ratio of new vs. old parameter values' log-posteriors
    # with Metropolis-Hastings adjustment for asymmetric proposal
    # Gelamn pg 279
    termA = newPost - termA.J
    termB = oldPost - termB.J
    
    # r > 1 is good
    ( r = exp( termA - termB ) )
    
    
    # probability of accepting this step
    ( p = ifelse(r<1, r, 1) )
    accept = rbinom(n=1, size=1, p=p)
    
    # save the acceptance indicator and ratio
    triedPars$accept[i] = accept
    triedPars$r[i] = r
    
    if (accept) {  # automatically accept improvement steps
      
      triedPars$scalars[i+1, "Z"] = candZ
      triedPars$scalars[i+1, "tau"] = candtau
      triedPars$scalars[i+1, "tauw"] = candtauw 
      
      triedPars$Zi[i+1,] = candZi
      triedPars$Ni[i+1,] = candNi
      
      oldPost = newPost
      
    } else {
      triedPars$scalars[i+1,] = as.numeric(lastScalarPars)
      triedPars$Zi[i+1,] =  triedPars$Zi[i,]
      triedPars$Ni[i+1,] = triedPars$Ni[i,]
    }
    
    # # for debugging only
    # triedPars$scalars[i:(i+1),]
    # triedPars$Zi[i:(i+1),]
    # triedPars$Ni[i:(i+1),]
  }
  
  
  return(triedPars)
}


init_sampler = function() {
  scalarPars <<- c("Z", "tau", "tauw")
  
  
  # for tracking parameters that we've tried
  triedPars <<- list( scalars = data.frame( matrix(NA,
                                                   nrow = B,
                                                   ncol = length(scalarPars) ) ),
                      
                      # vector params stored as df with a column for each study in meta-analysis
                      Zi = data.frame( matrix(NA,
                                              nrow = B,
                                              ncol = nrow(d) ) ),
                      
                      Ni = data.frame( matrix(NA,
                                              nrow = B,
                                              ncol = nrow(d) ) ) ) 
  
  names(triedPars$scalars) <<- scalarPars
  
  
  # track whether jump was accepted
  triedPars$accept <<- c()
  # and the ratio
  triedPars$r <<- c()
  
  # SD of jumping dist for parameters
  jump.SD <<- 0.5
  
  # number studies 
  k <<- nrow(d)
  
  # init pars
  #initPars = data.frame( Z = 0.5, tau = 0.5, tauw = 0.5 )
  triedPars$scalars[1,] <<- c(0.5, 0.5, 0.5)
  triedPars$Zi[1,] <<- 0.5
  triedPars$Ni[1,] <<- N.max
}




# LIKELIHOOD FNS FOR USE IN POSTERIOR --------------------------------------------------

# args with "i" subscript are vectors and the rest are scalars
log_joint_post = function(.Z,
                          .tau,
                          .tauw,
                          
                          # intermediate quantities that we don't observe
                          .Zi, # study means (intermediate parameter)
                          .Ni,
                          
                          # observed data
                          .Zhat,
                          
                          # fixed quantities
                          N.max,
                          hack){
  
  # # test only
  # .Z = 0
  # .Zi = c(1.2, 0, 0.5)
  # .Ni = c(1, 20, 5)
  # .Zhat = c(2.1, -0.2, 0.5)
  # .T2 = 0.5
  # .t2w = 0
  # hack = "affirm"
  # N.max = 20
  

  #browser()
  
  ### Priors
  # scalars
  Z.pdf = dnorm( .Z, mean = 0, sd = 10^2)
  # use Bai prior: 
  # "Half-Cauchy is a truncated standard t w/ 1 degree of freedom" from "RobustBayesianCopasModel.R"
  tau.pdf = dtt( .tau, df = 1, location = 0, scale = 1, left = 0 )
  tauw.pdf = dtt( .tauw, df = 1, location = 0, scale = 1, left = 0 )
  
  ### Likelihood term 1: .Zi | .Z, .T2
  # random intercept for each study 
  # but not for every draw within this study
  Zi.pdf = dnorm( as.numeric(.Zi), mean = .Z, sd = .tau )
  
  ### Likelihood term 2: .Ni | .Zi, .t2
  # statistical power term for geometric, given Zi
  SDi = sqrt(.tauw^2 + 1)  # marginal SD of the Z-statistics
  crit = qnorm(.975)
  Psig.pos = 1 - pnorm( (crit - .Zi) / SDi )
  Psig.neg = pnorm( (-crit - .Zi) / SDi )
  
  if ( hack == "signif" ) power_i = Psig.pos + Psig.neg
  if ( hack == "affirm" ) power_i = Psig.pos
  
  # truncated geometric (max = maxN)
  # dtrunc isn't vectorized
  Ni.pdf = vapply( X = seq( 1:length(.Ni) ),
                   FUN = function(i) dtrunc( .Ni[i], spec = "geom", prob = power_i[i], b = N.max ),
                   FUN.VALUE = -99)
  
  
  ### Likelihood term 3: .Zhat | .Zi, .Ni
  #Zhat.pdf = Zhat_pdf(...)
  Zhat.pdf = Zhat_pdf(.Zhat = .Zhat,
                      .Zi = .Zi,
                      .Ni = .Ni,
                      .t2w = .tauw^2,
                      N.max = N.max,
                      hack = hack)
  
  ### Full joint posterior
  return( log(Z.pdf) + 
            log(tau.pdf) + 
            log(tauw.pdf) + 
            sum( log(Zi.pdf) ) +
            sum( log(Ni.pdf) ) +
            sum( log(Zhat.pdf) ) )
  
}


# log_joint_post(   .Z = 0,
#                   .Zi = c(1.2, 0, 0.5),
#                   .Ni = c(1, 20, 20),
#                   .Zhat = c(2.1, -0.2, 0.5),
#                   .T2 = 0.5,
#                   .t2w = 0,
#                   hack = "affirm",
#                   N.max = 20 )


#@assumes uncorrelated draws for now
Zhat_pdf = Vectorize( function( .Zhat, 
                                .Zi, 
                                .Ni,
                                .t2w,
                                
                                N.max,
                                hack) {
  
  # is the final draw affirmative or not?
  if ( hack == "signif" ) favored = abs(.Zhat) > qnorm(0.975)
  if ( hack == "affirm" ) favored = .Zhat > qnorm(0.975)
  
  # impossible to have stopped before N.max without having gotten an affirmative result
  #@presence of zeroes in the lkl here might cause problems! 
  if ( .Ni < N.max & favored == FALSE ) return(0)
  
  # if we've run out of draws, then Zhat is just from the usual normal dist
  #  because we have to report whatever we get on the final draw (no truncation)
  #@ assumes uncorrelated draws
  SD = .t2w + 1  # marginal SD of the Z-statistics
  if ( .Ni == N.max) return( dnorm( .Zhat, mean = .Zi, sd = SD) )
  
  # but if we stopped before the max number of draws, then Zhat is from truncated normal
  #   because we know it was affirmative in order to stop
  #@ also assumes uncorrelated draws
  crit = qnorm(.975)
  upper.tail.pdf = dtrunc( .Zhat,
                           spec = "norm",
                           mean = .Zi,
                           sd = SD,
                           a = crit )
  lower.tail.pdf = dtrunc( .Zhat,
                           spec = "norm",
                           mean = .Zi,
                           sd = SD,
                           b = -crit )
  
  if ( hack == "signif" ) return( upper.tail.pdf + lower.tail.pdf )
  if ( hack == "affirm" ) return( upper.tail.pdf )
  
}, vectorize.args = c(".Zhat", ".Zi", ".Ni", ".t2w") )


# expect_equal( Zhat_pdf(.Zhat = 2.1,
#                        .Zi = 0,
#                        .Ni = 10,
#                        .t2w = 0,
#                        N.max = 10,
#                        hack = "signif"),
#               dnorm( 2.1, mean = 0, sd = 1 + 0) )
# 
# expect_equal( Zhat_pdf(.Zhat = 0.1,
#                        .Zi = 0,
#                        .Ni = 2,
#                        .t2w = 0,
#                        N.max = 10,
#                        hack = "signif"),
#               0 )
# 
# expect_equal( Zhat_pdf(.Zhat = -2.1,
#                        .Zi = 0,
#                        .Ni = 2,
#                        .t2w = 0,
#                        N.max = 10,
#                        hack = "affirm"),
#               0 )


# FOR SIMULATION --------------------------------------------------

# adapted from the earlier helper.R file

# simulate a single study from potentially heterogeneous meta-analysis distribution 
#  and from its own heterogeneous distribution

# works with Z-scores rather than effect sizes

library(Hmisc)

phack_study = function(N.max = Inf,  # max draws before giving up
                       Z,  # across-study mean
                       T2,  # across-study heterogeneity
                    
                       t2,  # within-study heterogeneity
                       
                       NSpreference,
                       hack = "affirm") {  # hack until significant ("signif") or until affirmative ("affirm")?
  
  
  # are we done hacking yet?
  success = FALSE
  
  # number of draws
  N = 0
  
  # keep track of all candidate Zhats that were drawn
  # for use later if we use the max of these (best nonsignificant p-value)
  Zhat.cands = c()
  
  # mean for this study 
  Zi = rnorm( mean = Z, 
              sd = sqrt(T2),
              n = 1)

  while ( success == FALSE ) {
    
    # stop if we reach the max number of draws
    if ( N == N.max ) break
    N = N + 1
    
    # "candidate" Zhat
    # variance: sqrt(t2) is variance of true effect for this draw; "+1" is variance of Zhat given its true mean
    Zhat.cand = rnorm( mean = Zi,
                       sd = sqrt(t2) + 1,
                       n = 1)
    
    Zhat.cands = c(Zhat.cands, Zhat.cand)
    
    if (hack == "signif") success = ( abs(Zhat.cand) > qnorm(0.975) )
    if (hack == "affirm") success = ( Zhat.cand > qnorm(0.975) )
  }
  

  # decide which one to report
  # if investigators just report the last NS draw, then winner is always the last draw
  #  (whether significant or not)
  if ( NSpreference == "last" ) Zhat.reported = Zhat.cands[ length(Zhat.cands) ]
  # if investigators report the best NS draw, then winner is always the best Z-score
  if ( NSpreference == "best" & hack == "affirm" ) Zhat.reported = max(Zhat.cands)
  if ( NSpreference == "best" & hack == "signif" ) Zhat.reported = max( abs(Zhat.cands) )
  

  return( llist( N,
                Zi,
                Zhat.cands,
                Zhat.reported,
                labels = FALSE) )
}


# phack_study( N.max = 20,
#              Z = 0,
#              T2 = 0.5,
#              t2 = 0,
#              NSpreference = "best")

# simulate meta-analysis of hacked studies
# returns the observed Zhat, but also unobserved quantities N and Zi
sim_meta = function(k, ...) {
  

  for ( i in 1:k ){
    studyRes = phack_study(...)
    
    newRow = data.frame( Zhati = studyRes$Zhat.report,
                         Zi = studyRes$Zi,
                         Ni = studyRes$N )
    
    if ( i == 1 ) res = newRow
    if ( i > 1 ) res = rbind(res, newRow)
    
  }
  
  return(res)

}

# d = sim_meta(k = 1500,
#              N.max = 100,
#              Z = 0,
#              T2 = 0.5,
#              t2 = 0,
#              NSpreference = "last")
# hist(Zhats, breaks = 20)

