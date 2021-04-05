

# LIKELIHOOD FNS FOR USE IN POSTERIOR --------------------------------------------------

# args with "i" subscript are vectors and the rest are scalars
log_joint_post = function(.Z,
                          .T2,
                          .t2w,
                          
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
  T2.pdf = dtt( sqrt(.T2), df = 1, location = 0, scale = 1, left = 0 )
  t2w.pdf = dtt( sqrt(.t2w), df = 1, location = 0, scale = 1, left = 0 )
  
  ### Likelihood term 1: .Zi | .Z, .T2
  # random intercept for this study 
  # but not for every draw within this study
  Zi.pdf = dnorm( .Zi, mean = .Z, sd = sqrt(.T2) )
  
  ### Likelihood term 2: .Ni | .Zi, .t2
  # statistical power term for geometric, given Zi
  SDi = .t2w + 1  # marginal SD of the Z-statistics
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
  #@ignores within-study heterogeneity for now
  
  #Zhat.pdf = Zhat_pdf(...)
  
  Zhat.pdf = Zhat_pdf(.Zhat = .Zhat,
                      .Zi = .Zi,
                      .Ni = .Ni,
                      .t2w = .t2w,
                      N.max = N.max,
                      hack = hack)
  
  ### Full joint posterior
  return( log(Z.pdf) + 
            log(T2.pdf) + 
            log(t2w.pdf) + 
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

