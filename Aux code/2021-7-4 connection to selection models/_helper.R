
# modified from understand_vevea_woods.R
# this is the -LL for JUST the affirmative results

minusLL_LTN = function( muhat,
                    t2,
                    yi, 
                    vi ) {
  # # TEMP ONLY
  # muhat = 0.2
  # t2 = 0.2
  # yi = d$yi
  # vi = d$vi

  
  # # take inverse of weights if needed
  # if ( any(weight > 1) ) weight = 1/weight
  
  # prob of normal RV being in the interval
  cuts = c(0.025, 1)  # cut point
  
  # max nonsignificant effect size for each study
  bij = sqrt(vi) * -qnorm(cuts[1])  
  
  # Z-score of the max NS effect sizes
  term1 = ( bij - muhat ) / sqrt( vi + t2 )
  
  # prob that a normal RV with mean mu and variance vi + t2 will
  #  be in its actually observed p-value interval
  #browser()
  p_nonsignif = pnorm(term1)
  # Bij = rep( NA, length(yi) )
  # Bij[ weight > 1 ] = p_nonsignif[ weight > 1 ]  # P(nonsignif) for NS studies
  # Bij[ weight == 1 ] = 1 - p_nonsignif[ weight == 1 ]  # P(signif) for S studies
  # table(Bij)
  
  Bi1 = p_nonsignif
  Bi2 = 1 - p_nonsignif
  
  # the "usual" part of the log-likelihood
  LL.term1 = -0.5 * sum( ( yi - muhat )^2 / ( vi + t2 ) ) - 
    0.5 * sum( log( vi + t2 ) )
  
  # # the selection model part - VERSION FROM REPLICATING VEVEA/WOODS
  # LL.term2A = ( Bi1 * min(weight) ) + ( Bi2 * 1 )
  # LL.term2 = sum( log( LL.term2A ) )
  # # when all weights are 1, LL.term2A + LL.term2B = number of studies
  # # and so LL.term2 is independent of muhat and tau2, so it doesn't affect estimation
  # #  at all, which makes sense :)
  
  # the selection model part - NEW VERSION
  # will indeed exactly coincide with the above when min(weight)=0
  LL.term2A = Bi2
  LL.term2 = sum( log( LL.term2A ) )
  # when all weights are 1, LL.term2A + LL.term2B = number of studies
  # and so LL.term2 is independent of muhat and tau2, so it doesn't affect estimation
  #  at all, which makes sense :)
  
  # plain version
  #return( -LL.term1 )
  
  # selection model version
  return( -(LL.term1 - LL.term2) )
}  
