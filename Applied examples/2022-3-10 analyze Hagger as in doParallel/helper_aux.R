

# verbatim from McShane supplement ("selection.meta.functions.R")
# Master heterogeneous negative log likelihood function:
# We use a z approximation to the likelihood because using normal mixture of non-central t distributions
# requires two numerical integrations per datapoint which is too computationally expensive but this 
# integration is analytic with the z approximation. Plus, a normal mixture of non-central t distributions
# is relatively normal.
onestep.heterogeneous.nll <- function(theta, z.obs, n1, n2, alpha=0.025){
  delta <- theta[1]						# True population average effect size
  tau <- theta[2]							# Heterogeneity
  w <- theta[3]							# Relative likelihood non-stat sig study is reported
  z.cut <- qnorm(1-alpha)
  d.obs <- z.obs / sqrt((n1*n2)/(n1+n2))
  d.cut <- z.cut / sqrt((n1*n2)/(n1+n2))
  k <- sum(z.obs<z.cut)
  
  s <- sqrt(tau^2 + 1/n1 + 1/n2)
  ll <- ifelse(k==0, 0, k*log(w))			# Equivalent to defining 0*log(0)=0 as is common
  ll <- ll + sum(dnorm(d.obs, delta, s, log=TRUE))
  ll <- ll - sum(log( w*pnorm(d.cut,delta,s) + (1 - pnorm(d.cut,delta,s)) ))
  -ll
}


# my own implementation of Vevea-Woods (hold weights constant)
# verbatim from "understand_vevea_woods.R"
minusLL = function( muhat,
                    t2,
                    yi = ds$yi, 
                    vi = ds$vi,
                    weight = ds$weight ) {
  # # TEMP ONLY
  # mu = 0.2
  # t2 = 0.2  
  # yi = ds$yi
  # vi = ds$vi
  # weight = ds$weight
  
  # take inverse of weights if needed
  if ( any(weight > 1) ) weight = 1/weight
  
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
  
  # the selection model part
  LL.term2A = ( Bi1 * min(weight) ) + ( Bi2 * 1 )
  LL.term2 = sum( log( LL.term2A ) )
  # when all weights are 1, LL.term2A + LL.term2B = number of studies
  # and so LL.term2 is independent of muhat and tau2, so it doesn't affect estimation
  #  at all, which makes sense :)
  
  # plain version
  #return( -LL.term1 )
  
  # selection model version
  return( -(LL.term1 - LL.term2) )
}  




# small helper fn

show_rep_res = function() {
  
  if( "optimx.Mhat.winner" %in% names(rep.res) ) {
    cat("\n")
    print( rep.res %>% select(method, Mhat, MLo, MHi,
                              Shat,
                              optim.converged,
                              optimx.Mhat.winner,
                              optimx.Nconvergers,
                              optimx.Pagree.of.convergers.Mhat.winner) %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  } else {
    cat("\n")
    print( rep.res %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  }
}