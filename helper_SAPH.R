
# NOTES ---------------------------------------------------------------

# keeping this script in general Code dir because it's a living work in progress  

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study

# If Nmax is small, rhoEmp (empirical autocorrelation of muin's) will be smaller
#  than rho. That's okay because it reflects small-sample bias in autocorrelation
# estimate itself, not a problem with the simulation code
# For more about the small-sample bias: # https://www.jstor.org/stable/2332719?seq=1#metadata_info_tab_contents



# FNS OF RTMA JEFFREYS PRIOR ---------------------------------------------------------------


# IMPORTANT NOTATION FOR THESE FNS:
# .Mu: mean of effect sizes, not Z-scores
# .T2t: total heterogeneity of effect sizes (=T2 + t2w in older notation, or t2a + t2w in newer notation)
# "T2t" is synonymous with "V" in parameters and results; "Tt" synonymous with "S" in params and results

# ~ 2021-11-22 Jeffreys prior fns ---------------


# agrees with weightr per "Repurpose TNE code.R"
# RTMA log-likelihood - now uses TNE version
# carefully structured for use with Deriv()
joint_nll_2 = function(.yi,
                       .sei,
                       .Mu,
                       .Tt = NULL,  # allow either parameterization
                       .T2t = NULL,
                       .tcrit = rep( qnorm(.975), length(.yi) ) ) {
  
  
  if ( is.null(.T2t) ) .T2t = .Tt^2
  if ( is.null(.Tt) ) .Tt = sqrt(.T2t)
  
  
  # as in TNE::nll()
  .dat = data.frame(yi = .yi, sei = .sei, crit = .tcrit)
  
  .dat = .dat %>% rowwise() %>%
    mutate( term1 = dmvnorm(x = as.matrix(yi, nrow = 1),
                            mean = as.matrix(.Mu, nrow = 1),
                            # deal with different SEs by just incorporating them into the total variance
                            sigma = as.matrix(.T2t + sei^2, nrow=1),
                            log = TRUE),
            
            term2 = log( pmvnorm( lower = -99,
                                  upper = crit * sei,
                                  mean = .Mu,
                                  sigma = .T2t + sei^2 ) ),
            
            nll.i = -term1 + term2 )
  
  
  nll = sum(.dat$nll.i)
  
  # # sanity checks
  # term1.new = log( dnorm( x = .dat$yi[1],
  #                         mean = .Mu,
  #                         sd = sqrt(.T2t + .dat$sei[1]^2) ) )
  # 
  # term2.new = log( pnorm( q = .dat$crit[1] * .dat$sei[1],
  #                         mean = .Mu,
  #                         sd = sqrt(.T2t + .dat$sei[1]^2) ) ) 
  # 
  # expect_equal( as.numeric(.dat$term1[1]), as.numeric(term1.new), tol = 0.001 )
  # expect_equal( as.numeric(.dat$term2[1]), as.numeric(term2.new), tol = 0.001 )
  # 
  # # another sanity check
  # library(truncnorm)
  # nll.new = -log( dtruncnorm( x = .dat$yi[1],
  #                             mean = .Mu,
  #                             sd = sqrt(.T2t + .dat$sei[1]^2),
  #                             a = -99,
  #                             b = .dat$sei[1] * .dat$crit[1] ) )
  # 
  # expect_equal( nll.new, as.numeric(.dat$nll.i[1]) , tol = 0.001)
  
  # return it
  return(nll)
  
  # from TNE's nll:
  # term1 = dnorm(x = .x,
  #               mean = .mu,
  #               sd = .sigma,  
  #               log = TRUE)
  # 
  # term2 = length(.x) * log( pmvnorm(lower = .a,
  #                                   upper = .b,
  #                                   mean = .mu,
  #                                   # note use of sigma^2 here because of pmvnorm's different parameterization:
  #                                   sigma = .sigma^2 ) ) 
}



# verbatim from TNE
E_fisher_TNE = function(.mu, .sigma, .n, .a, .b) {
  
  # prevent infinite cutpoints
  # if either cutpoint is infinite, there are numerical issues because the alpha*Z terms
  #  below are 0*Inf
  Za = max( -99, (.a - .mu) / .sigma )
  Zb = min( 99, (.b - .mu) / .sigma )
  
  alpha.a = dnorm(Za) / ( pnorm(Zb) - pnorm(Za) )
  alpha.b = dnorm(Zb) / ( pnorm(Zb) - pnorm(Za) )
  
  k11 = -(.n/.sigma^2) + (.n/.sigma^2)*( (alpha.b - alpha.a)^2 + (alpha.b*Zb - alpha.a*Za) )
  
  k12 = -( 2*.n*(alpha.a - alpha.b) / .sigma^2 ) +
    (.n/.sigma^2)*( alpha.a - alpha.b + alpha.b*Zb^2 - alpha.a*Za^2 +
                      (alpha.a - alpha.b)*(alpha.a*Za - alpha.b*Zb) )
  
  k22 = (.n/.sigma^2) - (3*.n*(1 + alpha.a*Za - alpha.b*Zb) / .sigma^2) +
    (.n/.sigma^2)*( Zb*alpha.b*(Zb^2 - 2) - Za*alpha.a*(Za^2 - 2) +
                      (alpha.b*Zb - alpha.a*Za)^2 )
  
  return( matrix( c(-k11, -k12, -k12, -k22),
                  nrow = 2,
                  byrow = TRUE ) )
}


# important: note that in this fn, critical value is on t/z scale, NOT raw scale
#  vs. in E_fisher_TNE, .b is on raw scale
E_fisher_RTMA = function( .sei, .Mu, .Tt, .tcrit = qnorm(0.975) ) {
  # get expected Fisher info for each observation separately, based on its unique SE
  # each observation is RTN, so can just use TNE result!!
  
  Efish.list = lapply( X = as.list(.sei),
                       FUN = function(.s) {
                         E_fisher_TNE( .mu = .Mu,
                                       .sigma = sqrt(.Tt^2 + .s^2), 
                                       .n = 1,
                                       .a = -99,
                                       .b = .tcrit*.s )
                       })
  
  # add all the matrices entrywise
  # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
  Efish.all = Reduce('+', Efish.list) 
  
  # cat("\nFirst observation Efish:")
  # print(Efish.list[[1]])
  
  return(Efish.all)
}

# ### Sanity check: Should agree with E_fisher_TNE if all SEs equal
# n = 10
# se = 2
# Mu = 0.1
# Tt = 2
# zcrit = 1.96
# ( Efish1 = E_fisher_RTMA( .sei = rep(se, n),
#                           .Mu = Mu,
#                           .Tt = Tt,
#                           .tcrit = rep(zcrit, n) ) )
# 
# ( Efish2 = E_fisher_TNE( .mu = Mu,
#                          .sigma = sqrt(Tt^2 + se^2),
#                          .n = n,
#                          .a = -99,
#                          .b = zcrit*se ) )
# expect_equal(Efish1, Efish2)
# 
# # c.f.: different SEs but with same mean across studies
# E_fisher_RTMA( .sei = runif(n = n, min = se - 1, max = se + 1),
#                .Mu = Mu,
#                .Tt = Tt,
#                .tcrit = rep(zcrit, n) )



lprior = function(.sei, .Mu, .Tt, .tcrit) {
  Efish = E_fisher_RTMA( .sei = .sei, .Mu = .Mu, .Tt = .Tt, .tcrit = .tcrit )
  log( sqrt( det(Efish) ) )
}


# with usePrior = FALSE, agrees with weightr per "Repurpose TNE code.R"
# .pars: (.Mu, .Tt) or (.Mu, .Tt2)
nlpost_jeffreys_RTMA = function( .pars,
                                 .par2is = "Tt",  # "Tt" or "T2t"
                                 .yi,
                                 .sei,
                                 .tcrit = qnorm(.975),
                                 
                                 # if .usePrior = FALSE, will just be the MLE
                                 .usePrior = TRUE) {
  
  
  if ( .pars[2] < 0 ) return(.Machine$integer.max)
  
  # variance parameterization
  if (.par2is == "T2t") {
    Mu = .pars[1]
    Tt = sqrt(.pars[2])
  }
  
  
  if (.par2is == "Tt") {
    Mu = .pars[1]
    Tt = .pars[2]
  }
  
  
  
  # negative log-likelihood
  # joint_nll_2 uses the TNE version
  nll.value = joint_nll_2( .yi = .yi,
                           .sei = .sei,
                           .Mu = Mu,
                           .Tt = Tt,
                           .tcrit = .tcrit )
  
  # log-prior
  # lprior uses the TNE expected Fisher and then just sums over observations
  if ( .usePrior == TRUE ) {
    prior.value = lprior(.sei = .sei,
                         .Mu = Mu,
                         .Tt = Tt,
                         .tcrit = .tcrit)
  } else {
    prior.value = 0
  }
  
  # negative log-posterior
  nlp.value = sum(nll.value) + prior.value
  
  if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  
  return(nlp.value)
}

# ESTIMATION METHOD FNS ----------------------------------------------

# ~~ Estimation Methods That ARE Structured for Use Inside run_method_safe --------


# Because these fns are run inside run_method_safe, the latter will handle editing rep.res
#  All of these fns should take get.CIs as an argument and return CIs as c(NA, NA) if not wanted
# Fns in this category need to return a dataframe with the below structure, although it's okay if they don't return all of these names since run_method_safe will handle that. Note that this is a LIST with a dataframe called "stats", not just the dataframe; this allows easy extension in case you want to return other objects, like model objects.

# return( list( stats = data.frame( method = "jeffreys-mcmc",
#                                   
#                                   Mhat = Mhat,
#                                   Shat = Shat,
#                                   
#                                   MhatSE = MhatSE,
#                                   ShatSE = ShatSE,
#                                   
#                                   # this will use same CI limits for all 3 pt estimates
#                                   MLo = M.CI[1],
#                                   MHi = M.CI[2],
#                                   
#                                   SLo = S.CI[1],
#                                   SHi = S.CI[2],
#                                   
#                                   stan.warned = stan.warned,
#                                   stan.warning = stan.warning,
#                                   MhatRhat = postSumm["mu", "Rhat"],
#                                   ShatRhat = postSumm["tau", "Rhat"],
#                                   
#                                   overall.error = error ) ) )


# 2021-9-2: MM audited fn by reading
# mu.start, sigma.start: start values for optimatization
# as illustrated in a sanity check after nlpost_simple, this fn's MAPs agree with
#  using mle() directly on nlpost_Jeffreys
# confirmed that this agrees with weightr when usePrior = FALSE; see "2021-11-23 repurpose TNE code"
estimate_jeffreys_RTMA = function( yi, 
                                   sei,
                                   par2is = "Tt",
                                   Mu.start,
                                   par2.start,
                                   tcrit,  # SCALAR?
                                   
                                   usePrior = TRUE,
                                   get.CIs,
                                   CI.method = "wald" ) {
  
  # #bm
  # #@TEST ONLY
  # dpn = dp[ dp$affirm == FALSE, ]
  # yi = dpn$yi
  # sei = sqrt(dpn$vi)
  # par2is = "Tt"
  # Mu.start = p$Mu
  # par2.start = sqrt(p$t2a + p$t2w)
  # tcrit = qnorm(.975)
  # usePrior = FALSE
  # get.CIs = TRUE
  # CI.method = "wald"
  
  
  ### Get MAP by Calling mle() ###
  # IMPORTANT: This fn cannot be moved outside the scope of estimate_jeffreys
  #  because mle() is too dumb to allow extra args (e.g., x) to be passed,
  #  so it's forced to rely on global vars
  #  and that's a problem with a doParallel loop
  #  if this fn is outside estimate_jeffreys, different parallel iterations will use each other's global vars
  
  ##### Get MLE with Main Optimizer #####
  #  expects yi, sei, and tcrit to be global vars (wrt this inner fn)
  # second parameter could be Tt or T2t depending on par2is argument above
  nlpost_simple_RTMA = function(.Mu, .par2) {
    
    nlpost_jeffreys_RTMA( .pars = c(.Mu, .par2),
                          .par2is = par2is,
                          .yi = yi,
                          .sei = sei,
                          .tcrit = tcrit,
                          .usePrior = usePrior )
  }
  
  if ( usePrior == FALSE ) main.optimizer = "BFGS" else main.optimizer = "Nelder-Mead"
  
  #**important: force use of Nelder-Mead optimization, which works better for Jeffreys
  #  (even though BFGS works better for MLE)
  # for more on this issue, see "2021-9-23 SD vs. var redux with Jeffreys.R"
  mle.obj = mle( minuslogl = nlpost_simple_RTMA,
                 start = list( .Mu = Mu.start, .par2 = par2.start),
                 method = main.optimizer )
  
  
  # not actually MLEs, of course, but rather MAPs
  mles = as.numeric(coef(mle.obj))
  
  # here could pull in TNE code to run optimx
  # see call to get_optimx_dataframe() inside estimate_mle() 
  
  # THIS BEHAVES WELL
  if ( par2is == "Tt" ) {
    # need this structure for run_method_safe to understand
    MuHat = mles[1]
    TtHat = mles[2]
  }
  
  # from TNE
  # THIS BEHAVES BADLY
  if ( par2is == "T2t" ) {
    # need this structure for run_method_safe to understand
    MuHat = mles[1]
    TtHat = sqrt(mles[2])
  }
  
  # recode convergence more intuitively
  # optim uses "0" to mean successful convergence
  optim.converged = attr(mle.obj, "details")$convergence == 0 
  
  ##### Inference #####
  # in case someone passes a set of params that aren't handled
  SEs = los = his = c(NA, NA)
  
  #profile.CI.error = NA
  
  if ( get.CIs == TRUE & CI.method == "wald" ) {
    # get Wald CI 
    # SEs for both parameters
    # these are from the observed Fisher info,
    #  as confirmed in "2021-8-9 Investigate inference.R"
    SEs = as.numeric( attr( summary(mle.obj), "coef" )[, "Std. Error" ] )
    
    # if needed, overwrite SE for second parameter so that it's 
    #  on Tt scale rather than T2t scale
    # delta method:
    # let g(y) = y^(1/2), where y=Tt
    #  so g'(y) = 0.5 * y^(-0.5)
    if ( par2is == "var" ) {
      TtHatSE = SEs[2] * 0.5*TtHat^(-0.5)
      SEs[2] = TtHatSE
    }
    
    los = mles - SEs * qnorm(0.975)
    his = mles + SEs * qnorm(0.975)
    
  } else if ( get.CIs == TRUE & CI.method != "wald" ) {
    warning("That CI.method isn't handled yet")
  }
  
  return( list( stats = data.frame( Mhat = MuHat, 
                                    Shat = TtHat,
                                    
                                    MhatSE = SEs[1],
                                    ShatSE = SEs[2],
                                    
                                    MLo = as.numeric(los[1]),
                                    MHi = as.numeric(his[1]),
                                    
                                    SLo = as.numeric(los[2]),
                                    SHi = as.numeric(his[2]),
                                    
                                    optim.converged = optim.converged ) ) )
}



# modified from TNE (2022-2-21)

# Fns in this category need to return a list with these elements:
#   Mhat, Vhat, M.SE, M.CI (2-vector), V.SE, V.CI (2-vector)
# Because these fns are run inside run_method_safe, the latter will handle editing rep.res
#  All of these fns should take get.CIs as an argument and return CIs as c(NA, NA) if not wanted

estimate_jeffreys_mcmc_RTMA = function(.yi,
                                       .sei,
                                       .tcrit, #VECTOR
                                       .Mu.start,
                                       .Tt.start,
                                       .stan.adapt_delta = 0.8,
                                       .stan.maxtreedepth = 10 ) {
  
  # LL and UU: cutpoints on RAW scale, not Z-scores
  # tau: SD, not variance
  # to do:
  #  - can get rid of mustarL and LL; they are just placeholders now
  # this model.text is from "2022-2-22 stan for Sherlock"
  model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	
	  // these will be overwritten for EACH observation
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;
    real sigma;
		real LL; 
		real UU;
		// will just be set to 1
		int n; 
    

		// this will be the TOTALS for all observations
		matrix[2,2] fishinfototal;
		fishinfototal[1,1] = 0;
  	fishinfototal[1,2] = 0;
  	fishinfototal[2,1] = 0;
  	fishinfototal[2,2] = 0;
		
		
		// MM: build a Fisher info matrix for EACH observation
		for (i in 1:k) {
		
		  // marginal SD for this one observation
		  sigma = sqrt(tau^2 + sei[i]^2);
		  
		  // upper truncation limit (i.e., affirmative threshold)
		  //   for THIS study, given its SE
		  UU = tcrit[i] * sei[i];
		
		  // standardized truncation limits
  		mustarL = -999;
  		mustarU = (UU - mu) / sigma;
  		
  		// because EACH fisher info below has n=1 only
  		n = 1; 
  		
  		// beginning of stuff that is not modified at all from TNE
  		// note that normal_lpdf, etc., parameterize in terms of SD, not var
  		//  the (0,1) below are *not* start values for MCMC
  		alphaL = exp( normal_lpdf(mustarL | 0, 1) - 
  	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
  	                normal_lcdf(mustarL | 0, 1) ) ); 
  	                
  		alphaU = exp( normal_lpdf(mustarU | 0, 1) - 
   	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
   	                normal_lcdf(mustarL | 0, 1) ) );
   	                
  		// second derivatives for Fisher info			
  		kmm = -n/sigma^2 + n/sigma^2 * ((alphaU-alphaL)^2 + alphaU*mustarU- alphaL*mustarL);
  		kms = -2*n/sigma^2 * (alphaL - alphaU) + 
  	   		  n/sigma^2 * (alphaL - alphaU + (alphaU*mustarU^2 - alphaL*mustarL^2) +
  			  				(alphaL-alphaU) * (alphaL*mustarL - alphaU*mustarU));
  		kss = n/sigma^2 - 3*n/sigma^2 * (1 + mustarL*alphaL - mustarU*alphaU) +
  	   			n/sigma^2 * (mustarU*alphaU*(mustarU^2 - 2) - mustarL*alphaL*(mustarL^2 - 2) +
  								(alphaU*mustarU - alphaL*mustarL)^2);
  		
  		fishinfo[1,1] = -kmm;
      fishinfo[1,2] = -kms;
      fishinfo[2,1] = -kms;
      fishinfo[2,2] = -kss;
  		// end stuff that is not modified at all from TNE
  		
  		// MM: add the new fisher info to the total one
  		fishinfototal = fishinfototal + fishinfo;
		}
		return sqrt(determinant(fishinfototal));
	}
}

data{
	int<lower=0> k;
  real sei[k];
  real tcrit[k];
	real y[k];
}

parameters{
  real mu;
	real<lower=0> tau;
}


model{
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit) );
	for(i in 1:k)
	      y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
}

// this chunk doesn't actually affect the model that's being fit to the data
//  (I think); it's just re-calculating the prior, lkl, and post to return to user
generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit) );
  real log_post;
  // this is just an intermediate quantity for log_lik
  real UU;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
      
      UU = tcrit[i] * sei[i];
      
      // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
      // see Truncation with upper bounds in Stan
	      if ( y[i] > UU )
          log_lik += negative_infinity();
        else
          log_lik += -1 * normal_lcdf(UU | mu, sqrt(tau^2 + sei[i]^2) ); 
  }
  // easy way to see what is essentially MLE: just remove prior from this
  log_post = log_prior + log_lik;
}
"


stan.model <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")



# #@TEST ONLY
# .yi = dn$yi
# .sei = dn$sei
# .tcrit = dn$tcrit
# .Mu.start = Mu
# .Tt.start = sqrt(T2)
# p = data.frame(stan.maxtreedepth = 10,
#                stan.iter = 2000,
#                  stan.adapt_delta = 0.8)

# handle scalar tcrit
if ( length(.tcrit) < length(.yi) ) .tcrit = rep( .tcrit, length(.yi) )

# prepare to capture warnings from Stan
stan.warned = 0
stan.warning = NA

# set start values for sampler
init.fcn <- function(o){ list(mu = .Mu.start,
                              tau = .Tt.start ) }

# like tryCatch, but captures warnings without stopping the function from
#  returning its results
withCallingHandlers({
  
  cat( paste("\n estimate_jeffreys_mcmc flag 1: about to call stan_model") )
  
  # necessary to prevent ReadRDS errors in which cores try to work with other cores' intermediate results
  # https://groups.google.com/g/stan-users/c/8snqQTTfWVs?pli=1
  options(mc.cores = parallel::detectCores())
  
  # "isystem" arg is just a placeholder to avoid Stan's not understanding special characters
  #  in getwd(), even though we don't actually use the dir at all
  # note: removing the isystem arg does NOT fix the very sporadic "error reading from connection" on cluster
  stan.model <- stan_model(model_code = model.text,
                           isystem = "~/Desktop")
  
  
  cat( paste("\n estimate_jeffreys_mcmc flag 2: about to call sampling") )
  
  post = sampling(stan.model,
                  cores = 1,
                  refresh = 0,
                  data = list( k = length(.yi),
                               sei = .sei,
                               tcrit = .tcrit,
                               y = .yi ),
                  
                  #iter = p$stan.iter,   
                  control = list(max_treedepth = p$stan.maxtreedepth,
                                 adapt_delta = p$stan.adapt_delta),
                  
                  init = init.fcn)
  
  
}, warning = function(condition){
  stan.warned <<- 1
  stan.warning <<- condition$message
} )

cat( paste("\n estimate_jeffreys_mcmc flag 3: about to call postSumm") )
postSumm = summary(post)$summary


# pull out best iterate to pass to MAP optimization later
ext = extract(post) # a vector of all post-WU iterates across all chains
best.ind = which.max(ext$lp__)  # single iterate with best log-posterior should be very close to MAP


# posterior means, posterior medians, and max-LP iterate
Mhat = c( postSumm["mu", "mean"],
          median( rstan::extract(post, "mu")[[1]] ),
          ext$mu[best.ind] )

Shat = c( postSumm["tau", "mean"],
          median( rstan::extract(post, "tau")[[1]] ),
          ext$tau[best.ind] )

# sanity check
expect_equal( Mhat[1], mean( rstan::extract(post, "mu")[[1]] ) )


# SEs
MhatSE = postSumm["mu", "se_mean"]
ShatSE = postSumm["tau", "se_mean"]
# how Stan estimates the SE: https://discourse.mc-stan.org/t/se-mean-in-print-stanfit/2869
expect_equal( postSumm["mu", "sd"],
              sd( rstan::extract(post, "mu")[[1]] ) )
expect_equal( MhatSE,
              postSumm["mu", "sd"] / sqrt( postSumm["mu", "n_eff"] ) )

# CI limits
S.CI = c( postSumm["tau", "2.5%"], postSumm["tau", "97.5%"] )
V.CI = S.CI^2
M.CI = c( postSumm["mu", "2.5%"], postSumm["mu", "97.5%"] )
# sanity check:
myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], 0.025 ),
                          quantile( rstan::extract(post, "mu")[[1]], 0.975 ) ) )
expect_equal(M.CI, myMhatCI)



# n.ests = length(Mhat)
# # OLD return structure (list)
# return( list( Mhat = Mhat,
#               Shat = Shat,
#               
#               MhatSE = rep(MhatSE, length(n.ests)),
#               ShatSE = rep(ShatSE, length(n.ests)),
#               
#               M.CI = M.CI,
#               S.CI = S.CI,
#               
#               stan.warned = stan.warned,
#               stan.warning = stan.warning,
#               
#               MhatRhat = postSumm["mu", "Rhat"],
#               ShatRhat = postSumm["tau", "Rhat"]
# ) )

# the point estimates are length 2 (post means, then medians),
#  but the inference is the same for each type of point estimate
return( list( stats = data.frame( 
  
  Mhat = Mhat,
  Shat = Shat,
  
  MhatSE = MhatSE,
  ShatSE = ShatSE,
  
  # this will use same CI limits for all 3 pt estimates
  MLo = M.CI[1],
  MHi = M.CI[2],
  
  SLo = S.CI[1],
  SHi = S.CI[2],
  
  stan.warned = stan.warned,
  stan.warning = stan.warning,
  MhatRhat = postSumm["mu", "Rhat"],
  ShatRhat = postSumm["tau", "Rhat"] ) ) )

}


# 2021-9-2: MM audited fn by reading through
# x: data vector
# doesn't handle the case CI.method = "profile" and par2is = "sd" (will just return NAs)
estimate_mle = function( x,
                         p,  # scenario parameters
                         par2is = "sd",
                         mu.start = 0,
                         sigma.start = 1,
                         get.CIs,  # this is kept separate from p for use inside boot fn (where we'll want to suppress getting CIs)
                         CI.method = "wald"
) {
  
  #@TEMP
  #write.csv( c(mu.start, sigma.start), "inside_mle_start_values.csv")
  
  ##### Get MLE with Main Optimizer (NM) #####
  # fn needs to be formatted exactly like this (no additional args)
  #  in order for mle() to understand
  nll_simple = function(.mu, .sigma) {
    nll(.pars = c(.mu, .sigma),
        par2is = par2is,
        .x = x, .a = p$a, .b = p$b)
  }
  
  myMLE = mle( minuslogl = nll_simple,
               method = "BFGS",
               start = list( .mu=mu.start, .sigma=sigma.start) )
  
  mles = as.numeric( coef(myMLE) )
  
  
  # this parameterization behaves badly!
  if ( par2is == "sd" ) {
    # need this structure for run_method_safe to understand
    Mhat = mles[1]
    Vhat = mles[2]^2
    Shat = mles[2]
  }
  
  # this parameterization behaves well
  if ( par2is == "var" ) {
    # need this structure for run_method_safe to understand
    Mhat = mles[1]
    Vhat = mles[2]
    Shat = sqrt(mles[2])
  }
  
  mles = c(Mhat, Vhat, Shat)
  
  # recode convergence more intuitively
  # optim uses "0" to mean successful convergence
  optim.converged = attr(myMLE, "details")$convergence == 0
  
  
  ##### Try Other Optimizers #####
  w = get_optimx_dataframe(.method = "mle",
                           x = x,
                           p = p,  # scenario parameters
                           par2is = par2is,
                           mu.start = mu.start,
                           sigma.start = sigma.start)
  
  
  ##### Inference #####
  # in case someone passes a set of params that aren't handled
  SEs = los = his = c(NA, NA)
  
  profile.CI.error = NA
  
  if ( get.CIs == TRUE & CI.method == "wald" & par2is == "sd" ) {
    # get Wald CI 
    # SEs for both parameters
    # these are from the observed Fisher info,
    #  as confirmed in "2021-8-9 Investigate inference.R"
    SEs = as.numeric( attr( summary(myMLE), "coef" )[, "Std. Error" ] )
    
    # fill in VhatSE using delta method
    # let g(y) = y^2, where y=Shat
    VhatSE = SEs[2] * 2 * Shat
    
    # include Vhat
    SEs = c(SEs[1], VhatSE, SEs[2])
    
    los = mles - SEs * qnorm(0.975)
    his = mles + SEs * qnorm(0.975)
    
  } else if ( get.CIs == TRUE & CI.method == "wald" & par2is == "var" ) {
    # get Wald CI 
    # SEs for both parameters
    # these are from the observed Fisher info,
    #  as confirmed in "2021-8-9 Investigate inference.R"
    SEs = as.numeric( attr( summary(myMLE), "coef" )[, "Std. Error" ] )
    
    # fill in ShatSE using delta method
    # let g(y) = y^(1/2), where y=Shat
    #  so g'(y) = 0.5 * y^(-0.5)
    ShatSE = SEs[2] * 0.5*Vhat^(-0.5)
    
    # include Vhat
    SEs = c(SEs[1], SEs[2], ShatSE)
    
    los = mles - SEs * qnorm(0.975)
    his = mles + SEs * qnorm(0.975)
    
  } else if ( get.CIs == TRUE & CI.method == "profile" & par2is == "var" ) {
    
    tryCatch({
      CIs = confint(myMLE)
      
      los = c( as.numeric( CIs[1,1] ),
               as.numeric( CIs[2,1] ),
               sqrt( as.numeric( CIs[2,1] ) ) )
      his = c( as.numeric( CIs[1,2] ), 
               as.numeric( CIs[2,2] ),
               sqrt( as.numeric( CIs[2,2] ) ) )
      
      ( res.SEs = as.numeric( attr( summary(myMLE), "coef" )[, "Std. Error" ] ) )
      
      # fill in ShatSE using delta method
      # let g(y) = y^(1/2), where y=Shat
      #  so g'(y) = 0.5 * y^(-0.5)
      SEs = c(res.SEs[1], res.SEs[2], res.SEs[2] * 0.5*Vhat^(-0.5) )
      
    }, error = function(err) {
      SEs <<- los <<- his <<- rep(NA, 3)
      profile.CI.error <<- err$message
    })
    
  }   
  
  ##### Organize Results #####
  return( list( Mhat = Mhat, 
                Vhat = Vhat,
                Shat = Shat,
                
                MhatSE = SEs[1],
                VhatSE = SEs[2],
                ShatSE = SEs[3],
                
                M.CI = as.numeric( c(los[1], his[1]) ),
                V.CI = as.numeric( c(los[2], his[2]) ),
                S.CI = as.numeric( c(los[3], his[3]) ),
                optim.converged = optim.converged,
                profile.CI.error = profile.CI.error,
                
                optimx.dataframe = w
  ) )
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



# older version (before TNE)
# use only the nonaffirms to get trunc MLE
#  and throws away the affirms
# **note that the returned Vhat is an estimate of T2 + t2w, not T2 itself
correct_meta_phack1 = function( .dp,  # published studies
                                .p  # parameters as dataframe
) { 
  
  
  # published nonaffirmatives only
  dpn = .dp[ .dp$affirm == FALSE, ]
  
  
  #@ assumes they all have same se:
  crit = unique(dpn$tcrit)
  
  
  ### MLEs from trunc normal ###
  # these are the MLEs of the *t-stats*
  #@ IMPORTANT: for convenience, this is using the normal distribution, 
  #  so won't work well for small m
  mle.fit = mle.tmvnorm( X = as.matrix(dpn$tstat, ncol = 1),
                         lower = -Inf,
                         upper = crit)
  #summary(mle.fit)
  mles = coef(mle.fit)
  
  #@can't get CIs yet because of package issue
  #  emailed Wilhelm about this on 2021-5-7
  # CIs.tstats = confint( profile(mles, dpn$tstat, method="BFGS", trace=TRUE) )
  
  # get Wald CI a different way
  tstat.mu.SE = attr( summary(mle.fit), "coef" )[ "mu_1", "Std. Error" ]
  tstat.mu.CI = c( mles[1] - tstat.mu.SE * qnorm(0.975),
                   mles[1] + tstat.mu.SE * qnorm(0.975) )
  
  # # pretty good :)
  # mles[1]; p$Mu/p$se
  # mles[2]; (1/p$se^2) * (p$T2 + p$t2w + p$se^2)
  
  # APPROXIMATELY rescale MLEs to represent effect sizes rather than tstats
  Mhat = mles[1] * .p$se
  # **use Vhat to represent MARGINAL heterogeneity (i.e., T2 + t2w)
  Vhat = ( mles[2] * .p$se^2 ) - .p$se^2
  
  # and rescale CI limits
  MhatCI = tstat.mu.CI * .p$se
  
  ### Sanity checks: Moments of published nonaffirms vs. theory ###
  # check that moments are what we expect
  # without delta method:
  
  theoryExpTstat = extrunc(spec = "norm",
                           mean =.p$Mu /.p$se,
                           #@doesn't use the delta-method thing
                           sd = sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ),
                           b = crit )
  
  theoryVarTstat = vartrunc(spec = "norm",
                            mean =.p$Mu /.p$se,
                            #@doesn't use the delta-method thing
                            sd = sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ),
                            b = crit )
  
  # delta-method version (not checked and seems not to work):
  # library(msm)
  # # https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
  # sd.y = .p$se * sqrt(.p$m)
  # viSE = sqrt( 2 * sd.y^4 / (.p$m-1) )
  # correctedSE = deltamethod( g = ~ x1/sqrt(x2),  # the t-stat
  #                            mean = c(.p$Mu, .p$se^2),
  #                            cov = matrix( c( .p$T2 + .p$t2w + .p$se^2, 0, 0, viSE^2 ),
  #                                          nrow = 2 ) )
  # 
  # extrunc(spec = "norm",
  #         mean =.p$Mu /.p$se,
  #         sd = correctedSE,
  #         b = crit )
  
  
  ### Return all the things ###
  return( list( metaCorr = data.frame( MhatCorr = Mhat,
                                       MhatLoCorr = MhatCI[1],
                                       MhatHiCorr = MhatCI[2],
                                       MhatCoverCorr = ( MhatCI[1] <=.p$Mu ) & ( MhatCI[2] >=.p$Mu ),
                                       VhatCorr = Vhat),
                
                
                # **note that all of these stats pertain to only published nonaffirmatives
                sanityChecks = data.frame( kNonaffirmPub = nrow(dpn),
                                           kNonaffirmPubUnhacked = sum( dpn$hack == "no" ),
                                           kNonaffirmPubHacked = sum( dpn$hack ==.p$hack ),
                                           
                                           ### check that moments of published nonaffirms are what we expect
                                           TheoryExpTstat = theoryExpTstat,
                                           # below should match TheoryExpTstat IF all nonaffirms are from unhacked
                                           #  study sets, but otherwise may not match:
                                           MeanTstat = mean(dpn$tstat), 
                                           # mean of t-stats from published nonaffirms from
                                           #  unhacked study sets (should match TheoryExpTstat)
                                           MeanTstatUnhacked = mean( dpn$tstat[dpn$hack == "no" ] ),
                                           # mean of t-stats from published nonaffirms from
                                           #  hacked study sets (may not match TheoryExpTstat)
                                           MeanTstatHacked = mean( dpn$tstat[dpn$hack ==.p$hack ] ),
                                           
                                           TheoryVarTstat = theoryVarTstat,
                                           EstVarTstat = var(dpn$tstat),
                                           # should match theory:
                                           EstVarTstatUnhacked = var( dpn$tstat[dpn$hack == "no" ] ),
                                           EstVarTstatHacked = var( dpn$tstat[dpn$hack ==.p$hack ] ),
                                           
                                           # MLEs of the t-stats themselves (before rescaling using the SE)
                                           # *marginal* t-stats (underlying distribution rather than truncated one)
                                           TheoryExpTstatMarg = .p$Mu/.p$se,
                                           tstatMeanMLE = mles[1],
                                           tstatMeanMLELo = tstat.mu.CI[1],
                                           tstatMeanMLEHi = tstat.mu.CI[2],
                                           tstatMeanMLECover = (tstat.mu.CI[1] <= .p$Mu/.p$se) & (tstat.mu.CI[2] >= .p$Mu/.p$se),
                                           
                                           # other stats
                                           Mean.yi = mean(dpn$yi),
                                           Mean.yi.Unhacked = mean( dpn$yi[dpn$hack == "no" ] ),
                                           Mean.yi.Hacked = mean( dpn$yi[dpn$hack ==.p$hack ] ),
                                           
                                           Mean.vi = mean(dpn$vi),
                                           Mean.vi.Unhacked = mean( dpn$vi[dpn$hack == "no" ] ),
                                           Mean.vi.Hacked = mean( dpn$vi[dpn$hack ==.p$hack ] ),
                                           
                                           # stats about all published studies
                                           dp.k = nrow(.dp),
                                           dp.kAffirm = sum(.dp$affirm == TRUE),
                                           dp.kNonaffirm = sum(.dp$affirm == FALSE),
                                           dp.Nrealized = mean(.dp$N)
                )
  ) )
  
}







# Notes from TNE:
# In order to catch errors from individual estimation methods safely and informatively,
#  in general the estimation method fns are structured st they can be run within the
#  fn run_method_safe, which automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error.

# ~~ Wrapper Fn to Safely Run a Method -------

# See note at the beginning of this script
#  this fn automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error

# Important: this fn works if method.fn() returns multiple rows
# BUT in that case, it assumes that the CIs are shared for all rows of that method

# expects global vars: all.errors, rep.res
# directly edits res via superassignment
run_method_safe = function( method.label,
                            method.fn,
                            .rep.res ) {
  
  cat( paste("\n run_method_safe flag 1: about to try running method", method.label) )
  
  
  tryCatch({
    
    new.rows = method.fn()$stats
    
    cat( paste("\n run_method_safe flag 2: done calling method.fn() for", method.label) )
    
    error = NA
    
  }, error = function(err) {
    # needs to be superassignment because inside the "error" fn
    error <<- err$message
    
    # only need one variable in the blank dataframe since bind_rows below
    #  will fill in the rest
    new.rows = data.frame( method = method.label )
    
  })
  
  new.rows = new.rows %>% add_column( method = method.label, .before = 1 )
  new.rows$overall.error = error
  
  # optimx.dataframe is itself a df, so needs to be handled differently
  # if ( !is.null(optimx.dataframe) ) new.row = bind_cols(new.row, optimx.dataframe)
  
  if ( nrow(.rep.res) == 0 ) .rep.res = new.rows else .rep.res = bind_rows(.rep.res, new.rows)
  return(.rep.res) 
  
}

# example of how to call it when method.fn takes args
# all.errors = c()
# if  exists("rep.res") ) r("rep.re("rep.re
# run_method_safe( method = "mle",
#                  method.fn = function() estimate_mles(x = x, get.CIs = TRUE ) )

# #### Sanity checks ####
# # fake method for estimating the moments, but it breaks if x<0 or x>5
# crappy_method = function(x) {
#   if ( x > 0 & x < 5 ) return( list(Mhat = x+1,
#                                     Vhat = x-1,
#                                     M.CI = c(NA, NA),
#                                     V.CI = c(NA, NA) ) )
#   if ( x <= 0 ) stop("Fake error A generated by method.fn!")
#   if ( x >= 5 ) stop("Fake error B generated by method.fn!")
# }
# 
# all.errors = c()
# if( exists("rep.res") ) rm(rep.res)
# 
# run_method_safe( "mle",
#                  method.fn = function() { crappy_method(-1) } )
# 
# # no error on this one
# run_method_safe( "mle",
#                  method.fn = function() { crappy_method(4) } )
# 
# 
# # this one will have a different error
# # no error on this one
# run_method_safe( "mle", method.fn = function() { crappy_method(40) } )
# 
# expect_equal( all.errors, c( "mle: Fake error A generated by method.fn!",
#                             "mle: Fake error B generated by method.fn!" ) )
# 
# expect_equal( rep.res,
#               data.frame( method = rep("mle", 3),
#                           Mhat = c(NA, 5, NA),
#                           Vhat = c(NA, 3, NA),
#                           MLo = rep(NA, 3),
#                           MHi = rep(NA, 3),
#                           VLo = rep(NA, 3),
#                           VHi = rep(NA, 3) ) )
# #### end sanity checks



# 2022-2-24: NO LONGER IN USE; NOW USING RUN_METHOD_SAFE ABOVE
# use only the nonaffirms to get trunc MLE
#  and throws away the affirms
# **note that the returned Vhat is an estimate of T2 + t2w, not T2 itself
correct_meta_phack3 = function( .dp,  # published studies
                                .p,   # parameters as dataframe
                                .method, # see options in estimate_jeffreys_RTMA
                                .Mu.start,
                                .par2.start
                                
) { 
  
  
  # published nonaffirmatives only
  dpn = .dp[ .dp$affirm == FALSE, ]
  
  if ( .method %in% c("mle-sd", "mle-var") ) {
    usePrior = FALSE
  } else if ( .method %in% c( "jeffreys-mode-sd", "jeffreys-mode-var" ) ) {
    usePrior = TRUE
  } else {
    stop("method not handled")
  }
  
  if ( .method %in% c("mle-sd", "jeffreys-mode-sd") ) par2is = "Tt"
  if ( .method %in% c("mle-var", "jeffreys-mode-var") ) par2is = "T2t"
  
  res = estimate_jeffreys_RTMA( yi = dpn$yi,
                                sei = sqrt(dpn$vi),
                                par2is = par2is,
                                Mu.start = .Mu.start,
                                par2.start = .par2.start,  
                                tcrit = dpn$tcrit,
                                
                                usePrior = usePrior,
                                get.CIs = p$get.CIs,
                                CI.method = "wald" )
  
  
  ### Return all the things ###
  return( list(metaCorr = data.frame( Mhat = res$MuHat,
                                      MhatLo = res$Mu.CI[1],
                                      MhatHi = res$Mu.CI[2],
                                      MhatCover = ( res$Mu.CI[1] <=.p$Mu ) & ( res$Mu.CI[2] >=.p$Mu ),
                                      
                                      Shat = res$TtHat,
                                      ShatLo = res$Tt.CI[1],
                                      ShatHi = res$Tt.CI[2],
                                      ShatCover = ( res$Tt.CI[1] <= .p$S ) & ( res$Tt.CI[2] >= .p$S ) )
  ))
  
  
}


# nicely report a metafor object with optional suffix to denote which model it is
# includes coverage
report_rma = function(.mod,
                      .Mu,  # true mean (to get coverage)
                      .suffix = "") {
  
  if ( !is.null(.mod) ) {
    .res = data.frame( .mod$b,
                       .mod$ci.lb,
                       .mod$ci.ub,
                       (.mod$ci.lb <= .Mu) & (.mod$ci.ub >= .Mu), 
                       .mod$pval )
  } else {
    .res = data.frame( NA,
                       NA,
                       NA,
                       NA, 
                       NA )
  }
  
  
  names(.res) = paste( c("Mhat", "MhatLo", "MhatHi", "MhatCover", "MhatPval"), .suffix, sep = "" )
  return(.res)
}


# ANALYSIS FNS: APPLIED EXAMPLES -----------------------------------------------

# just like correct_meta_phack1, but designed for applied examples instead of sim study
# use only the nonaffirms to get trunc MLE
#  and throws away the affirms
# **note that the returned Vhat is an estimate of T2 + t2w, not T2 itself
# crit: value at which to truncate the dataset; using a non-default value is only for sanity checks
correct_meta_phack2 = function( yi,
                                vi,
                                crit = qnorm(.975) ) { 
  
  
  d = data.frame(yi = yi, vi = vi)
  d$tstat = d$yi / sqrt(d$vi)
  d$affirm = d$tstat > crit
  
  
  # published affirmatives only
  dpn = d[ d$affirm == FALSE, ]
  
  ### MLEs from trunc normal ###
  # these are the MLEs of the *t-stats*
  #@ IMPORTANT: for convenience, this is using the normal distribution, 
  #  so won't work well for small m
  mle.fit = mle.tmvnorm( X = as.matrix(dpn$tstat, ncol = 1),
                         lower = -Inf,
                         upper = crit)
  mles = coef(mle.fit)
  
  # get Wald CI a different way
  tstat.mu.SE = attr( summary(mle.fit), "coef" )[ "mu_1", "Std. Error" ]
  tstat.mu.CI = c( mles[1] - tstat.mu.SE * qnorm(0.975),
                   mles[1] + tstat.mu.SE * qnorm(0.975) )
  
  
  # rescale MLEs to represent effect sizes rather than tstats
  # SEs could differ across studies, so use the tstat MLE to rescale for each study's vi
  #@ is that right?
  # this is using the approximation E[tstat] = E[yi/SE] \approx E[yi]/E[SE]
  #  for now because otherwise we'd have to treat the yi's as non-iid in the likelihood
  # the approximation holds if yi and SEi are independent
  #   (in the underlying population, I think)
  Mhat = mles[1] * mean( sqrt(dpn$vi) )
  # **use Vhat to represent MARGINAL heterogeneity (i.e., T2 + t2w)
  Vhat = mean( ( mles[2] * dpn$vi ) - dpn$vi )
  
  # and rescale CI limits
  #@think about this again in light of SEs differing across studies
  #MhatCI = tstat.mu.CI * .p$se
  MhatCI = c(NA, NA)
  
  ### Naive Meta-Analysis ###
  
  metaNaive = rma.uni( yi = yi,
                       vi = vi,
                       method = "REML", 
                       knha = TRUE )
  
  ### Return all the things ###
  return( list( metaCorr = data.frame( MhatCorr = Mhat,
                                       MhatLoCorr = MhatCI[1],
                                       MhatHiCorr = MhatCI[2],
                                       VhatCorr = Vhat),
                
                metaNaive = report_rma(.mod = metaNaive,
                                       .Mu = NA,
                                       .suffix = "Naive"),
                
                # **note that all of these stats pertain to only published nonaffirmatives
                sanityChecks = data.frame( kNonaffirmPub = nrow(dpn),
                                           
                                           # MLEs of the t-stats themselves (before rescaling using the SE)
                                           # marginal t-stats (underlying distribution rather than truncated one)
                                           
                                           
                                           tstatMeanMLE = mles[1],
                                           tstatMeanMLELo = tstat.mu.CI[1],
                                           tstatMeanMLEHi = tstat.mu.CI[2],
                                           
                                           tstatVarMLE = mles[2],
                                           
                                           
                                           # other stats
                                           Mean.yi = mean(dpn$yi),
                                           
                                           Mean.vi = mean(dpn$vi),
                                           
                                           # stats about all published studies
                                           dp.k = nrow(d),
                                           dp.kAffirm = sum(d$affirm == TRUE),
                                           dp.kNonaffirm = sum(d$affirm == FALSE)
                                           
                ),
                
                data = d,
                crit = crit
  ) )
  
}


# plot empirical data

# .obj: object returned by correct_meta_phack2
# showAffirms: should it show all studies, even affirms?
# black line = LOESS density of nonaffirms
# red line = MLE from RTMA (parametric counterpart to the above)
# blue line = LOESS density of all tstats (including affirms)
plot_trunc_densities = function(.obj,
                                showAffirms = FALSE) {
  
  # already has affirmative indicator
  d = .obj$data
  dn = d[d$affirm == FALSE,]
  
  tstatMeanMLE = .obj$sanityChecks$tstatMeanMLE
  tstatVarMLE = .obj$sanityChecks$tstatVarMLE
  
  xmin = floor(min(dn$tstat))
  xmax = ceiling(max(dn$tstat))
  
  p = ggplot(data = data.frame(x = c(xmin, 3)),
             aes(x)) +
    
    geom_vline(xintercept = 0,
               lwd = 1,
               color = "gray") +
    
    # geom_vline(xintercept = tstatMeanMLE,
    #            lty = 2,
    #            lwd = 1,
    #            color = "red") +
    
    
    # estimated density of estimates
    geom_density( data = dn,
                  aes(x = tstat),
                  adjust = .3 ) + 
    
    
    
    # estimated density from meta-analysis
    stat_function( fun = dtrunc,
                   n = 101,
                   args = list( spec = "norm",
                                mean = tstatMeanMLE,
                                sd = sqrt(tstatVarMLE),
                                b = .obj$crit),
                   #aes(y = .25 * ..count..),  # doesn't work
                   lwd = 1.2,
                   color = "red") +
    
    #bm: got histogram to work with density in terms of scaling, but not yet the stat_function
    
    
    ylab("") +
    #scale_x_continuous( breaks = seq(xmin, 3, 0.5)) +
    xlab("t-stat") +
    theme_minimal() +
    scale_y_continuous(breaks = NULL) +
    theme(text = element_text(size=16),
          axis.text.x = element_text(size=16))
  
  
  # also show density of all t-stats, not just the nonaffirms
  if ( showAffirms == TRUE ) {
    p = p + geom_density( data = d,
                          aes(x = tstat),
                          color = "blue",
                          adjust = .3 )
  }
  
  
  return(p)
  
}




# get MLE for ONE nonaffirmative study
# can either estimate the total within-study variance (t2w + SE^2)
#  by providing both .se and .t2w, or can allow it to be estimated by not providing 
# BUT in the latter case, could end up being smaller than SE^2,
#  which is illogical
one_nonaffirm_mle = function(.tstat,
                             .se = NA, # of yi, not tstat
                             .t2w = NA,
                             .crit
) {
  
  if (length(.tstat) > 1) stop("fn not intended for more than 1 observation")
  
  # if we're given the components of the WITHIN-study variance,
  #  treat as fixed rather than estimating it
  if ( !is.na(.se) & !is.na(.t2w) ) Vw = (1/.se^2) * (.t2w + .se^2)
  
  if ( exists("Vw") ) {
    mle.fit = mle.tmvnorm( X = as.matrix(.tstat, ncol = 1),
                           lower = -Inf,
                           upper = .crit,
                           # fixed parameters must be named exactly as
                           #  coef(mle.fit) returns them
                           fixed = list(sigma_1.1 = Vw) )
  } else{
    mle.fit = mle.tmvnorm( X = as.matrix(.tstat, ncol = 1),
                           lower = -Inf,
                           upper = .crit)
  }
  
  
  mles = coef(mle.fit)
  
  return( data.frame( tstat.mu.hat = as.numeric(mles[1]),
                      muiHat = as.numeric(mles[1]) * .se,
                      VwHat = as.numeric(mles[2]) ) )
  
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
  
  # version in correct_meta_phack1 for nonaffirms:
  # for large m, var(d$vi) above is tiny, so sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ) is 
  # basically the same as correctedSE above
  # theoryExpTstat = extrunc(spec = "norm",
  #                          mean =.p$Mu /.p$se,
  #                          #@doesn't use the delta-method thing
  #                          sd = sqrt( (1/.p$se^2) * (.p$T2 +.p$t2w +.p$se^2) ),
  #                          b = crit )
  
  
  res = data.frame( matrix( NA, nrow = 2, ncol = 2) )
  names(res) = c("affirms", "nonaffirms")
  row.names(res) = c("theoryExp", "empiricalExp")
  
  
  ## Expectation of affirmatives
  res[ "theoryExp", "affirms" ] = extrunc( spec = "norm",
                                           mean = Mu / se,
                                           sd = correctedSE,
                                           a = crit )
  
  res[ "empiricalExp", "affirms" ] = mean( dph$tstat[ dph$affirm == TRUE ] )
  
  
  ## Variance of affirmatives
  res[ "theoryVar", "affirms" ] = vartrunc( spec = "norm",
                                            mean = Mu / se,
                                            sd = correctedSE,
                                            a = crit )
  
  res[ "empiricalVar", "affirms" ] = var( dph$tstat[ dph$affirm == TRUE ] )
  
  
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



# changes from sim_meta:
# sei.expr
# k.pub.nonaffirm
sim_meta_2 = function(Nmax,  # max draws to try
                      Mu,  # overall mean for meta-analysis
                      t2a,  # across-study heterogeneity (NOT total heterogeneity)
                      
                      # study parameters, assumed same for all studies:
                      m,  # sample size for this study
                      t2w,  # within-study heterogeneity
                      true.sei.expr,  # TRUE SE string to evaluate
                      rho = 0,  # autocorrelation of muin's
                      
                      hack,  # mechanism of hacking for studies that DO hack (so not "no")
                      
                      k.pub.nonaffirm,  # number of published nonaffirmatives
                      prob.hacked,
                      
                      return.only.published = FALSE
) {
  
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k.pub.nonaffirm", "prob.hacked", "true.sei.expr")]
  
  
  if ( hack == "no" ) stop("hack should only be 'affirm' or 'signif' for this fn")
  
  k.pub.nonaffirm.achieved = 0
  i = 1
  
  while( k.pub.nonaffirm.achieved < k.pub.nonaffirm ) {
    
    is.hacked = rbinom(n = 1, size = 1, prob = prob.hacked)
    true.se = eval( parse( text = true.sei.expr) )
    
    # do we still need this??
    #if ( exists(".dat") ) startInd = max(.dat$study) + 1 else startInd = 1
    
    if ( is.hacked == 0 ) {
      # to generate unhacked studies, need to change argument "hack"
      .argsUnhacked = .args
      .argsUnhacked[ names(.args) == "hack" ] = "no"
      .argsUnhacked$se = true.se
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsUnhacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
      
    } else if ( is.hacked == 1 ) {
      
      # for unhacked studies, no need to change argument "hack"
      .argsHacked = .args
      .argsHacked$se = true.se
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsHacked )
      
      # add study ID
      newRows = newRows %>% add_column( .before = 1,
                                        study = i )
      
      if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    }
    
    i = i + 1
    k.pub.nonaffirm.achieved = sum( .dat$affirm == FALSE & .dat$Di == 1 ) 
    
  }
  
  # add more info to dataset
  .dat$k.underlying = length(unique(.dat$study))
  .dat$k.nonaffirm.underlying = length( unique( .dat$study[ .dat$affirm == FALSE ] ) )
  
  if ( return.only.published == TRUE ) .dat = .dat[ .dat$Di == 1, ]
  
  return(.dat)
}

# d = sim_meta_2(  # test only
#   Nmax = 20,
#   Mu = 1,
#   t2a = 0.1,
#   m = 50,
#   t2w = .5,
#   true.sei.expr = "runif( n = 1, min = 0.5, max = 2 )",
#   hack = "affirm",
#   return.only.published = FALSE,
#   rho=0,
#   k.pub.nonaffirm = 30,
#   prob.hacked = 0.4 )
# 
# d$k.nonaffirm.underlying[1]
# d$k.underlying[1]
# table(d$Di, d$affirm)
# #@KEEP/PULL IN THE SANITY CHECKS FROM THE VERSION BELOW


# KEEP THIS VERSION FOR REVERSE-COMPATIBILITY AND SANITY CHECKS
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
                             t2a,  # across-study heterogeneity (NOT total heterogeneity)
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
  # t2a = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "affirm"
  # rho=0
  
  # mean for this study set
  # doesn't have t2w because that applies to results within this study set
  mui = Mu + rnorm(mean = 0,
                   sd = sqrt(t2a),
                   n = 1)
  
  # TRUE SD (not estimated)
  #@SHOULD BE M-1, I THINK
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
#                       t2a = 0.1,
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
#                          t2a = 0.1,
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
#                       t2a = 0.1,
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
#                         t2a = 0.1,
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
# t2a = 0.5
# t2w = 0.3
# se = 1
# 
# d = d %>% rowwise() %>%
#   mutate( sim_one_study_set( Nmax = Inf,
#                              Mu = Mu,
#                              t2a = t2a,
#                              m = 50,
#                              t2w = t2w,
#                              se = se,
#                              hack = "affirm",
#                              return.only.published = TRUE ) )
# 
# summary(d$N)
# 
# # calculate noncentrality parameter
# ncp = Mu / sqrt( t2a + t2w + se^2 )
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
                         m,  # sample size
                         
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
  
  # newRow = bind_cols( corrObject$metaCorr,
  #                 corrObject$sanityChecks )
  #@TEMP: DON'T KEEP THE SANITY CHECKS BECAUSE CORRECT_META_PHACK2 doesn't have it
  newRow = corrObject$metaCorr
  
  newRow = newRow %>% add_column(.before = 1,
                                 methName = methName )
  
  
  # "if" condition is hacky way to deal with repRes = NA case
  if ( is.null( nrow(repRes) ) ) repRes = newRow else repRes = bind_rows(repRes, newRow)
  return(repRes)
}


my_ggsave = function(name,
                     width,
                     height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  setwd(.results.dir)
  ggsave( name,
          width = width, 
          height = height)
  
  setwd(.overleaf.dir)
  ggsave( name,
          width = width, 
          height = height)
}




# 2021-5-8 CLUSTER FNS ---------------------------------------------------------------

# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
sbatch_skeleton <- function() {
  return(
    "#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=JOBNAME
#################  
#a file for job output, you can check job progress
#SBATCH --output=OUTFILE
#################
# a file for errors from the job
#SBATCH --error=ERRORFILE
#################
#time you think you need; default is one hour
#SBATCH --time=JOBTIME
#################
#quality of service; think of it as job priority
#SBATCH --qos=QUALITY
#################
#submit to both owners and normal partition
#SBATCH -p normal,owners
#################
#number of nodes you are requesting
#SBATCH --nodes=NODENUMBER
#################
#memory per node; default is 4000 MB
#SBATCH --mem=MEMPERNODE
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=MAILTYPE
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=USER_EMAIL
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=TASKS_PER_NODE
#################
#SBATCH --cpus-per-task=CPUS_PER_TASK
#now run normal batch commands

ml load R/4.0.2
R -f PATH_TO_R_SCRIPT ARGS_TO_R_SCRIPT")
}



generateSbatch <- function(sbatch_params,
                           runfile_path = NA,
                           run_now = F) {
  
  #sbatch_params is a data frame with the following columns
  #jobname: string, specifies name associated with job in SLURM queue
  #outfile: string, specifies the name of the output file generated by job
  #errorfile: string, specifies the name of the error file generated by job
  #jobtime: string in hh:mm:ss format, max (maybe soft) is 48:00:00 
  #specifies the amoung of time job resources should be allocated
  #jobs still running after this amount of time will be aborted
  #quality: kind of like priority, normal works
  #node_number, integer: the number of nodes (computers w/16 cpus each) to allocate 
  #mem_per_node, integer: RAM, in MB, to allocate to each node
  #mailtype, string: ALL, BEGIN, END, FAIL: what types of events should you be notified about via email
  #user_email string: email address: email address to send notifications
  #tasks_per_node: integer, number of tasks, you should probably use 1
  #cpus_per_task: integer, 1-16, number of cpus to use, corresponds to number of available cores per task
  #path_to_r_script: path to r script on sherlock
  #args_to_r_script: arguments to pass to r script on command line
  #write_path: where to write the sbatch file
  #server_sbatch_path: where sbatch files will be stored on sherlock
  #runfile_path is a string containing a path at which to write an R script that can be used to run
  #the batch files generated by this function. 
  #if NA, no runfile will be written
  #run_now is a boolean specifying whether batch files should be run as they are generated
  
  sbatches <- list()
  if (!is.na(runfile_path)) {
    outfile_lines <- c(paste0("# Generated on ",  Sys.time()))
  }
  for (sbatch in 1:nrow(sbatch_params) ) {
    gen_batch <- sbatch_skeleton()
    #set job name
    if (is.null(sbatch_params$jobname[sbatch])) { 
      gen_batch <- gsub("JOBNAME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBNAME", sbatch_params$jobname[sbatch], gen_batch) 
    }
    #set outfile name
    if (is.null(sbatch_params$outfile[sbatch])) { 
      gen_batch <- gsub("OUTFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("OUTFILE", sbatch_params$outfile[sbatch], gen_batch) 
    }
    #set errorfile name
    if (is.null(sbatch_params$errorfile[sbatch])) { 
      gen_batch <- gsub("ERRORFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ERRORFILE", sbatch_params$errorfile[sbatch], gen_batch) 
    }
    #set jobtime
    if (is.null(sbatch_params$jobtime[sbatch])) { 
      gen_batch <- gsub("JOBTIME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBTIME", sbatch_params$jobtime[sbatch], gen_batch) 
    }
    #set quality
    if (is.null(sbatch_params$quality[sbatch])) { 
      gen_batch <- gsub("QUALITY", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("QUALITY", sbatch_params$quality[sbatch], gen_batch) 
    }
    #set number of nodes
    if (is.null(sbatch_params$node_number[sbatch])) { 
      gen_batch <- gsub("NODENUMBER", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("NODENUMBER", sbatch_params$node_number[sbatch], gen_batch) 
    }
    #set memory per node
    if (is.null(sbatch_params$mem_per_node[sbatch])) { 
      gen_batch <- gsub("MEMPERNODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MEMPERNODE", sbatch_params$mem_per_node[sbatch], gen_batch) 
    }
    #set requested mail message types
    if (is.null(sbatch_params$mailtype[sbatch])) { 
      gen_batch <- gsub("MAILTYPE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MAILTYPE", sbatch_params$mailtype[sbatch], gen_batch) 
    }
    #set email at which to receive messages
    if (is.null(sbatch_params$user_email[sbatch])) { 
      gen_batch <- gsub("USER_EMAIL", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("USER_EMAIL", sbatch_params$user_email[sbatch], gen_batch) 
    }
    #set tasks per node
    if (is.null(sbatch_params$tasks_per_node[sbatch])) { 
      gen_batch <- gsub("TASKS_PER_NODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("TASKS_PER_NODE", sbatch_params$tasks_per_node[sbatch], gen_batch) 
    }
    #set cpus per task
    if (is.null(sbatch_params$cpus_per_task[sbatch])) { 
      gen_batch <- gsub("CPUS_PER_TASK", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("CPUS_PER_TASK", sbatch_params$cpus_per_task[sbatch], gen_batch) 
    }
    #set path to r script
    if (is.null(sbatch_params$path_to_r_script[sbatch])) { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", sbatch_params$path_to_r_script[sbatch], gen_batch) 
    }
    #set args to r script
    if (is.null(sbatch_params$args_to_r_script[sbatch])) { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", sbatch_params$args_to_r_script[sbatch], gen_batch) 
    }
    
    #write batch file
    if (is.null(sbatch_params$write_path[sbatch])) { 
      cat(gen_batch, file = paste0("~/sbatch_generated_at_", gsub(" |:|-", "_", Sys.time()) ), append = F)
    } else { 
      cat(gen_batch, file = sbatch_params$write_path[sbatch], append = F)
    }
    
    if (!is.na(sbatch_params$server_sbatch_path[sbatch])) {
      outfile_lines <- c(outfile_lines, paste0("system(\"sbatch ", sbatch_params$server_sbatch_path[sbatch], "\")"))
    } 
    sbatches[[sbatch]] <- gen_batch
  }
  if (!is.na(runfile_path)) {
    cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
  }
  if(run_now) { system(paste0("R -f ", runfile_path)) } 
  
  return(sbatches)
}


########################### FN: STITCH RESULTS FILES ###########################

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # .results.singles.path = "/home/groups/manishad/MRM/sim_results/long"
  # .results.stitched.write.path = "/home/groups/manishad/MRM/sim_results/overall_stitched"
  # .name.prefix = "long_results"
  # .stitch.file.name="stitched.csv"
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  names = names( read.csv(keepers[1] )[-1] )
  
  # read in and rbind the keepers
  tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )
  s <- do.call(rbind, tables)
  
  names(s) = names( read.csv(keepers[1], header= TRUE) )
  
  if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
  return(s)
}

