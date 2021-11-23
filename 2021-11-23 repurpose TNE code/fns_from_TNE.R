

# verbatim from helper_TNE.R on 2021-11-23

# ~ NOTES ----------------------------------------------

# In order to catch errors from individual estimation methods safely and informatively,
#  in general the estimation method fns are structured st they can be run within the
#  fn run_method_safe, which automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error.

# But some methods that are especially delicate (e.g., mle-boot) have their own standalone
#  fns that do all of this (including writing to rep.res) because it seems likely that 
#  they could fail at different points (e.g., mle-boot might get a point estimate but
#  not BCa CIs), in which case we want to catch errors at each step individually to 
#  avoid discarding information.

#@Eventually remove any vestigial fns from MRM at the end
#@Also see if there are vestigial other fns in helper_TNE, but remember some are used by auxiliary scripts

# ~ PACKAGES USED FOR SANITY CHECKS THROUGHOUT THIS CODE ----------------------------------------------

# # commented out for use on cluster
# library(crayon)
# library(dplyr)
# library(foreach)
# library(doParallel)
# library(boot)
# library(metafor) 
# library(robumeta)
# library(data.table)
# library(purrr)
# library(metRology)
# library(fansi)
# library(MetaUtility)
# library(ICC)
# library(cfdecomp)
# library(tidyr)
# library(here)
# library(foreach)
# library(doParallel)
# library(dplyr)
# library(boot)
# library(purrr)
# library(robumeta)
# library(MetaUtility)
# library(truncdist)
# library(tibble)
# library(tmvtnorm)
# library(testthat)
# library(truncSP)
# library(truncnorm)

# ~ ESTIMATION METHOD FNS ----------------------------------------------


# ~~ Estimation Methods That Are Structured for Standalone Use --------

# Fns in this category need to modify and return .rep.res (passed as argument)
#  themselves, just like run_method_safe
# For all of these, .rep.res should be the existing results df from any previously run methods;
#  otherwise should be data.frame() so it's recognized as having 0 rows

# see note at beginning for why this is a standalone fn
# 2021-9-2: MM audited fn by reading
estimate_anything_boot = function(x,
                                  p,
                                  
                                  method.to.boot,
                                  .rep.res,
                                  parallel = "multicore",
                                  
                                  ... ) {  # here you can pass additional args that are method-specific
  
  
  # collect extra args to be passed to specific methods
  extra.args = list(...)
  
  if ( !method.to.boot %in% c("mle-wald", "jeffreys-wald-map") ) stop( "method not implemented in estimate_anything_boot" )
  
  # get estimates from the original sample; these are the ones we will bias-correct
  Mhat.old = .rep.res$Mhat[ .rep.res$method == method.to.boot ]
  Shat.old = .rep.res$Shat[ .rep.res$method == method.to.boot ]
  Vhat.old = .rep.res$Vhat[ .rep.res$method == method.to.boot ]
  
  # check that extra.args contains what we need to run this method
  if ( ! "mu.start" %in% names(extra.args) ) stop("Must pass mu.start")
  
  if ( ! "sigma.start" %in% names(extra.args) ) stop("Must pass sigma.start")
  
  ### Draw Resamples ####
  tryCatch({
    
    # # fake error to test outer tryCatch behavior
    # stop("Haha! Outer tryCatch error in bootstrap!")
    
    # this is just the resampling part, not the CI estimation
    boot.res = my_boot( data = x,
                        parallel = parallel,
                        R = p$boot.reps,
                        statistic = function(original, indices) {
                          
                          # plain nonparametric resample
                          x.bt = original[indices]
                          
                          tryCatch({
                            
                            if ( method.to.boot == "mle-wald" ) {
                              # get point estimates for bootstrap sample
                              ests.bt = estimate_mle( x = x.bt,
                                                      p = p,
                                                      mu.start = extra.args$mu.start,
                                                      sigma.start = extra.args$sigma.start,
                                                      # never need CIs here
                                                      get.CIs = FALSE )
                            }
                            
                            
                            if ( method.to.boot == "jeffreys-wald-map" ) {
                              
                              # get point estimates for bootstrap sample
                              ests.bt = estimate_jeffreys(x = x.bt,
                                                          p = p,
                                                          
                                                          mu.start = extra.args$mu.start,
                                                          sigma.start = extra.args$sigma.start,
                                                          # never need CIs here
                                                          get.CIs = FALSE )
                              
                            }
                            
                            
                            # return the stats of interest
                            # order of stats has to match indices in CI tryCatch loops below
                            # and in returned results because of bt.means and bt.sds
                            c( as.numeric( ests.bt$Mhat ),
                               as.numeric( ests.bt$Vhat ),
                               as.numeric( ests.bt$Shat ) ) 
                            
                          }, error = function(err){
                            # could catch the message here, but note that it will be one message for
                            #  EVERY resample that fails! (lots of output)
                            return( rep(NA, 3) )
                          } )  # end innermost tryCatch
                          
                        } )
    
    # boot diagnostics
    bt.pfails =  as.numeric( colMeans( is.na(boot.res$t) ) )  # proportion of boot reps that failed (NAs)
    bt.means = as.numeric( colMeans(boot.res$t, na.rm = TRUE) )
    bt.bias = bt.means - c( Mhat.old, Vhat.old, Shat.old )  # agrees with boot's own bias calculation
    bt.sds = apply( boot.res$t, 2, function(x) sd(x, na.rm = TRUE) )
    
    error = NA
    
    ### Inference ###
    # get CIs for each estimand individually in case some work and others don't
    
    # CIs for Mhat
    if ( p$get.CIs == TRUE ) {
      tryCatch({
        
        CI = boot.ci(boot.res, type = "bca", index = 1)
        # hacky way around the fact that boot()'s "error" about all values of t being equal is structurally
        #  just a printed string, so the tryCatch won't catch it
        if ( is.null(CI) ) stop("boot.ci returned NULL, probably because all values of t were the same")
        
        # put in nice vector format
        MhatBootCIs = c( CI[[4]][4], CI[[4]][5] )
      }, error = function(err){
        MhatBootCIs <<- c(NA, NA)
      } )
    } else {  # i.e., CIs not wanted
      MhatBootCIs = c(NA, NA)
    }
    
    # CIs for Vhat
    if ( p$get.CIs == TRUE ) {
      tryCatch({
        CI = boot.ci(boot.res, type = "bca", index = 2)
        if ( is.null(CI) ) stop("boot.ci returned NULL, probably because all values of t were the same")
        VhatBootCIs = c( CI[[4]][4], CI[[4]][5] )
        VhatBootSD = sd( boot.res$t[,2] )
      }, error = function(err){
        VhatBootCIs <<- c(NA, NA)
      } )
    } else {  # i.e., CIs not wanted
      VhatBootCIs = c(NA, NA)
    }
    
    # CIs for Shat
    if ( p$get.CIs == TRUE ) {
      tryCatch({
        CI = boot.ci(boot.res, type = "bca", index = 3)
        if ( is.null(CI) ) stop("boot.ci returned NULL, probably because all values of t were the same")
        ShatBootCIs = c( CI[[4]][4], CI[[4]][5] )
        ShatBootSD = sd( boot.res$t[,2] )
      }, error = function(err){
        ShatBootCIs <<- c(NA, NA)
      } )
    } else {  # i.e., CIs not wanted
      ShatBootCIs = c(NA, NA)
    }
    
    # this part happens only if bootstrapping fails completely
    # not just CIs
  }, error = function(err){
    # one list item for each stat of interest (3),
    #  and one sub-entry for lower/upper CI limit
    n.ests = 3  # needs to match number of params in boot.res
    MhatBootCIs <<- c(NA, NA)
    VhatBootCIs <<- c(NA, NA)
    ShatBootCIs <<- c(NA, NA)
    
    bt.means <<- rep(NA, n.ests)
    bt.sds <<- rep(NA, n.ests)
    bt.pfails <<- rep(NA, n.ests)
    bt.bias <<- rep(NA, n.ests)
    
    error <<- err$message
    
  } )  # end of the big tryCatch loop for the whole boot() cal
  
  
  ### Save Boot Results ###
  # add biased-corrected MLE results to iterate results
  new.row = data.frame( method = paste("boot-", method.to.boot, sep=""),
                        
                        Mhat = Mhat.old - bt.bias[1],
                        Vhat = Vhat.old - bt.bias[2],
                        Shat = Shat.old - bt.bias[3],
                        
                        # SE estimates are SDs of boot distribution
                        #  but note that these don't incorporate the BCa correction
                        MhatSE = bt.sds[1],
                        VhatSE = bt.sds[2],
                        ShatSE = bt.sds[3],
                        
                        MLo = MhatBootCIs[1],
                        MHi = MhatBootCIs[2],
                        
                        VLo = VhatBootCIs[1],
                        VHi = VhatBootCIs[2],
                        
                        SLo = ShatBootCIs[1],
                        SHi = ShatBootCIs[2],
                        
                        error = error,
                        # mean boot failures of individual iterates between the two parameters
                        bt.prop.resamples.failed = mean(bt.pfails) ) 
  # currently not keeping track of the other things
  
  if ( nrow(.rep.res) == 0 ) .rep.res = new.row else .rep.res = bind_rows(.rep.res, new.row)
  return(.rep.res)
  
}


# ~~ Estimation Methods That ARE Structured for Use Inside run_method_safe --------

# Fns in this category need to return a list with these elements:
#   Mhat, Vhat, M.SE, M.CI (2-vector), V.SE, V.CI (2-vector)
# Because these fns are run inside run_method_safe, the latter will handle editing rep.res
#  All of these fns should take get.CIs as an argument and return CIs as c(NA, NA) if not wanted

# 2021-9-14: MM audited fn by reading through
estimate_jeffreys_mcmc = function(x,
                                  p,
                                  mu.start,
                                  sigma.start) {
  
  # LL and UU: cutpoints on RAW scale, not Z-scores
  # sigma: SD, not variance
  model.text <- "
functions{
	real jeffreys_prior(real mu, real sigma, real LL, real UU, int n){
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;
		
		mustarL = (LL - mu) / sigma;
		mustarU = (UU - mu) / sigma;
		
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
		
		return sqrt(determinant(fishinfo));
	}
}

data{
	int<lower=0> n;
    real LL;
	real UU;
	real<lower=LL,upper=UU> y[n];
}

parameters{
    real mu;
	real<lower=0> sigma;
}

model{
	target += log( jeffreys_prior(mu, sigma, LL, UU, n) );
	for(i in 1:n)
        y[i] ~ normal(mu, sigma)T[LL,UU];
}

generated quantities{
  real log_lik;
  real log_prior = log(jeffreys_prior(mu, sigma, LL, UU, n));
  real log_post;
  log_lik = normal_lpdf(y | mu, sigma);
  log_lik += -n * log_diff_exp( normal_lcdf(UU | mu, sigma), normal_lcdf(LL | mu, sigma) );  							 
  log_post = log_lik + log_prior;
}
"

# prepare to capture warnings from Stan
stan.warned = 0
stan.warning = NA

# set start values for sampler
init.fcn <- function(o){ list(mu=mu.start, sigma=sigma.start) }

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
  
  # as in E_fisher(), prevent numerical issues due to infinite cutpoints
  .a = max( -99, p$a )
  .b = min( 99, p$b )
  
  cat( paste("\n estimate_jeffreys_mcmc flag 2: about to call sampling") )
  post = sampling(stan.model,
                  cores = 1,
                  refresh = 0,
                  data = list( n = p$n, LL = .a, UU = .b, y = x ),
                  
                  iter = p$stan.iter,   
                  control = list(max_treedepth = p$stan.maxtreedepth,
                                 adapt_delta = p$stan.adapt_delta),
                  
                  init = init.fcn)
  
  
}, warning = function(condition){
  stan.warned <<- 1
  stan.warning <<- condition$message
} )


cat( paste("\n estimate_jeffreys_mcmc flag 3: about to call postSumm") )
postSumm = summary(post)$summary


# posterior means, then medians
Mhat = c( postSumm["mu", "mean"], median( rstan::extract(post, "mu")[[1]] ) )
Shat = c( postSumm["sigma", "mean"], median( rstan::extract(post, "sigma")[[1]] ) )
Vhat = Shat^2
# sanity check
expect_equal( Mhat[1], mean( rstan::extract(post, "mu")[[1]] ) )


# SEs
MhatSE = postSumm["mu", "se_mean"]
ShatSE = postSumm["sigma", "se_mean"]
# because VhatSE uses delta method, VhatSE will be length 2 because Shat is length 2
VhatSE = ShatSE * 2 * Shat  
# how Stan estimates the SE: https://discourse.mc-stan.org/t/se-mean-in-print-stanfit/2869
expect_equal( postSumm["mu", "sd"], sd( rstan::extract(post, "mu")[[1]] ) )
expect_equal( MhatSE,
              postSumm["mu", "sd"] / sqrt( postSumm["mu", "n_eff"] ) )

# CI limits
S.CI = c( postSumm["sigma", "2.5%"], postSumm["sigma", "97.5%"] )
V.CI = S.CI^2
M.CI = c( postSumm["mu", "2.5%"], postSumm["mu", "97.5%"] )
# sanity check:
myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], 0.025 ),
                          quantile( rstan::extract(post, "mu")[[1]], 0.975 ) ) )
expect_equal(M.CI, myMhatCI)


# the point estimates are length 2 (post means, then medians),
#  but the inference is the same for each type of point estimate
return( list( Mhat = Mhat,
              Vhat = Vhat,
              Shat = Shat,
              
              MhatSE = rep(MhatSE, 2),
              VhatSE = VhatSE,  # already length 2 (see above)
              ShatSE = rep(ShatSE, 2),
              
              M.CI = M.CI,
              V.CI = V.CI,
              S.CI = S.CI,
              
              stan.warned = stan.warned,
              stan.warning = stan.warning,
              
              MhatRhat = postSumm["mu", "Rhat"],
              ShatRhat = postSumm["sigma", "Rhat"]
) )
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
  
  # 2021-11-19: try another optimizer (Nelder-Mead)
  myMLE.nm = mle( minuslogl = nll_simple,
                  method = "Nelder-Mead",
                  start = list( .mu=mu.start, .sigma=sigma.start) )
  mles.nm = as.numeric( coef(myMLE.nm) )
  
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
                
                Mhat.opt.diff = Mhat - mles.nm[1]
  ) )
}


# 2021-9-2: MM audited fn by reading
# mu.start, sigma.start: start values for optimatization
# as illustrated in a sanity check after nlpost_simple, this fn's MAPs agree with
#  using mle() directly on nlpost_Jeffreys
estimate_jeffreys = function( x,
                              p,
                              par2is = "sd",
                              mu.start,
                              sigma.start,
                              get.CIs,
                              CI.method = "wald" ) {
  
  ### Get MAP by Calling mle() ###
  # IMPORTANT: This fn cannot be moved outside the scope of estimate_jeffreys
  #  because mle() is too dumb to allow extra args (e.g., x) to be passed,
  #  so it's forced to rely on global vars
  #  and that's a problem with a doParallel loop
  #  if this fn is outside estimate_jeffreys, different parallel iterations will use each other's global vars
  nlpost_simple = function(.mu, .sigma) {
    nlpost.value = nlpost_Jeffreys(.pars = c(.mu, .sigma),
                                   par2is = par2is,
                                   .x = x, .a = p$a, .b = p$b) 
    return(nlpost.value)
  }
  
  #**important: force use of Nelder-Mead optimization, which works better for Jeffreys
  #  (even though BFGS works better for MLE)
  # for more on this issue, see "2021-9-23 SD vs. var reduc with Jeffreys.R"
  res = mle( minuslogl = nlpost_simple,
             start = list( .mu=mu.start, .sigma=sigma.start),
             method = "Nelder-Mead" )
  
  # not actually MLEs, of course, but rather MAPs
  mles = as.numeric(coef(res))
  
  # 2021-11-19: try another optimizer (BFGS)
  myMLE.bfgs = mle( minuslogl = nlpost_simple,
                    start = list( .mu=mu.start, .sigma=sigma.start),
                    method = "BFGS" )
  mles.bfgs = as.numeric( coef(myMLE.bfgs) )
  
  # THIS BEHAVES WELL
  if ( par2is == "sd" ) {
    # need this structure for run_method_safe to understand
    Mhat = mles[1]
    Vhat = mles[2]^2
    Shat = mles[2]
  }
  
  # THIS BEHAVES BADLY
  if ( par2is == "var" ) {
    # need this structure for run_method_safe to understand
    Mhat = mles[1]
    Vhat = mles[2]
    Shat = sqrt(mles[2])
  }
  
  # recode convergence more intuitively
  # optim uses "0" to mean successful convergence
  optim.converged = attr(res, "details")$convergence == 0 
  
  profile.CI.error = NA
  
  ests = c(Mhat, Vhat, Shat)
  
  ### Get CIs ###
  if ( get.CIs == TRUE & CI.method == "wald" & par2is == "sd" ) {
    # get Wald CI 
    # SEs for both parameters
    # this has its own tryCatch in case point estimation was possible, but not inference
    tryCatch({
      
      # IMPORTANT: Despite the name "Ofish", this is actually the Hessian of the nlpost, which incorporates the prior (see Bayesian Data Analysis, page 84)
      # not the Fisher, which would be from the lkl only
      Ofish = attr(res, "details")$hessian
      invFisher = solve(Ofish)
      
      # SEs for Mhat and Shat, leaving blank space for Vhat
      SEs = sqrt( c( invFisher[1,1], NA, invFisher[2,2] ) )
      # fill in VhatSE using delta method
      # let g(y) = y^2, where y=Shat
      SEs[2] = SEs[3] * 2 * Shat
      
      los = ests - SEs * qnorm(0.975)
      his = ests + SEs * qnorm(0.975)
      
    }, error = function(err) {
      SEs <<- los <<- his <<- rep(NA, 3)
      profile.CI.error <<- err$message
    })
    
  } else if ( get.CIs == TRUE & CI.method == "profile" & par2is == "sd" ) {
    
    tryCatch({
      # as confirmed in "2021-8-19 Investigate profile penalized LRT inference",
      #  these are indeed profile CIs
      CIs = confint(res)
      # NA's here represent Vhat
      los = c( as.numeric( CIs[1,1] ), as.numeric( CIs[2,1] )^2, as.numeric( CIs[2,1] ) )
      his = c( as.numeric( CIs[1,2] ), as.numeric( CIs[2,2] )^2, as.numeric( CIs[2,2] ) )
      
      ( res.SEs = as.numeric( attr( summary(res), "coef" )[, "Std. Error" ] ) )
      SEs = c(res.SEs[1], NA, res.SEs[2])
      
    }, error = function(err) {
      SEs <<- los <<- his <<- rep(NA, 3)
      profile.CI.error <<- err$message
    })
    
    
  } else if ( get.CIs == TRUE & CI.method == "wald" & par2is == "var" ) {
    # get Wald CI 
    # SEs for both parameters
    
    # this has its own tryCatch in case point estimation was possible, but not inference
    tryCatch({
      # as noted above, this is actually the Hessian of the nlpost, not the Fisher info
      Ofish = attr(res, "details")$hessian
      invFisher = solve(Ofish)
      
      # SEs for Mhat and Shat, leaving blank space for Shat
      SEs = sqrt( c( invFisher[1,1], invFisher[2,2], NA ) )
      # fill in ShatSE using delta method
      # let g(y) = y^(1/2), where y=Shat
      #  so g'(y) = 0.5 * y^(-0.5)
      SEs[3] = SEs[2] * 0.5*Vhat^(-0.5)
      
      los = ests - SEs * qnorm(0.975)
      his = ests + SEs * qnorm(0.975)
      
    }, error = function(err) {
      SEs <<- los <<- his <<- rep(NA, 3)
    })
    
  } else if ( get.CIs == TRUE & CI.method == "profile" & par2is == "var" ) {
    
    tryCatch({
      # as confirmed in "2021-8-19 Investigate profile penalized LRT inference",
      #  these are indeed profile CIs
      CIs = confint(res)
      # NA's here represent Shat
      los = c( as.numeric( CIs[1,1] ), as.numeric( CIs[2,1] ), sqrt( as.numeric( CIs[2,1] ) ) )
      his = c( as.numeric( CIs[1,2] ), as.numeric( CIs[2,2] ), sqrt( as.numeric( CIs[2,2] ) ) )
      
      
      ( res.SEs = as.numeric( attr( summary(res), "coef" )[, "Std. Error" ] ) )
      SEs = c(res.SEs[1], res.SEs[2], NA)
      SEs[3] = SEs[2] * 0.5*Vhat^(-0.5)
      
    }, error = function(err) {
      SEs <<- los <<- his <<- rep(NA, 3)
    })
    
  } else {  # i.e., if get.CIs == FALSE 
    SEs = los = his = c(NA, NA, NA)
  }
  
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
                # to match output of estimate_mle, return BFGS - NM
                Mhat.opt.diff = mles.bfgs[1] - Mhat, 
                profile.CI.error = profile.CI.error ) )
}


# ~~ Helpers for Above Estimation Methods -------

# only used in auxiliary code (so don't delete)
jeffreys_deviance_test_mu = function(.mu0, .resHA, sigma.start, x, p) {
  
  # test the hypothesis that mu = .mu0
  
  # get supremum of posterior under H0 (i.e., constrain mu = .mu0 and only estimate sigma) 
  resH0 = optim( par = sigma.start,
                 method = "Brent",  # can't use Nelder-Mead for 1D optimization
                 lower = -999, # Brent requires providing lower and upper limits
                 upper = 999,
                 # here .pars is just sigma
                 fn = function(.sigma, .x, .a, .b) {
                   # hold constant mu0
                   nlpost_Jeffreys( .pars = c(.mu0, .sigma), .x = .x, .a = .a, .b = .b )
                 },
                 .x = x,
                 .a = p$a,
                 .b = p$b )
  
  # should be LARGER (lower likelihood) than res$value because we're constraining parameters
  #resH0$value
  
  # deviance = 2( logL_H0 - logL_HA )
  # signs below are because I have negative lposts
  return( 2*(resH0$value - .resHA$value) )
}




# expected Fisher info for all n observations
#  from Leon's theory; checked against his code
E_fisher = function(.mu, .sigma, .n, .a, .b) {
  
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
# # sanity check against Leon's and Blake's
# setwd("~/Dropbox/Personal computer/Independent studies/2021/Truncated normal estimation (TNE)/Linked to OSF (TNE)/Git/Local code")
# source("helper.R")
# mine = E_fisher(.mu = 1, .sigma = 1, .n = 20, .a = 0, .b = 1)
# leon = find_K(m = 1, s = 1, n = 20, ul = 0, uh = 1)
# blake = fisher_info(mu = 1, sigma = 1, LL = 0, UU = 1, n = 20)
# expect_equal(mine, leon, blake) 


# 2 x 2 matrix of second derivatives of log-lkl, evaluated at .x
# aka observed Fisher information, aka negative Hessian
# 
# For MLEs, and presumably also Bayes, better to use observed Fisher vs.
#  expected for asymptotic variance estimation:
#  Efron, B.; Hinkley, D.V. (1978). "Assessing the accuracy of the maximum likelihood estimator: Observed versus expected Fisher Information"
#
# .x: vector all n observations
# .a: lower truncation point in RAW units (not a z-score)
# .b: upper truncation point in RAW units
O_fisher = function(.mu, .sigma, .x, .a, .b) {
  
  n = length(.x)
  
  Za = (.a - .mu)/.sigma
  Zb = (.b - .mu)/.sigma
  
  alpha.a = dnorm(Za) / ( pnorm(Zb) - pnorm(Za) )
  alpha.b = dnorm(Zb) / ( pnorm(Zb) - pnorm(Za) )
  
  # entry 1,1: dl/dmu^2
  H11 = (-n/.sigma^2) + (n/.sigma^2) * ( (alpha.b - alpha.a)^2 + Zb*alpha.b - Za*alpha.a )
  
  # entry 2,2: dl/d.sigma^2
  H22 = (n/.sigma^2) - ( 3 * sum( (.x - .mu)^2 ) / .sigma^4 ) +
    (n/.sigma^2)*( Zb*alpha.b*(Zb^2 - 2) - Za*alpha.a*(Za^2 - 2) +
                     (Zb*alpha.b - Za*alpha.a)^2 )
  
  # entry 1,2 and 2,1: dl / d.mu d.sigma
  H12 = ( -2 * sum(.x - .mu) / .sigma^3 ) +
    (n/.sigma^2)*( alpha.a - alpha.b + alpha.b*Zb^2 - alpha.a*Za^2 +
                     (alpha.a - alpha.b)*(alpha.a*Za - alpha.b*Zb) )
  
  return( matrix( c(-H11, -H12,
                    -H12, -H22), nrow = 2 ) )
}


# 2021-9-2: MM audited fn by reading
# negative log-posterior (for all n observations)
# .pars: (mu, sigma)
nlpost_Jeffreys = function(.pars, par2is = "sd", .x, .a, .b) {
  
  # variance parameterization
  if (par2is == "var") {
    
    .mu = .pars[1]
    .var = .pars[2]
    
    if ( .var < 0 ) return(.Machine$integer.max)
    
    # as in nll()
    term1 = dmvnorm(x = as.matrix(.x, nrow = 1),
                    mean = as.matrix(.mu, nrow = 1),
                    # sigma here is covariance matrix
                    sigma = as.matrix(.var, nrow=1),
                    log = TRUE)
    
    
    term2 = length(.x) * log( pmvnorm(lower = .a,
                                      upper = .b,
                                      mean = .mu,
                                      # remember sigma here is covariance matrix, not the SD
                                      sigma = .var ) ) 
    
    # prior
    term3 = log( sqrt( det( E_fisher(.mu = .mu, .sigma = sqrt(.var), .n = length(.x), .a = .a, .b = .b) ) ) )
    
    nlp.value = -( sum(term1) - term2 + term3 )
    
    if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  }
  
  # SD parameterization
  if (par2is == "sd") {
    
    .mu = .pars[1]
    .sigma = .pars[2]
    
    if ( .sigma < 0 ) return(.Machine$integer.max)
    
    # as in nll()
    term1 = dmvnorm(x = as.matrix(.x, nrow = 1),
                    mean = as.matrix(.mu, nrow = 1),
                    # sigma here is covariance matrix,
                    sigma = as.matrix(.sigma^2, nrow=1),
                    log = TRUE)
    
    
    term2 = length(.x) * log( pmvnorm(lower = .a,
                                      upper = .b,
                                      mean = .mu,
                                      # remember sigma here is covariance matrix, not the SD
                                      sigma = .sigma^2 ) ) 
    
    term3 = log( sqrt( det( E_fisher(.mu = .mu, .sigma = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) )
    
    nlp.value = -( sum(term1) - term2 + term3 )
    
    if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  }
  
  nlp.value
  
}

# same as nlpost_Jeffreys, but formatted for use with mle()
# fn needs to be formatted exactly like this (no additional args)
#  in order for mle() to understand
nlpost_simple = function(.mu, .sigma) {
  nlpost_Jeffreys(.pars = c(.mu, .sigma),
                  .x = x, .a = p$a, .b = p$b)
}


# ### Example and sanity check 
# x = c(0.385258397941442, 1.68066267127739, 0.729227742032434, 0.479120432291688, 
#       0.897279068695914, 0.0575356881970433, 0.165652783807015, 0.875647820475464, 
#       0.380168104717168, 0.825551957468494, 0.589842597791253, 0.402854395794205, 
#       2.7668857263465, 0.649703054651576, 0.83074621395, 0.407612065235468, 
#       1.49884534475478, 0.420770708293162, 0.699931456061883, 1.06188921169992
# )
# p = structure(list(rep.methods = "mle ; jeffreys ; jeffreys-lrt", 
#                    boot.reps = 1000, get.CIs = TRUE, mu = 0, V = 1, Za = 0, 
#                    Zb = 9, n = 20, a = 0, b = 9), row.names = c(NA, -1L), class = c("tbl_df", 
#                                                                                     "tbl", "data.frame"))
# a = p$Za * sqrt(p$V)
# b = p$Zb * sqrt(p$V)
# nlpost_Jeffreys(.pars = c(-1, 3), .x=x, .a = a, .b= b)
# 
# penMLE = mle( minuslogl = nlpost_simple,
#               start = list( .mu=1, .sigma=1) )
# 
# # c.f. estimate_jeffreys
# myMAP = estimate_jeffreys(x = x,
#                           p = p,
#                           mu.start = 3,
#                           sigma.start = 3,
#                           get.CIs = TRUE )
# 
# expect_equal( myMAP$Mhat, as.numeric( coef(penMLE)[1] ), tol = 0.0001 )
# expect_equal( myMAP$Shat, as.numeric( coef(penMLE)[2] ), tol = 0.0001 )


# log-prior for all n observations
# only used for illustrating prior in manuscript
# .pars: (mu, sigma)
lprior_Jeffreys = function(.pars, par2is = "sd", .n, .a, .b) {
  
  .mu = .pars[1]
  
  if (par2is == "var") .var = .pars[2]
  if (par2is == "sd") .var = .pars[2]^2
  
  # prevent numerical issues
  .a = max(.a, -99)
  .b = min(.b, 99)
  
  # as in nlpost_Jeffreys
  lprior = log( sqrt( det( E_fisher(.mu = .mu, .sigma = sqrt(.var), .n = .n, .a = .a, .b = .b) ) ) )
  lprior
  
  # if ( is.infinite(lprior) | is.na(lprior) ) {
  #   return(-.Machine$integer.max)
  # } else {
  #   return(lprior)
  # }
}



# 2021-9-2: MM audited fn by reading
# negative log-lkl (for all n observations)
# for use getting MLEs
# .pars: (mu, sigma) as a 2-vector
# par2is: says which parameterization to use
nll = function(.pars, .x, .a, .b, par2is = "sd") {
  
  # this parameterization behaves badly in sims!
  # see Project Log entry from 2021-8-25
  if ( par2is == "sd" ) {
    .mu = .pars[1]
    .sigma = .pars[2]
    
    # as in mle.tmvnorm, return a huge value if sigma < 0
    # see file "2021-8-22 Debug MLE" for why this is necessary
    if ( .sigma < 0 ) return(.Machine$integer.max)
    
    # write the nll without using dtruncnorm
    # this is taken closely from tmvnorm package's returned nll function
    term1 = dnorm(x = .x,
                  mean = .mu,
                  sd = .sigma,  
                  log = TRUE)
    
    term2 = length(.x) * log( pmvnorm(lower = .a,
                                      upper = .b,
                                      mean = .mu,
                                      # note use of sigma^2 here because of pmvnorm's different parameterization:
                                      sigma = .sigma^2 ) ) 
    
    # begin sanity checks against the other parameterization below
    # the point of this is to show that the terms of the nll agree exactly between
    #  the parameterizations, as they should, of course
    term1.check = dmvnorm(x = as.matrix(.x, nrow = 1),
                          mean = as.matrix(.mu, nrow = 1),
                          sigma = as.matrix(.sigma^2, nrow=1),
                          log = TRUE)
    
    term2.check = length(.x) * log( pmvnorm(lower = .a,
                                            upper = .b,
                                            mean = .mu,
                                            sigma = .sigma^2 ) )
    
    if ( !is.infinite(term1) & !is.infinite(term1.check) & any( abs(term1 - term1.check) > 0.001 ) ) stop("term1 and term1.check don't agree!")
    if ( !is.infinite(term2) & !is.infinite(term2.check) & any( abs(term2 - term2.check) > 0.001 ) ) stop("term2 and term2.check don't agree!")
    # end sanity checks
  }
  
  
  # I think this version always agrees with the tmvnorm package
  if ( par2is == "var" ) {
    .mu = .pars[1]
    .var = .pars[2]
    
    # as in mle.tmvnorm, return a huge value if sigma < 0
    # see file "2021-8-22 Debug MLE" for why this is necessary
    if ( .var < 0 ) return(.Machine$integer.max)
    
    # try rewriting without using dtruncnorm
    # this is taken closely from tmvnorm package's returned nll function
    term1 = dmvnorm(x = as.matrix(.x, nrow = 1),
                    mean = as.matrix(.mu, nrow = 1),
                    sigma = as.matrix(.var, nrow=1),
                    log = TRUE)
    
    # begin sqnity checks against the other parameterization
    # MLEs do NOT necessarily agree with pkg when using this one
    term1.check = dnorm(x = .x,
                        mean = .mu,
                        sd = sqrt(.var),
                        log = TRUE)
    
    summary( abs(term1 - term1.check) )
    if ( any( abs(term1 - term1.check) > 0.001 ) ) stop("term1 and term1.check don't agree!")
    
    
    term2 = length(.x) * log( pmvnorm(lower = .a,
                                      upper = .b,
                                      mean = .mu,
                                      sigma = .var ) ) 
    # end sanity checks
    
  }
  
  
  nll.value = -( sum(term1) - term2 )
  
  if ( is.infinite(nll.value) | is.na(nll.value) ) return(.Machine$integer.max)
  
  return(nll.value)
}

#nll(.pars = c(1,1), .x = x, .a = 0, .b = 9)

# ### Sanity check against "negloglik" as defined inside mle.tmvnorm:
# p = structure(list(rep.methods = "mle", 
#                    boot.reps = 50, get.CIs = TRUE, mu = 0, V = 1, Za = 0, Zb = 0.5, 
#                    n = 1000, a = 0, b = 0.5), row.names = c(NA, -1L), class = c("tbl_df", 
#                                                                                 "tbl", "data.frame"))
# p$a = p$mu + (p$Za * sqrt(p$V))
# p$b = p$mu + (p$Zb * sqrt(p$V))
# 
# x = rtrunc( n = p$n,
#             spec = "norm",
#             mean = p$mu,
#             sd = sqrt(p$V),
#             a = p$a,
#             b = p$b )
# 
# 
# pkgMLE = mle.tmvnorm( X = as.matrix(x, ncol = 1),
#                       lower = p$a,
#                       upper = p$b )
# ( pkgEsts = c( coef(pkgMLE)[1], sqrt( coef(pkgMLE)[2] ) ) ) # Mhat, Shat instead of Vhat
# 
# 
# pkg.nll = -as.numeric( logLik(pkgMLE) )
# 
# # compare to my fn
# my.nll = nll(.pars = c(pkgEsts[1], pkgEsts[2]),
#     .x = x,
#     .a = p$a,
#     .b = p$b )
# 
# expect_equal(pkg.nll, my.nll, tol = 0.001)
# ### end sanity check


# # sanity check against Blake's
# fisher_info <- function(mu, sigma, LL, UU, n){		
#   mustarL <- (LL - mu) / sigma
#   mustarU <- (UU - mu) / sigma
#   alphaL <- dnorm(mustarL) / (pnorm(mustarU) - pnorm(mustarL))
#   alphaU <- dnorm(mustarU) / (pnorm(mustarU) - pnorm(mustarL))
#   
#   kmm <- -n/sigma^2 + n/sigma^2 * ((alphaU-alphaL)^2 + alphaU*mustarU- alphaL*mustarL)
#   kms <- -2*n/sigma^2 * (alphaL - alphaU) + 
#     n/sigma^2 * (alphaL - alphaU + (alphaU*mustarU^2 - alphaL*mustarL^2) +
#                    (alphaL-alphaU) * (alphaL*mustarL - alphaU*mustarU))
#   kss <- n/sigma^2 - 3*n/sigma^2 * (1 + mustarL*alphaL - mustarU*alphaU) +
#     n/sigma^2 * (mustarU*alphaU*(mustarU^2 - 2) - mustarL*alphaL*(mustarL^2 - 2) +
#                    (alphaU*mustarU - alphaL*mustarL)^2)
#   
#   fishinfo <- matrix(NA, 2, 2)
#   fishinfo[1,1] <- -kmm
#   fishinfo[1,2] <- -kms
#   fishinfo[2,1] <- -kms
#   fishinfo[2,2] <- -kss
#   fishinfo				
# }
# nlp <- function(theta, y, LL, UU){
#   -sum(log(dtruncnorm(y, LL, UU, theta[1], theta[2]))) -
#     log(sqrt(det( fisher_info(theta[1], theta[2], LL, UU, length(y)) )))
# }
# mine = nlpost_Jeffreys(.pars = c(1,1), .x = c(.1, .75), .a = 0, .b = 1)
# blake = nlp(theta = c(1,1), y = c(.1, .75), LL = 0, UU = 1)
# expect_equal(mine, blake)



# ~~ Wrapper Fn to Safely Run a Method -------

# See note at the beginning of this script
#  this fn automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error

# Important: this fn works if method is length 2 because a single method.fn returns multiple rows
# BUT in that case, it assumes that the CIs are shared for all rows of that method

# expects global vars: all.errors, rep.res
# directly edits res via superassignment
run_method_safe = function( method,
                            method.fn,
                            .rep.res ) {
  
  cat( paste("\n run_method_safe flag 1: about to try running method", method) )
  
  
  tryCatch({
    
    stats = method.fn()
    
    cat( paste("\n run_method_safe flag 2: done calling method.fn() for", method) )
    
    # fill in things that not every method returns
    #  without this step, data.frame() below will say that 
    #  args imply differing numbers of rows because some things will be NULL instead of NA
    might.be.null = c("optim.converged",
                      "stan.warned",
                      "stan.warning",
                      "profile.CI.error",
                      "MhatRhat",
                      "ShatRhat",
                      "Mhat.opt.diff")
    
    for ( .var in might.be.null ) {
      if ( !.var %in% names(stats) ) stats[[.var]] = NA
    }
    
    Mhat = stats$Mhat
    Vhat = stats$Vhat
    Shat = stats$Shat
    
    MhatSE = stats$MhatSE
    VhatSE = stats$VhatSE
    ShatSE = stats$ShatSE
    
    M.CI = stats$M.CI
    V.CI = stats$V.CI
    S.CI = stats$S.CI
    
    optim.converged = stats$optim.converged
    Mhat.opt.diff = stats$Mhat.opt.diff
    stan.warned = stats$stan.warned
    stan.warning = stats$stan.warning
    profile.CI.error = stats$profile.CI.error
    MhatRhat = stats$MhatRhat
    ShatRhat = stats$ShatRhat
    
    error = NA
    
  }, error = function(err) {
    # needs to be superassignment because inside the "error" fn
    error <<- err$message
    Mhat <<- NA
    Vhat <<- NA
    Shat <<- NA
    
    MhatSE <<- NA
    VhatSE <<- NA
    ShatSE <<- NA
    
    M.CI <<- c(NA, NA)
    V.CI <<- c(NA, NA)
    S.CI <<- c(NA, NA)
    optim.converged <<- NA
    Mhat.opt.diff <<- NA
    profile.CI.error <<- NA
    stan.warned <<- NA
    stan.warning <<- NA
    MhatRhat <<- NA
    ShatRhat <<- NA
  })
  
  cat( paste("\n run_method_safe flag 3: about to make new.row for method", method) )
  
  # add the result to results dataframe
  # IMPORTANT: if you add anything here, need to update the above tryCatch 
  #  structure with another xxx <<- NA accordingly
  new.row = data.frame( method,
                        
                        Mhat = Mhat,
                        Vhat = Vhat,
                        Shat = Shat,
                        
                        MhatSE = MhatSE,
                        VhatSE = VhatSE,
                        ShatSE = ShatSE,
                        
                        
                        MLo = M.CI[1],
                        MHi = M.CI[2],
                        
                        VLo = V.CI[1],
                        VHi = V.CI[2],
                        
                        SLo = S.CI[1],
                        SHi = S.CI[2],
                        
                        optim.converged = optim.converged,
                        Mhat.opt.diff = Mhat.opt.diff,
                        profile.CI.error = profile.CI.error,
                        stan.warned = stan.warned,
                        stan.warning = stan.warning,
                        MhatRhat = MhatRhat,
                        ShatRhat = ShatRhat,
                        
                        overall.error = error) 
  
  
  if ( nrow(.rep.res) == 0 ) .rep.res = new.row else .rep.res = bind_rows(.rep.res, new.row)
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


# ~ FNS FOR SUMMARIZING SIMULATION DATA ----------------------------------------------


# fn for aggregating so we can look at different
#  iterate-level filtering rules
# .s: the iterate-level stitched data (not yet aggregated in any way)
# averagefn: fn to use when aggregating results across scenarios
# expected.sim.reps: only used for sanity checks
make_agg_data = function( .s,
                          .averagefn = "median",
                          badCoverageCutoff = 0.85,
                          expected.sim.reps = NA ){
  
  
  # make unique scenario variable, defined as scen.name AND method
  if ( !"unique.scen" %in% names(.s) ) .s$unique.scen = paste(.s$scen.name, .s$method)
  
  ##### Outcome and Parameter Variables #####
  # "outcome" variables used in analysis
  analysis.vars = c( 
    "Mhat",
    "Vhat",
    "Shat",
    
    "MLo",
    "MHi",
    
    "VLo",
    "VHi",
    
    "SLo",
    "SHi",
    
    
    ##### variables to be created in mutate below:
    
    "MhatBias",
    "VhatBias",
    "ShatBias",
    
    # "MhatRelBias",
    # "VhatRelBias",
    
    "MhatCover",
    "VhatCover",
    "ShatCover",
    
    "MhatWidth",
    "VhatWidth",
    "ShatWidth",
    
    "MhatRMSE",
    "VhatRMSE",
    "ShatRMSE",
    
    "MhatEstSE",
    "VhatEstSE",
    "ShatEstSE",
    
    "MhatEmpSE",
    "VhatEmpSE",
    "ShatEmpSE",
    
    # diagnostics regarding point estimation and CIs
    "MhatEstFail",
    "MhatCIFail",
    "ShatEstFail",
    "ShatCIFail",
    
    # diagnostics regarding bootstraps
    "BtPropResamplesFail",
    "BtMhatCIFail",
    "BtVhatCIFail",
    "BtShatCIFail",
    
    # other diagnostics
    "OptimConverged",
    "MhatOptimDiff",
    "MhatOptimDisagree",
    "StanWarned"
  )
  
  
  # variables that define the scenarios
  param.vars = c("unique.scen",  
                 "method",
                 "boot.reps",
                 "stan.iter",
                 "stan.adapt_delta",
                 "stan.maxtreedepth",
                 "trunc.type",
                 "prop.retained",
                 "mu",
                 "V",
                 "n")
  
  
  # sanity check to make sure we've listed all param vars
  t = .s %>% group_by_at(param.vars) %>% summarise(n())
  if ( !is.na(expected.sim.reps) ) {
    if ( max(t$`n()`) > expected.sim.reps ) stop("param.vars in make_agg_data might be missing something because grouping that way indicated some scenarios had more than expected.sim.reps")
  }
  
  
  ##### Overwrite Analysis Variables As Their Within-Scenario Means #####
  # organize variables into 3 mutually exclusive sets: 
  # - param.vars: parameter variables for grouping
  # - toDrop: variables to drop completely
  # - firstOnly: variables that are static within a scenario, for which we
  #   should just take the first one
  # - takeMean: variables for which we should take the mean within scenarios
  
  names(.s)[ !names(.s) %in% param.vars ]  # look at names of vars that need categorizing
  toDrop = c("rep.methods",
             "get.CIs",
             "error",
             "bt.prop.resamples.failed",
             "sim.reps",  # this is the INTENDED sim.reps, so confusing to retain it
             "rep.name",
             "doParallel.seconds",
             "optim.converged",
             "Mhat.opt.diff",
             "stan.warned")
  firstOnly = c("scen.name",
                "unique.scen",
                "Za",  # calculated from theory, so fixed within scen params
                "Zb"
  )
  
  ##### Add New Variables Calculated at Scenario Level #####
  
  # prevent errors for non-bootstrap methods
  if ( ! "bt.prop.resamples.failed" %in% names(.s) ) .s$bt.prop.resamples.failed = NA
  
  # if you have 10K iterates, script breaks from here forward if running locally
  # "vector memory limits"
  s2 = .s %>%
    rename(
      # static within scenario
      # just renaming for clarity
      MhatEstSE = MhatSE,
      VhatEstSE = VhatSE,
      ShatEstSE = ShatSE ) %>%
    
    # take just first entry of non-parameter variables that are static within scenarios
    group_by_at(param.vars) %>%
    mutate_at( firstOnly, 
               function(x) x[1] ) %>%
    
    # make certain ad hoc variables that don't conform to below rules
    # this step creates variables that are repeated for every rep within 
    #  a scenario, which is intentional
    
    # make variables that are calculated within scenarios
    # some of the vars are defined at the iterate level (i.e., they still vary within scen), 
    #  while others are calculated at the scen level (i.e., they are static within scen)
    # after this step, we take the means within scens of all these vars, which is immaterial
    #   for the ones that are already static within scenario
    group_by_at(param.vars) %>%
    
    mutate( sim.reps = n(),
            
            # varies within scenario
            MhatBias = Mhat - mu,
            VhatBias = Vhat - V,
            ShatBias = Shat - sqrt(V),
            
            # varies within scenario
            MhatCover = covers(truth = mu, lo = MLo, hi = MHi),
            VhatCover = covers(truth = V, lo = VLo, hi = VHi),
            ShatCover = covers(truth = sqrt(V), lo = SLo, hi = SHi),
            
            # varies within scenario
            MhatWidth = MHi - MLo,
            VhatWidth = VHi - VLo,
            ShatWidth = SHi - SLo,
            
            # static within scenario
            MhatRMSE = sqrt( meanNA( (Mhat - mu)^2 ) ),
            VhatRMSE = sqrt( meanNA( (Vhat - V)^2 ) ),
            ShatRMSE = sqrt( meanNA( ( Shat - sqrt(V) )^2 ) ),
            
            # static within scenario
            MhatEstFail = mean(is.na(Mhat)),
            MhatCIFail = mean(is.na(MLo)),
            ShatEstFail = mean(is.na(Shat)),
            ShatCIFail = mean(is.na(SLo)),
            
            # static within scenario
            BtPropResamplesFail = mean(bt.prop.resamples.failed),
            BtMhatCIFail = mean( is.na(MLo) ),
            BtVhatCIFail = mean( is.na(VLo) ),
            BtShatCIFail = mean( is.na(SLo) ),
            
            # static within scenario
            MhatEmpSE = sd(Mhat, na.rm = TRUE),
            VhatEmpSE = sd(Vhat, na.rm = TRUE),
            ShatEmpSE = sd(Shat, na.rm = TRUE),
            
            
            
            # varies within scenario
            # how much smaller is estimated SE compared to empirical one?
            MhatSEBias = MhatEstSE - MhatEmpSE,
            VhatSEBias = VhatEstSE - VhatEmpSE,
            ShatSEBias = ShatEstSE - ShatEmpSE,
            
            # varies within scenario
            MhatSERelBias = (MhatEstSE - MhatEmpSE) / MhatEmpSE, 
            VhatSERelBias = (VhatEstSE - VhatEmpSE) / VhatEmpSE,
            ShatSERelBias = (ShatEstSE - ShatEmpSE) / ShatEmpSE,
            
            # static within scenario
            OptimConverged = meanNA(optim.converged),
            MhatOptimDiff = meanNA(Mhat.opt.diff),
            MhatOptimDisagree = meanNA( abs( Mhat.opt.diff > 0.001 ) ),
            StanWarned = meanNA(stan.warned)
    ) 
  
  
  # now look for which variables should have their means taken
  # this step must happen here, after we've started making s2, 
  #  so that the takeMean vars are actually in s2
  ( takeMean = names(s2)[ !names(s2) %in% c(param.vars, toDrop, firstOnly) ] )
  # sanity check: have all variables been sorted into these categories?
  expect_equal( TRUE,
                all( names(s2) %in% c(param.vars, toDrop, firstOnly, takeMean) ) )
  
  
  ##### Aggregate to Scenario Level #####
  
  # calculate scenario-level averages, but keep dataset at the rep level
  #  for now to facilitate sanity checks
  # IMPORTANT: this uses meanNA regardless of the passed .avgfun 
  #  because right now we are only aggregating WITHIN scens
  #  so we should never use median 
  
  # don't try to drop vars that don't exist
  toDrop = toDrop[ toDrop %in% names(s2) ]
  
  s3 = s2 %>%
    # take averages of numeric variables
    group_by_at(param.vars) %>%
    mutate_at( takeMean,
               function(x) meanNA(x) ) %>%
    select( -all_of(toDrop) )
  
  
  # sanity check: name mismatches
  if ( length( analysis.vars[ !analysis.vars %in% names(s2) ] ) > 0 ) {
    stop("Might have name mismatches; edit analysis.vars in make_agg_data")
  }
  
  # sanity check: SDs of all analysis variables should be 0 within unique scenarios
  t = data.frame( s3 %>% group_by(unique.scen) %>%
                    summarise_at( analysis.vars, sd ) )
  
  t = t %>% select(-unique.scen)
  expect_equal( FALSE,
                any( !as.matrix( t[, 2:(ncol(t)) ] ) %in% c(0, NA, NaN) ) )
  # end sanity checks
  
  
  ##### Aggregate Data at Scenario Level #####
  # make aggregated data by keeping only first row for each 
  #  combination of scenario name and calib.method
  agg = s3[ !duplicated(s3$unique.scen), ]
  
  ##### create Variables That Are Defined At Scenario Rather Than Iterate Level #####
  agg = agg %>% mutate( BadMhatCover = MhatCover < badCoverageCutoff,
                        BadShatCover = ShatCover < badCoverageCutoff )
  
  # # absolute bias is now just the absolute value of bias
  # #  and this is only used for the regressions in Supplement
  # agg$PhatAbsBias = abs(agg$PhatBias)
  # agg$DiffAbsBias = abs(agg$DiffBias)
  
  ##### Make New Variables At Scenario Level ##### 
  # label methods more intelligently for use in plots
  agg$method.pretty.est = NA
  agg$method.pretty.est[ agg$method %in% c("mle-wald", "mle-profile") ] = "MLE"
  agg$method.pretty.est[ agg$method %in% c("boot-mle-wald") ] = "MLE + boot"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-wald-map", "jeffreys-profile-map") ] = "Jeffreys mode"
  agg$method.pretty.est[ agg$method %in% c("boot-jeffreys-wald-map") ] = "Jeffreys mode + boot"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-pmed") ] = "Jeffreys median"
  agg$method.pretty.est[ agg$method %in% c("jeffreys-mcmc-pmean") ] = "Jeffreys mean"
  table(agg$method, agg$method.pretty.est)
  
  agg$method.pretty.inf = NA
  agg$method.pretty.inf[ agg$method %in% c("mle-wald") ] = "MLE Wald"
  agg$method.pretty.inf[ agg$method %in% c("mle-profile") ] = "MLE profile"
  agg$method.pretty.inf[ agg$method %in% c("boot-mle-wald") ] = "MLE BCa"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-wald-map") ] = "Jeffreys mode Wald"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-profile-map") ] = "Jeffreys profile"
  agg$method.pretty.inf[ agg$method %in% c("boot-jeffreys-wald-map") ] = "Jeffreys mode BCa"
  agg$method.pretty.inf[ agg$method %in% c("jeffreys-mcmc-pmed", "jeffreys-mcmc-pmean") ] = "Jeffreys posterior quantiles"
  table(agg$method, agg$method.pretty.inf)
  
  agg$`Truncation type` = agg$trunc.type
  agg$`Truncation type`[ agg$`Truncation type` == "single" ] = "Single"
  agg$`Truncation type`[ agg$`Truncation type` == "double-symm" ] = "Symmetric double"
  agg$`Truncation type`[ agg$`Truncation type` == "double-asymm" ] = "Asymmetric double"
  
  agg$MethodUsesBoot = grepl( "boot", agg$method )
  
  
  return(agg %>% ungroup() )
}


# summarize performance metrics given a dataset (dat) that is already scenario-aggregated
#  looks for all variables with "Bias" or "Cover" in their names and takes their means

# IMPORTANT: this does NOT group by method intentionally

# description: description of the row to make a nice table
# .selectVars: "Mhat" or "Vhat" (by default looks for a global var called selectVars)
# .outcomeType: "est" or "inf" for which set of outcomes to include (also affects how methods are pooled together)
# .methodsToKeep: which methods to show in the table (names should be as in methods.pretty.est or 
#   methods.pretty.inf)
# .meanVars: which performance metrics to show, and in what order (can leave unspecified to get all of them)

my_summarise = function(dat,
                        description = NA,
                        .selectVars = selectVars, # "Mhat" or "Shat"
                        .outcomeType,  # "est" or "inf"
                        .methodsToKeep = NULL,  
                        .meanVars = NULL,
                        #badCoverageCutoff = 0.85,
                        #badWidthCutoff = 2,
                        averagefn = "mean"
){
  
  # # test only
  # dat = agg
  # badCoverageCutoff = 0.85
  # badWidthCutoff = 0.90
  # averagefn = "median.pctiles"
  # .selectVars = c("Mhat")
  
  
  # decide which method variable to use
  if ( .outcomeType == "est" ) dat$method = dat$method.pretty.est
  if ( .outcomeType == "inf" ) dat$method = dat$method.pretty.inf
  
  # subset methods if needed
  if ( !is.null(.methodsToKeep) ) dat = dat %>% filter( method %in% .methodsToKeep )
  
  
  # variables to use as columns (will take their averages)
  # if user specified the variables
  if ( !is.null(.meanVars) ) {
    meanVars = .meanVars
  } else {
    # if user didn't specify which variables to use
    meanVars = c( namesWith(pattern = "Bias", dat = dat), 
                  namesWith(pattern = "Cover", dat = dat),
                  namesWith(pattern = "Width", dat = dat),
                  namesWith(pattern = "RMSE", dat = dat),
                  namesWith(pattern = "Fail", dat = dat),
                  namesWith(pattern = "EmpSE", dat = dat),
                  namesWith(pattern = "Rhat", dat = dat) )
    
    if (.selectVars == "Mhat") meanVars = meanVars[ !grepl(pattern = "Vhat", x = meanVars) &
                                                      !grepl(pattern = "Shat", x = meanVars) ]
    
    if (.selectVars == "Vhat") meanVars = meanVars[ !grepl(pattern = "Mhat", x = meanVars) & 
                                                      !grepl(pattern = "Shat", x = meanVars) ]
    
    if (.selectVars == "Shat") meanVars = meanVars[ !grepl(pattern = "Mhat", x = meanVars) & 
                                                      !grepl(pattern = "Vhat", x = meanVars) ]
    
    meanVars = c(meanVars, "StanWarned")
  }
  
  
  if (averagefn == "mean") avgfun = function(x) meanNA(x)
  if (averagefn == "median") avgfun = function(x) medNA(x)
  if (averagefn == "median.pctiles") avgfun = function(x) medNA_pctiles(x)
  
  # make summary table
  tab = dat %>% 
    group_by(`Truncation type`, method) %>%
    summarise_at( .vars = meanVars, 
                  function(x) avgfun(x) )
  
  
  # nicely round the results, but only numeric columns
  var.classes = lapply(tab, class)
  varNumerics = sapply(var.classes, function(VEC) {
    any(VEC %in% c("numeric", "integer"))
  })
  num.vars = names(varNumerics)[varNumerics]
  
  tab = tab %>% mutate_at( .vars = num.vars,
                           .funs = function(x) round(x, 2) )
  
  # nicely order the truncation types
  tab = tab %>% arrange( factor(`Truncation type`,
                                levels = c("Single", "Symmetric double", "Asymmetric double") ) )
  
  
  #@OLDER VERSION (save?):
  # # if using median.pctiles, have to round only certain cols because others are strings
  # if ( averagefn == "median.pctiles" ){
  #   
  #   if (.selectVars == "Mhat"){
  #     tab$BadMhatCover = round(tab$BadMhatCover, 2)
  #     tab$BadMhatWidth = round(tab$BadMhatWidth, 2)
  #   }
  #   
  #   if (.selectVars == "Vhat"){
  #     tab$BadVhatCover = round(tab$BadVhatCover, 2)
  #     tab$BadVhatWidth = round(tab$BadVhatWidth, 2)
  #   } 
  #   
  # } else {  # otherwise round all columns
  #   tab = round( tab, 2 )
  # }
  
  tab = tab %>% add_column(n.scens = nrow(dat), .before = 1 )
  
  if ( !is.na(description) ) tab = tab %>% add_column(Scenarios = description, .before = 1)
  return(tab)
}


# for prettifying table of regression results
my_recode = function(x) {
  x = sub( pattern = "muN", replacement = "E[N]", x = x )
  x = sub( pattern = "true.effect.distnormal", replacement = "normal effects", x = x )
  x = sub( pattern = "clusteredTRUE", replacement = "clustered", x = x )
  x = sub( pattern = "contrast.extremeTRUE", replacement = "BC-rare", x = x )
  
  return(x)
}



# ~ ANALYSIS FNS ----------------------------------------------


# TAKEN FROM SAPH
# .obj: object returned by correct_meta_phack2
# showAffirms: should it show all studies, even affirms?
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




# .estName: "Mhat" or "Vhat"
# .outcomeName: "Bias", etc.
plot_by_n = function(.dp,
                     .estName,
                     .outcomeName,
                     .trunc.type,
                     
                     .y.breaks = NULL,
                     .include.ggtitle = TRUE, # by default uses .outcomeName, but can omit for manuscript figs 
                     .ggtitle.prefix = "",  # gets added to the automatic title for, e.g., misspecification scenarios
                     .showAllMethods = TRUE,  # TRUE = for supplement
                     
                     .writePlot = FALSE,
                     
                     # below needed only if .writePlot = TRUE
                     .name,
                     .results.dir,
                     .overleaf.dir.general
) {
  
  
  # ~~ Set PDF name ----
  if (.showAllMethods == TRUE) {
    name = paste( .trunc.type,
                  "_",
                  tolower(.estName),
                  "_n_",
                  tolower(.outcomeName),
                  "_plot.pdf",
                  sep = "" )
  } else {
    name = paste( .trunc.type,
                  "_",
                  tolower(.estName),
                  "_n_",
                  tolower(.outcomeName),
                  "_plot_mainText.pdf",
                  sep = "" )
  }
  
  
  
  # ~~ Set method variable for use in plot ----
  # will depend on whether the outcome is a point estimate or inference
  estimationOutcomes = c("Bias", "RMSE", "EmpSE", "EstFail")
  infOutcomes = c("Cover", "Width", "CIFail")
  
  if ( .outcomeName %in% estimationOutcomes ){
    .dp$method = .dp$method.pretty.est
    
    # omit methods for main text
    if ( .showAllMethods == FALSE ) .dp = .dp %>% filter( method.pretty.est %in% mainTextEstMethods )
    
  } else if ( .outcomeName %in% infOutcomes ) {
    .dp$method = .dp$method.pretty.inf
    
    if ( .showAllMethods == FALSE ) .dp = .dp %>% filter( method.pretty.inf %in% mainTextInfMethods )
  } 
  
  ### Force ordering of facets on trunc type
  if ( .trunc.type == "all" ) .dp$`Truncation type` = factor( .dp$`Truncation type`,
                                                              levels = c("Single", "Symmetric double", "Asymmetric double") )
  
  
  # ~~ Set ggplot color palette ----
  # to see all palettes:
  # par(mar=c(3,4,2,2))
  # display.brewer.all()
  n.colors.needed = length(unique(.dp$method))
  
  # NOT IN USE ANYMORE?
  .colors = brewer.pal(n = n.colors.needed, name = "Dark2")
  if( length(.colors) > n.colors.needed ) .colors = .colors[1:n.colors.needed]
  # this maps the colors onto levels of the factor
  names(.colors) = levels( factor(.dp$method) )
  
  # hex color picker:
  # https://www.webfx.com/web-design/color-picker/
  
  # make own color scale to force more striking colors for the key methods
  # first answer: https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
  if ( .outcomeName %in% estimationOutcomes ){
    
    # .colors[ grepl( pattern = "MLE", names(.colors) ) ] = "red"
    # .colors[ grepl( pattern = "Jeffreys", names(.colors) ) ] = "black"
    .colors[ names(.colors) %in% c("MLE", "MLE + boot") ] = "red"
    .colors[ names(.colors) %in% c("Jeffreys mode", "Jeffreys mode + boot") ] = "black"
    .colors[ names(.colors) %in% c("Jeffreys median") ] = "#0CB1E4"
    .colors[ names(.colors) %in% c("Jeffreys mean") ] = "#0FD284"
    
  } else if ( .outcomeName %in% infOutcomes ) {
    
    # .colors[ grepl( pattern = "MLE", names(.colors) ) ] = "red"
    # .colors[ grepl( pattern = "Jeffreys", names(.colors) ) ] = "black"
    
    .colors[ names(.colors) == "MLE profile" ] = "red"
    .colors[ names(.colors) == "MLE Wald" ] = "#D08C04"
    .colors[ names(.colors) == "MLE BCa" ] = "#D004C9"
    .colors[ names(.colors) == "Jeffreys mode Wald" ] = "#0CB1E4"
    .colors[ names(.colors) == "Jeffreys mode BCa" ] = "#0FD284"
    .colors[ names(.colors) == "Jeffreys posterior quantiles" ] = "black"
  } 
  
  myColorScale = scale_colour_manual(values = .colors)
  
  # ~~ Set ggplot linetype scale ----
  # by default, dotted lines
  .lty = rep("dotdash", nuni(.dp$method))
  names(.lty) = names(.colors)
  
  
  if ( .outcomeName %in% estimationOutcomes ){
    .lty[ names(.lty) == "MLE" ] = "solid"
    .lty[ names(.lty) == "Jeffreys mode" ] = "solid"
    .lty[ grepl( pattern = "boot", names(.lty) ) ] = "dotted"
    
  } else if ( .outcomeName %in% infOutcomes ) {
    .lty[ names(.lty) == "MLE profile" ] = "solid"
    .lty[ names(.lty) == "Jeffreys posterior quantiles" ] = "solid"
    .lty[ grepl( pattern = "BCa", names(.lty) ) ] = "dotted"
  } 
  myLtyScale = scale_linetype_manual(values = .lty)
  
  
  # ~~ Set presence of legend for multi-panel prettiness ----
  # for manuscript prettiness, only include legend for last subfigure
  #  of each panel
  include.legend = FALSE 
  # main text
  if ( .showAllMethods == FALSE & .outcomeName %in% c("RMSE", "EmpSE", "Width") ) include.legend = TRUE
  # Supplement
  if ( .showAllMethods == TRUE ) include.legend = TRUE
  
  
  if ( include.legend == TRUE ) legend.position = "bottom" else legend.position = "none"
  
  #browser()
  # ~~ Set plot dimensions that depend on main text vs. Supplement ----
  # dimensions should depend on number of facets, which in turn depends on 
  #  which trunc type
  width.coef = ifelse( .showAllMethods == FALSE, 6, 6*1.2 )
  height.coef = ifelse( .showAllMethods == FALSE, 5.5, 5.5*1.2 )
  
  if ( .trunc.type %in% c("single", "double-symm") ) {
    width = width.coef*3  # 3 * (number of cols in facet_wrap) works well
    height = height.coef*2  # 4 * (number of rows in facet_wrap) works well
    
    if ( .showAllMethods == TRUE ) height = height.coef*3  # because of legends
    
  } else if ( .trunc.type == "double-asymm" ) {
    width = width.coef*3
    height = height.coef*2.5
  } else if ( .trunc.type == "all" ) {
    width = width.coef*3
    height = height.coef*3
  }
  
  # for the plots with legends, give them a little extra vertical height
  #if (include.legend == TRUE) height.coef = height.coef*1.2
  
  # ~~ Make pretty "% retained" variable for facets ----
  .dp$cutpoints.pretty = paste(round( 100*.dp$prop.retained ),  # this is actually PROPORTION retained
                               "% retained",
                               " (",
                               round(.dp$Za, 2),
                               ", ",
                               round(.dp$Zb, 2),
                               ")",
                               
                               sep = "" )
  
  .dp$mu.pretty = paste( "mu = ", .dp$mu )
  
  yName = paste( .estName, .outcomeName, sep="" )
  .dp$Y = .dp[[ yName ]]
  
  
  # ~~ Make base plot ----
  # include titles for all Supplement figures
  title.string = ""
  if ( .showAllMethods == TRUE ) .include.ggtitle = TRUE
  if ( .include.ggtitle == TRUE ) {
    
    if ( .trunc.type == "single" ) trunc.type.pretty = "Single truncation, "
    if ( .trunc.type == "double-symm" ) trunc.type.pretty = "Symmetric double truncation, "
    if ( .trunc.type == "double-asymm" ) trunc.type.pretty = "Asymmetric double truncation, "
    if ( .trunc.type == "all" ) trunc.type.pretty = "All truncation types, "
    
    # e.g., add prefix that these are misspecified scenarios
    if ( .ggtitle.prefix != "" ) trunc.type.pretty = paste( .ggtitle.prefix,
                                                            tolower(trunc.type.pretty) )
    
    
    if ( .estName == "Mhat" ) title.string = bquote( .(trunc.type.pretty) ~ hat(mu) )
    if ( .estName == "Shat" ) title.string = bquote( .(trunc.type.pretty) ~ hat(sigma) )
  }
  
  
  p = ggplot( data = .dp,
              aes( x = n,
                   y = Y,
                   color = method,
                   linetype = method ) ) +
    
    #geom_point() +
    geom_line(lwd = 1.2) +
    
    # manually provided colors
    myColorScale +
    
    # manually provided linetypes
    myLtyScale +
    
    # base_size controls all text sizes; default is 11
    # https://ggplot2.tidyverse.org/reference/ggtheme.html
    theme_bw(base_size = 25) +
    ggtitle(title.string) +
    ylab(.outcomeName) +
    
    guides(color = guide_legend(title = "Method"),
           linetype = guide_legend(title = "Method") ) +
    
    theme( legend.position = legend.position,
           legend.key.size = unit(2, 'cm'),  # increase legend size
           text = element_text(face = "bold") ) +
    
    # use all values of
    #scale_x_log10( breaks = unique(.dp$n) )
    # use only some values
    scale_x_log10( breaks = c(10, 20, 50, 100, 200, 500, 1000) )
  
  
  
  # ~~ Add facetting depending on trunc type ----
  # add facetting if appropriate based on number of levels of facetting variables
  if ( nuni(.dp$`Truncation type`) > 1 ){
    p = p + facet_wrap(`Truncation type` ~ cutpoints.pretty)
  } else {
    p = p + facet_wrap( ~ cutpoints.pretty)
  }
  
  #@PREVIOUS
  # if (  .trunc.type == "all" & nuni(.dp$`Truncation type`)  ) p = p + facet_wrap(`Truncation type` ~ cutpoints.pretty)
  # if (  .trunc.type != "all" ) p = p + facet_wrap( ~ cutpoints.pretty)
  
  
  # ~~ Add horizontal reference line depending on outcome ----
  # choose intelligently based on outcome
  if ( .outcomeName == "Cover" ) refHeight = 0.95 else refHeight = 0
  p = p + geom_hline( yintercept = refHeight,
                      lty = 2,
                      color = "black" ) 
  
  
  
  # ~~ Set axis breaks that depend on main text vs. Supplement ----
  # will need to adjust for specific sims
  
  # ad hoc axis scalings for certain variables
  # zoom in more for the main text 
  # y.breaks: the one we'll actually use for plotting
  # .y.breaks: user-provided one
  y.breaks = NULL  
  if ( is.null(.y.breaks) ) {
    
    
    if ( .showAllMethods == FALSE ) {
      if ( .estName == "Mhat" & .outcomeName == "Bias" ) y.breaks = seq(-1.5, 1.5, 0.5)
      if ( .estName == "Mhat" & .outcomeName == "RMSE" ) y.breaks = seq(0, 5, 0.5)
      if ( .estName == "Mhat" & .outcomeName == "EmpSE" ) y.breaks = seq(0, 4, 0.5)
      if ( .estName == "Mhat" & .outcomeName == "Width" ) y.breaks = seq(0, 16, 2)
      
      if ( .estName == "Shat" & .outcomeName == "Bias" ) y.breaks = seq(-.8, .8, .2)
      if ( .estName == "Shat" & .outcomeName == "RMSE" ) y.breaks = seq(0, 1.8, 0.2)
      if ( .estName == "Shat" & .outcomeName == "EmpSE" ) y.breaks = seq(0, 1.4, 0.2)
      if ( .estName == "Shat" & .outcomeName == "Width" ) y.breaks = seq(0, 3, 0.5)
      
    } else {
      # zoom in less for Supplement
      if ( .estName == "Mhat" & .outcomeName == "Bias" ) y.breaks = seq(-3, 3, 1)
      if ( .estName == "Mhat" & .outcomeName == "RMSE" ) y.breaks = seq(0, 3, 0.5)
      if ( .estName == "Mhat" & .outcomeName == "EmpSE" ) y.breaks = seq(0, 10, 2)
      if ( .estName == "Mhat" & .outcomeName == "Width" ) y.breaks = seq(0, 16, 2)
      
      if ( .estName == "Shat" & .outcomeName == "Bias" ) y.breaks = seq(-.8, .8, .2)
      if ( .estName == "Shat" & .outcomeName == "RMSE" ) y.breaks = seq(0, 1.25, 0.25)
      if ( .estName == "Shat" & .outcomeName == "EmpSE" ) y.breaks = seq(0, 2, 0.5)
      if ( .estName == "Shat" & .outcomeName == "Width" ) y.breaks = seq(0, 3, 0.5)
    }
    
    # other outcomes follow rules or can just use default axis breaks
    # y.breaks are only still null if none of the above applied
    if ( is.null(y.breaks) ) {
      # set default breaks
      if ( .outcomeName %in% c("Cover", "StanWarned") ){
        y.breaks = seq(0, 1, .1)
        
      } else {
        # otherwise keep the default limits from ggplot
        y.breaks = ggplot_build(p)$layout$panel_params[[1]]$y$breaks
      }
    }
    
  } 
  
  # if user provided their own y.breaks
  if ( !is.null(.y.breaks) ) {
    y.breaks = .y.breaks
  }
  
  
  
  # use coord_cartesian so that lines/points that go outside limits look cut off
  #  rather than completely omitted
  p = p + coord_cartesian( ylim = c( min(y.breaks), max(y.breaks) ) ) +
    scale_y_continuous( breaks = y.breaks )
  
  # ~~ Rename y axis if needed ----
  if ( .outcomeName == "EmpSE" ) p = p + ylab("Empirical SE")
  if ( .outcomeName == "Cover" ) p = p + ylab("95% CI coverage")
  if ( .outcomeName == "Width" ) p = p + ylab("95% CI width")
  
  
  # ~~ Write plot ----
  if ( .writePlot == TRUE ) {
    my_ggsave( name = name,
               .width = width,
               .height = height,
               .results.dir = .results.dir,
               .overleaf.dir.general = .overleaf.dir.general )
  } 
  
  p
}


# make long-format dataset for plotting
# there will be 4 of these ([Mhat vs. Shat] x [est vs. inf])
# .outcomeType: "est" or "inf"
make_dpl = function(.dp,
                    .estName,
                    .outcomeType) {
  
  # set variables for each column of plots
  # put these in the desired order for factor leveling (and hence display in plots)
  if ( .outcomeType == "est" ) {
    .dp$method = .dp$method.pretty.est
    colVars = paste( .estName,
                     c("Bias", "EmpSE", "RMSE"), sep = "" )
  }
  
  
  if ( .outcomeType == "inf" ) {
    .dp$method = .dp$method.pretty.inf
    colVars = paste( .estName, c("Cover", "Width"), sep = "" )
  }
  
  
  # make summary table
  dpw = .dp %>% 
    group_by(`Truncation type`, method) %>%
    summarise_at( .vars = colVars, 
                  list( median = function(x) median(x, na.rm = TRUE),
                        lo95 = function(x) quantile(x, probs = .025, na.rm = TRUE),
                        hi95 = function(x) quantile(x, probs = .975, na.rm = TRUE),
                        
                        lo50 = function(x) quantile(x, probs = .25, na.rm = TRUE),
                        hi50 = function(x) quantile(x, probs = .75, na.rm = TRUE)
                  ) )
  
  
  # warnings are about setting row names on a tibble
  dpl = suppressWarnings( reshape(dpw,
                                  idvar=c("Truncation type", "method"),
                                  direction="long", 
                                  
                                  timevar = "outcome",
                                  
                                  varying=list(median = names_with(.dat = dpw, "median"),
                                               lo95 = names_with(.dat = dpw, "lo95"),
                                               hi95 = names_with(.dat = dpw, "hi95"),
                                               lo50 = names_with(.dat = dpw, "lo50"),
                                               hi50 = names_with(.dat = dpw, "hi50") ),
                                  
                                  v.names = c("median", "lo95", "hi95", "lo50", "hi50") ) )
  
  
  # hacky :(
  for ( i in unique(dpl$outcome) ) {
    dpl$outcome[ dpl$outcome == i ] = colVars[i]
  }
  
  
  
  ### Force ordering of factors for plotting joy
  # order truncation types
  dpl$`Truncation type` = factor( dpl$`Truncation type`, levels = c("Single", "Symmetric double", "Asymmetric double") )
  
  # order methods (y axis)
  if ( .outcomeType == "est" ) {
    dpl$method = factor( dpl$method,
                         levels = rev( c("MLE", "MLE + boot", "Jeffreys mode", "Jeffreys mean", "Jeffreys median", "Jeffreys mode + boot") ) )
    
  } else if ( .outcomeType == "inf" ) {
    dpl$method = factor( dpl$method,
                         levels = rev( c("MLE Wald", "MLE profile", "MLE BCa", "Jeffreys posterior quantiles", "Jeffreys mode Wald", "Jeffreys profile", "Jeffreys mode BCa" ) ) )
  }
  # order outcomes
  dpl$outcome = factor( dpl$outcome,
                        levels = colVars)
  
  return(dpl)
  
}


# .dpl: a wide plotting df created by make_dpl
plot_summary = function(.dpl,
                        .results.dir,
                        .overleaf.dir.general,
                        .showPlot = FALSE) {
  
  
  
  # figure out which estimate we're dealing with
  if ( grepl( pattern = "Mhat", x = .dpl$outcome[1] ) ) {
    .estName = "Mhat"
  } else if ( grepl( pattern = "Shat", x = .dpl$outcome[1] ) ) {
    .estName = "Shat"
  }
  
  # strip the "Mhat"/"Shat" prefix from outcomes
  outcomes = str_replace( string = unique(.dpl$outcome),
                          pattern = .estName,
                          replacement = "" )
  
  # focal methods will be emphasized visually in plot
  if ( "Bias" %in% outcomes ){
    .outcomeType = "est"
    focalMethods = c("MLE", "Jeffreys mode")
  } else if ( "Cover" %in% outcomes ) {
    .outcomeType = "inf"
    focalMethods = c("MLE profile", "Jeffreys posterior quantiles")
  }
  
  
  # ~~ Make plot for each outcome in dpl ------------
  # make a plot for each column (i.e., outcome)
  plotList = list()
  nOutcomes = length(outcomes)
  
  for ( i in 1:nOutcomes ) {
    
    .outcomeName = outcomes[i]
    xName = paste(.estName, .outcomeName, sep = "")
    
    
    temp = droplevels(.dpl %>% filter(outcome == xName) )
    
    
    # ~~ Axis scaling -------
    if (xName == "MhatBias") x.breaks = seq(-2, 1.5, .5)
    if (xName == "MhatEmpSE") x.breaks = seq(0, 6, 1)
    if (xName == "MhatRMSE") x.breaks = seq(0, 7, 1)
    
    if (xName == "ShatBias") x.breaks = seq(-0.8, 0.8, .2)
    if (xName == "ShatEmpSE") x.breaks = seq(0, 2, 0.5)
    if (xName == "ShatRMSE") x.breaks = seq(0, 1.25, 0.25)
    
    if (xName == "MhatWidth") x.breaks = seq(0, 10, 2)
    
    if (xName == "ShatWidth") x.breaks = seq(0, 3, 0.5)
    
    if ( xName %in% c("MhatCover",
                      "ShatCover",
                      "MhatEstFail",
                      "ShatEstFail",
                      "MhatCIFail",
                      "ShatCIFail") ) x.breaks = seq(0, 1, .2)
    
    
    # ~~ Base plot -------
    p = ggplot( data = temp, 
                aes(x = median,
                    y = method,
                    shape = method %in% focalMethods,
                    alpha = method %in% focalMethods)
    ) +
      
      geom_errorbarh( aes( xmin = lo95, xmax = hi95 ),
                      height = 0) +
      
      geom_errorbarh( aes( xmin = lo50, xmax = hi50 ),
                      height = 0,
                      lwd = 1.2,
                      color = "red") +
      
      geom_point(size = 2.6) +
      
      # emphasize the focal methods
      scale_shape_manual( values = c(1, 19) ) +
      # no using alpha right now, but it works
      scale_alpha_manual( values = c(1, 1) ) +
      
      xlab(.outcomeName) +
      ylab("") +
      
      coord_cartesian( xlim = c( min(x.breaks), max(x.breaks) ) ) +
      scale_x_continuous( breaks = x.breaks ) +
      
      facet_wrap( ~ `Truncation type`, ncol = 1, scales = "fixed" ) +
      
      theme_bw(base_size = 15) +
      theme(text = element_text(face = "bold"),
            legend.position = "none")
    
    # ~~~ Customize plot  -------
    # only include y tick marks on what will be the leftmost plot
    if ( !( .outcomeName %in% c("Bias", "Cover") ) ){
      p = p + theme( axis.text.y = element_blank() )
    }
    
    # ~~ Add horizontal reference line depending on outcome ----
    # choose intelligently based on outcome
    if ( .outcomeName == "Cover" ) ref.x = 0.95 else ref.x = 0
    p = p + geom_vline( xintercept = ref.x,
                        lty = 2,
                        color = "black" ) 
    
    # ~~ Rename y axis if needed ----
    if ( .outcomeName == "EmpSE" ) p = p + xlab("Empirical SE")
    if ( .outcomeName == "Cover" ) p = p + xlab("95% CI coverage")
    if ( .outcomeName == "Width" ) p = p + xlab("95% CI width")
    if ( .outcomeName %in% c("EstFail", "CIFail") ) p = p + xlab("Estimation failed")
    
    plotList[[i]] = p
  }
  
  # ~~ Nicely arrange plots as columns ------------
  
  # give extra space to leftmost one to accommodate y-axis labels
  # to check the width of these, export pCombined at dimensions 7 x 12.3
  rel.width = 1
  if ( nOutcomes == 4 ) rel.width = 1.6
  if ( nOutcomes == 3 ) rel.width = 1.5
  if ( nOutcomes == 2 ) rel.width = 1.8
  
  pCombined = cowplot::plot_grid(plotlist = plotList,
                                 #align = "v",
                                 nrow = 1,
                                 
                                 rel_widths = c(rel.width, rep(1, nOutcomes - 1) ) )
  
  pCombined
  
  # ~~ Save plot ------------
  
  name = paste( tolower(.estName),
                "_",
                .outcomeType,
                "_summaryplot.pdf",
                sep = "" )
  
  widthCoef = .6  # controls aspect ratio (larger = more flat; smaller = more tall)
  height = 7
  
  my_ggsave( name = name,
             .width = nOutcomes*height*widthCoef,
             .height = height,
             .results.dir = .results.dir,
             .overleaf.dir.general = .overleaf.dir.general )
  
  if (.showPlot == TRUE) pCombined
  
  
}

# make a quick summary table showing percent (if .varIsProp == TRUE) or mean (o.w.) of .var for each method
# .prefix: included in csv file name
quick_method_table = function(.agg = agg,
                              .var,
                              .varIsProp = FALSE,
                              .prefix = "",
                              .writeTable = FALSE ) {
  
  .agg$tempVar = .agg[[.var]]
  
  if ( .varIsProp == TRUE ) {
    t = .agg %>% group_by(method) %>%
      # "curly-curly" operator
      summarise( perc = round_prop_to_perc( mean(tempVar) ) )
    
  } else {
    t = .agg %>% group_by(method) %>%
      # "curly-curly" operator
      summarise( mean = my_round( mean(tempVar)  ) )
  }
  
  if ( .writeTable == TRUE ) {
    setwd(results.dir)
    setwd("Tables from R")
    fwrite(t, paste( .var, .prefix, "by_method.csv", sep = "_" ) )
  }
  
  return(t)
}


# makes wide table for a given outcome to see which method won in each scenario
quick_method_wide_table = function(.agg = agg,
                                   .var,
                                   .prefix = "",
                                   .writeTable = FALSE) {
  
  .agg$tempVar = .agg[[.var]]
  
  t = .agg %>% pivot_wider( id_cols = scen.name,
                            names_from = method,
                            values_from = tempVar )
  
  if ( .writeTable == TRUE ) {
    setwd(results.dir)
    setwd("Tables from R")
    fwrite(t, paste( .var, .prefix, "wide_by_method.csv", sep = "_" ) )
  }
  
  return(t)
  
}


# ~ SMALL MISC FNS ----------------------------------------------

# stands for "wipe results"
wr = function(){
  setwd(results.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  setwd(overleaf.dir.general)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
}

# stands for "view results"
vr = function(){
  setwd(results.dir)
  View( read.csv("stats_for_paper.csv") )
}



# round while keeping trailing zeroes
my_round = function(x, digits = 2) {
  formatC( round( x, digits ), format='f', digits=digits )
}

# take a proportion and present as a rounded percent
round_prop_to_perc = function(x, digits = 0) {
  formatC( round( 100*x, digits ), format='f', digits=digits )
}



# format a p-value with scientific notation stars for cutoffs
# star.cutoffs: cutoffs for *, **, ***, etc., provided in any order
format_pval = function( p,
                        digits = 3,
                        star.cutoffs = NA ) {
  
  if (p >= 0.01) string = as.character( my_round( p, digits ) )
  if (p < 0.01 & p > 10^-5 ) string = formatC( p, format = "E", digits = 2 )
  if ( p < 10^-5 ) string = "< 1E-05"
  
  if ( ! is.na(star.cutoffs[1]) ) {
    
    # put in descending order
    star.cutoffs = sort(star.cutoffs, decreasing = TRUE)
    
    for ( i in 1 : length(star.cutoffs) ) {
      if ( p < star.cutoffs[i] ) string = paste( string, "*", sep="" )
    }
  }
  
  return(string)
}

# example
# p = seq( 0, .2, 0.001 )
# vapply( p, format_pval, "asdf" )
# vapply( p, function(x) format_pval( x, star.cutoffs = c( 0.01, 0.05) ), "asdf" )

# make a string for estimate and CI
stat_CI = function(est, lo, hi, digits=2){
  paste( round(est, digits), " [", round(lo, digits), ", ", round(hi, digits), "]", sep = "" )
}
stat_CI( c(.5, -.1), c(.3, -.2), c(.7, .0) )

# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
# expects "study" to be a global var
update_result_csv = function( name,
                              .section = section,
                              value = NA,
                              print = FALSE ) {
  setwd(results.dir)
  
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(.section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    res.overleaf <<- read.csv( "stats_for_paper.csv",
                               stringsAsFactors = FALSE,
                               colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% res.overleaf$name) ) res.overleaf[ res.overleaf$name %in% name, ] <<- new.rows
    else res.overleaf <<- rbind(res.overleaf, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    res.overleaf <<- new.rows
  }
  
  # write to results directory
  # MUST use quote=FALSE or else you will get very cryptic errors from rlgetnum in Overleaf, even though stats_for_paper.csv will look fine
  write.csv( res.overleaf, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE)
  
  # also write to Overleaf
  setwd(overleaf.dir.general)
  write.csv( res.overleaf, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  if ( print == TRUE ) {
    View(res.overleaf)
  }
}



# one or both dirs can be NA
my_ggsave = function(name,
                     .width,
                     .height,
                     .results.dir = figures.results.dir,
                     .overleaf.dir.general = overleaf.dir.general) {
  
  dirs = c(.results.dir, .overleaf.dir.general)
  dirIsNA = sapply(dirs, is.na)
  validDirs = dirs[ !dirIsNA ]
  
  
  for ( dir in validDirs ) {
    setwd(dir)
    ggsave( name,
            width = .width,
            height = .height,
            device = "pdf" )
  }
}

#bm
# generate figure strings for whatever is the current working dir
# expect global vars: estName, outcomeName
overleaf_figure_strings = function() {
  #setwd(figures.results.dir)
  
  # # order numbers correctly so that, e.g., 20 comes before 100
  # # https://stackoverflow.com/questions/10777367/how-can-i-read-the-files-in-a-directory-in-sorted-order-using-r
  # library(gtools)
  # files = mixedsort(list.files())
  
  # estimation plots
  for ( .t in c("single", "double-symm", "double-asymm") ) {
    for ( .e in estNames ) {
      for ( .o in c("Bias", "EmpSE", "RMSE") ) {
        
        fig.string = paste( tolower(.t),
                            "_",
                            tolower(.e),
                            "_n_",
                            tolower(.o),
                            "_plot.pdf",
                            sep = "" )
        
        # \figcoef should be def'd in the Overleaf document to control all figures
        width.string = "\\figcoef\\textwidth"
        caption.string = paste("\\caption{\\label{sfig:",
                               .t,
                               "_",
                               tolower(.e),
                               "_",
                               tolower(.o),
                               "}}",
                               sep = "" )
        
        cat( "\\begin{figure}[H] \\centering \\includegraphics[width=",
             width.string,
             "]{R_objects/figures/",
             fig.string,
             "}",
             caption.string,
             "\\end{figure}",
             sep = "" )
        
        cat("\n\n")
        
      }
    }
  }
  
  
  # inference plots
  for ( .t in c("single", "double-symm", "double-asymm") ) {
    for ( .e in estNames ) {
      for ( .o in c("Cover", "Width") ) {
        
        fig.string = paste( tolower(.t),
                            "_",
                            tolower(.e),
                            "_n_",
                            tolower(.o),
                            "_plot.pdf",
                            sep = "" )
        
        # \figcoef should be def'd in the Overleaf document to control all figures
        width.string = "\\figcoef\\textwidth"
        caption.string = paste("\\caption{\\label{sfig:",
                               .t,
                               "_",
                               tolower(.e),
                               "_",
                               tolower(.o),
                               "}}",
                               sep = "" )
        
        cat( "\\begin{figure}[H] \\centering \\includegraphics[width=",
             width.string,
             "]{R_objects/figures/",
             fig.string,
             "}",
             caption.string,
             "\\end{figure}",
             sep = "" )
        
        cat("\n\n")
        
      }
    }
  }
  
  
  
  
  
}


# xtable docs:
# "p{3cm}" etc. for a LaTeX column of the specified width.
prettify_my_summarise = function(.tab,
                                 .label,
                                 .caption = "") {
  
  
  
  # omit the number of scenarios, but record it in the results csv
  if ( nuni(.tab$n.scens) == 1 ) {
    update_result_csv( name = paste( .label, "scens per row" ),
                       .section = section,
                       value = unique(.tab$n.scens) )
    
    .tab = .tab[ , names(.tab) != "n.scens" ]
  }
  
  
  # decide where to put horizontal rules
  # must be done before editing the trunc type variable below
  n.methods = nuni(.tab$method)
  hline.after = c(-1, 0, n.methods*seq( 1 : nuni(.tab$`Truncation type`) ) )
  
  # don't show the trunc type repeatedly
  .tab$`Truncation type`[ duplicated(.tab$`Truncation type`) ] = ""
  
  
  temp = xtable(.tab,
                booktabs = TRUE,
                caption = .caption)
  
  
  # prettify trunc types
  bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
  
  print(temp,
        include.rownames = FALSE,
        table.placement = "H",
        hline.after = hline.after,
        sanitize.colnames.function=bold )
  
}




# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}

# quick median with NAs removed
medNA = function(x){
  median(x, na.rm = TRUE)
}


# quick median with NAs removed and 10th and 90th percentiles
medNA_pctiles = function(x){
  paste( round( median(x, na.rm = TRUE), 2 ),
         " (",
         round( quantile(x, probs = 0.10, na.rm = TRUE), 2 ),
         ", ",
         round( quantile(x, probs = 0.90, na.rm = TRUE), 2 ),
         ")",
         sep = "" )
}


# take the logit of a probability, but truncate
#  to avoid infinities
truncLogit <- function(p) {
  p[p==0] = 0.001
  p[p==1] = 0.999
  log(p/(1-p))
}

expit = function(x) {
  exp(x) / (1 + exp(x))
}

# calculate I^2 from t^2 and N
I2 = function(t2, N) {
  t2 / (t2 + 4/N)
}

# check CI coverage
covers = function( truth, lo, hi ) {
  return( (lo <= truth) & (hi >= truth) )
}

# get names of dataframe containing a string
namesWith = function(pattern, dat){
  names(dat)[ grepl(pattern = pattern, x = names(dat) ) ]
}

# quick length(unique)
nuni = function(x) {
  length(unique(x))
}


# given the truncation type ("single", "double-symm", or "double-asymm"),
#  and the percentage of distribution that is to be retained,
#  calculate (Za, Zb)
calculate_cutpoints = function( .trunc.type,
                                .prop.retained,
                                .dist = "norm",
                                .dist.t.df = NULL ){
  
  n.cuts = length(.prop.retained)
  
  if ( .dist == "t" & is.null(.dist.t.df) ) stop("Need to provide .dist.t.df if choosing .dist='t'")
  
  ### single truncation
  if (.trunc.type == "single") {
    if (.dist == "norm") Za = qnorm(p = .prop.retained, lower.tail = FALSE)
    if (.dist == "t") Za = qt(p = .prop.retained, df = .dist.t.df, lower.tail = FALSE)
    Zb = rep(Inf, n.cuts)
    
    ### symmetric double truncation
  } else if (.trunc.type == "double-symm") {
    if (.dist == "norm") Zb = qnorm( p = (1 - .prop.retained )/2, lower.tail = FALSE )
    if (.dist == "t") Zb = qt( p = (1 - .prop.retained )/2, df = .dist.t.df, lower.tail = FALSE )
    Za = -Zb
    
    ### asymmetric double truncation
  } else if (.trunc.type == "double-asymm") {
    Za = -0.50 # fixed to a "middle" value compared to other scenarios
    
    # check if it's even possible to get the desired prop.retained 
    if (.dist == "norm") above.Za = 1 - pnorm(Za)
    if (.dist == "t") above.Za = 1 - pt(Za, df = .dist.t.df)
    
    if ( any(above.Za - .prop.retained < 0) ) {
      warning("Some .prop.retained values are so large that it's impossible to have asymmetric double-truncation with the hard-coded Za for that truncation type")
    }
    if (.dist == "norm") Zb = suppressWarnings( qnorm( above.Za - .prop.retained,
                                                       lower.tail = FALSE ) )
    if (.dist == "t") Zb = suppressWarnings( qt( above.Za - .prop.retained,
                                                 df = .dist.t.df,
                                                 lower.tail = FALSE ) )
  }
  
  .dat = data.frame( Za, Zb )
  
  # sanity check: are we actually retained the desired proportion?
  if (.dist == "norm") {
    .dat = .dat %>% rowwise() %>%
      mutate( actual.prop.retained = pnorm(Zb) - pnorm(Za) )
  } else if (.dist == "t") {
    .dat = .dat %>% rowwise() %>%
      mutate( actual.prop.retained = pt(Zb, df = .dist.t.df) - pt(Za, df = .dist.t.df) )
  }
  discrep = .dat$actual.prop.retained - .prop.retained
  discrep = discrep[ !is.na(discrep) ]  # handle NAs from double-asymm
  expect_equal( discrep, rep( 0, length(discrep) ), tol = 0.001 )
  
  # return results
  return(.dat)
  
}


# calculate_cutpoints( .trunc.type = "single",
#                      .prop.retained = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9) )
# 
# # should be about the same as above (t with huge df)
# calculate_cutpoints( .trunc.type = "single",
#                      .prop.retained = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9),
#                      .dist = "t",
#                      .dist.t.df = 1000)
# # will be more extreme than normal (t with small df)
# calculate_cutpoints( .trunc.type = "single",
#                      .prop.retained = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9),
#                      .dist = "t",
#                      .dist.t.df = 5)
# 
# calculate_cutpoints( .trunc.type = "double-symm",
#                      .prop.retained = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9) )
# 
# calculate_cutpoints( .trunc.type = "double-asymm",
#                      .prop.retained = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9) )



names_with = function(.dat, .pattern) {
  names(.dat)[ grepl(pattern = .pattern, x = names(.dat) ) ]
}


# generate from truncated noncentral t by specifying underlying mean and variance
draw_truncated_t = function(.n, 
                            .df,
                            .mean,
                            .sd,
                            .a,
                            .b) {
  
  # init truncated and untruncated data vector
  xt = c(NA)
  xu = c(NA)
  n.draws = 0
  
  # draw each observation
  while ( length(xt) < (.n+1) ) {
    
    # first generate one observation from untruncated central t
    new.xc = rt( n = 1, df = .df )
    n.draws = n.draws + 1
    
    # change the observation's location and scale in untruncated data
    # https://stackoverflow.com/questions/17843497/sampling-from-a-t-distribution-in-r
    new.xu = new.xc * sqrt( .sd^2 * (.df-2)/.df ) + .mean
    xu = c(xu, new.xu)
    
    # check if it falls within truncation points
    if ( new.xu > .a & new.xu < .b ) xt = c(xt, new.xu)
    
  }
  xt = xt[-1]  # remove the NA at the beginning
  xu = xu[-1]  
  return( list(xt = xt,
               empirical.prop.retained = length(xt)/n.draws,
               empirical.mean.untruncated = mean(xu),
               empirical.sd.untruncated = sd(xu) ) )
}



# res = draw_truncated_t(.n = 1000,
#                        .df = 5,
#                        .mean = 0,
#                        .sd = 1,
#                        .a = 1.143215,
#                        .b = Inf)
# res$empirical.prop.retained  # expect close to 10% because I got this cutpoint from doParallel
# res$empirical.mean.untruncated
# res$empirical.sd.untruncated
# 
# # sanity check
# # yes, untruncated moments match regardless of df
# res = draw_truncated_t(.n = 1000,
#                       .df = 5,
#                       .mean = 1.1,
#                       .sd = 2.3,
#                       .a = 2.250481,
#                       .b = Inf)
# res$empirical.prop.retained
# res$empirical.mean.untruncated
# res$empirical.sd.untruncated                 


# # sanity check
# # if huge n and not truncated, should have the desired mean and sd
# xt = draw_truncated_t(.n = 1000,
#                       .df = 5,
#                       .mean = 0,
#                       .sd = 1,
#                       .a = -Inf,
#                       .b = Inf)$xt
# mean(xt); 0
# sd(xt); 1

# ~ GENERIC FNS FOR INTERACTING WITH CLUSTER --------------------------------------------

# looks at results files to identify sbatches that didn't write a file
# .max.sbatch.num: If not passed, defaults to largest number in actually run jobs.

sbatch_not_run = function(.results.singles.path,
                          .results.write.path,
                          .name.prefix,
                          .max.sbatch.num = NA ) {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # extract job numbers
  sbatch.nums = as.numeric( unlist( lapply( strsplit( keepers, split = "_"), FUN = function(x) x[5] ) ) )
  
  # check for missed jobs before the max one
  if ( is.na(.max.sbatch.num) ) .max.sbatch.num = max(sbatch.nums)
  all.nums = 1 : .max.sbatch.num
  missed.nums = all.nums[ !all.nums %in% sbatch.nums ]
  
  # give info
  print( paste("The max job number is: ", max(sbatch.nums) ) )
  print( paste( "Number of jobs that weren't run: ",
                ifelse( length(missed.nums) > 0, length(missed.nums), "none" ) ) )
  
  if( length(missed.nums) > 0 ) {
    setwd(.results.write.path)
    write.csv(missed.nums, "missed_job_nums.csv")
  }
  
  return(missed.nums)
  
}

# missed.nums = sbatch_not_run( "/home/groups/manishad/multTest/sim_results/short",
#                 "/home/groups/manishad/multTest/sim_results",
#                 .name.prefix = "short_results" )
# scp mmathur@sherlock:/share/PI/manishad/multTest/sim_results/missed_job_nums.csv ~/Desktop



# These just generate the sbatch files
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

ml load v8
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




# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path,
                        .results.stitched.write.path=.results.singles.path,
                        .name.prefix,
                        .stitch.file.name="stitched_model_fit_results.csv") {
  
  # .results.singles.path = "/home/groups/manishad/TNE/sim_results/long"
  # .results.stitched.write.path = "/home/groups/manishad/TNE/sim_results/overall_stitched"
  # .name.prefix = "long_results"
  # .stitch.file.name="stitched.csv"
  
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # read in and rbind the keepers
  tables = lapply( keepers, function(x) fread(x, header=TRUE)[ , -1] )
  
  # check for name mismatches
  allNames = lapply( tables, function(x) names(x) )
  firstNames = allNames[[1]]
  namesMatchFirst = unlist( lapply( allNames, function(x) all(x == firstNames) ) )
  
  # look for bad variable names
  bad = which( namesMatchFirst == FALSE )
  
  if ( length(bad) > 0 ){
    browser()
    badExampleNames = allNames[[ bad[1] ]]
    firstNames[ !(firstNames %in% badExampleNames) ]
  }
  
  # bind them
  # rbindlist is much faster than rbind or even bind_rows:
  #  https://rstudio-pubs-static.s3.amazonaws.com/406521_7fc7b6c1dc374e9b8860e15a699d8bb0.html
  s = rbindlist(tables)
  
  names(s) = firstNames
  
  if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
  fwrite(s, paste(.results.stitched.write.path,
                  .stitch.file.name,
                  sep="/") )
  return(s)
}



# ~ GENERIC TWEAK TO BOOT'S BOOTSTRAPPING FN ----------------------------------- 

# see section called "MM additions"
# I minimally modified the function so that it can proceed even if some of the bootstrap iterates run into errors
# (in this case, Fisher convergence issues) because the boot package version gets confused about dimension mismatches

# source internal boot package functions
# source("bootfuns.R")

my_boot = function (data, statistic, R, sim = "ordinary", stype = c("i", 
                                                                    "f", "w"), strata = rep(1, n), L = NULL, m = 0, weights = NULL, 
                    ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ..., 
                    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
                                                                               1L), cl = NULL) 
{
  
  
  call <- match.call()
  stype <- match.arg(stype)
  if (missing(parallel)) 
    parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
    warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0', so ignored")
    simple <- FALSE
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  n <- NROW(data)
  if ((n == 0) || is.null(n)) 
    stop("no data in call to 'boot'")
  temp.str <- strata
  strata <- tapply(seq_len(n), as.numeric(strata))
  t0 <- if (sim != "parametric") {
    if ((sim == "antithetic") && is.null(L)) 
      L <- empinf(data = data, statistic = statistic, stype = stype, 
                  strata = strata, ...)
    if (sim != "ordinary") 
      m <- 0
    else if (any(m < 0)) 
      stop("negative value of 'm' supplied")
    if ((length(m) != 1L) && (length(m) != length(table(strata)))) 
      stop("length of 'm' incompatible with 'strata'")
    if ((sim == "ordinary") || (sim == "balanced")) {
      if (isMatrix(weights) && (nrow(weights) != length(R))) 
        stop("dimensions of 'R' and 'weights' do not match")
    }
    else weights <- NULL
    if (!is.null(weights)) 
      weights <- t(apply(matrix(weights, n, length(R), 
                                byrow = TRUE), 2L, normalize, strata))
    if (!simple) 
      i <- index.array(n, R, sim, strata, m, L, weights)
    original <- if (stype == "f") 
      rep(1, n)
    else if (stype == "w") {
      ns <- tabulate(strata)[strata]
      1/ns
    }
    else seq_len(n)
    t0 <- if (sum(m) > 0L) 
      statistic(data, original, rep(1, sum(m)), ...)
    else statistic(data, original, ...)
    rm(original)
    t0
  }
  else statistic(data, ...)
  pred.i <- NULL
  fn <- if (sim == "parametric") {
    ran.gen
    data
    mle
    function(r) {
      dd <- ran.gen(data, mle)
      statistic(dd, ...)
    }
  }
  else {
    if (!simple && ncol(i) > n) {
      pred.i <- as.matrix(i[, (n + 1L):ncol(i)])
      i <- i[, seq_len(n)]
    }
    if (stype %in% c("f", "w")) {
      f <- freq.array(i)
      rm(i)
      if (stype == "w") 
        f <- f/ns
      if (sum(m) == 0L) 
        function(r) statistic(data, f[r, ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, 
      ], ...)
    }
    else if (sum(m) > 0L) 
      function(r) statistic(data, i[r, ], pred.i[r, ], 
                            ...)
    else if (simple) 
      function(r) statistic(data, index.array(n, 1, sim, 
                                              strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    }
    else if (have_snow) {
      list(...)
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                             ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else lapply(seq_len(RR), fn)
  #t.star <- matrix(, RR, length(t0))  # ~~~ MM commented out
  
  
  # ~~~~~ MM added
  # number of non-NULL elements of the results vector
  #RR = length(unlist(res))
  # changed to accommodate multi-argument returns:
  RR = length(res)
  nulls = sapply( res, is.null)
  res = res[ !nulls ]
  t.star <- matrix(, RR, length(t0))
  # if failed reps return NAs, then R = original RR still
  # without this, boot.CI gets confused about number of replicates
  R = RR
  
  #mean(is.na(t.star))
  # ~~~~~ end of MM additions
  
  
  for (r in seq_len(RR)) t.star[r, ] <- res[[r]]
  if (is.null(weights)) 
    weights <- 1/tabulate(strata)[strata]
  boot.return(sim, t0, t.star, temp.str, R, data, statistic, 
              stype, call, seed, L, m, pred.i, weights, ran.gen, mle)
}




