

# RTMA log-likelihood - now uses TNE version
# carefully structured for use with Deriv()
joint_nll_2 = function(.yi, .sei, .Mu, .Tt, .crit = qnorm(.975)) {
  
  .T2t = .Tt^2
  # as in TNE::nll() instead
  .dat = data.frame(yi = .yi, sei = .sei)
  
  .dat = .dat %>% rowwise() %>%
    mutate( term1 = dmvnorm(x = as.matrix(yi, nrow = 1),
                            mean = as.matrix(.Mu, nrow = 1),
                            sigma = as.matrix(.T2t + sei^2, nrow=1),
                            log = TRUE),
            
            term2 = log( pmvnorm( lower = -99,
                                  upper = .crit * sei,
                                  mean = .Mu,
                                  sigma = .T2t + sei^2 ) ) )
  
  
  -sum(.dat$term1) - sum(.dat$term2)
  
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


E_fisher_RTMA = function( .sei, .Mu, .Tt, .crit = qnorm(0.975) ) {
  # get expected Fisher info for each observation separately, based on its unique SE
  # each observation is RTN, so can just use TNE result!!
  Efish.list = lapply( X = as.list(.sei),
                       FUN = function(.s) {
                         E_fisher_TNE( .mu = .Mu,
                                       .sigma = sqrt(.Tt^2 + .s^2), 
                                       .n = 1,
                                       .a = -99,
                                       .b = .crit*.s )
                       })
  
  # add all the matrices entrywise
  # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
  Efish.all = Reduce('+', Efish.list) 
  
  return(Efish.all)
}


lprior = function(.sei, .Mu, .Tt, .crit) {
  Efish = E_fisher_RTMA( .sei = .sei, .Mu = .Mu, .Tt = .Tt, .crit = .crit )
  log( sqrt( det(Efish) ) )
}


# .pars: (.Mu, .Tt) or (.Mu, .Tt2)
nlpost_jeffreys_RTMA = function( .pars,
                                 .par2is = "Tt",  # "Tt" or "Tt2"
                                 .yi,
                                 .sei,
                                 .crit = qnorm(.975),
                                 
                                 # if .usePrior = FALSE, will just be the MLE
                                 .usePrior = TRUE) {
  
  # variance parameterization
  if (.par2is == "Tt2") stop("Var parameterization not handled yet")
  
  if (.par2is == "Tt") {
    Mu = .pars[1]
    Tt = .pars[2]
    
    if ( Tt < 0 ) return(.Machine$integer.max)
    
    # negative log-likelihood
    # joint_nll_2 uses the TNE version
    nll.value = joint_nll_2( .yi = .yi,
                             .sei = .sei,
                             .Mu = Mu,
                             .Tt = Tt,
                             .crit = .crit )
    
    # log-prior
    # lprior uses the TNE expected Fisher and then just sums over observations
    if ( .usePrior == TRUE ) {
      prior.value = lprior(.sei = .sei,
                           .Mu = Mu,
                           .Tt = Tt,
                           .crit = .crit)
    } else {
      prior.value = 0
    }

    # negative log-posterior
    nlp.value = sum(nll.value) + prior.value
    
    if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  }
  
  return(nlp.value)
}


# 2021-9-2: MM audited fn by reading
# mu.start, sigma.start: start values for optimatization
# as illustrated in a sanity check after nlpost_simple, this fn's MAPs agree with
#  using mle() directly on nlpost_Jeffreys
estimate_jeffreys_RTMA = function( yi,
                                   sei,
                                   par2is = "Tt",
                                   Mu.start,
                                   Tt.start,
                                   crit,
                                   
                                   usePrior = TRUE,
                                   get.CIs,
                                   CI.method = "wald" ) {
  
  ### Get MAP by Calling mle() ###
  # IMPORTANT: This fn cannot be moved outside the scope of estimate_jeffreys
  #  because mle() is too dumb to allow extra args (e.g., x) to be passed,
  #  so it's forced to rely on global vars
  #  and that's a problem with a doParallel loop
  #  if this fn is outside estimate_jeffreys, different parallel iterations will use each other's global vars
  
  #  expects yi, sei, and crit to be global vars
  nlpost_simple_RTMA = function(.Mu, .Tt) {
    
    nlpost_jeffreys_RTMA( .pars = c(.Mu, .Tt),
                          .par2is = "Tt",
                          .yi = yi,
                          .sei = sei,
                          .crit = crit,
                          .usePrior = usePrior )
  }

  #**important: force use of Nelder-Mead optimization, which works better for Jeffreys
  #  (even though BFGS works better for MLE)
  # for more on this issue, see "2021-9-23 SD vs. var reduc with Jeffreys.R"
  res = mle( minuslogl = nlpost_simple_RTMA,
             start = list( .Mu = Mu.start, .Tt = Tt.start),
             method = "Nelder-Mead" )
  
  
  # not actually MLEs, of course, but rather MAPs
  mles = as.numeric(coef(res))
  
  # 2021-11-19: try another optimizer (BFGS)
  myMLE.bfgs = mle( minuslogl = nlpost_simple_RTMA,
                    start = list( .Mu = Mu.start, .Tt = Tt.start),
                    method = "BFGS" )
  mles.bfgs = as.numeric( coef(myMLE.bfgs) )
  
  # THIS BEHAVES WELL
  if ( par2is == "Tt" ) {
    # need this structure for run_method_safe to understand
    Mu.hat = mles[1]
    Tt.hat = mles[2]
  }
  
  # from TNE
  # # THIS BEHAVES BADLY
  # if ( par2is == "var" ) {
  #   # need this structure for run_method_safe to understand
  #   Mhat = mles[1]
  #   Vhat = mles[2]
  #   Shat = sqrt(mles[2])
  # }
  
  # recode convergence more intuitively
  # optim uses "0" to mean successful convergence
  optim.converged = attr(res, "details")$convergence == 0 
  
  profile.CI.error = NA
  
  
  # ### Get CIs ###
  # if ( get.CIs == TRUE & CI.method == "wald" & par2is == "Tt" ) {
  #   # get Wald CI 
  #   # SEs for both parameters
  #   # this has its own tryCatch in case point estimation was possible, but not inference
  #   tryCatch({
  #     
  #     # IMPORTANT: Despite the name "Ofish", this is actually the Hessian of the nlpost, which incorporates the prior (see Bayesian Data Analysis, page 84)
  #     # not the Fisher, which would be from the lkl only
  #     Ofish = attr(res, "details")$hessian
  #     invFisher = solve(Ofish)
  #     
  #     # SEs for Mhat and Shat, leaving blank space for Vhat
  #     SEs = sqrt( c( invFisher[1,1], NA, invFisher[2,2] ) )
  #     # fill in VhatSE using delta method
  #     # let g(y) = y^2, where y=Shat
  #     SEs[2] = SEs[3] * 2 * Shat
  #     
  #     los = ests - SEs * qnorm(0.975)
  #     his = ests + SEs * qnorm(0.975)
  #     
  #   }, error = function(err) {
  #     SEs <<- los <<- his <<- rep(NA, 3)
  #     profile.CI.error <<- err$message
  #   })
  #   
  # } else if ( get.CIs == TRUE & CI.method == "profile" & par2is == "sd" ) {
  #   
  #   tryCatch({
  #     # as confirmed in "2021-8-19 Investigate profile penalized LRT inference",
  #     #  these are indeed profile CIs
  #     CIs = confint(res)
  #     # NA's here represent Vhat
  #     los = c( as.numeric( CIs[1,1] ), as.numeric( CIs[2,1] )^2, as.numeric( CIs[2,1] ) )
  #     his = c( as.numeric( CIs[1,2] ), as.numeric( CIs[2,2] )^2, as.numeric( CIs[2,2] ) )
  #     
  #     ( res.SEs = as.numeric( attr( summary(res), "coef" )[, "Std. Error" ] ) )
  #     SEs = c(res.SEs[1], NA, res.SEs[2])
  #     
  #   }, error = function(err) {
  #     SEs <<- los <<- his <<- rep(NA, 3)
  #     profile.CI.error <<- err$message
  #   })
  #   
  #   
  # } else if ( get.CIs == TRUE & CI.method == "wald" & par2is == "var" ) {
  #   # get Wald CI 
  #   # SEs for both parameters
  #   
  #   # this has its own tryCatch in case point estimation was possible, but not inference
  #   tryCatch({
  #     # as noted above, this is actually the Hessian of the nlpost, not the Fisher info
  #     Ofish = attr(res, "details")$hessian
  #     invFisher = solve(Ofish)
  #     
  #     # SEs for Mhat and Shat, leaving blank space for Shat
  #     SEs = sqrt( c( invFisher[1,1], invFisher[2,2], NA ) )
  #     # fill in ShatSE using delta method
  #     # let g(y) = y^(1/2), where y=Shat
  #     #  so g'(y) = 0.5 * y^(-0.5)
  #     SEs[3] = SEs[2] * 0.5*Vhat^(-0.5)
  #     
  #     los = ests - SEs * qnorm(0.975)
  #     his = ests + SEs * qnorm(0.975)
  #     
  #   }, error = function(err) {
  #     SEs <<- los <<- his <<- rep(NA, 3)
  #   })
  #   
  # } else if ( get.CIs == TRUE & CI.method == "profile" & par2is == "var" ) {
  #   
  #   tryCatch({
  #     # as confirmed in "2021-8-19 Investigate profile penalized LRT inference",
  #     #  these are indeed profile CIs
  #     CIs = confint(res)
  #     # NA's here represent Shat
  #     los = c( as.numeric( CIs[1,1] ), as.numeric( CIs[2,1] ), sqrt( as.numeric( CIs[2,1] ) ) )
  #     his = c( as.numeric( CIs[1,2] ), as.numeric( CIs[2,2] ), sqrt( as.numeric( CIs[2,2] ) ) )
  #     
  #     
  #     ( res.SEs = as.numeric( attr( summary(res), "coef" )[, "Std. Error" ] ) )
  #     SEs = c(res.SEs[1], res.SEs[2], NA)
  #     SEs[3] = SEs[2] * 0.5*Vhat^(-0.5)
  #     
  #   }, error = function(err) {
  #     SEs <<- los <<- his <<- rep(NA, 3)
  #   })
  #   
  # } else {  # i.e., if get.CIs == FALSE 
  #   SEs = los = his = c(NA, NA, NA)
  # }
  
  #@temp only
  SEs = los = his = c(NA, NA, NA)
  profile.CI.error = NA
  
  return( list( MuHat = Mu.hat, 
                TtHat = Tt.hat,
                
                MuHatSE = SEs[1],
                TtHatSE = SEs[2],
                
                Mu.CI = as.numeric( c(los[1], his[1]) ),
                Tt.CI = as.numeric( c(los[2], his[2]) ),
                
                optim.converged = optim.converged,
                # to match output of estimate_mle, return BFGS - NM
                Mhat.opt.diff = mles.bfgs[1] - Mu.hat, 
                profile.CI.error = profile.CI.error ) )
}


