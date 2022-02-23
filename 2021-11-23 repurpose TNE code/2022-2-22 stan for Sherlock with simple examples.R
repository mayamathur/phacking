

# to update helper_SAPH.R, etc., to cluster
scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/helper_SAPH.R mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH


library(rstan)
library(tidyverse)
library(dplyr)
library(tibble)
# data-wrangling packages
library(here)
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(xtable)
library(testthat)
library(Deriv)
library(mosaic)
library(hpa)
library(pracma)
library(truncnorm)
library(tmvtnorm)
library(Hmisc)
library(truncdist)
library(weightr)

path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

# to do:
#  - can get rid of mustarL and LL; they are just placeholders now
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


options(mc.cores = parallel::detectCores())



# ~ EXAMPLES TO RUN ON CLUSTER -------------------------------------------

# ~~ Example 1: Simple fake data (not from sim_meta) -------------------------

.Mu.start = 0
.Tt.start = 1
k = 100
.sei = runif( n = k, 0.1, 2 )
.yi = rnorm( n = k,
             mean = .Mu.start,
             sd = sqrt( .Tt.start^2 + .sei^2 ) )

.tcrit = rep(1.96, k)

# if you want to truncate
keep = .yi < .tcrit * .sei
.yi = .yi[ keep == TRUE ]
.tcrit = .tcrit[ keep == TRUE ]
.sei = .sei[ keep == TRUE ]

# set start values for sampler
init.fcn <- function(o){ list(mu = .Mu.start,
                              tau = .Tt.start ) }

post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(.yi),
                             sei = .sei,
                             tcrit = .tcrit,
                             y = .yi ),
                iter=2000,  
                
                init = init.fcn)


# ~~ Example 2: Data from sim_meta -------------------------


# all methods agree closely for this example! :)

Mu = .1
T2 = .25  # ACROSS-study 
t2w = .25  # WITHIN-STUDY
k = 1000
true.sei.expr = "runif(n = 1, min = 0.1, max = 1)"
m = 500


# SAVE: this is how data were generated
d = sim_meta_2(Nmax = 1,
               Mu = Mu,
               T2 = T2,
               m = m,
               t2w = t2w,
               true.sei.expr = true.sei.expr,
               hack = "affirm",
               
               rho = 0,
               
               k.pub.nonaffirm = 500,
               prob.hacked = 0,
               
               return.only.published = FALSE )


dn = d %>% filter(Di == 1 & affirm == FALSE)
dn$sei = sqrt(dn$vi)
# sanity check
all( dn$yi <= dn$tcrit * dn$sei )

yi = dn$yi
sei = dn$sei
tcrit = dn$tcrit

#@ for now, start at truth (later should start at best MCMC iterate as in TNE)
Mu.start = Mu
Tt.start = sqrt(T2)



# set start values for sampler
init.fcn <- function(o){ list(mu = Mu.start,
                              tau = Tt.start ) }

### Version 0: MCMC ###
# post means: (0.13, 0.75)
# post medians: (0.13, 0.74)
post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(yi),
                             sei = sei,
                             tcrit = tcrit,
                             y = yi ),
                iter=2000,  
                
                init = init.fcn)


### Version 1: MLE ###
# (0.13, 0.74)
res.MLE.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    par2.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = FALSE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MLE.1$MuHat
res.MLE.1$TtHat

### Version 2: MAP (SD parameterization) ###
# (0.14, 0.75)
res.MAP.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    par2.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MAP.1$MuHat
res.MAP.1$TtHat


# ~~ Example 3: Hagger -------------------------


setwd(path)
setwd("misc")
dm = fread("prepped_hagger_meta_data.csv")


# SUMMARY OF THIS SECTION:
# MLE = 3.31 (seems very big, but agrees with weightr)
# MAP = 305!!! (horrible)

# scenario: different SEs as in Hagger meta, and we have to use the estimated SEs
#  rather than the true ones

tcrit = qnorm(.975)

dn = dm %>% filter( d/sqrt(var) < zcrit )

# set vars for all methods
yi = dn$d
sei = sqrt(dn$var)
kn = nrow(dn)
Mu.start = 0
Tt.start = 1

# sanity check 
all( dn$yi <= tcrit * sei )

# for sampling call happiness
tcrit = rep(tcrit, length(yi))


# set start values for sampler
init.fcn <- function(o){ list(mu = Mu.start,
                              tau = Tt.start ) }

### Version 0: MCMC ###

# With default stan control parameters:
# post means: (4.08, 0.73)
# post medians: (2.49, 0.62)
# Rhats: 1.06, 1.09

# With max_treedepth = 20, adapt_delta = 0.98:
# post means: (4.20, 0.74); 95% CI: [0.97, 19.37]
# post medians: (2.51, 0.62); 95% CI: [0.28, 1.97]
# Rhats: 1.01, 1.01
post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(yi),
                             sei = sei,
                             tcrit = tcrit,
                             y = yi ),
                iter=2000,  
                control = list(max_treedepth = 20,
                               adapt_delta = 0.98),
                
                init = init.fcn)


### Version 1: MLE ###
# (3.31, 0.73)
res.MLE.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    par2.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = FALSE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MLE.1$MuHat
res.MLE.1$TtHat

### Version 2: MAP (SD parameterization) ###
# (305, 8)!!
# but that's not starting at the best MCMC value
res.MAP.1 = estimate_jeffreys_RTMA( yi = yi,
                                    sei = sei,
                                    par2is = "Tt",
                                    Mu.start = Mu.start,
                                    par2.start = Tt.start,
                                    tcrit = tcrit,
                                    
                                    usePrior = TRUE,
                                    get.CIs = TRUE,
                                    CI.method = "wald" )

res.MAP.1$MuHat
res.MAP.1$TtHat



