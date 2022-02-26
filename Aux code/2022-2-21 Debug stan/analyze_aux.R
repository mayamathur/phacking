
# Goal: Why is the model not running at all?
# code is from helper_SAPH::estimate_jeffreys_RTMA on 2022-2-21


# for cluster:
# srun --mem=32G --time=4:00:00 --pty bash


library(rstan)

# first 30 observations of a dataframe
.yi = c(0.426251447318531, -0.274990385250978, -1.34416733163289, -0.132628196524813, 
        -1.07698887654361, 0.288016076714556, -0.280095486816409, 0.479115006596442, 
        0.178435851264946, -1.62618107209329, -0.514258399124606, 0.418936675881624, 
        0.341423767215658, -0.151980156337454, -0.700114988499077, -0.409716483485722, 
        -0.254066897447185, 1.2236451198118, 0.355352916095302, 0.470713258158004, 
        -1.29771677134009, 1.18371578091549, -0.59548594078466, 0.0301412541156913, 
        0.504749345722535, 0.0581004814175708, -0.0719312212494192, -0.450800849474051, 
        -1.72030207382822, -0.566558437014423)

.sei = c(0.653941663515478, 0.935502004400324, 0.657946935979067, 0.343265360845246, 
         0.698270521109975, 0.186932530494982, 0.306383087120182, 0.34112245234349, 
         0.293304997517225, 0.764446813023732, 0.136181335477904, 0.591601485861967, 
         0.353403324040768, 0.59888237542934, 0.497496505528266, 0.578369401560294, 
         0.274902197809735, 0.894906700567698, 0.505578796847308, 0.466748233209642, 
         0.619338959219028, 0.922918175225004, 0.641840610807037, 0.203034613543438, 
         0.297932447600233, 0.227229934342497, 0.626820173595342, 0.926976252855868, 
         0.742427708721177, 0.610892579685947)

.tcrit = c(1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769, 1.96472939098769, 1.96472939098769, 
           1.96472939098769, 1.96472939098769)

.Mu.start = 0.1
.Tt.start = 0.5



# FULL VERSION (BREAKS) --------------------------------------
model.text <- "
functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
    real sigma;
		real LL; 
		real UU;
		// will just be set to 1
		int n; 
    
    // this will now be the fisher info for EACH obs
		matrix[2,2] fishinfo;
		
		// this will be the TOTAL fisher info
		matrix[2,2] fishinfototal;
		

		// MM: build a Fisher info matrix for EACH observation
		for (i in 1:k) {
		
		  sigma = sqrt(tau^2 + sei[i]^2);
		  LL = -999;
		  UU = tcrit[i] * sei[i];
		
  		mustarL = (LL - mu) / sigma;
  		mustarU = (UU - mu) / sigma;
  		
  		// because EACH fisher info below has n=1 only
  		n = 1; 
  		
  		// MM: BELOW NEEDS NO EDITS FROM TNE
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
  		// MM: END STUFF THAT NEEDS NO EDITS FROM TNE
  		
  		// MM: add the new fisher info to the total one
  		fishinfototal = fishinfototal +  fishinfo;
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
        y[i] ~ normal(mu, sqrt(tau^2 + sei[i]^2))T[-999, tcrit[i] * sei[i]];
}

generated quantities{
  real log_lik = 0;
  // note that the log_prior is over all k observations already
  real log_prior = log(jeffreys_prior(mu, tau, k, sei, tcrit));
  real log_post;
  
  for ( i in 1:k ){
  
  	real LL = -999;
		real UU = tcrit[i] * sei[i];
  
    log_lik += normal_lpdf(y | mu, sqrt(tau^2 + sei[i]^2));
    log_lik += -1 * log_diff_exp( normal_lcdf(UU | mu, sqrt(tau^2 + sei[i]^2) ), normal_lcdf(LL | mu, sqrt(tau^2 + sei[i]^2) ) );  	
  }
  
  log_post = log_lik + log_prior;
}
"



# #3: SIMPLIFY PRIOR AND LPDF  (same error) --------------------------------------
model.text <- "
functions{

// now just a uniform prior
	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
		return 0;
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
	  y[i] ~ normal(mu, tau)T[-999, tcrit[i] * sei[i]];
        // y[i] ~ normal(mu, tau)T[-999, tcrit[i] * sei[i]];
}

generated quantities{
  real log_lik = 0;
  real log_prior = log(jeffreys_prior(mu, tau, k, sei, tcrit));
  real log_post;
  
  for ( i in 1:k ){
  
  	// real LL = -999;
		// real UU = tcrit[i] * sei[i];
  
    // log_lik += normal_lpdf(y | mu, sqrt(tau^2 + sei[i]^2));
    // log_lik += -1 * log_diff_exp( normal_lcdf(UU | mu, sqrt(tau^2 + sei[i]^2) ), normal_lcdf(LL | mu, sqrt(tau^2 + sei[i]^2) ) );
    
      log_lik += normal_lpdf(y | mu, sqrt(tau^2 + sei[i]^2));

  }
  
  log_post = log_lik + log_prior;
}
"


# #3: SIMPLIFY PRIOR AND LPDF; remove truncation  (same error)
model.text <- "
functions{

// now just a uniform prior
	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
		return 0;
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
	//CHANGED
        y[i] ~ normal(mu, tau);
}

generated quantities{
  real log_lik = 0;
  real log_prior = log(jeffreys_prior(mu, tau, k, sei, tcrit));
  real log_post;
  
  for ( i in 1:k ){
  
  	// real LL = -999;
		// real UU = tcrit[i] * sei[i];
  
    // log_lik += normal_lpdf(y | mu, sqrt(tau^2 + sei[i]^2));
    // log_lik += -1 * log_diff_exp( normal_lcdf(UU | mu, sqrt(tau^2 + sei[i]^2) ), normal_lcdf(LL | mu, sqrt(tau^2 + sei[i]^2) ) );
    
    //CHANGED
      log_lik += normal_lpdf(y | mu, tau);

  }
  
  log_post = log_lik + log_prior;
}
"










# #4: COMPLETELY REMOVE PRIOR - WORKS ON SHERLOCK --------------------------------------
model.text <- "

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
	target += 0;
	for(i in 1:k)
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
}

generated quantities{
  real log_lik = 0;
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y | mu, sqrt(tau^2 + sei[i]^2) );
  }
  
  log_post = log_lik;
}
"


stan.model <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")

.Mu.start = 0
.Tt.start = 1
k = 100
.sei = runif( n = k, 0.1, 2 )
.yi = rnorm( n = k,
             mean = .Mu.start,
             sd = sqrt( .Tt.start^2 + .sei^2 ) )

.tcrit = rep(0, k)


# set start values for sampler
init.fcn <- function(o){ list(mu = .Mu.start,
                              tau = .Tt.start ) }



options(mc.cores = parallel::detectCores())

# package example
# m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
# f <- sampling(m, iter = 1)  # I reduced iter from 100 to 1
# that DOES prevent crashing
# maybe try running it on cluster to see if that prevents crashing?
# even using iter=1 below doesn't work, so I do think it's a memory issue
# try on cluster
# given that #4 above seems ok for syntax but #3 doesn't, seems like the way the prior is written in #3 and earlier 
#  is wrong

# on Sherlock, this WORKS and gets reasonable estimates of (mu, tau)
post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(.yi),
                             sei = .sei,
                             tcrit = .tcrit,
                             y = .yi ),
                iter=2000,  
                
                init = init.fcn)






# #5 ADD PRIORS ONE AT A TIME --------------------------------------


# add uniform prior - WORKS on Sherlock
model.text <- "

// STEP 1 (WORKS): Add just the fn itself but don't call it
functions{
	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	// need to return 1 here since we'll be taking log
		return 1;
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
	//target += 0;
	// STEP 2: Replace target so that now we're adding in the fake prior that's always 0
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit) );
	for(i in 1:k)
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
}

generated quantities{
  real log_lik = 0;
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y | mu, sqrt(tau^2 + sei[i]^2) );
  }
  
  log_post = log_lik;
}
"



# start adding Jeffreys prior but don't actually calculate Fisher info - WORKS
model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
    real sigma;
		real LL; 
		real UU;
		// will just be set to 1
		int n; 
    
    // this will now be the fisher info for EACH obs
		matrix[2,2] fishinfo;
		
		// this will be the TOTAL fisher info
		matrix[2,2] fishinfototal;
		
		
		// WORKS: make fake fisherinfo whose determinant is 1
		fishinfo[1,1] = 1;
  	fishinfo[1,2] = 1;
  	fishinfo[2,1] = 1;
  	fishinfo[2,2] = 2;
  		
  	// doesn't print anything... :(
  	print(fishinfo[1,1]);
  		
  	return sqrt( determinant(fishinfo) );

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
	//target += 0;
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit) );
	for(i in 1:k)
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
}

generated quantities{
  real log_lik = 0;
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y | mu, sqrt(tau^2 + sei[i]^2) );
  }
  
  log_post = log_lik;
}
"




# this WORKS
model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
    real sigma;
		real LL; 
		real UU;
		// will just be set to 1
		int n; 
    
    // this will now be the fisher info for EACH obs
		matrix[2,2] fishinfo;
		
		// this will be the TOTAL fisher info
		matrix[2,2] fishinfototal;
		
		
		// MM: build a Fisher info matrix for EACH observation
		// this WORKS
		for (i in 1:k) {
		
		  sigma = sqrt(tau^2 + sei[i]^2);
		  LL = -999;
		  UU = tcrit[i] * sei[i];
		
  		mustarL = (LL - mu) / sigma;
  		mustarU = (UU - mu) / sigma;
  		
  		// because EACH fisher info below has n=1 only
  		n = 1; 
  		
  		// MM: BELOW NEEDS NO EDITS FROM TNE
  		// note that normal_lpdf, etc., parameterize in terms of SD, not var
  		//  the (0,1) below are *not* start values for MCMC
  		// THIS IS THE LINE THAT BREAKS! Error in sampler$call_sampler(args_list[[i]]) : Initialization failed
  		// I know that mustarU is FINE and mustarL is BAD
  		mustarL = -999; // fine if you use these!
  		alphaL = exp( normal_lpdf(mustarL | 0, 1) - 
  	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
  	                normal_lcdf(mustarL | 0, 1) ) ); 
  	                
  	                
  	 // fake version: WORKS
  	 // so something is wrong with mustarL, mustarU above
  	  // alphaL = exp( normal_lpdf(-9 | 0, 1) - 
  	  //              log_diff_exp( normal_lcdf(0 | 0, 1),
  	   //             normal_lcdf(0 | 0, 1) ) ); 
  	                
  		//  alphaU = exp( normal_lpdf(mustarU | 0, 1) - 
   	  //              log_diff_exp( normal_lcdf(mustarU | 0, 1),
   	  //              normal_lcdf(mustarL | 0, 1) ) );
  		
  		// second derivatives for Fisher info			
  		// kmm = -n/sigma^2 + n/sigma^2 * ((alphaU-alphaL)^2 + alphaU*mustarU- alphaL*mustarL);
  		// kms = -2*n/sigma^2 * (alphaL - alphaU) + 
  	  // 		  n/sigma^2 * (alphaL - alphaU + (alphaU*mustarU^2 - alphaL*mustarL^2) +
  		//	  				(alphaL-alphaU) * (alphaL*mustarL - alphaU*mustarU));
  		// kss = n/sigma^2 - 3*n/sigma^2 * (1 + mustarL*alphaL - mustarU*alphaU) +
  	 //  			n/sigma^2 * (mustarU*alphaU*(mustarU^2 - 2) - mustarL*alphaL*(mustarL^2 - 2) +
  		//						(alphaU*mustarU - alphaL*mustarL)^2);
  		
  		// fishinfo[1,1] = -kmm;
  		// fishinfo[1,2] = -kms;
  		// fishinfo[2,1] = -kms;
  		// fishinfo[2,2] = -kss;
  		// MM: END STUFF THAT NEEDS NO EDITS FROM TNE
  		
  		// MM: add the new fisher info to the total one
  		// fishinfototal = fishinfototal +  fishinfo;
		}
		
		
		// WORKS: make fake fishinfototal whose determinant is 1
		fishinfototal[1,1] = 1;
  	fishinfototal[1,2] = 1;
  	fishinfototal[2,1] = 1;
  	fishinfototal[2,2] = 2;
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
	//target += 0;
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit) );
	for(i in 1:k)
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
}

generated quantities{
  real log_lik = 0;
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y | mu, sqrt(tau^2 + sei[i]^2) );
  }
  
  log_post = log_lik;
}
"




# WORKS
model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
    real sigma;
		real LL; 
		real UU;
		// will just be set to 1
		int n; 
    
    // this will now be the fisher info for EACH obs
		matrix[2,2] fishinfo;
		
		// this will be the TOTAL fisher info
		matrix[2,2] fishinfototal;
		
		
		// MM: build a Fisher info matrix for EACH observation
		// this WORKS
		for (i in 1:k) {
		
		  sigma = sqrt(tau^2 + sei[i]^2);
		  UU = tcrit[i] * sei[i];
		
  		mustarL = -999;
  		mustarU = (UU - mu) / sigma;
  		
  		// because EACH fisher info below has n=1 only
  		n = 1; 
  		
  		// note that normal_lpdf, etc., parameterize in terms of SD, not var
  		//  the (0,1) below are *not* start values for MCMC
  		alphaL = exp( normal_lpdf(mustarL | 0, 1) - 
  	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
  	                normal_lcdf(mustarL | 0, 1) ) ); 
  	                
  		alphaU = exp( normal_lpdf(mustarU | 0, 1) - 
   	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
   	                normal_lcdf(mustarL | 0, 1) ) );
   	                
   	  // runs fine if everything below in for-loop is commented out
  		
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
  		
  		// MM: add the new fisher info to the total one
  		// 2022-2-22: BREAKS IF YOU UNCOMMENT THIS LINE (even just this line)
  		 // fishinfototal = fishinfototal +  fishinfo;
		}
		
		// WORKS: make fake fishinfototal whose determinant is 1
		fishinfototal[1,1] = 1;
  	fishinfototal[1,2] = 1;
  	fishinfototal[2,1] = 1;
  	fishinfototal[2,2] = 2;
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
	//target += 0;
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit) );
	for(i in 1:k)
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
}

generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit) );
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y | mu, sqrt(tau^2 + sei[i]^2) );
  }
  
  // note that prior isn't used yet
  log_post = log_lik;
}
"




# WORKS :)
# not yet implemented in this:
#  - truncation in lpdfs
#  - fishintototal is still fake
model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit){
	
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
		// this WORKS
		for (i in 1:k) {
		
		  sigma = sqrt(tau^2 + sei[i]^2);
		  UU = tcrit[i] * sei[i];
		
  		mustarL = -999;
  		mustarU = (UU - mu) / sigma;
  		
  		// because EACH fisher info below has n=1 only
  		n = 1; 
  		
  		// note that normal_lpdf, etc., parameterize in terms of SD, not var
  		//  the (0,1) below are *not* start values for MCMC
  		alphaL = exp( normal_lpdf(mustarL | 0, 1) - 
  	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
  	                normal_lcdf(mustarL | 0, 1) ) ); 
  	                
  		alphaU = exp( normal_lpdf(mustarU | 0, 1) - 
   	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
   	                normal_lcdf(mustarL | 0, 1) ) );
   	                
   	  // runs fine if everything below in for-loop is commented out
  		
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
  		
  		// MM: add the new fisher info to the total one
  		// 2022-2-22: BREAKS IF YOU UNCOMMENT THIS LINE (even just this line)
  		 fishinfototal = fishinfototal + fishinfo;
		}
		
		// WORKS: make fake fishinfototal whose determinant is 1
		//fishinfototal[1,1] = 1;
  	//fishinfototal[1,2] = 1;
  	//fishinfototal[2,1] = 1;
  	//fishinfototal[2,2] = 2;
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
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
}

generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit) );
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf( y | mu, sqrt(tau^2 + sei[i]^2) );
  }
  
  // note that prior isn't used yet
  log_post = log_prior + log_lik;
}
"









# WORKS :)
# first attempt at truncated
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

# .Mu.start = 0
# .Tt.start = 1
# k = 100
# .sei = runif( n = k, 0.1, 2 )
# .yi = rnorm( n = k,
#              mean = .Mu.start,
#              sd = sqrt( .Tt.start^2 + .sei^2 ) )
# 
# .tcrit = rep(1.96, k)

# # if you want to truncate
# keep = .yi < .tcrit * .sei
# .yi = .yi[ keep == TRUE ]
# .tcrit = .tcrit[ keep == TRUE ]
# .sei = .sei[ keep == TRUE ]
# 
# # set start values for sampler
# init.fcn <- function(o){ list(mu = .Mu.start,
#                               tau = .Tt.start ) }




options(mc.cores = parallel::detectCores())

post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(.yi),
                             sei = .sei,
                             tcrit = .tcrit,
                             y = .yi ),
                iter=2000,  
                
                init = init.fcn)








