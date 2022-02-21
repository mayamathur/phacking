
# Goal: Why is the model not running at all?
# code is from helper_SAPH::estimate_jeffreys_RTMA on 2022-2-21


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



# FULL VERSION (BREAKS)
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



# #3: SIMPLIFY PRIOR AND LPDF  (same error)
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


# #4: COMPLETELY REMOVE PRIOR  (same error)
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
        y[i] ~ normal(mu, tau);
}

generated quantities{
  real log_lik = 0;
  real log_post;
  
  for ( i in 1:k ){
      log_lik += normal_lpdf(y | mu, tau);
  }
  
  log_post = log_lik;
}
"


stan.model <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")


# set start values for sampler
init.fcn <- function(o){ list(mu = .Mu.start,
                              tau = .Tt.start ) }

# bm: this model seems to at least try to run, but crashed R



# bm: really vague error:
# "[1] "Error in sampler$call_sampler(args_list[[i]]) : Initialization failed."
#  error occurred during calling the sampler; sampling not done"
post = sampling(stan.model,
                cores = 1,
                refresh = 0,
                data = list( k = length(.yi),
                             sei = .sei,
                             tcrit = .tcrit,
                             y = .yi ),
                
                
                init = init.fcn)
