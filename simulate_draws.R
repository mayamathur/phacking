
# Simulate distribution of p-values if there are various amounts of p-hacking (or limiting case with N->infty)

################################ N DRAWS VS. PAFFIRM ################################ 

N.max = Inf

# meta-analytic distribution
Mu = 0
T2 = 0

# study parameters
t2 = 0  # within-study heterogeneity due to model choice
n = 50
se = 1
# k = 1  # heterogeneity ratio of STUDY, conditional on its mean
# prop = 0.95





sim.reps = 500


library(foreach)
library(doParallel)
registerDoParallel(cores=16)

rs = foreach( i = 1:sim.reps,
              .combine=rbind,
              .errorhandling = "stop"  # this shouldn't happen
) %dopar% {
  
  phack_study( N.max = N.max,
               Mu = Mu,
               T2 = T2,
               n = n,
               t2 = t2,
               se = se,
               hack = "affirm")
  
}  ### end parallelized loop


# distribution of only the significant p-values given hacking
# Unif(0,0.05) even with t2 > 0!!! (study heterogeneity)
summary(rs$N)
summary(rs$pval)

mean(rs$pval<0.05) / 0.05  # if no heterogeneity, this is similar to eta
# related to prob of getting a significant result (i.e., power) in N.max trials, 
#  and I don't think this depends on correlation of tests?
#  bm :)


hist(rs$pval, breaks = 30)
mean(rs$ybar)  # with no heterogeneity and mu=0, won't actually be biased unless we hack to get affirmative status rather than significance


# increasing t2 or T2 results in right skew, which makes sense

# for un-hacked distribution (i.e., N.max = 1)
# this is the distribution of p-values in a meta-analysis with no heterogeneity 
#  and eta = infinity (we only see the significant p-values)
mean(rs$pval<0.05)  # should be 5% if no heterogeneity and Mu=0
hist(rs$pval[rs$pval<0.05], breaks = 30)

# ** when do p-values bunch up under 0.05??


# need to look at Simonsohn 2014's "modeling p-hacking" section to understand
#  when p-hacking will create left skew

# I think that even when draws are correlated, this wouldn't affect E[N], only the
#  distribution of N? So could still do sensitivity analyses?

# look at bias in POINT ESTIMATES when hacking with no correlation



