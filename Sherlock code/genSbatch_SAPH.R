
# SET SIMULATION PARAMETERS MATRIX -----------------------------------------

# FOR CLUSTER USE
path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

allPackages = c("here",
                "magrittr",
                "dplyr",
                "data.table",
                "tidyverse",
                "tidyr",
                "metafor",
                "robumeta",
                "testthat",
                "truncdist",
                "gmm",
                "tmvtnorm",
                "doParallel",
                "foreach")


( packagesNeeded = allPackages[ !( allPackages %in% installed.packages()[,"Package"] ) ] )
if( length(packagesNeeded) > 0 ) install.packages(packagesNeeded)

# load all packages
lapply( allPackages,
        require,
        character.only = TRUE)

#**you need to see all "TRUE" printed by this in order for the package to actually be loaded

# set up sim params for cluster

# main scenarios of interest:
# 1. Nmax = 1, k.hacked = 0, rho = 0 (basically a sanity check)
# 2. Nmax > 1, k.hacked = 0, rho = 0.9 
# 3. Nmax > 1, k.hacked = 50, rho = 0 or 0.9 (conservative?)

# "full version"
scen.params = tidyr::expand_grid(
  rep.methods = "naive ; gold-std ; maon ; 2psm ; jeffreys-mcmc ; jeffreys-sd ; jeffreys-var ; mle-sd ; mle-var",
  #rep.methods = "jeffreys-mcmc ; jeffreys-sd ; jeffreys-var",
                                  # args from sim_meta_2
                                  Nmax = 10,
                                  Mu = 0.1,
                                  t2a = 0.25,
                                  t2w = 0.25,
                                  m = 50,
                                  true.sei.expr = "runif(n = 1, min = 0.1, max = 1)",
                                  hack = "affirm",
                                  rho = 0,
                                  k.pub.nonaffirm = 50,
                                  prob.hacked = 0.2,

                                  # Stan control args
                                  stan.maxtreedepth = 20,
                                  stan.adapt_delta = 0.98,

                                  get.CIs = TRUE,
                                  run.optimx = TRUE )

# # just for testing MCMC issues
# scen.params = tidyr::expand_grid(
#   rep.methods = "jeffreys-mcmc",
#   
#   # args from sim_meta_2
#   Nmax = 10,
#   Mu = 0.1,
#   t2a = 0.25,
#   t2w = 0.25,
#   m = 50,
#   true.sei.expr = "runif(n = 1, min = 0.1, max = 1)",
#   hack = "affirm",
#   rho = 0,
#   k.pub.nonaffirm = 50,
#   prob.hacked = 0,
#   
#   # Stan control args
#   stan.maxtreedepth = 10,
#   stan.adapt_delta = 0.8,
#   
#   get.CIs = TRUE,
#   run.optimx = FALSE )


# # OLD - Do I still need these?
# # hold constant the number of UNHACKED studies to 800
# scen.params$k[ scen.params$k.hacked == 800 ] = 800*2
# table(scen.params$k, scen.params$k.hacked)
# 
# # remove nonsense combinations
# # rho > 0 is pointless if there's only 1 draw
# scen.params = scen.params %>% dplyr::filter( !(rho > 0 & Nmax == 1) )

# add scen numbers
scen.params = scen.params %>% add_column( scen = 1:nrow(scen.params),
                                          .before = 1 )


( n.scen = nrow(scen.params) )
# look at it
head( as.data.frame(scen.params) )

# write the csv file of params (to Sherlock)
setwd(path)
write.csv( scen.params, "scen_params.csv", row.names = FALSE )


########################### GENERATE SBATCHES ###########################

# load functions for generating sbatch files
source("helper_SAPH.R")

# number of sbatches to generate (i.e., iterations within each scenario)
n.reps.per.scen = 1000  
n.reps.in.doParallel = 10  #@update these
( n.files = ( n.reps.per.scen / n.reps.in.doParallel ) * n.scen )


path = "/home/groups/manishad/SAPH"

scen.name = rep( scen.params$scen, each = ( n.files / n.scen ) )
jobname = paste("job", 1:n.files, sep="_")
outfile = paste("rm_", 1:n.files, ".out", sep="")
errorfile = paste("rm_", 1:n.files, ".err", sep="")
write_path = paste(path, "/sbatch_files/", 1:n.files, ".sbatch", sep="")
runfile_path = paste(path, "/testRunFile.R", sep="")


# 2022-2-27: timing benchmark: with all methods and sim.reps = 1, took 2.5 min
sbatch_params <- data.frame(jobname,
                            outfile,
                            errorfile,
                            jobtime = "01:00:00",  #@update this; was 1:00:00 with all methods and sim.reps = 10
                            quality = "normal",
                            node_number = 1,
                            mem_per_node = 64000,
                            mailtype =  "NONE",
                            user_email = "mmathur@stanford.edu",
                            tasks_per_node = 16,
                            cpus_per_task = 1,
                            path_to_r_script = paste(path, "/doParallel_SAPH.R", sep=""),
                            args_to_r_script = paste("--args", jobname, scen.name, sep=" "),
                            write_path,
                            stringsAsFactors = F,
                            server_sbatch_path = NA)

generateSbatch(sbatch_params, runfile_path)

n.files

# run just the first one
# sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/1.sbatch

path = "/home/groups/manishad/SAPH"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 1:100) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
}




######## If Running Only Some Jobs To Fill Gaps ########

# run in Sherlock ml load R
path = "/home/groups/manishad/SAPH"
setwd(path)
source("functions_SAPH.R")

missed.nums = sbatch_not_run( "/home/groups/manishad/SAPH/sim_results/long",
                              "/home/groups/manishad/SAPH/sim_results",
                              .name.prefix = "long_results",
                              .max.sbatch.num = 1440 )



setwd( paste(path, "/sbatch_files", sep="") )
for (i in missed.nums) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
}