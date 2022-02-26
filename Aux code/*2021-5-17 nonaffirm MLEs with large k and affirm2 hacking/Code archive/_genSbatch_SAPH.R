
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

scen.params = expand_grid( Mu = 0.1,
                           T2 = c(0, 0.25),
                           m = 500,
                           t2w = c(0, 0.25),
                           se = 0.5,
                           
                           Nmax = c(1, 10),
                           hack = "affirm2",
                           rho = c(0, 0.9),
                           
                           k = 800,
                           k.hacked = c(0, 800) )

# hold constant the number of UNHACKED studies to 800
scen.params$k[ scen.params$k.hacked == 800 ] = 800*2
table(scen.params$k, scen.params$k.hacked)

# remove nonsense combinations
# rho > 0 is pointless if there's only 1 draw
scen.params = scen.params %>% dplyr::filter( !(rho > 0 & Nmax == 1) )

scen.params = scen.params %>% add_column( scen = 1:nrow(scen.params),
                                          .before = 1 )


( n.scen = nrow(scen.params) )
# look at it
head( as.data.frame(scen.params) )

# write the csv file of params (to Sherlock)
write.csv( scen.params, "scen_params.csv", row.names = FALSE )


########################### GENERATE SBATCHES ###########################

# load functions for generating sbatch files
source("helper_SAPH.R")

# number of sbatches to generate (i.e., iterations within each scenario)
n.reps.per.scen = 500  
n.reps.in.doParallel = 100  #@update these
( n.files = ( n.reps.per.scen / n.reps.in.doParallel ) * n.scen )


path = "/home/groups/manishad/SAPH"

scen.name = rep( scen.params$scen, each = ( n.files / n.scen ) )
jobname = paste("job", 1:n.files, sep="_")
outfile = paste("rm_", 1:n.files, ".out", sep="")
errorfile = paste("rm_", 1:n.files, ".err", sep="")
write_path = paste(path, "/sbatch_files/", 1:n.files, ".sbatch", sep="")
runfile_path = paste(path, "/testRunFile.R", sep="")

# time with k = 100, reps.in.doParallel = 100:
# > summary(t$doParallelMin)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.03563  0.06377  0.68358  1.46078  0.99570 78.15405 

sbatch_params <- data.frame(jobname,
                            outfile,
                            errorfile,
                            jobtime = "2:00:00",  #@update this
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



# max hourly submissions seems to be 300, which is 12 seconds/job
path = "/home/groups/manishad/SAPH"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 1:n.files) {
  #system( paste("sbatch -p owners /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
  #Sys.sleep(2)  # delay in seconds
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