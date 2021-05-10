
# SET SIMULATION PARAMETERS MATRIX -----------------------------------------

# FOR CLUSTER USE
path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

library(here, lib.loc = "/home/groups/manishad/Rpackages/")
# data-wrangling packages
library(magrittr, lib.loc = "/home/groups/manishad/Rpackages/")
library(dplyr, lib.loc = "/home/groups/manishad/Rpackages/")
library(data.table, lib.loc = "/home/groups/manishad/Rpackages/")
library(tidyverse)
library(tidyr, lib.loc = "/home/groups/manishad/Rpackages/")
# meta-analysis packages
library(metafor, lib.loc = "/home/groups/manishad/Rpackages/")
library(robumeta, lib.loc = "/home/groups/manishad/Rpackages/")
# other
library(testthat)
# for this project
library(truncdist, lib.loc = "/home/groups/manishad/Rpackages/")
#library(ExtDist)
library(gmm, lib.loc = "/home/groups/manishad/Rpackages/")  # https://stackoverflow.com/questions/63511986/error-package-or-namespace-load-failed-for-gmm-in-dyn-loadfile-dllpath-dl
library(tmvtnorm, lib.loc = "/home/groups/manishad/Rpackages/")

# in case packages need to be installed
# install.packages("tidyr", lib = "/home/groups/manishad/Rpackages/")


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
                           hack = "affirm",
                           rho = c(0, 0.9),
                           
                           k = 100,
                           k.hacked = c(0, 50) )

# remove nonsense combinations
# rho > 0 is pointless if there's only 1 draw
scen.params = scen.params %>% dplyr::filter( !(rho > 0 & Nmax == 1) )



( n.scen = nrow(scen.params) )
# look at it
head( as.data.frame(scen.params) )

# write the csv file of params (to Sherlock)
write.csv( scen.params, "scen_params.csv", row.names = FALSE )


########################### GENERATE SBATCHES ###########################

# load functions for generating sbatch files
source("helper_SAPH.R")

# number of sbatches to generate (i.e., iterations within each scenario)
n.reps.per.scen = 5  # if you want to generate only 1 file, set this to 10
n.reps.in.doParallel = 5  #@update these
( n.files = ( n.reps.per.scen / n.reps.in.doParallel ) * n.scen )


path = "/home/groups/manishad/SAPH"

scen.name = rep( scen.params$scen.name, each = ( n.files / n.scen ) )
jobname = paste("job", 1:n.files, sep="_")
outfile = paste("rm_", 1:n.files, ".out", sep="")
errorfile = paste("rm_", 1:n.files, ".err", sep="")
write_path = paste(path, "/sbatch_files/", 1:n.files, ".sbatch", sep="")
runfile_path = paste(path, "/testRunFile.R", sep="")

sbatch_params <- data.frame(jobname,
                            outfile,
                            errorfile,
                            jobtime = "0:30:00",  #@update this
                            quality = "normal",
                            node_number = 1,
                            mem_per_node = 64000,
                            mailtype =  "NONE",
                            user_email = "mmathur@stanford.edu",
                            tasks_per_node = 16,
                            cpus_per_task = 1,
                            path_to_r_script = paste(path, "/doParallel.R", sep=""),
                            args_to_r_script = paste("--args", jobname, scen.name, sep=" "),
                            write_path,
                            stringsAsFactors = F,
                            server_sbatch_path = NA)

generateSbatch(sbatch_params, runfile_path)

n.files



# max hourly submissions seems to be 300, which is 12 seconds/job
path = "/home/groups/manishad/SAPH"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 1:1) {
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