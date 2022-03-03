

# to be run by stitch.sbatch or manually

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
start.num = as.numeric( args[1] )  # starting results number to stitch
stop.num = as.numeric( args[2] )  # stopping results number to stitch

path = "/home/groups/manishad/SAPH"
setwd(path)
source("helper_SAPH.R")

######## STITCH LONG FILES ########

library(data.table)
library(dplyr)
# s = stitch_files(.results.singles.path = "/home/groups/manishad/SAPH/sim_results/long",
#                  .results.stitched.write.path = "/home/groups/manishad/SAPH/sim_results/overall_stitched",
#                  .name.prefix = "long_results",
#                  .stitch.file.name="stitched.csv")

.results.singles.path = "/home/groups/manishad/SAPH/long_results"
.results.stitched.write.path = "/home/groups/manishad/SAPH/stitched_results"
.name.prefix = "long_results"
.stitch.file.name="stitched.csv"

# get list of all files in folder
all.files = list.files(.results.singles.path, full.names=TRUE)

# we only want the ones whose name includes .name.prefix
keepers = all.files[ grep( .name.prefix, all.files ) ]
length(keepers)

# grab variable names from first file
names = names( read.csv(keepers[1] ) )

# read in and rbind the keepers
tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )

# sanity check: do all files have the same names?
# if not, could be because some jobs were killed early so didn't get doParallelTime
#  variable added at the end
#  can be verified by looking at out-file for a job without name "doParallelTime"
allNames = lapply( tables, names )
# # find out which jobs had wrong number of names
# lapply( allNames, function(x) all.equal(x, names ) )
# allNames[[1]][ !allNames[[1]] %in% allNames[[111]] ]

# bind_rows works even if datasets have different names
#  will fill in NAs
s <- do.call(bind_rows, tables)



names(s) = names( read.csv(keepers[1], header= TRUE) )

if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )

# # are we there yet?
# nrow(s) / (1600*500)  # main sims: 1600*500, bias correction sims: 32*3*500
# length(unique(s$scen.name))  # main sims: 1600; bias correction sims: 26

##### Look for Missed Jobs #####
# look for missed jobs
missed.nums = sbatch_not_run( "/home/groups/manishad/SAPH/long_results",
                              "/home/groups/manishad/SAPH/long_results",
                              .name.prefix = "long",
                              .max.sbatch.num = 100)

path = "/home/groups/manishad/SAPH"

setwd( paste(path, "/sbatch_files", sep="") )
for (i in missed.nums) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/", i, ".sbatch", sep="") )
}

##### Move to Desktop #####
# Sherlock -> Desktop
scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/sim_results/overall_stitched/stitched.csv ~/Desktop




##### Quick Look at Results #####

library(cli, lib.loc = "/home/groups/manishad/Rpackages/")
library(fansi, lib.loc = "/home/groups/manishad/Rpackages/")
library(utf8, lib.loc = "/home/groups/manishad/Rpackages/")
library(rlang, lib.loc = "/home/groups/manishad/Rpackages/")
library(crayon, lib.loc = "/home/groups/manishad/Rpackages/")
library(dplyr, lib.loc = "/home/groups/manishad/Rpackages/")
library(foreach, lib.loc = "/home/groups/manishad/Rpackages/")
library(doParallel, lib.loc = "/home/groups/manishad/Rpackages/")
library(boot, lib.loc = "/home/groups/manishad/Rpackages/")
library(metafor, lib.loc = "/home/groups/manishad/Rpackages/")
library(robumeta, lib.loc = "/home/groups/manishad/Rpackages/")
library(data.table, lib.loc = "/home/groups/manishad/Rpackages/")
library(purrr, lib.loc = "/home/groups/manishad/Rpackages/")
library(metRology, lib.loc = "/home/groups/manishad/Rpackages/")


s = s[ ,!is.na(names(s)) ]

s %>% group_by(scen.name) %>%
  mutate(PhatRelBias = mean( abs(Phat - TheoryP[1])/TheoryP[1] ),
         EstVarRelBias = mean( abs(EstVar - V[1])/V[1] ) ) %>%
  group_by(calib.method) %>%
  summarise( n(),
             PhatRelBiasMn = mean(PhatRelBias),
             EstVarRelBiasMn = mean(EstVarRelBias) )
# expect n=16000 for each calib.method since there are 96 total rows (3 each per scenario)

table(is.na(s$PhatLo))
table(is.na(s$DiffLo))

# assumes a single scenario
mean(s$Phat)
table(s$TheoryP)
bias = mean(s$Phat) - s$TheoryP
mean(s$PhatBtMn, na.rm=TRUE); bias  # hope this is equal to the bias


mean(s$Diff)
table(s$TheoryDiff)
bias = mean(s$Diff) - s$TheoryDiff
mean(s$DiffBtMn, na.rm=TRUE); bias  # hope this is equal to the bias





# stitch on Sherlock
# sbatch -p qsu,normal,owners /home/groups/manishad/SAPH/stitch_sbatch_files/stitch_4.sbatch
# sacct --jobs=49474291 --format=User,JobID,account,Timelimit,elapsed,ReqMem,MaxRss,ExitCode

# # move it to Desktop
# scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/stitched_results/stitched.csv ~/Desktop

