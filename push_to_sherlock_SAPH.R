



####################### CHECK IN ####################### 
# see the sbatches
cd /home/groups/manishad/SAPH/sbatch_files

sbatch -p qsu,owners,normal /home/groups/manishad/SAPH/sbatch_files/1.sbatch

# check on my running or pending jobs
squeue -u mmathur -t RUNNING
squeue -u mmathur -t PENDING

# /home/groups/manishad/SAPH/sbatch_files/1.sbatch


# see the datasets
vim /home/groups/manishad/SAPH/sim_results/long/long_results_job_1_.csv
cd /home/groups/manishad/SAPH/sim_results/long
ls -l . | egrep -c '^-'

###### See the Errors #####
vim /home/groups/manishad/SAPH/sbatch_files/rm_1507.err

# see the scen parameters
nano /home/groups/manishad/SAPH/scen_params.csv

# see the stitched results
nano /home/groups/manishad/SAPH/sim_results/overall_stitched/sti*
  
  
  
# CODE -> SHERLOCK ----------------------------------

# push helper.SAPH
scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/helper_SAPH.R mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH

scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(SAPH\)/Code\ \(git\)/2021-5-9\ nonaffirm\ MLEs\ on\ Sherlock/doParallel.R mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH


# SHERLOCK -> DESKTOP ----------------------------------

scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/stitched_results/stitched.csv ~/Desktop


scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/* ~/Desktop


####################### SHERLOCK -> DESKTOP (DEBUGGING) ####################### 

# move error file to Desktop
scp -r mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/sbatch_files/rm_1250.err ~/Desktop

# move one sbatch file to Desktop
scp -r mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPH/sbatch_files/2296.sbatch ~/Desktop


####################### RUN SBATCH ####################### 

# run one of them
sbatch -p qsu,normal,owners /home/groups/manishad/SAPH/sbatch_files/1.sbatch




####################### RESULTS -> DESKTOP FOR ANALYSIS ####################### 

scp mmathur@login.sherlock.stanford.edu /home/groups/manishad/SAPH/results/overall_stitched/stitched.csv ~/Desktop

####################### CLEAN UP ####################### 

# clean up the directory
rm /home/groups/manishad/SAPH/sim_results/*
  
  rm /home/groups/manishad/SAPH/sim_results/overall_stitched/*
  
  
# clean up "rm" files
  rm /home/users/mmathur/rm*
  rm /home/groups/manishad/SAPH/rm*
  rm /home/groups/manishad/SAPH/sbatch_files/rm*
  rm /home/groups/manishad/SAPH/sbatch_files/slurm*
  
  
  # remove all sbatches
  #  rm -r /home/groups/manishad/SAPH/sbatch_files/*
  