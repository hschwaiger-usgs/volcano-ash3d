#!/bin/bash

# EmbarPar_Ash3d.sh is a script that drives a batch of single-processor
# Ash3d runs, waits for the batch to complete, then launches a subsequent
# batch, repeating until the full job is complete.
# This script uses the srun_Ash3d.sh script which is designed for slurm
# management and can be used in place of a slurm manager.
#
# If you are running on a system with a slurm manager, this script is not
# necessary. Simply edit the relavent bits of srun_Ash3d.sh and run via:
#   sbatch srun_Ash3d.sh

#####################   Set up script functions ###################################
RUN=F     # if T, then run Ash3d
PROC=T    # if T, then run processing of jobs

#####################   Starting value for runs ###################################
# Total number of runs performed = dirmax*(cyclemax+1)
# Runs are numbered from RunStartNumber to (RunStartNumber + dirmax*(cyclemax+1))
RunStartNumber=1     # Run number for first run in the series
dirmax=3             # Number of simultaneous runs (1 to 50)
cyclemax=15          # Number of cycles (0 to ????)
totruns=$(( $dirmax * ($cyclemax + 1) ))
echo "totruns=$totruns"

#  A bit of error-checking
NPROC=`nproc --all`
if (( $NPROC < $dirmax )); then
  echo "ERROR: More simultaneous jobs are requested than cpu's on this system."
  echo "       jobs = $dirmax"
  echo "       cpus = $NPROC"
  exit 1
fi

#####################   Set up and run models  ####################################
if [[ "$RUN" == "T" ]]; then
  for (( icycle=0;icycle<=$cyclemax;icycle++ )); do
    # Write table headers for input values
    echo "Setting up and running models for cycle $icycle of $cyclemax"
    # Start looping through directories
    for (( idir=1;idir<=$dirmax;idir++ )); do
      irun=`echo "$RunStartNumber - 1 + $icycle * $dirmax + $idir" | bc -l`
      ./srun_Ash3d.sh $irun > logfile_ash3d_$irun.txt 2>&1 &
    done
    echo "All done setting up and starting jobs. Waiting . . ."
    wait            #wait until background jobs have completed
  done
else
  echo "Skipping running Ash3d jobs"
fi

#####################   Now processing runs    ####################################
if [[ "$PROC" == "T" ]]; then
  for (( icycle=0;icycle<=$cyclemax;icycle++ )); do
    echo "Processing models for cycle $icycle of $cyclemax"
    # Start looping through directories
    for (( idir=1;idir<=$dirmax;idir++ )); do
      irun=`echo "$RunStartNumber - 1 + $icycle * $dirmax + $idir" | bc -l`
      ./srun_ProcessResults.sh $irun > logfile_proc_$irun.txt 2>&1 &
    done
    wait            #wait until background jobs have completed
  done
else
  echo "Skipping processing Ash3d jobs"
fi




