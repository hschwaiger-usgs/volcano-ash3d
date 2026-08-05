#!/bin/bash

# EmbarPar_Ash3d.sh is a script that drives a batch of single-processor
# Ash3d runs, waits for the batch to complete, then launches a subsequent
# batch, repeating until the full job is complete.
# This script uses the srun_Ash3d.sh script which is designed for slurm
# management and can be used in place of a slurm manager.
#
# If you are running on a system with a slurm manager, simple edit the
# relavent bits of srun_Ash3d.sh and run via:
#   sbatch srun_Ash3d.sh

##### STARTING VALUE FOR RUNS
# Total number of runs performed = dirmax*(cyclemax+1)
# Runs are numbered from RunStartNumber to (RunStartNumber + dirmax*(cyclemax+1))
RunStartNumber=1     # Run number for first run in the series
dirmax=3             # Number of simultaneous runs (1 to 50)
cyclemax=15          # Number of cycles (0 to ????)
totruns=$(( $dirmax * ($cyclemax + 1) ))
echo "totruns=$totruns"

#####################   Set up and run models  ####################################
for (( icycle=0;icycle<=$cyclemax;icycle++ )); do
  # Write table headers for input values
  echo "setting up and running models"
  # Start looping through directories
  for (( idir=1;idir<=$dirmax;idir++ )); do
    irun=`echo "$RunStartNumber - 1 + $icycle * $dirmax + $idir" | bc -l`
    ./srun_Ash3d.sh $irun > logfile.txt 2>&1 &
  done
  echo "All done setting up and starting jobs. Waiting . . ."
  wait            #wait until background jobs have completed
done

# Processing
for (( icycle=0;icycle<=$cyclemax;icycle++ )); do
  echo "Processing models"
  # Start looping through directories
  for (( idir=1;idir<=$dirmax;idir++ )); do
    irun=`echo "$RunStartNumber - 1 + $icycle * $dirmax + $idir" | bc -l`
    ./srun_ProcessResults.sh $irun
  done
  wait            #wait until background jobs have completed
done





