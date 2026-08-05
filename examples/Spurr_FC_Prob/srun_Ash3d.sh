#!/bin/bash

# srun_Ash3d.sh: a bash script that drives a series of Ash3d simulations in embarrassingly
# parallel mode (each Ash3d instance in serial, but lots of jobs in parallel) with
# parameters determined from an input table (${RUNTABLE}) which modifies select parameters
# in a template control file (${CTRTEMPLATE}). These tabular values can either be a
# systematic marching through parameter-space, or with stochastic variables with
# values drawn from a preferred distribution for Monte Carlo analysis. The run table,
# input_table.txt, can be generated with a script such as gen_table_param_march.sh,
# which steps through a 2-d parameter-space (time and plume height).
# This script with create a run directory SIMULATION_DIR/Workspace/Dir[jobID].
# All output will remain in those run directories with the expectation that a
# second slurm script, srun_ProcessResults.sh, with the same SBATCH settings,
# will be run to process the output, e.g. rename output files, zip log or ASCII files,
# remove netcdf output if too big, etc.

# Note: For batch jobs, this script includes slurm directives below and can 
#        be submitted to the slurm workload manager with the following command:
#   sbatch srun_Ash3d.sh
#
#        Alternatively, you can run this without a slurm manager by calling from a
#        wrapper script.
#
#################################################
# SLURM hints
#  To see the status of all jobs, type "squeue"
#  To see status of all jobs with user, type "squeue -u USERNAME"
#  To cancel a slurm-managed job, type "scancel job#"
#  To see HPC infrastructure (partitians, CPUs, time-limits, etc.), type "sinfo"
#  To see the resources, cpu's in usage etc., type "sstat job#"
#  To evaluate resources used for a completed job (e.g. to optimize future requests), type:
#    sacct --user=USERNAME --start=2025-09-4 --format=jobid,jobname,partition,alloccpus,elapsed,totalcpu,maxrss,reqmem,state,exitcode

# sbatch directives; detailed in https://slurm.schedmd.com/sbatch.html
#  Allocations		Specify an allocation account if you have multiple	-A	--account=account_no
#  Partitions		Specify a partition					-p	--partition=partition_name
#  Sending email	Receive email at beginning or end of job completion		--mail-type=type
#  Email address	Email address to receive the email				--mail-user=user
#  Number of nodes	The number of nodes needed to run the job		-N	--nodes=nodes
#  Number of tasks	The total number of cores needed to run the job		-n	--ntasks=processes
#  Minimum memory	Minimum memory (Mb) required per usable allocated CPU.		--mem-per-cpu=
#  Bind tasks		Bind tasks according to application hints.			--hint=compute_bound
#  Wall time		The max. amount of time your job will run for		-t	--time=wall-time
#  Job Name		Name your job so you can identify in queue			--job-name=<jobname>
#  Output file		Standard output of slurm script				-o	--output=filename

# Bits to edit/verify
#job-name=Ash3d-Spurr  # Job name
#account=vhp           # Account should be listed on ssh splash screen
#array=1-48%12         # This is the range of run IDs from the input_table.txt file: 1-100 with 20 at a time.
#time=0-06:00:00       # Total allowed time in D-HH:MM:SS format; 2 days is the max.
                       #  = expected max time for an individual run * safety fac (~ 1.2)

#################################################
### sbatch settings for a Hovenweep job-array
#################################################
#SBATCH --job-name=Ash3d_Spurr               ## Only a label, but nice to identify in squeue
#SBATCH --array=1-48%12                      ## List of job numbers to run, followed by % and # of simul.jobs
#SBATCH --time=0-03:00:00                    ## Anticipated time of individual job * 1.2
#SBATCH --output=%A-%a-Ash3d-Spurr.out       ## 'A' is the job ID, 'a' is the array #job number
#SBATCH --mail-user=USERNAME@DOMAIN.NAME     ## 
#SBATCH --account=vhp                        ## Job will not run without accounting
#SBATCH --partition=cpu                      ## 
#SBATCH --mail-type=ALL                      ## email user for all event types
#SBATCH --ntasks=1                           ## srun jobs managed by sbatch so this script is 1-job only
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread                 ## This prevents a 1-core job from hogging the whole cpu
##-----------------------------------------------
#### Additional bits if we are running Ash3d with OpenMP
##SBATCH --cpus-per-task=4                   ## num of threads for OpenMP executables
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK ## actually set the environment variable so as not to hog threads
#################################################
## then submit with sbatch 

###THINGS TO CHECK BEFORE STARTING A RUN
# 1) Check slurm settings
#      --job-name=Ash3d_Spurr    # Not critical, but useful identifier in job list
#      --array=1-48%12           # Specific runs to launch in range or comma format: 1-10 or 1,4,5,..
#                                #  The bit after the % indicates the number of simultaneous jobs
#      --time=0-03:00:00         # Make sure this is adaquate for each Ash3d instance
#
# 2) Make sure your Ash3d input files are complete. e.g. control file works and is in input_files/
#    also any other auxillary file such as Sites.txt as the Airport/POI file.
#    You should probably test a run with the template control file. Make sure all the output files that
#    you want are generated, but no more. In particular, if you write out netcdf files, be judicious about
#    writing the full concentraion array (e.g. Use 'yes 2   #Write out 3-D ash concentration at specified times?' in
#    block 4, line 15; with the '2' indicating only output products, not full concentration).
#    In block 3, line 2; consider if your source might need wind values above that provided by the NWP file.
#    In block 3, line 4; only stop the run at 99% mass deposited if you are not running a cryptotephra case.
# 
# 3) Make sure you have generated an input_table.txt file. You can use gen_table_param_march.sh as an example.
#    You should test this input_table.txt file with the template control file by running:
#      Ash3d_ASCII_GenCTR Ash3d_template.inp input_table.txt 1
#    This will generate a new control file using the parameters from run=1 of the input_table.txt file.
#
# 4) Edit the section below (Local environment, workspace, input files, etc.) to suit your specifications.
#    In particular, make sure ${WRKHOME} is correct for your job.
# 
# 5) Maybe run this script with --array=1%1  to verify that the directories, soft links and overall program-flow
#    behave as expected.
#
# 6) Lastly, make sure your post-processing script (srun_ProcessResults.sh) has the appropriate slurm settings
#    and moves files around, zips things, and otherwise processes your data to your liking.


rc=0         # Initialize return code to 0
# Parsing command-line arguments
NARGS=$#

#################################################
# module load everything needed for running: shouldn't need anything for the Ash3d part, but maybe post-processing
#################################################

#  Get the runID for this instance
irun=0
if [ $NARGS -gt 0 ]; then
  irun=$1
  echo "Running job #$irun via bash."
else
  irun=$SLURM_ARRAY_TASK_ID
  echo "Running job #$irun via slurm management."
fi
if [ $irun -eq 0 ]; then
  echo "Run number not set. Exiting."
  exit
fi

###############################################################################
##### Local environment, workspace, input files, etc. #########################
#####    This is the user-editable section            #########################
###############################################################################
runmax=48    # Maximum runID for this script
             # This should be <= the runs in input_table.txt and should be consistent with #SBATCH --array=?????
##### Names of directories where runs are performed
WRKHOME=~/work/USGS/Software/GIT/volcano-ash3d/examples/Spurr_FC_Prob
RUNDIRS=${WRKHOME}/Workspace                              # directory where runs are performed
##### Names of the template Ash3d control file and the input_table.txt listing run modifications
CTRTEMPLATE=${WRKHOME}/input_files/Ash3d_template.inp
RUNTABLE=${WRKHOME}/input_files/input_table.txt
LOCATIONFILE=${WRKHOME}/input_files/Sites.txt
##### Names of directories that contain programs, utilities, shared data
USGSROOT=/opt/USGS
WINDROOT=/data/WindFiles
TOPOROOT=/data/TOPO/GEBCO/GEBCO_23
ASH3DHOME=${USGSROOT}/Ash3d
###############################################################################
###############################################################################
###############################################################################


#####################   Make the directory  ###################################
# Make run numbers into five-digit numbers
if [[ ${irun} -lt 10 ]]; then
  RunNumber="000${irun}"
elif [[ ${irun} -lt 100 ]]; then
  RunNumber="00${irun}"
elif [[ ${irun} -lt 1000 ]]; then
  RunNumber="0${irun}"
else
  RunNumber="${irun}"
fi
if test -r  ${RUNDIRS}/Dir${RunNumber}; then
  rm -f ${RUNDIRS}/Dir${RunNumber}/*
else
  mkdir -p ${RUNDIRS}/Dir${RunNumber}
fi
cd    ${RUNDIRS}/Dir${RunNumber}
ln -s ${WINDROOT} Wind
#ln -s ${LOCATIONFILE} .

# make the new input file
#echo "${ASH3DHOME}/bin/tools/Ash3d_ASCII_GenCTR ${CTRTEMPLATE}  ${RUNTABLE} ${RunNumber}"
if [ -n "$SLURM_JOB_ID" ]; then
  srun ${ASH3DHOME}/bin/tools/Ash3d_ASCII_GenCTR ${CTRTEMPLATE}  ${RUNTABLE} ${RunNumber}
else
  ${ASH3DHOME}/bin/tools/Ash3d_ASCII_GenCTR ${CTRTEMPLATE}  ${RUNTABLE} ${RunNumber}
fi
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
    echo "Error running Ash3d_ASCII_GenCTR for irun=${irun}: rc=$rc"
    echo "Check output file for more information; maybe runID is too high for the table."
    exit 1
fi

#run the model
export ASH3DVERB=7                                      # quash all but essential stdout (only stderr and logfile)
#${ASH3DHOME}/bin/Ash3d ash3d_input.inp                 # used for testing
if [ -n "$SLURM_JOB_ID" ]; then
  srun ${ASH3DHOME}/bin/Ash3d ash3d_input.inp
else
  ${ASH3DHOME}/bin/Ash3d ash3d_input.inp
fi

