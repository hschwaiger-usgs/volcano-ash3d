#!/bin/bash

# srun_ProcessResults.sh: a bash script that processes a series of Ash3d simulations in embarrassingly
# parallel mode (each instance in serial, but lots of jobs in parallel).
# which steps through a 2-d parameter-space (time and plume height).
# This script is meant to process the jobs run with srun_Ash3d.sh and will
# process data in directory SIMULATION_DIR/Workspace/Dir[jobID].
# In this example, the output final deposit ASCII file is copied to an output
# directory, but other processing can be included such as generating figures,
# zipping output files or removing netcdf files, etc.

# Note: For batch jobs, this script includes slurm directives below and can 
#        be submitted to the slurm workload manager with the following command:
#   sbatch srun_ProcessResults.sh
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
#time=0-01:00:00       # Total allowed time in D-HH:MM:SS format; 2 days is the max.
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
  echo "Processing job #$irun via bash."
else
  irun=$SLURM_ARRAY_TASK_ID
  echo "Processing job #$irun via slurm management."
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
OUTPUTDIR=${WRKHOME}/run_output                         # directory containing output
FileDate=`date "+%Y%b%d"`                               # date, to be appended to file names
##### Names of the template Ash3d control file and the input_table.txt listing run modifications
CTRTEMPLATE=${WRKHOME}/input_files/Ash3d_template.inp
RUNTABLE=${WRKHOME}/input_files/input_table.txt
##### Names of directories that contain programs, utilities, shared data
USGSROOT=/opt/USGS
WINDROOT=/data/WindFiles
TOPOROOT=/data/TOPO/GEBCO/GEBCO_23
ASH3DHOME=${USGSROOT}/Ash3d

#### Output subdirectories.
## If the runs are done in a group of 1,000 or so at a time, each group of 1,000 will
## be stored in a directory named according to the run date.
ZIPFILEDIR=${OUTPUTDIR}/zip_files                          # location of zip files
DEPOSITDIR_ESRI=${OUTPUTDIR}/DepositFiles                  # location of DepositFiles
ARRIVALTIMESDIR=${OUTPUTDIR}/ArrivalTimes                  # location of ArrivalTimes
INFILEDIR=${OUTPUTDIR}/input_files                         # location of reformatted deposit files
MAPDIR=${OUTPUTDIR}/MapFiles                               # location of gif maps
PROFDIR=${OUTPUTDIR}/ProfileFiles

# See if the output directory exists.
echo "Checking for existence of output directory"
if test -r ${OUTPUTDIR}; then   
  echo "    ${OUTPUTDIR} exists."
else
  mkdir -p ${OUTPUTDIR}
  #echo "    Error.  ${OUTPUTDIR} does not exist.  Stopping"
  #exit 1
fi

#Check for existence of subdirectories
mkdir -p ${DEPOSITDIR_ESRI}
#mkdir -p ${ARRIVALTIMESDIR}
mkdir -p ${INFILEDIR}
#mkdir -p ${ZIPFILEDIR}
mkdir -p ${MAPDIR}
#mkdir -p ${OUTPUTDIR}/summary

#####################   Set up and run models  ####################################

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
cd ${RUNDIRS}/Dir${RunNumber}

echo "        Copying deposit file"
cp DepositFile_____final.dat ${DEPOSITDIR_ESRI}/Run${RunNumber}.txt
#sed 's/$/\r/' DepositFile_____final.dat > DepositFile_ESRI_ASCII.txt
#cp DepositFile_ESRI_ASCII.txt ${DEPOSITDIR_ESRI}/Run${RunNumber}.txt

if [ -n "$SLURM_JOB_ID" ]; then
  srun ASH3DPLOT=6 ${ASH3DHOME}/bin/Ash3d_PostProc 3d_tephra_fall.nc 5 3
else
  ASH3DPLOT=6 ${ASH3DHOME}/bin/Ash3d_PostProc 3d_tephra_fall.nc 5 3
fi
cp Ash3d_Deposit____final.png ${MAPDIR}/Ash3d_Deposit____final_Run${RunNumber}.png

#flip and rename input file
echo "        flipping and renaming input file"
cp ash3d_input.inp ${INFILEDIR}/Run${RunNumber}.inp

