#!/bin/bash

# Shell script to run all the different species 
# while setting the various parameters

# Now set PBS options
# set job name
#PBS -N peru

# keep job output
#PBS -k o

# set to the bash shell for PBS
#PBS -S /bin/bash

## Request the required number of cpus per job
#PBS -l nodes=1:ppn=48

# combine standard error and output files
#PBS -j  oe

# send email at end of job and if aborted
#PBS -m ae
# to email address
#PBS -M  robert.bagchi@uconn.edu

echo "Running  Peru spatial analyses" $PBS_JOBID

# set global parameters for analysis
# scl is size (length and breadth) of grid squares in km
## For the full analyses try using scl = 5
## or set scl to 0 which uses code in R to set scl = 5
## for species with mean dispersal < 5 and 10 othewise
export nclust=$PBS_NUM_PPN
export arrayid=$PBS_ARRAYID
export nsim=999
echo 'number of cpus requested = ' $nclust

~/programs/R/R-3.2.3/bin/R CMD BATCH ~/Peru/PeruDispersal_v4_bbcsrv3.R ~/Peru/perudisp$PBS_ARRAYID\_job$PBS_JOBID.Rout

