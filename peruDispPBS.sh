#!/bin/bash

# Shell script to run all the different species 
# while setting the various parameters

## use current working directory
cd $PBS_O_WORKDIR


# Now set PBS options
# set job name
#PBS -N perudisp

# keep job output
#PBS -k o

# set to the bash shell for PBS
#PBS -S /bin/bash

## Request the required number of cpus per job
#PBS -l nodes=1:ppn=16

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
##export arrayid=$PBS_ARRAYID

## set default values for various variables
## these can be overriden from the command line
## To run the code after
## resetting nsim and rmax do
## qsub -v nsim=999,rmax=10 peruDispersal.sh

nsim=${nsim:-99}
rmax=${rmax:-15}

echo 'number of cpus requested = ' $nclust

~/programs/R/R-3.2.3/bin/R CMD BATCH ~/Peru/PeruAnalysis/PeruDispersal_v7_bbcsrv3.R ~/Peru/progreports/perudisp\_v7\_niter$nsim\_rmax$rmax\_job$PBS_JOBID.Rout

# R CMD BATCH ~/Peru/PeruDispersal_v4_bbcsrv3.R ~/Peru/perudisp$PBS_ARRAYID\_job$PBS_JOBID.Rout
