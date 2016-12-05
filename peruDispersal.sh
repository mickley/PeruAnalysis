## set shell
#!/bin/bash

# Shell script to run all the different species 
# while setting the various parameters

# set job name
#$ -N peruDispersal

## use current working directory
#$ -cwd

# set to the bash shell 
#$ -S /bin/bash

# send to general queue
#$ -q highmem.q

## Request the required number of cpus per job
#$ -pe smp 8

# combine standard error and output files
#$ -j y

# send email at end of job and if aborted
#$ -m ae
# to email address
#$ -M  robert.bagchi@uconn.edu

## load R
module load R/3.2.2

echo "Running  Peru spatial analyses"

# set global parameters for analysis
export nclust=$NSLOTS

## set default values for various variables
## these can be overriden from the command line
## To run the code after
## resetting nsim and rmax do
## qsub -v nsim=999,rmax=10 peruDispersal.sh
nsim=${nsim:-99}
rmax=${rmax:-15}
abiotic=${abiotic:-1}

echo 'number of cpus requested = ' $nclust

R CMD BATCH ~/Peru/July2016/PeruAnalysis/PeruDispersal_v7_bbcsrv3.R ~/Peru/July2016/progreports/perudisp\_v7\_niter$nsim\_rmax$rmax\_abiotic$abiotic.Rout

