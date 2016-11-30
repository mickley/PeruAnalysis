## set shell
#!/bin/bash

# Shell script to run on mulitple cores
# while setting the various parameters

# Now set PBS options

# set job name

#PBS -N  peruAOVSaps

# keep job output
#PBS -k o

# set to the bash shell for PBS
#PBS -S /bin/bash

## Request the required number of cpus per job
#PBS -l nodes=1:ppn=16

# combine standard error and output files
#PBS -j  oe

# send email at beginning, end of job and if aborted
#PBS -m bae
# to email address
#PBS -M  robert.bagchi@uconn.edu

## Iterate through interaction levels
#PBS -t 1-2
echo "Running  Peru spatial analysis Anova " $PBS_JOBID


export nclust=$PBS_NUM_PPN
export arrayid=$PBS_ARRAYID

## use current working directory
cd $PBS_O_WORKDIR

## set default values for various variables
## these can be overriden from the command line
## To run the code after
## resetting nsim and rmax do
## qsub -v nsim=999,rmax=10 PeruAnalysisAnova_juvs_PBS.sh

nsim=${nsim:-99}
rmax=${rmax:-15}
abiotic=${rmax:-1}

echo 'number of cpus requested = ' $nclust
echo 'rmax = ' $rmax
echo 'abiotic = ' $abiotic

# execute the R commands in the R script
~/programs/R/R-3.2.3/bin/R CMD BATCH ~/Peru/PeruAnalysis/Peru_anovatestsSaps_v7.R ~/Peru/progreports/peru_anovaSaps\_v7\_niter$nsim\_rmax$rmax\_$arrayid.Rout
