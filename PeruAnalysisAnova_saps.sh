## set shell
#!/bin/bash

# Now set SGE options

# set job name
#$ -N peruAnovaSaps

# use current working directory
#$ -cwd

## set to bash shell
#$ -S /bin/bash

# send to high memory queue
#$ -q all.q

## Request the required number of cpus per job
#$ -pe smp 8

# combine standard error and output files
#$ -j y

# send email at end (e) of the job of if aborted (a)
#$ -m ae
# to email address
#$ -M robert.bagchi@uconn.edu

## Iterate through interaction levels
#$ -t 1-3

## load R
module load R/3.2.2

echo "Running  Peru spatial analyses Anova"

# set global parameters for analysis
export nclust=$NSLOTS
export arrayid=$SGE_TASK_ID

## set default values for various variables
## these can be overriden from the command line
## To run the code with default values of noqw and nopl but
## resetting nsim and rmax do
## qsub -v nsim=999,rmax=10 PeruAnalysisAnova_saps.sh
${noqw:0}
${nopl:0}
${nsim:99}
${rmax:15}

echo 'number of cpus requested = ' $nclust

# execute the R commands in the R script R_script
R CMD BATCH ~/Peru/July2016/PeruAnalysis/Peru_anovatestsSaps_v6.R ~/Peru/July2016/peru_anovaSaps\_v6\_pl$nopl\_qw$noqw\_niter$nsim\_rmax$rmax\_$arrayid.Rout
