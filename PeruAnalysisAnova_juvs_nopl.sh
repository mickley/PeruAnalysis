## set shell
#!/bin/bash

# Now set SGE options

# set job name
#$ -N peruAnovaJuvs

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
#$ -t 2-4

## load R
module load R/3.2.2

echo "Running  Peru spatial analyses Anova"

# set global parameters for analysis
export nclust=$NSLOTS
export arrayid=$SGE_TASK_ID
export nopl=1
export noqw=0
export nsim=99

echo 'number of cpus requested = ' $nclust

# execute the R commands in the R script R_script
R CMD BATCH ~/Peru/July2016/PeruAnalysis/Peru_anovatestsJuvs_v6.R ~/Peru/July2016/peru_anovaJuvs_v6\_pl$nopl\_$arrayid.Rout
