rm(list=ls())

##devtools::install_github('robertbagchi/ReplicatedPointPatterns')

library(parallel)
library(ReplicatedPointPatterns)

print(paste('code run on', Sys.Date()))
rseed <- 1234 ## set the seed for reproducibility

## nopl <- 0; noqw <- 0; intlevel <- 2; ncore <- 3; nboot <- 5

intlevel <- as.numeric(Sys.getenv('arrayid')) ## interaction level
rmax <- as.numeric(Sys.getenv('rmax')) ## maximum distance
ncore <-  as.numeric(Sys.getenv('nclust')) ## number of cpus
nboot <- as.numeric(Sys.getenv('nsim')) ## number of simulations

print(paste("interaction =", intlevel,
            "max distance =", rmax,
            "No. simulations =", nboot, 
            "No. cpus=", ncore))

type <- "sqrtNxNyweights"

## Load data
load(file='../data/data4peruanalysisv7.RData')

## Remove unknown dispersal  species from within and between cohort analyses
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel,  Unkwn !=1) 
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Unkwn !=1)

## alternative version - removing abiotic dispersed species
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c,  Abiotic !=1) 
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel.c, Abiotic !=1)


####################################################################
## Fitting the models
####################################################################

## Now setting up multiple anova tests, dependent on the value of intlevel which is defined in
## the shell script
covs  <- c('stage', 'huntpres', 'HSD')
fixform <- paste(covs, collapse='+')

if(intlevel > 1)
    fixform <- paste0("(", paste(covs, collapse='+'), ")^", intlevel)

## Define terms to assess
terms <- apply(combn(covs, intlevel), 2, paste, collapse=':')

intMod.bi <- lmeHyperframe(hyperdat.bi.sel.c, 0:rmax,
                                fixed=fixform,
                                random="1|site/Spp",
                                computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.sel.c, 0:rmax,
                                 fixed=fixform,
                                 random="1|site/Spp",
                                 computeK=FALSE)

tests.uni <- lapply(terms, function(term) {
    list(
         'dtable' =  bootstrap.compare.lme(mods = intMod.uni,
           term=term, dists=list(1:rmax, 1:5, 6:10, 11:15),
           nboot=nboot, ncore=ncore,
           iseed=rseed)
         )})

tests.uni <- list(anovas=tests.uni, model=intMod.uni, type=type)

tests.bi <- lapply(terms, function(term){
  list(
       'dtable' = bootstrap.compare.lme(mods = intMod.bi,
         term=term,
         dists=list(1:rmax, 1:5, 6:10, 11:15),
         nboot=nboot, ncore=ncore, iseed=rseed)
       )})

anovaJuvs <- list('uni' = tests.uni, 'bi' = tests.bi,
                  intlevel=intlevel, type=type)

objname <- paste0("anovaJuvs_noabiotic", intlevel)
assign(objname, anovaJuvs)


dir.name <- paste0('../results/results_', Sys.Date())
system(paste('mkdir -p', dir.name))

save(list=objname, file=paste0(dir.name, '/', 'Peru_v7_',
                     intlevel, '_wayAnovaJuvs', '_rmax', rmax,
                     '_nsim', nboot, '_NoAbiotic', '.RData'))
