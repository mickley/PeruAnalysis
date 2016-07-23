rm(list=ls())

##devtools::install_github('robertbagchi/ReplicatedPointPatterns')

library(parallel)
library(ReplicatedPointPatterns)

print(paste('code run on', Sys.Date()))
rseed <- 1234 ## set the seed for reproducibility

## nopl <- 0; noqw <- 0; intlevel <- 2; ncore <- 3; nboot <- 5
nopl <- as.numeric(Sys.getenv('nopl')) ## exclude Pseudolmedia laevis?
noqw <- as.numeric(Sys.getenv('noqw')) ## exclude Quarribea wittii?
intlevel <- as.numeric(Sys.getenv('arrayid')) ## interaction level
ncore <-  as.numeric(Sys.getenv('nclust')) ## number of cpus
nboot <- as.numeric(Sys.getenv('nsim')) ## number of simulations
rmax <- 15
print(paste("interaction =", intlevel, 
	"No. cpus=", ncore, 
	"No. simulations =", nboot, 
            "P. laevis excluded (1=yes, 0=no) =", nopl,
            "Q. wittii excluded (1=yes, 0=no) =", noqw))

type <- paste("No_Plaevis =", nopl, "sqrtNxNyweights")

## Load data
load(file='../data/data4peruanalysisv6.RData')
##load(file='~/Peru/data4peruanalysis3.RData')

## Remove unknown and abiotic dispersed species from within and between cohort analyses
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel,  Unkwn !=1) 
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Unkwn !=1)

hyperdat.bi.saps.c <- subset(hyperdat.bi.sel.c, stage == 'S')
hyperdat.uni.saps.c <- subset(hyperdat.uni.sel.c, stage == 'S')


## Removing Pseudolmedia laevis if required
if(nopl==1){
    hyperdat.bi.saps.c <- subset(hyperdat.bi.saps.c,
                                !(Spp == 'Pseudolmedia laevis'))
    hyperdat.uni.saps.c <- subset(hyperdat.uni.saps.c,
                                 !(Spp == 'Pseudolmedia laevis'))
}

## Removing Quarribea wittii if required
if(noqw==1){
    hyperdat.bi.saps.c <- subset(hyperdat.bi.saps.c,
                                !(Spp == 'Quararibea wittii'))
    hyperdat.uni.saps.c <- subset(hyperdat.uni.saps.c,
                                 !(Spp == 'Quararibea wittii'))
}
## Fitting the models

## Now setting up multiple anova tests, dependent on the value of intlevel which is defined in
## the shell script
covs  <- c('comp', 'huntpres', 'HSD')
fixform <- paste(covs, collapse='+')

if(intlevel > 1)
    fixform <- paste0("(", paste(covs, collapse='+'), ")^", intlevel)

## Define terms to assess
terms <- apply(combn(covs, intlevel), 2, paste, collapse=':')
terms <- as.list(terms[grep('comp', terms)])

intMod.bi <- lmeHyperframe(hyperdat.bi.saps.c, 0:rmax,
                                fixed=fixform,
                                random="1|site/Spp",
                                computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.saps.c, 0:rmax,
                                 fixed=fixform,
                                 random="1|site/Spp",
                                 computeK=FALSE)

tests.uni <- lapply(terms, function(term) {
  list(
       ## 'd10' = bootstrap.compare.lme(mods = intMod.uni, term=term,
       ##   dists=1:10, nboot=nboot,
       ##   ncore=ncore, iseed=rseed),
       'd15' =  bootstrap.compare.lme(mods = intMod.uni,
         term=term, dists=1:15,
         nboot=nboot, ncore=ncore,
         iseed=rseed)
       )})

tests.uni <- list(anovas=tests.uni, model=intMod.uni, type=type)

tests.bi <- lapply(terms, function(term){
    list(
         ## 'd10' = bootstrap.compare.lme(mods = intMod.bi,
         ##   term=term,
         ##   dists=1:10, nboot=nboot, ncore=ncore, iseed=rseed),

         'd15' = bootstrap.compare.lme(mods = intMod.bi,
           term=term,
           dists=1:15, nboot=nboot, ncore=ncore,iseed=rseed)
)})

tests.bi <- list(anovas=tests.bi, model=intMod.bi)

anovaSaps <- list('uni' = tests.uni, 'bi' = tests.bi,
                  intlevel=intlevel, type=type)

objname <- paste0("anovaSaps", intlevel)
assign(objname, anovaSaps)

dir.name <- paste0('../results/results_', Sys.Date())
system(paste('mkdir -p', dir.name))

save(list=objname, file=paste0(dir.name, 'Peru_v6_',
                     intlevel, '_wayAnovaSaps_nopl', nopl,
                     '_noqw', noqw, '.RData'))
    

