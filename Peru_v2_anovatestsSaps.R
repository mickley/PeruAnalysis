rm(list=ls())

##devtools::install_github('robertbagchi/ReplicatedPointPatterns')

library(parallel)
library(ReplicatedPointPatterns)

print(paste('code run on', Sys.Date()))
set.seed(1234)
## nopl <- 0; intlevel <- 2; ncore <- 3; nboot <- 99
nopl <- as.numeric(Sys.getenv('nopl')) ## include Pseudolmedia laevis?
intlevel <- as.numeric(Sys.getenv('arrayid')) ## interaction level
ncore <-  as.numeric(Sys.getenv('nclust')) ## number of cpus
nboot <- as.numeric(Sys.getenv('nsim')) ## number of simulations

print(paste("interaction =", intlevel, 
	"No. cpus=", ncore, 
	"No. simulations =", nboot, 
	"P. laevis excluded (1=yes, 0=no) =", nopl))



type <- paste("No_Plaevis =", nopl, "sqrtNxNyweights")

## Load data
##load(file='data4peruanalysis3.RData')
load(file='~/Peru/data4peruanalysis3.RData')

## Remove unknown and abiotic dispersed species from within and between cohort analyses
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel,  Abiotic !=1 & Unkwn !=1) 
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Abiotic!=1  & Unkwn !=1)

hyperdat.bi.saps.c <- subset(hyperdat.bi.sel.c, stage == 'S')
hyperdat.uni.saps.c <- subset(hyperdat.uni.sel.c, stage == 'S')


## Removing Pseudolmedia laevis if required
if(nopl==1){
    hyperdat.bi.saps.c <- subset(hyperdat.bi.saps.c,
                                !(Spp == 'Pseudolmedia laevis'))
    hyperdat.uni.saps.c <- subset(hyperdat.uni.saps.c,
                                 !(Spp == 'Pseudolmedia laevis'))
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

intMod.bi <- lmeHyperframe(hyperdat.bi.saps.c, 0:25,
                                fixed=fixform,
                                random="1|site/Spp",
                                computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.saps.c, 0:25,
                                 fixed=fixform,
                                 random="1|site/Spp",
                                 computeK=FALSE)

tests.uni <- lapply(terms, function(term) {
    list(
        ## 'd10' = bootstrap.compare.lme(mods = intMod.uni,
        ##                               term=term,
        ##                               dists=1:10, nboot=nboot, ncore=ncore),
        'd25' =  bootstrap.compare.lme(mods = intMod.uni,
                                       term=term,
                                       dists=1:25, nboot=nboot, ncore=ncore),
        'type' = type
    )})
tests.uni <- list(anovas=tests.uni, model=intMod.uni)

tests.bi <- lapply(terms, function(term){
    list(
        ## 'd10' = bootstrap.compare.lme(mods = intMod.bi,
        ##                               term=term,
        ##                               dists=1:10, nboot=nboot, ncore=ncore),
        'd25' = bootstrap.compare.lme(mods = intMod.bi,
                                      term=term,
                                      dists=1:25, nboot=nboot, ncore=ncore),
        'type' = type)
})
tests.bi <- list(anovas=tests.bi, model=intMod.bi)

anovaSaps <- list('uni' = tests.uni, 'bi' = tests.bi, intlevel=intlevel)

objname <- paste0("anovaSaps", intlevel)
assign(objname, anovaSaps)

save(list=objname, file=paste0('~/Peru/Peru_', intlevel, '_wayAnovaSaps', '_nopl', nopl,'.RData'))
    

