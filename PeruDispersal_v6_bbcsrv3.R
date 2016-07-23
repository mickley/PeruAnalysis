################################################################################
## Fitting models to the bivariate and univariate data
################################################################################
rm(list=ls())
library(parallel)
## library(Rmpi)
## library(snow) ## the snow library integrates mpi capabilities into lapply etc
library(ReplicatedPointPatterns)

print(paste('code run on', Sys.Date()))
rseed <- 1234 ## set the seed for reproducibility

##nboot <- 19; ncore <- 3;nopl <- 0;  ## for testing
rmax <- 15
nopl <- as.numeric(Sys.getenv('nopl')) ## exclude Pseudolmedia laevis?
ncore <- as.numeric(Sys.getenv("nclust"))
nboot <- as.numeric(Sys.getenv('nsim')) ## number of simulations

################################################################################
## Load pre-processed data file
################################################################################

load(file='../data/data4peruanalysisv6.RData')

################################################################################
### Subsetting data
################################################################################
## abiotic species may be a bit odd and poorly represented. Have to remove from analysis
##hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Abiotic!=1 & Unkwn !=1)  ## biivariate
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel,  Unkwn !=1)  ## biivariate
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel,  Unkwn !=1) ## univariate


## Removing Pseudolmedia laevis if required
if(nopl==1){
    hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c,
                                !(Spp == 'Pseudolmedia laevis'))
    hyperdat.uni.sel.c <- subset(hyperdat.uni.sel.c,
                                 !(Spp == 'Pseudolmedia laevis'))
}

## Fit models
intMod.bi <- lmeHyperframe(hyperdat.bi.sel.c, 0:rmax,         # Bivariate
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.sel.c, 0:rmax,     # Univariate
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

## construct data to make predictions with envelopes for figures
preddat.int <-  expand.grid(stage=c('S', 'J'),
                            comp=c('con', 'het'),
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
## RHS of model formula
allform <- update(formula(intMod.uni[[1]]), NULL~.)

## Make model matrix
modmat.int <-  model.matrix(allform, data=preddat.int)

## pull out conspecific - heterospecific distance
modmat.int <-  modmat.int[preddat.int$comp=='con',]-modmat.int[preddat.int$comp=='het',]

## Bootstrapping
intMod.bi.boot <-  bootstrap.t.CI.lme(intMod.bi, lin.comb.Ct=modmat.int, nboot=nboot,
                                      alpha=0.05, ncore=ncore, iseed=rseed, cltype="PSOCK")

intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni, lin.comb.Ct=modmat.int, nboot=nboot,
                                       alpha=0.05, ncore=ncore, iseed=rseed, cltype="PSOCK")

results <- list('bi' = list('model' = intMod.bi, 'boot' = intMod.bi.boot),
                'uni' = list('model '= intMod.uni, 'boot' = intMod.uni.boot),
                'preddat' = preddat.int)

objname <- paste0('results.nopl', nopl)
assign(results, objname)

save(list=objname, file="PeruDispersal.RData")
