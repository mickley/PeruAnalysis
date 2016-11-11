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

##nboot <- 19; ncore <- 3; rmax <- 15;  ## for testing
rmax <- as.numeric(Sys.getenv('rmax')) ## maximum distance
# nopl <- as.numeric(Sys.getenv('nopl')) ## exclude Pseudolmedia laevis?
# noqw <- as.numeric(Sys.getenv('noqw')) ## exclude Quarraribea laevis?
ncore <- as.numeric(Sys.getenv("nclust"))
nboot <- as.numeric(Sys.getenv('nsim')) ## number of simulations

print(paste("No. cpus=", ncore, 
            "No. simulations =", nboot))


################################################################################
## Load pre-processed data file
################################################################################

load(file='../data/data4peruanalysisv6.RData')

################################################################################
### Subsetting data
################################################################################
## abiotic species may be a bit odd and poorly represented. Have to remove from analysis
##hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Abiotic!=1 & Unkwn !=1) ## bivariate
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel,  Unkwn !=1)
 ## univariate
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel,  Unkwn !=1)


## Exploring model fits suggests some species are big outliers in the  bivariate analysis
## We are removing them here.
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c, !((sp.id == 108 & pname=='BM') |
                                                   (sp.id==314 & pname=='CC2') |
                                                   (sp.id==38 & pname=='LA') |
                                                   (sp.id==144 & pname=='LA'))
)


## Fit models
## Bivariate
intMod.bi <- lmeHyperframe(hyperdat.bi.sel.c, 0:rmax,         
                         fixed="stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)
## Univariate
intMod.uni <- lmeHyperframe(hyperdat.uni.sel.c, 0:rmax,    
                         fixed="stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)


preddat.int <-  expand.grid(stage=c('S', 'J'), 
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
## RHS of model formula
allform <- update(formula(intMod.uni[[1]]), NULL~.)
## Make model matrix
modmat.int <-  model.matrix(allform, data=preddat.int)

preddat.int <-  expand.grid(stage=c('S', 'J'), 
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
## RHS of model formula
allform <- update(formula(intMod.uni[[1]]), NULL~.)
## Make model matrix
modmat.int <-  model.matrix(allform, data=preddat.int)

## Bootstrapping
system.time(intMod.bi.boot <-  bootstrap.t.CI.lme(intMod.bi, lin.comb.Ct=modmat.int, nboot=nboot, alpha=0.05,
                                                  ncore=7))

system.time(intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni, lin.comb.Ct=modmat.int, nboot=nsim,
                                                   alpha=0.05, ncore=7))


## Bootstrapping
intMod.bi.boot <-  bootstrap.t.CI.lme(intMod.bi,
                                      lin.comb.Ct=modmat.int,
                                      nboot=nboot, alpha=0.05,
                                      ncore=ncore, iseed=rseed,
                                      cltype="PSOCK")

intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni,
                                       lin.comb.Ct=modmat.int,
                                       nboot=nboot, alpha=0.05,
                                       ncore=ncore, iseed=rseed,
                                       cltype="PSOCK")

results <- list('bi' = list('model' = intMod.bi,
                  'boot' = intMod.bi.boot),
                'uni' = list('model '= intMod.uni,
                  'boot' = intMod.uni.boot),
                'preddat' = preddat.int)
## Save object
objname <- 'results'

assign(objname, results)
dir.name <- paste0('../results/results_', Sys.Date())

system(paste('mkdir -p', dir.name))

save(list=objname, file=paste0(dir.name, '/',
                     paste0("PeruDispersal",
                             '_rmax', rmax,
                            ".RData")
                     ))
