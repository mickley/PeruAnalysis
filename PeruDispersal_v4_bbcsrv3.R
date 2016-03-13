################################################################################
## Fitting models to the bivariate and univariate data
################################################################################
rm(list=ls())
load(file='data4peruanalysis3.RData')
##nsim <- 19; ncore <- 3 ## for testing
nsim <- as.numeric(Sys.getenv("nsim"))
ncore <- as.numeric(Sys.getenv("PBS_NUM_PPN"))

################################################################################
### Subsetting data
################################################################################
## abiotic species may be a bit odd and poorly represented. Have to remove from analysis
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Abiotic==0 & Unkwn !=1)  ## biivariate
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Abiotic==0  & Unkwn !=1) ## univariate

## For analysis excluding P. laevis need to uncomment these lines
## hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c,
##                             !(Spp == 'Pseudolmedia laevis'))

## hyperdat.uni.sel.c <- subset(hyperdat.uni.sel.c,
##                             !(Spp == 'Pseudolmedia laevis'))


## Fit models
intMod.bi <- lmeHyperframe(hyperdat.bi.sel.c, 0:25,         # Bivariate
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.sel.c, 0:25,     # Univariate
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

## construct data to make predictions with envelopes for figures
preddat.int <-  expand.grid(stage=c('S', 'J'), comp=c('con', 'het'),
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
## RHS of model formula
allform <- update(formula(intMod.uni[[1]]), NULL~.)
## Make model matrix
modmat.int <-  model.matrix(allform, data=preddat.int)
## pull out conspecific - heterospecific distance
modmat.int <-  modmat.int[preddat.int$comp=='con',]-modmat.int[preddat.int$comp=='het',]

## construct data to make predictions for figures with
preddat.int <-  expand.grid(stage=c('S', 'J'), comp=c('con', 'het'),
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
## RHS of model formula
allform <- update(formula(intMod.uni[[1]]), NULL~.)
## Make model matrix
modmat.int <-  model.matrix(allform, data=preddat.int)
## pull out conspecific - heterospecific distance
modmat.int <-  modmat.int[preddat.int$comp=='con',]-modmat.int[preddat.int$comp=='het',]

## Bootstrapping
intMod.bi.boot <-  bootstrap.t.CI.lme(intMod.bi, lin.comb.Ct=modmat.int, nboot=nsim,
                                      alpha=0.05, ncore=ncore)

intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni, lin.comb.Ct=modmat.int, nboot=nsim,
                                       alpha=0.05, ncore=ncore)
