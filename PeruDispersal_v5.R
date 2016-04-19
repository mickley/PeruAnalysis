################################################################################
### File: PeruDispersal.R
### Version: 5
### Author: Robert Bagchi
### Email: robert.bagchi@uconn.edu
### Date last modified: 15 April 2016
### Description: Code to analyse spatial patterns of saplings and juvenile
### trees in six Peruvian forests plots subjected to varying levels of
### defaunation. This is the version used immediately
### prior to submission of the paper.
################################################################################

################################################################################
### Preliminaries
################################################################################
rm(list=ls())


library(spatstat)
library(ggplot2)
library(cowplot)
library(grid)
library(gtable)
library(reshape)
## if you need to install
##devtools::install_github('robertbagchi/ReplicatedPointPatterns')
library(ReplicatedPointPatterns)

################################################################################
## Define some useful functions
################################################################################

## function to convert K functions to L functions
K2L <- function(k,r) {
    L <- sqrt(k/(pi))
    return(cbind(r=r, L=L))}

## a function that takes the output from the model and puts it in
## a convenient format for plotting
kfunclme2plot <- function(res, preddat)
{
    do.call('rbind', mapply( function(kx, lx, ux, preddat, nmx)
        {
            cbind(distance=as.numeric(nmx), K=as.vector(kx), lcl=lx, ucl=ux, preddat)
        },

                            kx=res$estimator,
                            lx=res$lower,
                            ux=res$upper,
                            nmx=as.list(names(res$estimator)),
                            MoreArgs=list(preddat=preddat),
                            SIMPLIFY=FALSE) )
}
## A function that plots the BLUEs against distance
plot.kfunctionlme <- function(x){
  require(ggplot2)
  dists <- as.numeric(names(x))
  dat <- rbind(dists, sapply(x, fixef))
  parnames <- rownames(dat)[-1]
  dat <- melt(dat, ~distance)
  names(dat)[c(1, 3)] <-  c("Parameter", "K")
  dat <- dat[dat$Parameter !='dists',]
  dat$Parameter <- factor(dat$Parameter, levels=parnames)
  dat$K[dat$Parameter == "(Intercept)"] <- dat$K[dat$Parameter == "(Intercept)"] -
    pi*dat$distance[dat$Parameter == "(Intercept)"]^2
  ggplot(dat, aes(x=distance, y=K)) + geom_line() +
    facet_wrap(~Parameter, scales="free_y")  + theme_bw()
}


## Function that standardises a vector
standardise <-  function(x) (x - mean(x))/sd(x)

## Function to extract random effect BLUPs to allow plotting
extractModRanefs <-  function(mod, level=1){
    require(reshape2)
    ranefs <- do.call('cbind', lapply(mod, function(mod) ranef(mod)[[level]]))
    ranefs$id <- rownames(ranefs)
    names(ranefs)[1:length(mod)] <-  paste('distance', 1:length(mod), sep='.')
    ranefs <-  melt(ranefs, id='id')
    ranefs$distance <- as.numeric(do.call('rbind',
                                          strsplit(as.character(ranefs$variable), split='.', fixed=T))[,2])
    return(ranefs)
}


################################################################################
## ggplot theme for paper
################################################################################
theme_cust <- theme(axis.text = element_text(size=7),
                    axis.title.y=element_text(size=10, face='bold'),
                    plot.margin=unit(c(0.5, 0.3, 0.5, 0 ),  'lines'),
                    strip.background=element_rect(fill='white'),
                    strip.text=element_text(size=8),
                    legend.text=element_text(size=8),
                    legend.title=element_text(size=10),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank())

#########################################################################
### Read in and process data
#########################################################################
## read in the data
trees <- read.csv("../../data/mdd6plots_data_4jul2013.csv") ## tree locations

dispersal <-  read.csv("../../data/MDDspp_dispersal_syndromes_CURR.csv") ## dispersal data

## put into 1/0 format
dispersal[, c('LV', 'SV', 'SB', 'Bat', 'TR', 'Wind', 'Expl', 'Unkwn', 'HSD', 'HID')] <-
  apply(dispersal[, c('LV', 'SV', 'SB', 'Bat', 'TR', 'Wind', 'Expl', 'Unkwn', 'HSD', 'HID')],
        2, function(x) ifelse(is.na(x), 0, x))

## combine wind and explosive dispersal into a single variable
dispersal$Abiotic <-  dispersal$Wind + dispersal$Expl

## site information
sitedat <- data.frame(site=c("BM",  "CashuTr12", "CashuTr3",  "LA", "RA","TRC"),
                             forest=c('BM', 'Cashu', 'Cashu', 'LA', 'RA', 'TRC'),
                             pname = factor(c('BM', 'CC2', 'CC1', 'LA', 'RA', 'TRC'),
                                            levels=c('CC1', 'CC2', 'TRC', 'LA', 'BM', 'RA')),
                             hunted=c('hunted', 'intact', 'intact',
                                     'hunted', 'hunted', 'intact'),
                             lg.primates=c(0, 37.5, 37.5, 3.9, 0, 15.9),
                             lg.birds=c(2.1, 29.2, 29.2, 18.2, 8, 44))
## construct hunting pressure variable
sitedat$huntpres <- with(sitedat, (standardise(lg.primates) + standardise(lg.birds))*(-1))

sitedat[order(sitedat$pname),] ## quick look

## First deal with duplicated points (same site, coordinates, code and Type) -
## sometimes these have different tags, but for our purposes, 2 points at the
## same location cause problems, so we add a little bit of noise to them
## (less than the resolution of the original data)
dups <- duplicated(trees[,c('Site', 'E', 'N', 'code', 'Type')])
sum(dups)## a total of 102 such points out of 67000, so not many
trees$E[dups]  <-  trees$E[dups] + runif(sum(dups), -0.1, 0.1)
trees$N[dups]  <-  trees$N[dups] + runif(sum(dups), -0.1, 0.1)

trees.old <- trees ## make a copy for comparison
trees <-trees[-c(grep('Astrocaryum', trees$Spp),
                 grep('Attalea', trees$Spp),
                 grep('Iriartea', trees$Spp)),] ## removing palms
## warnings can be ignored - some string issues, but they still get removed
## I double checked
table(trees$Type)

## now remove adults that are really too small to be adults. 
## Adults must have a dbh > 10 and > min(median dbh or >30 cm)
tree.dbhcut <- aggregate(diam~Spp, data=trees, subset=Type=='A',
                         function(x) {
                             cutoff <- max(10, median(x, na.rm=T))
                             return(min(30, cutoff))
                         })
summary(tree.dbhcut)
names(tree.dbhcut)[2] <- 'diam.cutoff'
hist(tree.dbhcut$diam.cutoff)

trees <- merge(trees, tree.dbhcut, all=T) ## merge in the cutoff point
names(trees)
## now remove all "adults" below the cut-off size
trees <- trees[!(trees$Type == 'A' & trees$diam < trees$diam.cutoff),]

## remove any species that don't have at least 2 individuals in at least 1 site
splist <- apply(table(trees$Spp, trees$Site), 1, function(x) any(x>2))
trees <- trees[trees$Spp %in% names(splist)[splist],]
nlevels(trees$Spp)
trees <- droplevels(trees) ## get rid of unused levels that might cause problems
                                        #later
nlevels(trees$Spp) ## leaves us with 906 species

## split the data by site
trees <- split(trees, f=trees$Site)

## Pull out some summary stats
sapply(trees, nrow)
sapply(trees, function(x) length(unique(x$Spp)))

## defining the window at CC12 which is not quite rectangular
cc12  <- owin(poly=list(x=c(225, 75, 75, 24.9, 25, 0, 0, 225),
                y= c(200, 200, 125, 125, 75, 75, 0, 0)))

plot(cc12, axes=T) ## looks good
##plot(rotate(winlist$CashuTr12$J, 3*pi/180))
## defining a list of windows for each site
## Each comprises of a window for each size class.
winlist <- list(BM=list(A=owin(c(0, 200), c(0, 200)),
                  J=owin(c(50, 150), c(50, 150)),
                  S=owin(c(50, 150), c(50, 150))),
                CashuTr12=list(A=cc12,
                  J=owin(c(75, 175), c(50, 150)),
                  S=owin(c(75, 175), c(50, 150))),
                CashuTr3=list(A=owin(c(0, 200), c(0, 200)),
                  J=owin(c(55, 145), c(55, 145)),
                  S=owin(c(55, 145), c(55, 145))),
                LA=list(A=owin(c(0, 200), c(0, 200)),
                  J=owin(c(50, 150), c(50, 150)),
                  S=owin(c(50, 150), c(50, 150))),
                RA=list(A=owin(c(0, 200), c(0, 200)),
                  J=owin(c(50, 150), c(50, 150)),
                  S=owin(c(50, 150), c(50, 150))),
                TRC=list(A=owin(c(0, 200), c(0, 200)),
                  J=owin(c(50, 150), c(50, 150)),
                  S=owin(c(50, 150), c(50, 150))))

par(mfrow=c(2,3))
lapply(winlist, function(x) { plot(x$A, axes=T, asp=1)
                              plot(x$J, add=T, border=2)}) ## looks ok

## A function that works out the list of species at a site that
## have at least one pair in all classes.
getSpeciesList <- function(dat){
    counts <- table(dat$code, dat$Type)
    counts <- counts[counts[,'A']>1 & counts[,'J']>1 & counts[, 'S']>1,]
    counts <- counts[rownames(counts)!='-999',]
    counts <- counts[rev(order(counts[, 'A'])),]
    rownames(counts)}


## make a species list for each plot (based on the species that have sufficient
## stems in each size class (i.e. at least 1 pair or two stems)
splists <- lapply(trees, function(dat) getSpeciesList(dat))

## turn the data from each species into a point pattern, split up by
## age class
## names(ppp.sps$CashuTr12)
## a few adult points in CC12 are right on the boundary.
## this leads to the warnings below, that can be
## ignored
ppp.sps <- mapply(function(spid, treedat, win){
    ppp.sp <- lapply(as.list(spid), function(id, treedat, win){
        condat <- treedat[treedat$code==id,] ## subset data from species
        ppp.con <- mapply(function(Tx, Wx){ ## turn into ppp object
            pppx <- ppp(x=Tx$E, y=Tx$N, window=Wx)
            return(pppx)}, Tx=split(condat, f=condat$Type), Wx=win, SIMPLIFY=FALSE)
        ## note that species are split by size-class (Type) when input to the function
        ## now heterospecifics
        hetdat <- treedat[treedat$code!=id,] ## subset out all heterspecifics
        ppp.het <- mapply(function(Tx, Wx){
            Tx <- Tx[!duplicated(Tx[,c('E', 'N')]),] ## remove colocated stems
            pppx <- ppp(x=Tx$E, y=Tx$N, window=Wx)

            return(pppx)}, Tx=split(hetdat, f=hetdat$Type), Wx=win, SIMPLIFY=FALSE)
        names(ppp.het) <- paste(names(ppp.het), 'h', sep='') ## add 'h' suffix to het names
        ppp.sp <- c(ppp.con, ppp.het) ## combine cons and hets
        return(ppp.sp)},
                     treedat=treedat, win=win)
    names(ppp.sp) <- spid
    return(ppp.sp)
  }, spid=splists, treedat=trees, win=winlist, SIMPLIFY=FALSE)

## Univariate K functions
hyperdat.uni  <- mapply(function(dat, win)
    {

        hyper <- hyperframe(
            sp.id = rep(names(dat), 4), ## 6 times for S+J x con+het
            stage = factor(rep(c('S', 'J'), each=2*length(dat)), levels=c('S', 'J')),
            comp = rep(rep(c('con', 'het'), each=length(dat)), 2),
            pppx=c(rep(lapply(dat, function(x) x$S), 2),
                rep(lapply(dat, function(x) x$J),2)),
            pppy=c(lapply(dat, function(x) x$S),
                lapply(dat, function(x) x$Sh),
                lapply(dat, function(x) x$J),
                lapply(dat, function(x) x$Jh))
            )

        smallwindow <- win$S ## smallwindow <- winlist[[1]]$S
        hyper$N1 <- sapply(hyper$pppx, npoints)
        hyper$N2 <- sapply(hyper$pppy, npoints)

        hyper <- cbind.hyperframe(hyper,
                                  dispersal[match(hyper$sp.id, dispersal$code),
                                            c('Spp', 'LV', 'SV', 'SB', 'Bat',
                                              'TR', 'Wind', 'Expl',
                                              'Unkwn', 'Abiotic', 'HSD', 'HID')])

        ppp.s <- with(hyper, superimpose(i=pppx, j=pppy, W=pppy$window))

        hyper$K <- lapply(ppp.s, function(ppp.s)
            Kmulti(ppp.s, I=marks(ppp.s)=='i', J= marks(ppp.s)=='j',
                   r=0:25, correction="border", ratio=TRUE))
        hyper$K[hyper$comp=='con'] <-
            lapply(hyper$pppx[hyper$comp=='con'], function(px)
                Kmulti(px, I=inside.owin(px, w=smallwindow),
                       J=inside.owin(px, w=Window(px)), r=0:25, correction="border",
                       ratio=TRUE))

        hyper$wts <-  mapply(function(pppx, pppy)
            kfunc.weights.calc(pppx=pppx, pppy=pppy, r=0:25,
                               correction='border', type='sqrtnxny_A'),
                             pppx=hyper$pppx, pppy=hyper$pppy, SIMPLIFY=FALSE)
        return(hyper)
    }, dat = ppp.sps, win=winlist, SIMPLIFY=FALSE)

## Add site names
hyperdat.uni <- mapply(function(dat, site) {
                       dat$site <- site
                       return(dat)
                   }, dat=hyperdat.uni, site=as.list(names(hyperdat.uni)),
                   SIMPLIFY=FALSE)

## work out which species in which sites have insufficient data
## This might be because too many individuals are on the border and
## so are removed when border corrections are applied.
univarKmods.bysite <-  lapply(hyperdat.uni, function(dat){
    lmeHyperframe(dat, 0:25, fixed="stage*comp", random="1|site/Spp",
                  computeK=FALSE, weights.type='nxny_A', minsamp=2)})
##  These are returned by the function - remove them from data set
removedsp <- sapply(univarKmods.bysite, function(x) attr(x, 'removed.species'))


## remove species with insufficient data
hyperdat.uni.sel <- mapply(function(dat, sp) {
    sp.rm <- unique(dat[sp, 'sp.id']$sp.id)
    print(length(sp.rm))
    dat <- dat[!(dat$sp.id %in% sp.rm),]
    return(dat)},  hyperdat.uni, removedsp, SIMPLIFY=FALSE)

## make one hyperframe with all the data
hyperdat.uni.sel <- do.call('rbind', hyperdat.uni.sel)

## add site information
hyperdat.uni.sel <- cbind.hyperframe(
    hyperdat.uni.sel,
    sitedat[match(hyperdat.uni.sel$site, sitedat$site),
                   c('forest', 'pname', 'hunted', 'lg.primates', 'lg.birds', 'huntpres')])

## Convert to a factor from character
hyperdat.uni.sel$site <- factor(hyperdat.uni.sel$site)
summary(hyperdat.uni.sel)


## Bivariate K-function data frame
hyperdat.bi  <- lapply(ppp.sps, function(dat)
    {
        hyper <- hyperframe(
            sp.id = rep(names(dat), 4), ## 4 times for S+J x con+het
            stage = factor(rep(c('S', 'J'), each=2*length(dat)), levels=c('S', 'J')),
            comp = rep(rep(c('con', 'het'), each=length(dat)), 2),
            pppx=c(rep(lapply(dat, function(x) x$S), 2),
                rep(lapply(dat, function(x) x$J),2)),
            pppy = rep(c(lapply(dat, function(x) x$A),
                lapply(dat, function(x) x$Ah)), 2)
            )
        hyper$N1 <- sapply(hyper$pppx, npoints)
        hyper$N2 <- sapply(hyper$pppy, npoints)

        hyper <- cbind.hyperframe(hyper,
                                  dispersal[match(hyper$sp.id, dispersal$code),
                                            c('Spp', 'LV', 'SV', 'SB', 'Bat', 'TR',
                                              'Wind', 'Expl',
                                              'Unkwn', 'Abiotic','HSD', 'HID')])
        ppp.s <- with(hyper, superimpose(i=pppx, j=pppy, W=pppy$window))
        hyper$K <- lapply(ppp.s, function(x)
            Kmulti(x, I=x$marks =='i', J=x$marks=='j',
                   r=0:25, corr='border', ratio=TRUE))

        hyper$wts <-  mapply(function(pppx, pppy)
            kfunc.weights.calc(pppx=pppx, pppy=pppy, r=0:25,
                               correction='border', type='sqrtnxny_A'),
                             pppx=hyper$pppx, pppy=hyper$pppy, SIMPLIFY=FALSE)

        return(hyper)
    })

## add site names
hyperdat.bi <- mapply(function(dat, site) {
                       dat$site <- site
                       return(dat)
                   }, dat=hyperdat.bi, site=as.list(names(hyperdat.bi)),
                   SIMPLIFY=FALSE)

## find species and site combinations with too few data - technically
## don't need to do this (could set minsamp=0), but that would mean we
## analysed different species for the univariate and bivariate case.
## This way both are the same.
bivarKmods.bysite <-  lapply(hyperdat.bi, function(dat){
    lmeHyperframe(dat, 0:25, fixed="stage*comp", random="1|site/Spp",
                  computeK=FALSE, weights.type='sqrtnxny_A', minsamp=2)})

removedsp.bi <- sapply(bivarKmods.bysite, function(x) attr(x, 'removed.species'))

## mapply(function(sp, dat)
##     unique(dat[sp,]$Spp), removedsp, hyperdat.uni, SIMPLIFY=FALSE)
## Remove species with too few individuals
hyperdat.bi.sel <- mapply(function(dat, sp) {
    sp.rm <- unique(dat[sp, 'sp.id']$sp.id)
    print(length(sp.rm))
    dat <- dat[!(dat$sp.id %in% sp.rm),]
    return(dat)},  hyperdat.bi, removedsp.bi, SIMPLIFY=FALSE)


## data set with all data
hyperdat.bi.sel <- do.call('rbind', hyperdat.bi.sel)

hyperdat.bi.sel <- cbind.hyperframe(
    hyperdat.bi.sel,
    sitedat[match(hyperdat.bi.sel$site, sitedat$site),
                   c('forest', 'pname', 'hunted', 'lg.primates', 'lg.birds', 'huntpres')])

hyperdat.bi.sel$site <- factor(hyperdat.bi.sel$site)

## Save the objects required for analysis so we don't
## always have to repeat processing (and for transfer to cluster)
save('dispersal', "extractModRanefs", "getSpeciesList",
          'hyperdat.uni.sel', 'hyperdat.bi.sel', 'K2L', 'kfunclme2plot',
          'plot.kfunctionlme', 'ppp.sps', 'theme_cust', 'trees', 'winlist',
          'sitedat', 'standardise', file='data4peruanalysis4.RData')



################################################################################
## Fitting models to the bivariate and univariate data
################################################################################
##rm(list=ls())
##load(file='data4peruanalysis4.RData')
## Can start from here and jump to models to save time

################################################################################
### Subsetting data
################################################################################

## abiotic species may be a bit odd and poorly represented. Have to remove from analysis
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Abiotic!=1 & Unkwn !=1) 

summary(hyperdat.bi.sel.c)

## Now univariate K funcs
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Abiotic!=1  & Unkwn !=1)


## Extracing some summaries of the data to report in the paper
## Some numbers on number of individuals/species etc etc

## mean reliance on different dispersers
disp.list <-  aggregate(cbind(LV, SV, SB, Bat, TR) ~Spp, subset(hyperdat.uni.sel.c,
                                                                stage=='S' & comp=='con'), mean)

apply(disp.list[, -1], 2, function(x) sum(x>0)) ## number dependent to some degree
apply(disp.list[, -1], 2, function(x) sum(x==1)) ## number totally dependent

aggregate(N1 ~ comp, hyperdat.uni.sel.c, sum) ## number of stems 
aggregate(N1 ~ stage*comp, hyperdat.uni.sel.c, sum) ## number of stems per cohort
aggregate(Spp ~ 1, hyperdat.uni.sel.c, function(x) length(unique(x))) ## number of species
aggregate(Spp ~ pname, hyperdat.uni.sel.c, function(x) length(unique(x))) # species per site
mean(aggregate(Spp ~ pname, hyperdat.uni.sel.c, function(x) length(unique(x)))$Spp) ## mean S
## total number of adults
aggregate(N2 ~ stage*comp, hyperdat.bi.sel.c,  sum)

## For analysis excluding P. laevis need to uncomment these lines
## hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c,
##                             !(Spp == 'Pseudolmedia laevis'))

## hyperdat.uni.sel.c <- subset(hyperdat.uni.sel.c,
##                             !(Spp == 'Pseudolmedia laevis'))


## Fitting the models

## set number of iterations for bootstrapping - keep small for trials (takes a long time)
## Code was run on a grid computer for paper
nsim <- 19

## First sapling only models
hyperdat.bi.sap.c <- subset(hyperdat.bi.sel.c, stage=='S')
hyperdat.uni.sap.c <- subset(hyperdat.uni.sel.c, stage=='S')

sapMod.bi <- lmeHyperframe(hyperdat.bi.sap.c, 0:25,
                         fixed="comp*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

sapMod.uni <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
                         fixed="comp*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)
## ## set up model matrix for predictions
preddat.sap <-  expand.grid(stage=c('S'), comp=c('con', 'het'),
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
sapform <- update(formula(sapMod.uni[[1]]), NULL~.)
modmat.sap <-  model.matrix(sapform, data=preddat.sap)

## ## Want to show differences between conspecifics and heterospecifics, so set up
## ## model matrix to do this
modmat.sap <-  modmat.sap[preddat.sap$comp=='con',]-modmat.sap[preddat.sap$comp=='het',]

## ## Bootstrapping
sapMod.bi.boot <-  bootstrap.t.CI.lme(sapMod.bi, lin.comb.Ct=modmat.sap, nboot=nsim, alpha=0.05,
                                    ncore=7)

sapMod.uni.boot <-  bootstrap.t.CI.lme(sapMod.uni, lin.comb.Ct=modmat.sap, nboot=nsim,
                                       alpha=0.05, ncore=7)

## ## Extract data for plot (redundant because later plot including juveniles
## includes both cohorts - however helpful for initial analysis
sap.bi.plotdat <- kfunclme2plot(sapMod.bi.boot, preddat=subset(preddat.sap, comp=='con'))
sap.uni.plotdat <- kfunclme2plot(sapMod.uni.boot, preddat=subset(preddat.sap, comp=='con'))
names(sap.bi.plotdat)[2:4] <- paste(names(sap.bi.plotdat)[2:4], 'bi', sep='.')
names(sap.uni.plotdat)[2:4] <- paste(names(sap.uni.plotdat)[2:4], 'uni', sep='.')

## ## COmbine data from within and between cohort analyses
sap.plotdat <- merge(sap.uni.plotdat, sap.bi.plotdat,
                     by=c('distance', 'stage', 'comp', 'huntpres', 'HSD'), all=T)
## ## Make labels more reader friendly
sap.plotdat$huntpres <-  factor(sap.plotdat$huntpres, labels=c('Low', 'High'))
sap.plotdat$HSD <-  factor(sap.plotdat$HSD, labels=c('Insensitive', 'Sensitive'))

## ## Make plots
pl.sap.bi <- ggplot(sap.plotdat, aes(x=distance, y=K.bi, ymin=lcl.bi, ymax=ucl.bi)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD~  huntpres) +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])))

pl.sap.uni <- ggplot(sap.plotdat, aes(x=distance, y=K.uni, ymin=lcl.uni, ymax=ucl.uni)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD ~ huntpres) +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
   scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])))

plot_grid(pl.sap.uni, pl.sap.bi, labels= c('Uni', 'bi'))


## Anova type analysis for table 2 - run on cluster. Repeated here with
## very few iterations to demonstrate code. 999 used on cluster.

## main effects
sapMod.uni1way <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
                         fixed="comp+huntpres+HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

anova1waysaps <-  bootstrap.compare.lme(sapMod.uni1way, term='comp',
                          dists=1:25, nboot=nsim, ncore=7)

## 2 way interactions
sapMod.uni2way <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
                         fixed="(comp+huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE)

anova2waysaps <- sapply(c('comp:huntpres', 'comp:HSD'), function(term)
{
    bootstrap.compare.lme(sapMod.uni2way, term=term,
                          dists=1:25, nboot=nsim, ncore=7)
}, simplify=FALSE)



## 3 way interaction
anova3waysaps <- bootstrap.compare.lme(sapMod.uni, term="comp:huntpres:HSD", dists=1:25, nboot=nsim, ncore=7)

anova1waysaps[1:2]
sapply(anova2waysaps, function(x) x[1:2])
anova3waysaps[1:2]

## load the results of the anovas - run on computing grid
## This will only run after the other analyses have been run and results saved.
## This analysis includes all except for wind dispersed species.

## function to extract the data we need from these results
extractAnovas <- function(Aobj, type){
    Aobj <- Aobj[[type]]$anovas
    sapply(Aobj, function(x) x$d25[c("term", "D", "p")])
}

### saplings only
## with P laevis
##system("scp bbcsrv3:~/Peru/March2016/results/Peru_v4_?_wayAnovaSaps_nopl0.RData ../../March16/results/")
## extract the relevant data
saps.Anova.withpl <- sapply(c('uni', 'bi'), function(pairing){
    saps.AnovaUni_withpl <- t(Reduce('cbind', sapply(1:3, function(i){
    cat(i)
    objname <- load(paste0("../../March16/results/Peru_v4_", i, "_wayAnovaSaps_nopl0.RData"))
    Aobj <- get(objname)
    print(Aobj[c('intlevel', 'type')]) ## check it's the correct data
    extractAnovas(Aobj, pairing)
})))}, simplify=FALSE)
saps.Anova.withpl

## without P laevis
##system("scp bbcsrv3:~/Peru/March2016/results/Peru_v4_?_wayAnovaSaps_nopl1.RData ../../March16/results/")
## check effect of excluding Q wittii
saps.Anova.nopl <- sapply(c('uni', 'bi'), function(pairing){
    saps.AnovaUni_withpl <- t(Reduce('cbind', sapply(1:3, function(i){
    cat(i)
    objname <- load(paste0("../../March16/results/Peru_v4_", i, "_wayAnovaSaps_nopl1.RData"))
    Aobj <- get(objname)
    print(Aobj[c('intlevel', 'type')]) ## check it's the correct data
    extractAnovas(Aobj, pairing)
})))}, simplify=FALSE)
saps.Anova.nopl

## Full models with Juveniles and saplings (pair2)
intMod.bi <- lmeHyperframe(hyperdat.bi.sel.c, 0:25,
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.sel.c, 0:25,
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

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
intMod.bi.boot <-  bootstrap.t.CI.lme(intMod.bi, lin.comb.Ct=modmat.int, nboot=nsim, alpha=0.05,
                                    ncore=7)

intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni, lin.comb.Ct=modmat.int, nboot=nsim,
                                       alpha=0.05, ncore=7)

## can load results from cluster
##system("scp bbcsrv3:~/Peru/March2016/results/PeruDispersal_v4_nopl_0.RData ../../March16/results/")

objname <- load('../../March16/results/PeruDispersal_v4_nopl_0.RData')
intMod.boot <- get(objname)
intMod.bi.boot <- intMod.boot$bi$boot
intMod.uni.boot <- intMod.boot$uni$boot
preddat.int <- intMod.boot$preddat

################################################################################
## Now plotting the results
################################################################################
## extract predictions and cis from model output
int.bi.plotdat <- kfunclme2plot(intMod.bi.boot, preddat=subset(preddat.int, comp=='con'))
int.uni.plotdat <- kfunclme2plot(intMod.uni.boot, preddat=subset(preddat.int, comp=='con'))
names(int.bi.plotdat)[2:4] <- paste(names(int.bi.plotdat)[2:4], 'bi', sep='.')
names(int.uni.plotdat)[2:4] <- paste(names(int.uni.plotdat)[2:4], 'uni', sep='.')
## Combine within and between cohort results
int.plotdat <- merge(int.uni.plotdat, int.bi.plotdat, by=c('distance', 'stage', 'comp', 'huntpres', 'HSD'), all=T)

## More reader friendly labels
int.plotdat$stage <- factor(int.plotdat$stage, labels=c('Saplings', 'Juveniles'))
int.plotdat$huntpres <-  factor(int.plotdat$huntpres, labels=c('Low', 'High'))
int.plotdat$HSD <-  factor(int.plotdat$HSD, labels=c('Insensitive', 'Sensitive'))

## stat = identity hardwires the order - fixing so that isn't a problem
int.plotdat <- int.plotdat[order(int.plotdat$distance),]

## Make plots

pl.int.bi <-
    ggplot(int.plotdat, aes(x=distance, y=K.bi, ymin=lcl.bi, ymax=ucl.bi,
                    colour=stage, fill=stage)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD~  huntpres, scales="free_y") +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)',
         y=expression(paste('[', (hat(K)(r)[con]-hat(K)(r)[het]), ']')),
        color='Cohort', fill='Cohort')
pl.int.bi


pl.int.uni <- ggplot(int.plotdat, aes(x=distance, y=K.uni, ymin=lcl.uni, ymax=ucl.uni,
                                      colour=stage, fill=stage)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD~  huntpres, scales="free_y") +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)',
         y=expression(paste('[', (hat(K)(r)[con]-hat(K)(r)[het]), ']')),
         color='Cohort', fill='Cohort')
pl.int.uni

## Format figures (1 and 3)
z.bi <- ggplot_gtable(ggplot_build(pl.int.bi))

##gtable_show_layout(z.bi)
## add a label at the top and side
z.bi <- gtable_add_rows(z.bi, height=unit(z.bi$heights[[3]], 'cm'), pos=2)
z.bi <- gtable_add_cols(z.bi, width=unit(z.bi$widths[[7]], 'cm'), pos=7)
z.bi <- gtable_add_cols(z.bi, width=z.bi$widths[[2]], pos=1)

z.bi <- gtable_add_grob(z.bi,
                     list(rectGrob(gp = gpar(fill=NA, col=NA, lwd=0)),
                          textGrob("Relative conspecific clustering", rot = 90,
                                   gp = gpar(fontsize=12))),
                     t=3, l=2, b=9, r=2, name = paste(runif(2)))

z.bi <- gtable_add_grob(z.bi,
                     list(rectGrob(gp=gpar(fill='grey', col=NA, lwd=0.2)),
                          textGrob("Hunting pressure",
                                   gp = gpar(fontsize=12))),
                     t=3, l=5, b=3, r=7, name = paste(runif(2)))
z.bi <- gtable_add_grob(z.bi,
                     list(rectGrob(gp = gpar(fill='grey', col=NA, lwd=0.2)),
                          textGrob("Disperser sensitivity", rot = -90,
                                   gp = gpar(fontsize=12))),
                     t=5, l=9, b=7, r=9, name = paste(runif(2)))
ggdraw(z.bi)


## Univariate plot
z.uni <- ggplot_gtable(ggplot_build(pl.int.uni+theme_cust))
##gtable_show_layout(z.uni)
## add a label at the top and side
## extra row
z.uni <- gtable_add_rows(z.uni, height=unit(z.uni$heights[[3]], 'cm'),  pos=2)
## extra column
z.uni <- gtable_add_cols(z.uni, width=unit(z.uni$widths[[7]], 'cm'), pos=7)
z.uni <- gtable_add_cols(z.uni, width=z.uni$widths[[2]], pos=1)

z.uni <- gtable_add_grob(z.uni,
                     list(rectGrob(gp = gpar(fill=NA, col=NA, lwd=0)),
                          textGrob("Relative  conspecific clustering", rot = 90,
                                   gp = gpar(fontsize=12))),
                     t=3, l=2, b=9, r=2, name = paste(runif(2)))

## add hunting pressure label
z.uni <- gtable_add_grob(z.uni,
                         list(rectGrob(gp=gpar(fill='grey', col=NA, lwd=0.2)),
                              textGrob("Hunting pressure",
                                       gp = gpar(fontsize=12))),
                         t=3, l=5, b=3, r=7, name = paste(runif(2)))
## add disperser sensitivity label
z.uni <- gtable_add_grob(z.uni,
                         list(rectGrob(gp = gpar(fill='grey', col=NA, lwd=0.2)),
                              textGrob("Disperser sensitivity", rot = -90,
                                       gp = gpar(fontsize=12))),
                         t=5, l=9, b=7, r=9, name = paste(runif(2)))

ggdraw(z.uni)

## Do not run this unless you want to change figures (and have done at least 999 iterations
## commented out to prevent mistakes
## pdf(file='../huntingpaper/figures/Figure1_uniwithpl.pdf', width=4.5, height=4)
## ggdraw(z.uni)
## dev.off()

## png(file='../huntingpaper/figures/Figure1_uniwithpl.png', width=4.5, height=4, units='in', res=600)
## ggdraw(z.uni)
## dev.off()

## pdf(file='../huntingpaper/figures/Figure3_biwithpl.pdf', width=4.5, height=4)
## ggdraw(z.bi)
## dev.off()

## png(file='../huntingpaper/figures/Figure3_biwithpl.png', width=4.5, height=4, units='in', res=600)
## ggdraw(z.bi)
## dev.off()

### anovas with all stems
## with P laevis
##system("scp bbcsrv3:~/Peru/March2016/results/Peru_v4_?_wayAnovaJuvs_nopl0.RData ../../March16/results/")
## extract the relevant data
all.Anova.withpl <- sapply(c('uni', 'bi'), function(pairing){
    t(Reduce('cbind', sapply(2:4, function(i){
        cat(i)
        objname <- load(paste0("../../March16/results/Peru_v4_", i, "_wayAnovaJuvs_nopl0.RData"))
        Aobj <- get(objname)
        print(Aobj[c('intlevel', 'type')]) ## check it's the correct data
        extractAnovas(Aobj, pairing)
    })))}, simplify=FALSE)
all.Anova.withpl

## without P laevis
##system("scp bbcsrv3:~/Peru/March2016/results/Peru_v4_?_wayAnovaJuvs_nopl1.RData ../../March16/results/")
all.Anova.nopl <- sapply(c('uni', 'bi'), function(pairing){
    t(Reduce('cbind', sapply(2:4, function(i){
        cat(i)
        objname <- load(paste0("../../March16/results/Peru_v4_", i, "_wayAnovaJuvs_nopl1.RData"))
        Aobj <- get(objname)
        print(Aobj[c('intlevel', 'type')]) ## check it's the correct data
        extractAnovas(Aobj, pairing)
    })))}, simplify=FALSE)
all.Anova.nopl


## Anova type analysis for table 3 - run on cluster for 999 iterations
## (see Peru_anovatestsJuvs.R, but code  provided here too
## main effects
juvMod.uni1way <- lmeHyperframe(hyperdat.uni.sel.c, 0:25,
                         fixed="(comp+stage + huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE
)
anova2wayjuvs <-  bootstrap.compare.lme(juvMod.uni1way, term='comp',
                          dists=1:25, nboot=nsim, ncore=7)

## 2 way interactions
sapMod.uni2way <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
                         fixed="(comp+huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE)

anova2waysaps <- sapply(c('comp:huntpres', 'comp:HSD'), function(term)
{
    bootstrap.compare.lme(sapMod.uni2way, term=term,
                          dists=1:25, nboot=nsim, ncore=7)
}, simplify=FALSE)



## 3 way interaction
anova3waysaps <- bootstrap.compare.lme(sapMod.uni, term="comp:huntpres:HSD", dists=1:25, nboot=nsim, ncore=3)

anova1waysaps[1:2]
sapply(anova2waysaps, function(x) x[1:2])
anova3waysaps[1:2]

################################################################################
## Analyses of individual species
## For figure 2 and supplement figure 2
################################################################################
summary(hyperdat.uni.sel)

hyperdat.uni.ind <- subset(hyperdat.uni.sel, Abiotic < 1 & Unkwn < 1)

## K functions based on < 15 individuals are meaningless on their own
## (although they can contribute with low weight to an overall analysis)
##so removing these.
hyperdat.uni.ind <- subset(hyperdat.uni.sel, N1 > 15)
hyperdat.uni.ind$Spp <- droplevels(hyperdat.uni.ind$Spp)
counts <- table(hyperdat.uni.ind$Spp, hyperdat.uni.ind$site)
counts <- as.matrix(counts)
counts <- counts[apply(counts, 1, function(x) any(x[c(1, 4, 5)] > 3) & any(x[c(2, 3, 6)] > 3)),]
sp.keep <- rownames(counts)
sp.keep

## Pulling out selected species
hyperdat.uni.ind <- subset(hyperdat.uni.sel, Spp %in% sp.keep)
hyperdat.bi.ind <- subset(hyperdat.bi.sel, Spp %in% sp.keep)
hyperdat.uni.ind <- hyperdat.uni.ind[order(hyperdat.uni.ind$N1),]
ntotals <- aggregate(N1 ~ Spp, data=subset(hyperdat.uni.ind, comp=='con'), sum)
ntotals <- ntotals[rev(order(ntotals$N1)),]
disptype <- aggregate(HSD ~ Spp, data=hyperdat.uni.ind, mean)

## Choose a few to plot based on abundance
disptype <- merge(disptype, ntotals)
disptype <- disptype[order(disptype$N1, decreasing=T),]
disptype <- disptype[order(disptype$HSD>0.5),]

## Function to pull out conspecific - heterospecific difference
## The subsetting in this function which will determine which species are included
## Including 1:4 so that 4 most abundant species are included
## If more are included we need to modify the plot code too
Kdiffcalc <- function(con, het) cbind(r=con$r, Kdiff=con$border-het$border)

ind.funcs <- do.call('rbind', lapply(sp.keep[sp.keep %in% ntotals$Spp[c(1:4)]], function(sp, hyp){
    dat <-hyp[hyp$Spp == sp,]
    dat.con <- subset(dat, comp=='con')

    dat.con$Kdiff <- mapply(Kdiffcalc,
                            con=dat.con$K, het=subset(dat, comp=='het')$K,
                            SIMPLIFY=FALSE)
    sp.dat <- do.call('rbind', lapply(1:nrow(dat.con), function(i){
        dat <- dat.con[i,]
        Kdiff <- data.frame(dat$Kdiff)
        names(Kdiff) <- c('r', 'Kdiff')
        dat <- cbind(as.data.frame(dat[, c('stage', 'N1', 'N2', 'Spp', 'site',
                                           'pname', 'hunted', 'huntpres')]), Kdiff)
        return(dat)}))
}, hyp=hyperdat.uni.sel))
warnings() ## can ignore (just formatting)
summary(ind.funcs)
ind.funcs$stage <- factor(ind.funcs$stage, labels=c('Saplings', 'Juveniles'))
##ind.funcs$Spp <-  factor(ind.funcs$Spp, levels=ntotals$Spp)
ind.funcs$Spp <-  factor(ind.funcs$Spp, levels=disptype$Spp)

ind.funcs$N <-  ntotals$N1[match(ind.funcs$Spp, ntotals$Spp)]
ind.funcs$HSD <- disptype$HSD[match(ind.funcs$Spp, disptype$Spp)]
ind.funcs$splab <- paste(ind.funcs$Spp, ' (n = ', ind.funcs$N, ')', sep='')
ind.funcs$hunted <- factor(ind.funcs$hunted, levels=c('intact', 'hunted'),
                           labels=c('Low', 'High'))

library(RColorBrewer)

sp.plot  <- ggplot(ind.funcs, aes(x=r, y=Kdiff, colour=hunted, group=stage:site), fill='white')+
    geom_line() + facet_grid(stage~ Spp, scales='free') +

    scale_colour_brewer(palette='Dark2') +
    ##  scale_colour_manual(values = rev(brewer.pal(3,"Set1"))[2:3]) +
    geom_abline(intercept=0, lty='dashed', col='grey50') + theme_bw() + theme_cust +
    labs(x='Distance, r, (m)', y=expression(hat(K)(r)[con]-hat(K)(r)[het]),
         color='Hunting\npressure')
sp.plot <- sp.plot + theme(strip.text.x=element_text(size=7, face='italic'))
sp.plot


z.sp <- ggplot_gtable(ggplot_build(sp.plot))
##gtable_show_layout(z.sp)
## add a label at the top and side
z.sp <- gtable_add_rows(z.sp, height=z.sp$heights[3], pos=2)
z.sp <- gtable_add_cols(z.sp, width=z.sp$widths[[11]], pos=11)
z.sp <- gtable_add_cols(z.sp, width=z.sp$widths[[2]], pos=1)
##sapply(z.test$grobs[grep('strip.absoluteGrob', sapply(z.test$grobs, function(x) x$name))], function(x) x$children[[1]]$gp$fill)

z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp = gpar(fill=NA, col=NA, lwd=0)),
                          textGrob("Relative conspecific clustering", rot = 90,
                                   gp = gpar(fontsize=12))),
                     t=4, l=2, b=8, r=2, name = paste(runif(2)))
z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp=gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Insensitive Dispersers",
                                   gp = gpar(fontsize=10))),
                     t=3, l=5, b=3, r=7, name = paste(runif(2)))
z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp=gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Sensitive Dispersers",
                                   gp = gpar(fontsize=10))),
                     t=3, l=9, b=3, r=11, name = paste(runif(2)))

z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp = gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Stage", rot = -90,
                                   gp = gpar(fontsize=10))),
                     t=5, l=13, b=7, r=13, name = paste(runif(2)))
ggdraw(z.sp)
## pdf('../huntingpaper/figures/Figure2_bysp.pdf', width=8, height=4)
## ggdraw(z.sp)
## dev.off()
## png('../huntingpaper/figures/Figure2_bysp.png', width=8, height=4, units='in', res=600)
## ggdraw(z.sp)
## dev.off()

## Repeat for all species to make supplemental figure

ind.funcs.all <- do.call('rbind', lapply(sp.keep[sp.keep %in% ntotals$Spp], function(sp, hyp){
    dat <-hyp[hyp$Spp == sp,]
    dat.con <- subset(dat, comp=='con')

    dat.con$Kdiff <- mapply(Kdiffcalc,
                            con=dat.con$K, het=subset(dat, comp=='het')$K,
                            SIMPLIFY=FALSE)
    sp.dat <- do.call('rbind', lapply(1:nrow(dat.con), function(i){
        dat <- dat.con[i,]
        Kdiff <- data.frame(dat$Kdiff)
        names(Kdiff) <- c('r', 'Kdiff')
        dat <- cbind(as.data.frame(dat[, c('stage', 'N1', 'N2', 'Spp', 'site',
                                           'pname', 'hunted', 'huntpres')]), Kdiff)
        return(dat)}))
}, hyp=hyperdat.uni.sel))


ind.funcs.all$stage <- factor(ind.funcs.all$stage, labels=c('Saplings', 'Juveniles'))
##ind.funcs$Spp <-  factor(ind.funcs$Spp, levels=ntotals$Spp)
## order by disperser sensitivity
ind.funcs.all$Spp <-  factor(ind.funcs.all$Spp, levels=disptype$Spp)

ind.funcs.all$N <-  ntotals$N1[match(ind.funcs.all$Spp, ntotals$Spp)]
ind.funcs.all$HSD <- disptype$HSD[match(ind.funcs.all$Spp, disptype$Spp)]
ind.funcs.all$splab <- paste(ind.funcs.all$Spp, ' (n = ', ind.funcs.all$N, ')', sep='')
ind.funcs.all$hunted <- factor(ind.funcs.all$hunted, levels=c('intact', 'hunted'),
                           labels=c('Low', 'High'))
head(ind.funcs.all)
sp.plot.all  <- ggplot(ind.funcs.all, aes(x=r, y=Kdiff, colour=hunted, group=stage:site),
                       fill='white')+
    geom_line() + facet_grid( splab + HSD  ~stage, scales='free') +

    scale_colour_brewer(palette='Dark2') +
    geom_abline(intercept=0, lty='dashed', col='grey50') + theme_bw() + theme_cust +
    labs(x='Distance, r, (m)', y=expression(hat(K)(r)[con]-hat(K)(r)[het]),
         color='Hunting\npressure')

sp.plot.all <- sp.plot.all + theme(strip.text.y=element_text(size=7, face='italic'))
sp.plot.all

z.sp.all <- ggplot_gtable(ggplot_build(sp.plot.all))
##gtable_show_layout(z.sp.all)
## add a label at the top and side
z.sp.all <- gtable_add_rows(z.sp.all, height=z.sp.all$heights[[3]], pos=2)
z.sp.all <- gtable_add_cols(z.sp.all, width=z.sp.all$widths[[8]], pos=8)
##ggdraw(z.sp.all)

##sapply(z.test$grobs[grep('strip.absoluteGrob', sapply(z.test$grobs, function(x) x$name))], function(x) x$children[[1]]$gp$fill)

z.sp.all <- gtable_add_grob(z.sp.all,
                     list(rectGrob(gp=gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Stage",
                                   gp = gpar(fontsize=10))),
                     t=3, l=4, b=3, r=6, name = paste(runif(2)))

z.sp.all <- gtable_add_grob(z.sp.all,
                     list(rectGrob(gp = gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Disperser sensitivity", rot = -90,
                                   gp = gpar(fontsize=10))),
                     t=5, l=9, b=41, r=9, name = paste(runif(2)))

z.sp.all <- gtable_add_cols(z.sp.all, width=z.sp.all$widths[[2]], pos=1)

z.sp.all <- gtable_add_grob(z.sp.all,
                     list(rectGrob(gp = gpar(fill=NA, col=NA, lwd=0)),
                          textGrob("Relative Conpecific clustering", rot = 90,
                                   gp = gpar(fontsize=10))),
                     t=5, l=2, b=40, r=2, name = paste(runif(2)))


ggdraw(z.sp.all)

## pdf('../huntingpaper/figures/FigS2_byspecies.pdf', width=10, height=40)
## ggdraw(z.sp.all)
## dev.off()


## Plot figure S1 (Univariate K functions without P. laevis)
##system("scp bbcsrv3:~/Peru/March2016/results/PeruDispersal_v4_nopl_1.RData ../../March16/results/")

objname <- load('../../March16/results/PeruDispersal_v4_nopl_1.RData')
intModnopl.boot <- get(objname)
intModnopl.uni.boot <- intModnopl.boot$uni$boot
preddatnopl.int <- intModnopl.boot$preddat

################################################################################
## Now plotting the results
################################################################################
## extract predictions and cis from model output
int.uninopl.plotdat <- kfunclme2plot(intModnopl.uni.boot,
                                     preddat=subset(preddatnopl.int, comp=='con'))
names(int.uninopl.plotdat)[2:4] <- paste(names(int.uninopl.plotdat)[2:4], 'uni', sep='.')
## Combine within and between cohort results

## More reader friendly labels
int.uninopl.plotdat$stage <- factor(int.uninopl.plotdat$stage,
                                    labels=c('Saplings', 'Juveniles'))
int.uninopl.plotdat$huntpres <-  factor(int.uninopl.plotdat$huntpres, labels=c('Low', 'High'))
int.uninopl.plotdat$HSD <-  factor(int.uninopl.plotdat$HSD,
                                   labels=c('Insensitive', 'Sensitive'))

## Make plots

pl.int.uninopl <- ggplot(int.uninopl.plotdat,
                         aes(x=distance, y=K.uni, ymin=lcl.uni, ymax=ucl.uni,
                                      colour=stage, fill=stage)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD~  huntpres, scales="free_y") +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])),
         color='Cohort', fill='Cohort')
pl.int.uninopl

## Format figure s1
z.uni.s1 <- ggplot_gtable(ggplot_build(pl.int.uninopl))
##gtable_show_layout(z.uni)
## add a label at the top and side
## extra row
z.uni.s1 <- gtable_add_rows(z.uni.s1, height=unit(z.uni.s1$heights[[3]], 'cm'),  pos=2)
## extra column
z.uni.s1 <- gtable_add_cols(z.uni.s1, width=unit(z.uni.s1$widths[[7]], 'cm'), pos=7)
## add hunting pressure label
z.uni.s1 <- gtable_add_grob(z.uni.s1,
                         list(rectGrob(gp=gpar(fill='grey', col=NA, lwd=0.2)),
                              textGrob("Hunting pressure",
                                       gp = gpar(fontsize=12, fontface='bold'))),
                         t=3, l=4, b=3, r=6, name = paste(runif(2)))
## add disperser sensitivity label
z.uni.s1 <- gtable_add_grob(z.uni.s1,
                         list(rectGrob(gp = gpar(fill='grey', col=NA, lwd=0.2)),
                              textGrob("Disperser sensitivity", rot = -90,
                                       gp = gpar(fontsize=12, fontface='bold'))),
                         t=5, l=8, b=7, r=8, name = paste(runif(2)))
ggdraw(z.uni.s1)
## Do not run this unless you want to change figures (and have done at least 999 iterations
## commented out to prevent mistakes
## pdf(file='../huntingpaper/figures/FigureS1_uniwitouthpl.pdf', width=4.5, height=4)
## ggdraw(z.uni.s1)
## dev.off()

## png(file='../huntingpaper/figures/FigureS1_uniwithoutpl.png', width=4.5, height=4, units='in', res=600)
## ggdraw(z.uni.s1)
## dev.off()
##

## End of file
