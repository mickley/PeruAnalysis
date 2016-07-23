################################################################################
### File: PeruDispersal.R
### Version: 4
### Author: Robert Bagchi
### Email: robert.bagchi@uconn.edu
### Date last modified: 28 February 2016
### Description: Code to analyse spatial patterns of saplings and juvenile
### trees in six Peruvian forests plots subjected to varying levels of
### defaunation.
################################################################################

################################################################################
### Preliminaries
################################################################################
rm(list=ls())


library(spatstat)
library(ggplot2)
library(cowplot)
library(grid)
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
trees <- read.csv('../data/mdd6plots_data_4jul2013.csv') ## tree locations

dispersal <-  read.csv('../data/MDDspp_dispersal_syndromes_CURR.csv') ## dispersal data

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


## remove Juveniles with diameter > 5 cm
trees <- trees[!(trees$Type=='J' & trees$diam >5),]
table(trees$Type, trees$Site)
table(trees.old$Type, trees.old$Site) ## seems like we lose relatively few trees
                                        # by removing the ones >5
## now remove adults that are really too small to be adults. Follwoing
## Varun's rules that adults must have a dbh > median dbh or >30 cm
tree.dbhcut <- aggregate(diam~Spp, data=trees, subset=Type=='A',
                         function(x) {
                             maxd <- max(x)
                             if(maxd < 25) return(10)
                             else return(min(30, median(x, na.rm=T)))
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
nlevels(trees$Spp) ## leaves us with 958 species

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
          'sitedat', 'standardise', file='data4peruanalysis3.RData')



################################################################################
## Fitting models to the bivariate and univariate data
################################################################################
##rm(list=ls())
##load(file='data4peruanalysis3.RData')
## Can start from here and jump to models to save time

################################################################################
### Subsetting data
################################################################################


## abiotic species may be a bit odd and poorly represented. Have to remove from analysis
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Abiotic==0 & Unkwn !=1) 

summary(hyperdat.bi.sel.c)

## Now univariate K funcs
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Abiotic==0  & Unkwn !=1)

nrow(hyperdat.uni.sel.c)/(2*2*2*2)

## Extracing some summaries of the data to report in the paper
## Some numbers on number of individuals/species etc etc

## mean reliance on different dispersers
disp.list <-  aggregate(cbind(LV, SV, SB, Bat, TR) ~Spp, subset(hyperdat.uni.sel.c,
                                                                stage=='S' & comp=='con'), mean)

apply(disp.list[, -1], 2, function(x) sum(x==1)) ## number totally dependent
apply(disp.list[, -1], 2, function(x) sum(x>0)) ## number dependent to some degree


aggregate(N1 ~ stage*comp, hyperdat.uni.sel.c, sum)
aggregate(Spp ~ site, hyperdat.uni.sel.c, function(x) length(unique(x)))
mean(aggregate(Spp ~ site, hyperdat.uni.sel.c, function(x) length(unique(x)))$Spp)
aggregate(Spp ~ 1, hyperdat.uni.sel.c, function(x) length(unique(x)))
aggregate(N2 ~ stage*comp, hyperdat.bi.sel.c,  sum)

## For analysis excluding P. laevis need to uncomment these lines
## hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c,
##                             !(Spp == 'Pseudolmedia laevis'))

## hyperdat.uni.sel.c <- subset(hyperdat.uni.sel.c,
##                             !(Spp == 'Pseudolmedia laevis'))


## Fitting the models

## set number of iterations for bootstrapping - keep small for trials (takes a long time)
## Code was run on a grid computer for paper
nsim <- 99

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
## set up model matrix for predictions
preddat.sap <-  expand.grid(stage=c('S'), comp=c('con', 'het'),
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))
sapform <- update(formula(sapMod.uni[[1]]), NULL~.)

modmat.sap <-  model.matrix(sapform, data=preddat.sap)
## Want to show differences between conspecifics and heterospecifics, so set up
## model matrix to do this
modmat.sap <-  modmat.sap[preddat.sap$comp=='con',]-modmat.sap[preddat.sap$comp=='het',]

## Bootstrapping
sapMod.bi.boot <-  bootstrap.t.CI.lme(sapMod.bi, lin.comb.Ct=modmat.sap, nboot=nsim, alpha=0.05,
                                    ncore=3)

sapMod.uni.boot <-  bootstrap.t.CI.lme(sapMod.uni, lin.comb.Ct=modmat.sap, nboot=nsim,
                                       alpha=0.05, ncore=3)

## Extract data for plot
sap.bi.plotdat <- kfunclme2plot(sapMod.bi.boot, preddat=subset(preddat.sap, comp=='con'))
sap.uni.plotdat <- kfunclme2plot(sapMod.uni.boot, preddat=subset(preddat.sap, comp=='con'))
names(sap.bi.plotdat)[2:4] <- paste(names(sap.bi.plotdat)[2:4], 'bi', sep='.')
names(sap.uni.plotdat)[2:4] <- paste(names(sap.uni.plotdat)[2:4], 'uni', sep='.')

## COmbine data from within and between cohort analyses
sap.plotdat <- merge(sap.uni.plotdat, sap.bi.plotdat,
                     by=c('distance', 'stage', 'comp', 'huntpres', 'HSD'), all=T)
## Make labels more reader friendly
sap.plotdat$huntpres <-  factor(sap.plotdat$huntpres, labels=c('Low', 'High'))
sap.plotdat$HSD <-  factor(sap.plotdat$HSD, labels=c('Insensitive', 'Sensitive'))

## Make plot
pl.sap.bi <- ggplot(sap.plotdat) +
    geom_smooth(aes(x=distance, y=K.bi, ymin=lcl.bi, ymax=ucl.bi),
                stat='identity', colour='black') +
    facet_grid(HSD~  huntpres) +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])))

pl.sap.uni <- ggplot(sap.plotdat) +
    geom_smooth(aes(x=distance, y=K.uni, ymin=lcl.uni, ymax=ucl.uni),
                stat='identity', col='black') +
    facet_grid(HSD ~ huntpres) +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
##   scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])))

plot_grid(pl.sap.uni, pl.sap.bi, labels= c('Uni', 'bi'))
dev.new(); pl.sqrtnowind
pl.sqrtnowind <- plot_grid(pl.sap.uni, pl.sap.bi, labels= c('Uni', 'bi'))
ggsave(file='Saps.sqrt.nowind.pdf', pl.sqrtnowind)

## Explore BLUPS - shows the outliers
ranef.sp <- extractModRanefs(sapMod.uni, level=2)
sp.cols <- do.call('rbind', strsplit(as.character(ranef.sp$id), split='/', fixed=TRUE))
colnames(sp.cols) <-  c('site', 'sp')
ranef.sp <- cbind(sp.cols, ranef.sp)
## Note the outliers
qplot(x=distance, y=value, data=ranef.sp, colour=sp, facets=~site, geom='line')

nsim <- 999
## Anova type analysis for table 2
anova3waysaps <- bootstrap.compare.lme(sapMod.uni, term="comp:huntpres:HSD", dists=1:25, nboot=nsim, ncore=3)
anova3waysaps[1:2]
names(anova3waysaps)
hist(anova3waysaps$D.boot); abline(v=anova3waysaps$D, col=2)
sapMod.uni2way <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
                         fixed="(comp+huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE)
fixef(sapMod.uni2way[[1]])
anova2waysaps <- sapply(c('comp:huntpres', 'comp:HSD'), function(term)
{
    bootstrap.compare.lme(sapMod.uni2way, term=term,
                          dists=1:25, nboot=nsim, ncore=3)
}, simplify=FALSE)

sapply(anova2waysaps, function(x) x[1:2])


sapMod.uni1way <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
                         fixed="comp+huntpres+HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

anova1waysaps <-  bootstrap.compare.lme(sapMod.uni1way, term='comp',
                          dists=1:25, nboot=nsim, ncore=3)

?RNGkind
?clusterSetRNGStream

## load the results of the anovas - run on computing grid
## This will only run after the other analyses have been run and results saved.
## all data (except for wind dispersed species)
## function to extract the data we need from these results
extractAnovas <- function(Aobj, type){
summary(Aobj[[type]]$anovas$d25)
sapply(Aobj[[type]]$anovas, function(x) x[['d25']])
       [c("term", "D", "p")])
}


### saplings only
##system('scp uchc:~/Peru/Peru_?_wayAnovaSaps_nopl0.RData results/')

load('results/Peru_1_wayAnovaSaps_nopl0.RData')
names(anovaSaps$uni$anovas)
saps.AnovaUni <- extractAnovas(anovaSaps, 'uni')
saps.AnovaBi <- extractAnovas(anovaSaps$anovas, 'bi')


load('results/Peru_2_wayAnovaSaps_nopl0.RData')
saps.AnovaUni <- cbind(saps.AnovaUni, extractAnovas(anovaSaps, 'uni'))
saps.AnovaBi <- cbind(saps.AnovaBi, extractAnovas(anovaSaps, 'bi'))
fixef(anovaSaps$uni[[3]])
hist(
    anovaSaps$uni[[1]]$d25$D.boot
)
    D.boot
load('results/Peru_3_wayAnovaSaps_nopl0.RData')
saps.AnovaUni <- cbind(saps.AnovaUni, extractAnovas(anovaSaps, 'uni'))
saps.AnovaBi <- cbind(saps.AnovaBi, extractAnovas(anovaSaps, 'bi'))

saps.AnovaUni
saps.AnovaBi

## ## after removing P laevis
load('Peru_1_wayAnovaSaps_nopl1.RData')

saps.AnovaUni_nopl <- extractAnovas(anovaSaps, 'uni')
saps.AnovaBi_nopl <- extractAnovas(anovaSaps, 'bi')


load('Peru_2_wayAnovaSaps_nopl1.RData')
saps.AnovaUni_nopl <- cbind(saps.AnovaUni_nopl,
                            extractAnovas(anovaSaps, 'uni'))
saps.AnovaBi_nopl <- cbind(saps.AnovaBi_nopl,
                           extractAnovas(anovaSaps, 'bi'))


load('Peru_3_wayAnovaSaps_nopl1.RData')
saps.AnovaUni_nopl <- cbind(saps.AnovaUni_nopl,
                            extractAnovas(anovaSaps, 'uni'))
saps.AnovaBi_nopl <- cbind(saps.AnovaBi_nopl,
                           extractAnovas(anovaSaps, 'bi'))

saps.AnovaUni_nopl
saps.AnovaBi_nopl




## ## Pull out some summaries
## sapply(anova.1way$uni1way, function(x)
##     sapply(x, function(d) d[1:2]), simplify=FALSE)
## sapply(anova.2way$uni2way, function(x)
##     sapply(x, function(d) d[1:2]), simplify=FALSE)

## sapply(anova.3way$uni3way[c('d10', 'd25')], function(x) x[c('D', 'p', 'term')])

## names(anova.3way)

## sapply(anova.1way$bi1way, function(x)
##     sapply(x, function(d) d[1:2]), simplify=FALSE)
## sapply(anova.2way$bi2way, function(x)
##     sapply(x, function(d) d[1:2]), simplify=FALSE)

## sapply(anova.3way$bi3way[c('d10', 'd25')], function(x) x[c('D', 'p', 'term')])


load('Peru_2wayAnovaJuvs_sqrtnxnynowind.RData')
load('Peru_3wayAnovaJuvs_sqrtnxnynowind.RData')
load('Peru_4wayAnovaJuvs_sqrtnxnynowind.RData')

sapply(anova.2wayJuvs$uni2way[c('d10', 'd25')], function(x) x[c('D', 'p', 'term')])

## sapply(anova.3wayJuvs$uni3way, function(x)
##     sapply(x[c('d10', 'd25')], function(d) d[c('D', 'p', 'term')]), simplify=FALSE)

## sapply(anova.4wayJuvs$uni4way, function(x) x[c('D', 'p', 'term')])


## sapply(anova.2wayJuvs$bi2way[c('d10', 'd25')], function(x) x[c('D', 'p', 'term')])

## sapply(anova.3wayJuvs$bi3way, function(x)
##     sapply(x[c('d10', 'd25')], function(d) d[c('D', 'p', 'term')]), simplify=FALSE)

## sapply(anova.4wayJuvs$bi4way, function(x) x[c('D', 'p', 'term')])

## table(hyperdat.uni.sel.c$Spp)
## Full models

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
                                    ncore=3)

intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni, lin.comb.Ct=modmat.int, nboot=nsim,
                                       alpha=0.05, ncore=3)

## can load results from cluster
rm(intMod.uni.boot, intMod.bi.boot)
load('dispersalBootstrap.RData')


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

## Make plots
pl.int.bi <- ggplot(int.plotdat) +
    geom_smooth(aes(x=distance, y=K.bi, ymin=lcl.bi, ymax=ucl.bi, colour=stage, fill=stage),
                stat='identity') +
    facet_grid(HSD~  huntpres, scales="free_y") +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])),
        color='Cohort', fill='Cohort')
pl.int.bi
## Note that on windows machines the legend gets a little messed up.

pl.int.uni <- ggplot(int.plotdat) +
    geom_smooth(aes(x=distance, y=K.uni, ymin=lcl.uni, ymax=ucl.uni, colour=stage, fill=stage),
                stat='identity') +
    facet_grid(HSD~  huntpres, scales="free_y") +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
    scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])),
         color='Cohort', fill='Cohort')
pl.int.uni


## Format figures
library(gtable)
z.bi <- ggplot_gtable(ggplot_build(pl.int.bi))
##gtable_show_layout(z.bi)
## add a label at the top and side
z.bi <- gtable_add_rows(z.bi, height=z.bi$heights[[3]], pos=2)
z.bi <- gtable_add_cols(z.bi, width=z.bi$widths[[7]], pos=7)
z.bi <- gtable_add_grob(z.bi,
                     list(rectGrob(gp=gpar(fill='grey', col=NA, lwd=0.2)),
                          textGrob("Hunting pressure",
                                   gp = gpar(fontsize=12, fontface='bold'))),
                     t=3, l=4, b=3, r=6, name = paste(runif(2)))
z.bi <- gtable_add_grob(z.bi,
                     list(rectGrob(gp = gpar(fill='grey', col=NA, lwd=0.2)),
                          textGrob("Disperser sensitivity", rot = -90,
                                   gp = gpar(fontsize=12, fontface='bold'))),
                     t=5, l=8, b=7, r=8, name = paste(runif(2)))
ggdraw(z.bi)
## Once again, some problems doing this on Windows - works fine on Ubuntu

## Univariate plot
z.uni <- ggplot_gtable(ggplot_build(pl.int.uni+theme_cust))
##gtable_show_layout(z.uni)
## add a label at the top and side
## extra row
z.uni <- gtable_add_rows(z.uni, height=z.uni$heights[[3]], pos=2)
## extra column
z.uni <- gtable_add_cols(z.uni, width=z.uni$widths[[7]], pos=7)
## add hunting pressure label
z.uni <- gtable_add_grob(z.uni,
                         list(rectGrob(gp=gpar(fill='grey', col=NA, lwd=0.2)),
                              textGrob("Hunting pressure",
                                       gp = gpar(fontsize=12, fontface='bold'))),
                         t=3, l=4, b=3, r=6, name = paste(runif(2)))
## add disperser sensitivity label
z.uni <- gtable_add_grob(z.uni,
                         list(rectGrob(gp = gpar(fill='grey', col=NA, lwd=0.2)),
                              textGrob("Disperser sensitivity", rot = -90,
                                       gp = gpar(fontsize=12, fontface='bold'))),
                         t=5, l=8, b=7, r=8, name = paste(runif(2)))

dir('../huntingpaper/figures')
## Do not run this unless you want to change figures (and have done at least 999 iterations
pdf(file='../huntingpaper/figures/Figure1_uniwithpl.pdf', width=4.5, height=4)
ggdraw(z.uni)
dev.off()

png(file='../huntingpaper/figures/Figure1_uniwithpl.png', width=4.5, height=4, units='in', res=600)
ggdraw(z.uni)
dev.off()

pdf(file='../huntingpaper/figures/Figure3_biwithpl.pdf', width=4.5, height=4)
ggdraw(z.bi)
dev.off()

png(file='../huntingpaper/figures/Figure3_biwithpl.png', width=4.5, height=4, units='in', res=600)
ggdraw(z.bi)
dev.off()

## load('Peru_1wayAnova_sqrtnxnynowindnoPl.RData')
## load('Peru_2wanyAnova_sqrtnxnynowindnopl.RData')
## load('Peru_3wayAnova_sqrtnxnynowindnopl.RData')

### Anova tye analysis including both stages
## with P laevis

load("results/Peru_1_wayAnovaJuvs_nopl0.RData")

juvs.AnovaUni <-  extractAnovas(anovaJuvs, 'uni')
juvs.AnovaBi <-  extractAnovas(anovaJuvs, 'bi')

load("results/Peru_2_wayAnovaJuvs_nopl0.RData")
juvs.AnovaUni <- cbind(juvs.AnovaUni, extractAnovas(anovaJuvs, 'uni'))
juvs.AnovaBi <- cbind(juvs.AnovaBi, extractAnovas(anovaJuvs, 'bi'))

load("results/Peru_3_nopl0_wayAnovaJuvs.RData")
juvs.AnovaUni <- cbind(juvs.AnovaUni, extractAnovas(anovaJuvs, 'uni'))
juvs.AnovaBi <- cbind(juvs.AnovaBi, extractAnovas(anovaJuvs, 'bi'))

load("results/Peru_4_wayAnovaJuvs_nopl0.RData")
juvs.AnovaUni <- cbind(juvs.AnovaUni, extractAnovas(anovaJuvs, 'uni'))
juvs.AnovaBi <- cbind(juvs.AnovaBi, extractAnovas(anovaJuvs, 'bi'))

juvs.AnovaUni
juvs.AnovaBi



## after removing P laevis
load("Peru_1_wayAnovaJuvs_nopl.RData")
juvs.AnovaUni <-  extractAnovas(anovaJuvs, 'uni')
juvs.AnovaBi <-  extractAnovas(anovaJuvs, 'bi')

load("Peru_2_wayAnovaJuvs_nopl.RData")
juvs.AnovaUni <- cbind(juvs.AnovaUni, extractAnovas(anovaJuvs, 'uni'))
juvs.AnovaBi <- cbind(juvs.AnovaBi, extractAnovas(anovaJuvs, 'bi'))

load("Peru_3_wayAnovaJuvs_nopl.RData")
juvs.AnovaUni <- cbind(juvs.AnovaUni, extractAnovas(anovaJuvs, 'uni'))
juvs.AnovaBi <- cbind(juvs.AnovaBi, extractAnovas(anovaJuvs, 'bi'))

load("Peru_4_wayAnovaJuvs_nopl.RData")
juvs.AnovaUni <- cbind(juvs.AnovaUni, extractAnovas(anovaJuvs, 'uni'))
juvs.AnovaBi <- cbind(juvs.AnovaBi, extractAnovas(anovaJuvs, 'bi'))

juvs.AnovaUni
juvs.AnovaBi



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

## Chose a few to plot
disptype <- merge(disptype, ntotals)
disptype <- disptype[order(disptype$N1, decreasing=T),]
disptype <- disptype[order(disptype$HSD>0.5),]

##subset(hyperdat.uni.ind, Spp == 'Quararibea wittii')
## Function to pull out conspecific - heterospecific difference
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
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])),
         color='Hunting\npressure')
sp.plot <- sp.plot + theme(strip.text.x=element_text(size=7, face='italic'))
sp.plot


library(gtable)

z.sp <- ggplot_gtable(ggplot_build(sp.plot))
##gtable_show_layout(z.sp)
## add a label at the top and side
z.sp <- gtable_add_rows(z.sp, height=z.sp$heights[[3]], pos=2)
z.sp <- gtable_add_cols(z.sp, width=z.sp$widths[[11]], pos=11)

##sapply(z.test$grobs[grep('strip.absoluteGrob', sapply(z.test$grobs, function(x) x$name))], function(x) x$children[[1]]$gp$fill)

z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp=gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Insensitive Dispersers",
                                   gp = gpar(fontsize=10, fontface='bold'))),
                     t=3, l=4, b=3, r=6, name = paste(runif(2)))
z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp=gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Sensitive Dispersers",
                                   gp = gpar(fontsize=10, fontface='bold'))),
                     t=3, l=8, b=3, r=10, name = paste(runif(2)))

z.sp <- gtable_add_grob(z.sp,
                     list(rectGrob(gp = gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Stage", rot = -90,
                                   gp = gpar(fontsize=10, fontface='bold'))),
                     t=5, l=12, b=7, r=12, name = paste(runif(2)))
ggdraw(z.sp)
pdf('../huntingpaper/figures/Fig2_bysp.pdf', width=8, height=4)
ggdraw(z.sp)
dev.off()
png('../huntingpaper/figures/Fig2_bysp.png', width=8, height=4, units='in', res=600)
ggdraw(z.sp)
dev.off()

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

summary(ind.funcs.all)
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
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])),
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
                                   gp = gpar(fontsize=10, fontface='bold'))),
                     t=3, l=4, b=3, r=6, name = paste(runif(2)))

z.sp.all <- gtable_add_grob(z.sp.all,
                     list(rectGrob(gp = gpar(fill='grey', col='grey', lwd=0.2)),
                          textGrob("Disperser sensitivity", rot = -90,
                                   gp = gpar(fontsize=10, fontface='bold'))),
                     t=5, l=9, b=41, r=9, name = paste(runif(2)))
ggdraw(z.sp.all)

pdf('../huntingpaper/figures/FigS2_byspecies.pdf', width=10, height=40)
ggdraw(z.sp.all)
dev.off()


