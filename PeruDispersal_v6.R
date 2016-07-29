################################################################################
### File: PeruDispersal.R
### Version: 6
### Author: Robert Bagchi
### Email: robert.bagchi@uconn.edu
### Date last modified: 22 July 2016
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
## if you need to install ReplicatedPointPatterns 
## (note you need devtools installed)
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
trees <- read.csv("../../data/mdd6plots_data_22jul2016.csv") ## tree locations

dispersal <-  read.csv("../../data/MDDspp_dispersal_syndromes_22jul2016.csv") ## dispersal data
sitedat <- read.csv("../../data/sitedata.csv") ## site information

## Define maximum range over which to consider spatial patterns
rmax <- 15
## put dispersal into 1/0 format
dispersal[, c('LV', 'SV', 'SB', 'Bat', 'TR', 'Wind', 'Expl', 'Unkwn', 'HSD', 'HID')] <-
  apply(dispersal[, c('LV', 'SV', 'SB', 'Bat', 'TR', 'Wind', 'Expl', 'Unkwn', 'HSD', 'HID')],
        2, function(x) ifelse(is.na(x), 0, x))

## combine wind and explosive dispersal into a single variable
dispersal$Abiotic <-  dispersal$Wind + dispersal$Expl



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
nrow(trees)
trees.old <- trees ## make a copy for comparison
trees <-trees[-c(grep('Astrocaryum', trees$Spp),
                 grep('Attalea', trees$Spp),
                 grep('Iriartea', trees$Spp),
                 grep("Chamaedorea", trees$Spp)),] ## removing palms

nrow(trees.old) - nrow(trees) ## lose 3727 palms
table(trees$Type)

## We need to remove species that

## 1. unidentified (spcode = -999)
## 2. have no dispersal data
## 3. Juveniles > 5 cm dbh
## 4. bivariate: have no adults or saps or juvs in the central ha of a site
##    univariate: have < 2 juvs or saplings in the central ha. of a site

## remove unidentified stems
sum(trees$code==-999) ## removes 7161 unidentified stems
trees <-  trees[trees$code != -999,] 
nrow(trees) ## leaves us with 56215 stems

## now remove juveniles between 5-10 cm
## there are a few stems with missing dbh data - remove them too
trees <- trees[!(trees$diam > 5 & trees$Type=='J'),]
trees <- trees[!is.na(trees$Type),]
summary(trees)
nrow(trees) ## lose 3008 stems. Seems relatively minor.

## Now get the abundance of each species x type x site combination
## and merge it into the dataset
abund <- aggregate(Tag ~ Spp + code + Type + Site, data=trees, length)
summary(abund)
names(abund)[5] <- 'abund'

library(tidyr)
abund <- spread(abund, Type, abund, fill=0)
names(abund)[4:6] <- paste('N', names(abund[4:6]), sep='.')
abund$hasdisp <- abund$code %in% dispersal$code
## now working out if there is enough replication at each 
## species and site combination for 

abund$Quni <- abund$N.S > 1 & abund$N.J > 1 & abund$hasdisp
abund$Qbi <- abund$N.A > 0 & abund$N.J > 0 & abund$N.S > 0 & abund$hasdisp


## Some preliminary numbers of how many stems and species qualify
sum(abund$Quni); sum(abund$Qbi) ## 601 & 624 species
apply(abund[abund$Quni, c('N.J', 'N.S')], 2, sum) ## 17470 juvs and 9709 saps 
apply(abund[abund$Qbi, c('N.A', 'N.J', 'N.S')], 2, sum) ## 9667 ads 18570 juvs and 8854 saps 


trees <- merge(trees, abund, by=c('Spp', 'code', 'Site'), all =T)
trees$Spp <- droplevels(trees$Spp)
nlevels(trees$Spp)
summary(trees)

### Analysis of number of abundance as a function of hunting pressure and 
## dispersal syndrome
# trees.all <- trees[trees$code %in% dispersal$code,]
# 
# trees.all <- cbind(trees.all, dispersal[match(trees.all$code, dispersal$code),
#                                     c('Spp', 'LV', 'SV', 'SB', 'Bat',
#                                       'TR', 'Wind', 'Expl',
#                                       'Unkwn', 'Abiotic', 'HSD', 'HID')])
# table(trees.all$Type)
# length(unique(trees.all$code))
# trees.all <- cbind(trees.all,
#   sitedat[match(trees.all$Site, sitedat$site),
#           c('forest', 'pname', 'hunted', 'lg.primates', 'lg.birds', 'huntpres')])
# 
# trees.all$Site <- factor(trees.all$Site)
# ## For the analysis of counts we need to aggregate the data so we have 
# ## number of stems within each species and their 
# library(tidyr)
# library(dplyr)
# trees.all <- trees.all[, c('Site', 'pname', 'hunted', 
#                            'lg.primates', 'lg.birds', 'huntpres',
#                            'Spp', 'code', 'Type', 'Tag', 'LV',
#                            'SV', 'SB', 'Bat', 'TR', 'Wind', 'Expl', 'Unkwn', 'Abiotic',
#                            'HSD')]
# trees.all$count <- rep(1, nrow(trees.all))
# trees.all.ag <-  aggregate(count ~ Site + pname + hunted + lg.primates + lg.birds +
#                             huntpres + Spp + code + Type + 
#                             LV + SV + SB + Bat + TR + Wind + Expl + 
#                             Unkwn + Abiotic + HSD, data=trees.all, sum)
# library(lme4)
# library(afex)
# library(car)
# summary(trees.all.ag)
# trees.all.ag <-  subset(trees.all.ag, Abiotic !=1  & Unkwn !=1)
# 
# trees.all.ag <- droplevels(trees.all.ag)
# 
# test.modS <-  glmer.nb(count ~ hunted*HSD + (1|Site) + (1|Spp), 
#                    data=trees.all.ag, subset=Type=='S')
# 
# 
# test.modJ <-  glmer.nb(count ~ hunted *HSD + (1|Site) + (1|Spp), 
#                        data=trees.all.ag, subset=Type=='J')
# 
# test.modA <-  glmer.nb(count ~ hunted * HSD + (1|Site) + (1|Spp), 
#                        data=trees.all.ag, subset=Type=='A')
# 
# summary(test.modS)
# summary(test.modJ)
# summary(test.modA)
# Anova(test.modA)
# 
# trees.all.ag$area <- ifelse(trees.all.ag$Type=='A', 4, 1)
# 
# test.mod <-  glmer.nb(count ~ Type*huntpres*HSD + (1|Site) + (1|Spp)  + offset(log(area)), 
#                         data=trees.all.ag)
# test.mod <- update(test.mod, start=pars)
# summary(test.mod)
# devfun <- update(test.mod, devFunOnly=TRUE)
# Anova(test.mod, type=3)
# summary(test.mod2)
# plot(test.mod)
# plot(test.mod, resid(., type='pearson') ~ fitted(.), type=c('p', 'smooth'))
# summary(test.mod)
# Anova(test.mod)
# names(trees.all.ag)
# trees.counts <- spread(trees.all.ag[, c('pname', 'hunted', 'huntpres', 'Spp', 'Type', 'count')], Spp, count, fill=0)
# nms <- aggregate(count~Spp, data=trees.all.ag, sum)
# nms <- nms[rev(order(nms$count)),]
# 
# library(vegan)
# ord <- metaMDS(trees.counts[, -c(1:4)])
# plot(ord)
# plot(ord, type='n')
# text(ord, display="sites", labels=paste(trees.counts$Type), 
#      col=ifelse(trees.counts$hunted=='hunted', 'red', 'black'))
# plot(ord, type='n')
# text(ord, display="sites", labels=trees.counts$pname, 
#      col=c(1, 4, 2)[as.numeric(trees.counts$Type)])
# ordiellipse(ord, trees.counts$hunted)
# 
# ord.fit <- envfit(ord ~ huntpres+Type, data=trees.counts)
# plot(ord.fit)
# ord.fit
# 
# plot(ord)
# text(ord, display='species', select=head(nms$Spp, 20), cex=0.5, col=)
# orditorp(ord, display='species')

## split the data by site
trees <- split(trees, f=trees$Site)

## Pull out some summary stats
sapply(trees, nrow)
sapply(trees, function(x) length(unique(x$Spp)))

## defining the window at CC12 which is not quite rectangular
cc12  <- owin(poly=list(x=c(225, 75, 75, 24.9, 25, 0, 0, 225),
                y= c(200, 200, 125, 125, 75, 75, 0, 0)))

plot(cc12, axes=T) ## looks good

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


## make a species list for each plot (based on the species that have sufficient
## stems in each size class 
## Univariate analysis
splists.uni <- lapply(trees, function(dat) unique(dat$code[dat$Quni]))

## bivariate analysis
splists.bi <- lapply(trees, function(dat) unique(dat$code[dat$Qbi]))

## turn the data from each species into a point pattern, split up by
## age class
## a few adult points in CC12 are right on the boundary.
## this leads to the warnings below, that can be
## ignored
ppp.sps.uni <- mapply(function(spid, treedat, win){
    ppp.sp <- lapply(as.list(spid), function(id, treedat, win){
        treedat <- droplevels(treedat[treedat$Type !='A',])
        win <- win[c('J', 'S')]
        condat <- treedat[treedat$code==id,] ## subset data from species
        ppp.con <- mapply(function(Tx, Wx){ ## turn into ppp object
            pppx <- ppp(x=Tx$E, y=Tx$N, window=Wx)
            return(pppx)}, Tx=split(condat, f=condat$Type), Wx=win, SIMPLIFY=FALSE)
        ## note that species are split by size-class (Type) when input to the function
        ## now heterospecifics
        hetdat <- treedat[treedat$code!=id,] ## subset out all heterspecifics
                                            ## Note that this will include species
                                            ## not included in elsewhere in the analysis
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
  }, spid=splists.uni, treedat=trees, win=winlist, SIMPLIFY=FALSE)

## Now for the bivariate case.
ppp.sps.bi <- mapply(function(spid, treedat, win){
  ppp.sp <- lapply(as.list(spid), function(id, treedat, win){
    condat <- treedat[treedat$code==id,] ## subset data from species
    ppp.con <- mapply(function(Tx, Wx){ ## turn into ppp object
      pppx <- ppp(x=Tx$E, y=Tx$N, window=Wx)
      return(pppx)}, Tx=split(condat, f=condat$Type), Wx=win, SIMPLIFY=FALSE)
    ## note that species are split by size-class (Type) when input to the function
    ## now heterospecifics
    hetdat <- treedat[treedat$code!=id,] ## subset out all heterspecifics
    ## Note that this will include species
    ## not included in elsewhere in the analysis
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
}, spid=splists.bi, treedat=trees, win=winlist, SIMPLIFY=FALSE)


## Remove species with adults only on the edge of the plot
ppp.sps.bi <-  
  lapply(ppp.sps.bi, function(site.dat){
    site.dat[sapply(site.dat, function(sp.dat){
    npoints(sp.dat$A[dilation(Window(sp.dat$J), 15)])}) > 0]})

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
                   r=0:rmax, correction="border", ratio=TRUE))
        hyper$K[hyper$comp=='con'] <-
            lapply(hyper$pppx[hyper$comp=='con'], function(px)
                Kmulti(px, I=inside.owin(px, w=smallwindow),
                       J=inside.owin(px, w=Window(px)), r=0:rmax, correction="border",
                       ratio=TRUE))

        hyper$wts <-  mapply(function(pppx, pppy)
            kfunc.weights.calc(pppx=pppx, pppy=pppy, r=0:rmax,
                               correction='border', type='sqrtnxny_A'),
                             pppx=hyper$pppx, pppy=hyper$pppy, SIMPLIFY=FALSE)
        return(hyper)
    }, dat = ppp.sps.uni, win=winlist, SIMPLIFY=FALSE)

## Add site names
hyperdat.uni <- mapply(function(dat, site) {
                       dat$site <- site
                       return(dat)
                   }, dat=hyperdat.uni, site=as.list(names(hyperdat.uni)),
                   SIMPLIFY=FALSE)

## work out which species in which sites have insufficient data
## This might be because too many individuals are on the border and
## so are removed when border corrections are applied - need at least one
## individual within the border region
univarKmods.bysite <-  lapply(hyperdat.uni, function(dat){
    lmeHyperframe(dat, 0:rmax, fixed="stage*comp", random="1|Spp",
                  computeK=FALSE, weights.type='nxny_A', minsamp=1)})
##  The row names of removed species are returned by the function - 
## remove them from data set
removedsp <- sapply(univarKmods.bysite, function(x) attr(x, 'removed.species'))

## Find the removed species names
rmsp.uni <- mapply(function(dat, sp){
  rmsp <- droplevels(unique(dat[sp,]$Spp))}, 
  hyperdat.uni, removedsp)


## remove species with insufficient data
hyperdat.uni.sel <- mapply(function(dat, sp) {
    sp.rm <- unique(dat[sp, 'sp.id']$sp.id) ## pulls out the rows
    print(length(sp.rm))
    dat <- dat[!(dat$sp.id %in% sp.rm),]
    return(dat)},  hyperdat.uni, removedsp, SIMPLIFY=FALSE)
## gives us the number of speciesxsite combinations and for each site.
sum(sapply(hyperdat.uni.sel, function(x) length(unique(x$Spp))))
sapply(hyperdat.uni.sel, function(x) length(unique(x$Spp)))

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
hyperdat.bi  <- lapply(ppp.sps.bi, function(dat)
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
                   r=0:rmax, corr='border', ratio=TRUE))

        hyper$wts <-  mapply(function(pppx, pppy)
            kfunc.weights.calc(pppx=pppx, pppy=pppy, r=0:rmax,
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

## data set with all data
hyperdat.bi.sel <- do.call('rbind', hyperdat.bi) ## note that we don't have 
                                                     # to worry about edge effects.

hyperdat.bi.sel <- cbind.hyperframe(
    hyperdat.bi.sel,
    sitedat[match(hyperdat.bi.sel$site, sitedat$site),
                   c('forest', 'pname', 'hunted', 'lg.primates', 'lg.birds', 'huntpres')])

hyperdat.bi.sel$site <- factor(hyperdat.bi.sel$site)

## One species, Unonopsis floribunda seems to cause problems(sp.id 729)
## at CC2. Removing that manually
hyperdat.bi.sel <- subset(hyperdat.bi.sel, !(sp.id=="729" & pname=='CC2'))

## Save the objects required for analysis so we don't
## always have to repeat processing (and for transfer to cluster)
save('dispersal', "sitedat", "extractModRanefs",
          'hyperdat.uni.sel', 'hyperdat.bi.sel', 'K2L', 'kfunclme2plot',
          'plot.kfunctionlme', 'ppp.sps.uni', "ppp.sps.bi", 'theme_cust', 'trees', 'winlist',
          'sitedat', 'standardise', file='data4peruanalysisv6.RData')

################################################################################
## Fitting models to the bivariate and univariate data
################################################################################
##rm(list=ls())
##setwd("D:/rbagchi/Documents/Dropbox/Varun/finalAnalysis/PeruAnalysis")
##load(file='../data/data4peruanalysisv6.RData')
## Can start from here and jump to models to save time

################################################################################
### Subsetting data
################################################################################

# ## abiotic species may be a bit odd and poorly represented. could remove from analysis
## but not doing so for now.
# hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Abiotic!=1 & Unkwn !=1) 
## Removing species with unknown dispersal
hyperdat.bi.sel.c <- subset(hyperdat.bi.sel, Unkwn !=1)

## Now univariate K funcs
##hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Abiotic!=1  & Unkwn !=1)
hyperdat.uni.sel.c <- subset(hyperdat.uni.sel, Unkwn !=1)

## Extracing some summaries of the data to report in the paper
## Some numbers on number of individuals/species etc etc

## mean reliance on different dispersers
disp.list <-  aggregate(cbind(LV, SV, SB, Bat, TR) ~Spp, subset(hyperdat.uni.sel.c,
                                                                stage=='S' & comp=='con'), mean)


## Need to revise this code because we have to take into account both univariate and bivariate
## analyses separately now. 
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
 # hyperdat.bi.sel.c <- subset(hyperdat.bi.sel.c,
 #                             !(Spp == 'Pseudolmedia laevis'))
 # 
 # hyperdat.uni.sel.c <- subset(hyperdat.uni.sel.c,
 #                             !(Spp == 'Pseudolmedia laevis'))
 # 

## Fitting the models

## set number of iterations for bootstrapping - keep small for trials (takes a long time)
## Code was run on a grid computer for paper
nsim <- 19
rmax <- 15
## First sapling only models
hyperdat.bi.sap.c <- subset(hyperdat.bi.sel.c, stage=='S')
hyperdat.uni.sap.c <- subset(hyperdat.uni.sel.c, stage=='S')

## Set intact forests as the reference level when checking the
## effect of using a categorical hunting treatment
hyperdat.uni.sap.c$hunted <- relevel(hyperdat.uni.sap.c$hunted, 'intact')
hyperdat.bi.sap.c$hunted <- relevel(hyperdat.bi.sap.c$hunted, 'intact')

sapMod.bi <- lmeHyperframe(hyperdat.bi.sap.c, 0:rmax,
                         fixed="comp*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

# sapMod.bi <- lmeHyperframe(hyperdat.bi.sap.c, 0:rmax,
#                            fixed="comp*hunted*HSD",
#                            random="1|site/Spp",
#                            computeK=FALSE)



sapMod.uni <- lmeHyperframe(hyperdat.uni.sap.c, 0:rmax,
                         fixed="comp*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

# sapMod.uni <- lmeHyperframe(hyperdat.uni.sap.c, 0:25,
#                             fixed="comp*hunted*HSD",
#                             random="1|site/Spp",
#                             computeK=FALSE)


## ## set up model matrix for predictions
preddat.sap <-  expand.grid(stage=c('S'), comp=c('con', 'het'),
                            huntpres=quantile(sitedat$huntpres, c(0.25, 0.75)),
                            HSD = c(0, 1))

# preddat.sap <-  expand.grid(stage=c('S'), comp=c('con', 'het'),
#                             hunted=c('intact', 'hunted'),
#                             HSD = c(0, 1))
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
## ## COmbine data from within and between cohort analyses

# sap.plotdat <- merge(sap.uni.plotdat, sap.bi.plotdat,
#                      by=c('distance', 'stage', 'comp', 'hunted', 'HSD'), all=T)


## ## Make labels more reader friendly
sap.plotdat$huntpres <-  factor(sap.plotdat$huntpres, labels=c('Low', 'High'))
##sap.plotdat$huntpres <-  factor(sap.plotdat$hunted, labels=c('Low', 'High'))
sap.plotdat$HSD <-  factor(sap.plotdat$HSD, labels=c('Insensitive', 'Sensitive'))

## ## Make plots
pl.sap.bi <- ggplot(sap.plotdat, aes(x=distance, y=K.bi, ymin=lcl.bi, ymax=ucl.bi)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD~  huntpres) +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
  ##  scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])))

pl.sap.uni <- ggplot(sap.plotdat, aes(x=distance, y=K.uni, ymin=lcl.uni, ymax=ucl.uni)) +
    geom_ribbon(alpha=0.5, colour=NA) +  geom_line() +
    facet_grid(HSD ~ huntpres) +
    geom_abline(intercept=0, slope=0, linetype='dashed', col='grey') +
    theme_bw() + theme_cust +
##   scale_color_brewer(palette='Set1')+ scale_fill_brewer(palette='Pastel1') +
    labs(x='Distance, r, (m)', y=expression(bold(hat(K)(r)[con]-hat(K)(r)[het])))

pl.saps <- plot_grid(pl.sap.uni, pl.sap.bi, labels= c('Uni', 'bi'))

##ggsave("saplingplot.pdf", pl.saps)
## Anova type analysis for table 2 - run on cluster. Repeated here with
## very few iterations to demonstrate code. 999 used on cluster.

## main effects
sapMod.uni1way <- lmeHyperframe(hyperdat.uni.sap.c, 0:rmax,
                         fixed="comp+huntpres+HSD",
                         random="1|site/Spp",
                         computeK=FALSE)


anova1waysaps <-  bootstrap.compare.lme(sapMod.uni1way, term='comp',
                          dists=1:rmax, nboot=nsim, ncore=7)

## 2 way interactions
sapMod.uni2way <- lmeHyperframe(hyperdat.uni.sap.c, 0:rmax,
                         fixed="(comp+huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE)

anova2waysaps <- sapply(c('comp:huntpres', 'comp:HSD'), function(term)
{
    bootstrap.compare.lme(sapMod.uni2way, term=term,
                          dists=1:rmax, nboot=nsim, ncore=7)
}, simplify=FALSE)



## 3 way interaction
anova3waysaps <- bootstrap.compare.lme(sapMod.uni, term="comp:huntpres:HSD", dists=1:rmax, nboot=nsim, ncore=7)

anova1waysaps[1:2]
sapply(anova2waysaps, function(x) x[1:2])
anova3waysaps[1:2]

## load the results of the anovas - run on computing grid
## This will only run after the other analyses have been run and results saved.
## This analysis includes all except for wind dispersed species.

## function to extract the data we need from these results
extractAnovas <- function(Aobj, type){
    Aobj <- Aobj[[type]]$anovas
    sapply(Aobj, function(x)  x[['dtable']][c("term", "D", "p")], simplify=TRUE)
}


### saplings only
## with P laevis
## system("scp bbcsrv3:~/Peru/July2016/results/consolidatedResults/Peru_v6_?_wayAnovaSaps_nopl0_noqw0_rmax15.RData ../results/")
## extract the relevant data
saps.Anova.withpl <- sapply(c('uni', 'bi'), function(pairing){
    saps.AnovaUni_withpl <- t(Reduce('cbind', sapply(1:3, function(i){
        cat(i)

        objname <- load(paste0("../results/Peru_v6_", i,
                               "_wayAnovaSaps_nopl0_noqw0_rmax10.RData"))
        Aobj <- get(objname)
        print(Aobj[c('intlevel', 'type')]) ## check it's the correct data
        extractAnovas(Aobj, pairing)
    })))}, simplify=FALSE)
(Aobj$uni$anovas[[1]]$dtable$dists)
summary(Aobj$uni$model)
saps.Anova.withpl

## without P laevis
## system("scp bbcsrv3:~/Peru/July2016/results/consolidatedResults/Peru_v6_?_wayAnovaSaps_nopl1_noqw0_rmax10.RData ../results/")

## check effect of excluding P laevis
saps.Anova.nopl <- sapply(c('uni', 'bi'), function(pairing){
    saps.AnovaUni_withoutpl <- t(Reduce('cbind', sapply(1:3, function(i){
    cat(i)
    objname <- load(paste0("../results/Peru_v6_", i, "_wayAnovaSaps_nopl1_noqw0_rmax10.RData"))
    Aobj <- get(objname)
    print(Aobj[c('intlevel', 'type')]) ## check it's the correct data
    extractAnovas(Aobj, pairing)
})))}, simplify=FALSE)
saps.Anova.nopl


## Full models with Juveniles and saplings (pair2)
intMod.bi <- lmeHyperframe(hyperdat.bi.sel.c, 0:rmax,
                         fixed="comp*stage*huntpres*HSD",
                         random="1|site/Spp",
                         computeK=FALSE)

intMod.uni <- lmeHyperframe(hyperdat.uni.sel.c, 0:rmax,
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
system.time(intMod.bi.boot <-  bootstrap.t.CI.lme(intMod.bi, lin.comb.Ct=modmat.int, nboot=nsim, alpha=0.05,
                                    ncore=5))

system.time(intMod.uni.boot <-  bootstrap.t.CI.lme(intMod.uni, lin.comb.Ct=modmat.int, nboot=nsim,
                                       alpha=0.05, ncore=8))

system("scp bbcsrv3:~/Peru/July2016/results/results_2016-07-24/PeruDispersal_nopl0.RData ../results/")## can load results from cluster


objname <- load('../results/PeruDispersal_nopl0.RData')
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
juvMod.uni1way <- lmeHyperframe(hyperdat.uni.sel.c, 0:rmax,
                         fixed="(comp+stage + huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE
)
anova2wayjuvs <-  bootstrap.compare.lme(juvMod.uni1way, term='comp',
                          dists=1:25, nboot=nsim, ncore=7)

## 2 way interactions
sapMod.uni2way <- lmeHyperframe(hyperdat.uni.sap.c, 0:rmax,
                         fixed="(comp+huntpres+HSD)^2",
                         random="1|site/Spp",
                         computeK=FALSE)

anova2waysaps <- sapply(c('comp:huntpres', 'comp:HSD'), function(term)
{
    bootstrap.compare.lme(sapMod.uni2way, term=term,
                          dists=1:rmax, nboot=nsim, ncore=7)
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

hyperdat.uni.ind <- subset(hyperdat.uni.sel, Unkwn < 1)

## K functions based on < 15 individuals are meaningless on their own
## (although they can contribute with low weight to an overall analysis)
##so removing these.
hyperdat.uni.ind <- subset(hyperdat.uni.sel, N1 > 15)
hyperdat.uni.ind$Spp <- droplevels(hyperdat.uni.ind$Spp)
counts <- table(hyperdat.uni.ind$Spp, hyperdat.uni.ind$site)
counts <- as.matrix(counts)
## only include species that are represented in at least 1 hunted and one intact site (note
# that if a site has a species there will be 4 rows (saps + juvs x con  + het)
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


## Top 4 species
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
library(ggrepel)
ind.funcs$site.lab <- ind.funcs$pname
ind.funcs$site.lab[ind.funcs$r != 23] <- NA

sp.plot  <- ggplot(ind.funcs, aes(x=r, y=Kdiff, colour=hunted, group=stage:site), fill='white')+
    geom_line() + 
  geom_text_repel(aes(label=site.lab), size=2, show.legend=FALSE) + ##,) segment.color=0) +
  facet_grid(stage~ Spp, scales='free') +
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
pdf('../huntingpaper/figures/Figure2_bysp.pdf', width=8, height=4)
ggdraw(z.sp)
dev.off()
png('../huntingpaper/figures/Figure2_bysp.png', width=8, height=4, units='in', res=600)
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


ind.funcs.all$stage <- factor(ind.funcs.all$stage, labels=c('Saplings', 'Juveniles'))

##ind.funcs$Spp <-  factor(ind.funcs$Spp, levels=ntotals$Spp)
## order by disperser sensitivity
ind.funcs.all$Spp <-  factor(ind.funcs.all$Spp, levels=disptype$Spp)

ind.funcs.all$N <-  ntotals$N1[match(ind.funcs.all$Spp, ntotals$Spp)]
ind.funcs.all$HSD <- disptype$HSD[match(ind.funcs.all$Spp, disptype$Spp)]
ind.funcs.all$splab <- paste(ind.funcs.all$Spp, '\n(n = ', ind.funcs.all$N, ')', sep='')
hsdbysp <- aggregate(HSD ~ splab, data=ind.funcs.all, unique) 
hsdbysp <- hsdbysp[order(hsdbysp$HSD),]
ind.funcs.all$splab <- factor(ind.funcs.all$splab, levels=hsdbysp$splab) 
ind.funcs.all$hunted <- factor(ind.funcs.all$hunted, levels=c('intact', 'hunted'),
                           labels=c('Low', 'High'))
ind.funcs.all$site.lab <- ind.funcs.all$pname
ind.funcs.all$site.lab[ind.funcs.all$r !=23] <- NA

plotgrps <- list(hsdbysp$splab[1:7], hsdbysp$splab[8:13], hsdbysp$splab[14:19]) 
sp.plot.all  <- lapply(plotgrps, function(grp) 
  {
  ggplot(subset(ind.funcs.all, splab %in% grp), aes(x=r, y=Kdiff, colour=hunted, group=stage:site),
         fill='white')+
    geom_line() + geom_text_repel(aes(label=site.lab), size=2, show.legend=FALSE) +
    facet_grid( splab + HSD  ~stage, scales='free') +
    
    scale_colour_brewer(palette='Dark2') +
    geom_abline(intercept=0, lty='dashed', col='grey50') + theme_bw() + theme_cust +
    labs(x='Distance, r, (m)', y=expression(hat(K)(r)[con]-hat(K)(r)[het]),
         color='Hunting\npressure') +
    theme(strip.text.y=element_text(size=7, face='italic'))
})  


z.sp.all <- lapply(sp.plot.all, function(pl) ggplot_gtable(ggplot_build(pl)))
##gtable_show_layout(z.sp.all[[1]])
## add a labels at the top and sides

z.sp.all <- lapply(z.sp.all, function(z)
{
  z <- gtable_add_rows(z, height=unit(z$heights[[3]], 'cm'), pos=2)
  z <- gtable_add_cols(z, width=unit(z$widths[[8]], 'cm'), pos=8)

  ## add header row with stage label
  z <- gtable_add_grob(z,list(rectGrob(gp=gpar(fill='grey', col='grey', lwd=0.2)),
                       textGrob("Stage",
                                gp = gpar(fontsize=10))),
                  t=3, l=4, b=3, r=6, name = paste(runif(2)))
  ## add column with disperser senstitivity label
  barlngth <-  length(z$heights) -3

  z <- gtable_add_grob(z,list(rectGrob(gp = gpar(fill='grey', col='grey', lwd=0.2)),
                                   textGrob("Disperser sensitivity", rot = -90,
                                            gp = gpar(fontsize=10))),
                              t=5, l=9, b=barlngth, r=9, name = paste(runif(2)))
  
  z <- gtable_add_cols(z, width=z$widths[[2]], pos=1)
  
  z <- gtable_add_grob(z, list(rectGrob(gp = gpar(fill=NA, col=NA, lwd=0)),
                                   textGrob("Relative Conpecific clustering", rot = 90,
                                            gp = gpar(fontsize=10))),
                              t=5, l=2, b=barlngth, r=2, name = paste(runif(2)))
  return(z)
  })
ggdraw(z.sp.all[[2]])  
  
pdf('../huntingpaper/figures/FigS2_byspecies_3pg.pdf', paper='A4', width=7, height=10)
lapply(z.sp.all, function(z) ggdraw(z))
dev.off()


## Plot figure S1 (Univariate K functions without P. laevis)
##system("scp bbcsrv3:~/Peru/July2016/results/results_2016-07-24/PeruDispersal_nopl1.RData ../results/")

objname <- load('../results/PeruDispersal_nopl1.RData')
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
dev.new()
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
