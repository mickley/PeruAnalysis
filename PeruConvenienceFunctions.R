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
  
  dat <- do.call('rbind', lapply(x, tidy, effects='fixed'))
  dists <- as.numeric(names(x))
  dat$distance <- rep(dists, each=length(fixef(x[[1]])))

  ggplot(dat, aes(x=distance, y=estimate)) + geom_line() +
    facet_wrap(~term, scales="free_y")  + theme_bw()
}

## Function that standardises a vector
standardise <-  function(x) (x - mean(x))/sd(x)

## Function to extract random effect BLUPs to allow plotting
extractModRanefs <-  function(mod, level=1){
  require(reshape2)
  ranefs <- lapply(mod, function(mod) ranef(mod))
  if(class(ranefs[[1]])[2] =='list') 
     ranefs <- do.call('cbind', lapply(ranefs, function(re) re[[level]]))
  else
  ranefs <- do.call('cbind', ranefs)  
  ranefs$id <- rownames(ranefs)
  names(ranefs)[1:length(mod)] <-  paste('distance', 1:length(mod), sep='.')
  ranefs <-  melt(ranefs, id='id')
  ranefs$distance <- as.numeric(do.call('rbind',
                                        strsplit(as.character(ranefs$variable), split='.', fixed=T))[,2])
  return(ranefs)
}


plotModResids <-  function(mod){
  require(tidyr)
  
  resids<- do.call('cbind', lapply(mod, function(mod) resid(mod)))
  nms <- names(attr(mod[[1]]$terms, 'dataClasses'))
  nms <- nms[nms!='K']
  dat <- mod[[1]]$data[,nms]
  
  colnames(resids)[1:length(mod)] <-  paste('distance', 1:length(mod), sep='.')
  resids <- tbl_df(cbind(dat, labs =rownames(resids), resids))
  resids <- gather_(resids, "distance", "resid", paste0('distance.', 1:15))
  
  resids <- separate(resids, distance, c('d', 'distance'), '[.]') %>%
    separate(labs, c('site', 'species'), '/')
  resids$distance <- as.numeric(resids$distance)
  return(resids)
  
}

## a function that plots the coefficient and confidence intervals
## from a bootstrap
effectplot.Kfunctionlme <- function(bootobj){
  dat <- bootobj$modelpars
  require(broom)
  require(tidyr)
  dat <-  tidy(t(as.data.frame(dat)))
  names(dat)[1] <- 'labels'
  dat <- separate(dat, labels, into=c("distance", "term"),sep="[.]")
  
  dat$distance <- as.numeric(as.character(dat$distance))
  dat$term <- factor(dat$term)
  dat$term <- factor(dat$term, levels=levels(dat$term)[order(
    sapply(levels(dat$term), function(x) 
      length(strsplit(x, split=":", fixed=T)[[1]])))])
  
  ggplot(dat, aes(x=distance, y=estimate, ymin=X2.5., ymax=X97.5.)) +
    geom_ribbon(colour=NA, alpha=0.3) + geom_line() + 
    geom_hline(yintercept=0, linetype='dotted') + facet_wrap(~term, scale='free')
}
