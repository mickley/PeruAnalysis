## A function that plots the BLUEs against distance
plot.klmerHyper <- function(x, ...){
  require(ggplot2)
  dat <- do.call('rbind', lapply(x, tidy, effects='fixed'))
  dists <- as.numeric(names(x))
  dat$distance <- rep(dists, each=length(fixef(x[[1]])))
  
  ggplot(dat, aes(x=distance, y=estimate)) + geom_line() +
    facet_wrap(~term, scales="free_y")  + theme_bw()
}

## Function that standardises a vector
standardise <-  function(x) (x - mean(x))/sd(x)


## a function that plots the coefficient and confidence intervals
## from a bootstrap
plot.klmerci <- function(x){
  dat <- x$pars_fixed
  require(tidyr)
  dat <-  tidy(t(as.data.frame(dat)))

  names(dat)[1] <- 'labels'
  
  dat <- separate(dat, labels, into=c("distance", "term"),sep="[.]")
  
  dat$distance <- as.numeric(as.character(dat$distance))
  dat$term <- factor(dat$term)
  dat$term <- factor(dat$term, 
                     levels=levels(dat$term)[order(
                       sapply(levels(dat$term), function(x) 
                         length(strsplit(x, split=":", fixed=T)[[1]])))])
  
  ggplot(dat, aes(x=distance, y=estimate, ymin=X2.5., ymax=X97.5.)) +
    geom_ribbon(colour=NA, alpha=0.3) + geom_line() + 
    geom_hline(yintercept=0, linetype='dotted') + facet_wrap(~term, scale='free')
}



## Function to extract random effect BLUPs to allow plotting
extractModRanefs <-  function(mod){
  mod_re <- lapply(mod, ranef)
  
  ranefs <- do.call('rbind', mapply(function(re, d)
  {
    do.call('rbind', mapply(function(re, level, d)
    {
      if(is.null(re))
        return(NULL)
      else
      {
        names(re) <- sapply(names(re), function(s) gsub("[[:punct:]]","", s))
        dat <- data.frame(level=level, group=rownames(re), distance=d, re)
        return(dat)
      }
    }, re=re, level=as.list(names(re)), MoreArgs=list(d=d), SIMPLIFY=FALSE))
  }, re=mod_re, d=as.list(names(mod_re)), SIMPLIFY=FALSE))
}

extractModResids <-  function(mod)
  {
  resids<- do.call('cbind', lapply(mod, function(mod) {
    if(is.null(mod))
      return(NULL)
    else
      return(resid(mod))
      return(NULL)
  }))
  distances <- as.numeric(names(mod)[!sapply(mod, is.null)])
  colnames(resids) <- distances
  coefdat <- mod[!sapply(mod, is.null)][[1]]@frame
  coefdat <- coefdat[, !(names(coefdat) %in% c("k", "(weights)"))]
  
  resids <- tbl_df(cbind(coefdat, resids))
  resids <-  gather(resids, key="distance", 
                    value="resid", as.character(distances))
  resids$distance <- as.numeric(resids$distance)
  return(resids)
}


## a function that takes the output from the model and puts it in
## a convenient format for plotting

klmerci2plot <- function(x, preddat)
{
  
  preds <- tbl_df(aperm(x$predictions, c(2, 1, 3)))
  preds <- cbind(preddat, preds)
  preds <- gather(preds, key="key", value="value",  
                  (ncol(preddat)+1):ncol(preds))  %>%
    separate(key, into=c("distance", "quant"), sep="[.]",
             convert=TRUE) %>% 
    spread(key = quant, value=value)
  return(preds)
}

## A function to calculate P-values from a klmerci object
pvalCalc <-  function(x, mods, dists=NULL)
{
  
  bootobj <- attr(x, "bootobj")
  
  if(is.null(dists)) 
    dists <- 1:length(bootobj[[1]])
  
  t_obs <- sapply(mods, function(m) lme4::fixef(m)/sqrt(Matrix::diag(vcov(m))))
  
  t_sim <- lapply(bootobj, function(sim)
  {
    sapply(sim, function(sim_d) sim_d$beta_r/sqrt(diag(sim_d$vcov_r)))
  })
  t_sim_mat <- do.call(abind::"abind", args=list(what = t_sim, along=3))
  
  t_sim_sd <- apply(t_sim_mat, c(1, 2), sd)
  
  T_obs <- apply(t_obs/t_sim_sd, 1, mean)
  
  T_sim <- sapply(t_sim, function(sim, s, d)
  {
    apply(sim[,d]/s[,d], 1,  mean)
  }, s=t_sim_sd, d=dists)
  
  pval <- apply(T_sim, 1, function(x) 1-2*abs((sum(x>0)+1)/(length(x) + 1) - 0.5))
  return(cbind(T=T_obs, p=pval))  
}

