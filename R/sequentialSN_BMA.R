## input data/initial values functions for Bayesian model averaging
#### 1. additive
stan.in_seq_add <- function(infile=tempdata, x="Chla", y="part_microcystin",
                            n.chains=nchains, priorYr_rv=prior.rv, priorSn_rv=priorSn.rv[[3]]){
  cens <- infile$part_microcystin<=0.1
  ncens <- sum(cens)
  logx_bar <- mean(log(infile[,x]), na.rm=T)
  n <- dim(infile)[1]-ncens
  xrange <- range(log(infile[,x]))
  theta <- 0.01*diff(xrange)
  if (ncens>0 & n>0){
    xtemp <- infile[,x]
    ytemp <- infile[,y]
    ycens <- log(ytemp[cens])
    xcens <- log(xtemp[cens])
    y <- log(ytemp[!cens])
    x <- log(xtemp[!cens])
  } else {
    xtemp <- infile[,x]
    ytemp <- infile[,y]
    y <- log(ytemp)
    x <- log(xtemp)
  } 
  prior_yr0 <- priorN(priorYr_rv$B0)
  prior_yrD <- priorN(priorYr_rv$De)
  prior_yrP <- priorN(priorYr_rv$Ph)
  prior_snd0 <- priorN(priorSn_rv$d0)
  prior_sndp <- priorN(priorSn_rv$dp)
  prior_sndd <- priorN(priorSn_rv$dd)
  if (ncens>0 & n>0){
    bugs.data <- list(N=n, Ncens=ncens, y=array(y, dim=n),
                      ycens=array(ycens, dim=ncens),
                      x=array(x, dim=n), xcens=array(xcens,dim=ncens),
                      theta=theta, #beta1=0,
                      mu_dp = prior_sndp$m0,
                      sigma_dp=prior_sndp$sd,
                      mu_d0 = prior_snd0$m0,
                      sigma_d0=prior_snd0$sd,
                      mu_dd = prior_sndd$m0,
                      sigma_dd=prior_sndd$sd,
                      mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                      muD = prior_yrD$m0, sigmad = prior_yrD$sd,
                      muP = prior_yrP$m0, sigmap = prior_yrP$sd)
  } else {
    bugs.data <- list(N=ncens+n, y=y, x=x, 
                      theta=theta, #beta1=0,
                      mu_dp = prior_sndp$m0,
                      sigma_dp=prior_sndp$sd,
                      mu_d0 = prior_snd0$m0,
                      sigma_d0=prior_snd0$sd,
                      mu_dd = prior_sndd$m0,
                      sigma_dd=prior_sndd$sd,
                      mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                      muD = prior_yrD$m0, sigmad = prior_yrD$sd,
                      muP = prior_yrP$m0, sigmap = prior_yrP$sd)
  }
  inits <- list()
  for (i in 1:n.chains){
    inits[[i]] <- list(beta0=rnorm(1, prior_yr0$m0, prior_yr0$sd),
                       delta=abs(rnorm(1,prior_yrD$m0,prior_yrD$sd)),
                       phi=runif(1, range(x)[1], range(x)[2]),
                       d0=runif(1), dd=runif(1), dp=runif(1),
                       sigma=runif(1))
  }
  para <- c("beta0", "delta", "phi", "sigma", "d0", "dd", "dp","log_lik")
  return(list(para=para, data=bugs.data,
              inits=inits,n.chains=n.chains, theta=theta,
              cens=c(n, ncens)))
}

#### 2. nested
stan.in_seq_nst <- function(infile=tempdata, x="Chla", y="part_microcystin",
                            n.chains=nchains, priorYr_rv=priorYr.rvnstd,
                            priorSn_rv=priorSn.rvnstd){
  cens <- infile$part_microcystin<=0.1
  ncens <- sum(cens)
  logx_bar <- mean(log(infile[,x]), na.rm=T)
  n <- dim(infile)[1]-ncens
  xrange <- range(log(infile[,x]))
  theta <- 0.01*diff(xrange)
  if (ncens>0 & n>0){
    xtemp <- infile[,x]
    ytemp <- infile[,y]
    ycens <- log(ytemp[cens])
    xcens <- log(xtemp[cens])
    y <- log(ytemp[!cens])
    x <- log(xtemp[!cens])
  } else {
    xtemp <- infile[,x]
    ytemp <- infile[,y]
    y <- log(ytemp)
    x <- log(xtemp)
  } 
  prior_yr0 <- priorN(priorYr_rv$B0)
  prior_yrD <- priorN(priorYr_rv$De)
  prior_yrP <- priorN(priorYr_rv$Ph)
  prior_snd0 <- priorIG(priorSn_rv$sigmaS0)
  prior_sndd <- priorIG(priorSn_rv$sigmaSd)
  prior_sndp <- priorIG(priorSn_rv$sigmaSp)
  if (ncens>0 & n>0){
    bugs.data <- list(N=n, Ncens=ncens, y=array(y, dim=n),
                      ycens=array(ycens, dim=ncens), x=array(x, dim=n),
                      xcens=array(xcens, dim=ncens),
                      theta=theta, #beta1=0,
                      al0=prior_snd0$alp, bt0=prior_snd0$bet,
                      alP=prior_sndp$alp, btP=prior_sndp$bet, 
                      alD=prior_sndd$alp, btD=prior_sndd$bet,
                      mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                      muD = prior_yrD$m0, sigmad = prior_yrD$sd,
                      muP = prior_yrP$m0, sigmap = prior_yrP$sd)
  } else {
    bugs.data <- list(N=ncens+n, y=y, x=x, 
                      theta=theta, #beta1=0,
                      al0=prior_snd0$alp, bt0=prior_snd0$bet,
                      alP=prior_sndp$alp, btP=prior_sndp$bet,
                      alD=prior_sndd$alp, btD=prior_sndd$bet,
                      mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                      muD = prior_yrD$m0, sigmad = prior_yrD$sd,
                      muP = prior_yrP$m0, sigmap = prior_yrP$sd)
  }
  inits <- list()
  for (i in 1:n.chains){
    inits[[i]] <- list(beta0=rnorm(1,  prior_yr0$m0, prior_yr0$sd),
                       delta=abs(rnorm(1,prior_yrD$m0,prior_yrD$sd)),
                       phi=runif(1, range(x)[1], range(x)[2]),
                       d0=runif(1), dd=runif(1), dp=runif(1),
                       sigma=runif(1), sigma_ddsq=runif(1),
                       sigma_d0sq=runif(1), sigma_dpsq=runif(1))
  }
  para <- c("beta0", "delta", "phi","sigma", "d0", "dd","dp",
            "sigma_dp", "sigma_d0","sigma_dd","log_lik")
  return(list(para=para, data=bugs.data, inits=inits,
              n.chains=n.chains, theta=theta, cens=c(n, ncens)))
}

#### 3. no season
stan.in_seq_nsn <- function(infile=tempdata, x="Chla", y="part_microcystin",
                            n.chains=nchains, priorYr_rv=prior.rv){
  cens <- infile$part_microcystin<=0.1
  ncens <- sum(cens)
  logx_bar <- mean(log(infile[,x]), na.rm=T)
  n <- dim(infile)[1]-ncens
  xrange <- range(log(infile[,x]))
  theta <- 0.01*diff(xrange)
  if (ncens>0 & n>0){
    xtemp <- infile[,x]
    ytemp <- infile[,y]
    ycens <- log(ytemp[cens])
    xcens <- log(xtemp[cens])
    y <- log(ytemp[!cens])
    x <- log(xtemp[!cens])
  } else {
    xtemp <- infile[,x]
    ytemp <- infile[,y]
    y <- log(ytemp)
    x <- log(xtemp)
  } 
  prior_yr0 <- priorNIG(priorYr_rv$B0, priorYr_rv$sigma0)
  prior_yrD <- priorNIG(priorYr_rv$De, priorYr_rv$sigmad)
  prior_yrP <- priorNIG(priorYr_rv$Ph, priorYr_rv$sigmap)
  m0 <- prior_yr0$m0
  mD <- prior_yrD$m0
  mP <- prior_yrP$m0
  lmbd0=prior_yr0$lmbd
  lmbdD=prior_yrD$lmbd
  lmbdP=prior_yrP$lmbd
  al0=prior_yr0$alp
  alP=prior_yrP$alp
  alD=prior_yrD$alp
  bt0=prior_yr0$bet
  btP=prior_yrP$bet
  btD=prior_yrD$bet
  s0 <- sqrt(bt0/(al0-1))
  sD <- sqrt(btD/(alD-1))
  sP <- sqrt(btP/(alP-1))
  if (ncens>0 & n>0){
    bugs.data <- list(N=n, Ncens=ncens, y=array(y, dim=n),
                      ycens=array(ycens, dim=ncens),
                      x=array(x, dim=n), xcens=array(xcens,dim=ncens),
                      theta=theta, #beta1=0,
                      m0 = m0,  mD = mD, mP = mP,
                      lmbd0=lmbd0, lmbdD=lmbdD, lmbdP=lmbdP,
                      al0=al0, alP=alP, alD=alD,
                      bt0=bt0, btP=btP, btD=btD)
  } else {
    bugs.data <- list(N=ncens+n, y=y, x=x, 
                      theta=theta, #beta1=0,
                      m0 = m0,  mD = mD, mP = mP,
                      lmbd0=lmbd0, lmbdD=lmbdD, lmbdP=lmbdP,
                      al0=al0, alP=alP, alD=alD,
                      bt0=bt0, btP=btP, btD=btD)
  }
  inits <- list()
  for (i in 1:n.chains){
    inits[[i]] <- list(beta0=rnorm(1, m0, s0),
                       delta=abs(rnorm(1,m0,sD)),
                       phi=runif(1, range(x)[1], range(x)[2]),
                       mu0=rnorm(1, m0, s0), 
                       muD=abs(rnorm(1, mD,sD)),
                       muP=rnorm(1, mP, sP),
                       sigma=runif(1), sigmaPsq=runif(1), sigmaDsq=runif(1),
                       sigma0sq=runif(1))
  }
  para <- c("beta0", "delta", "phi","sigma",
            "mu0", "muD","muP", "sigmap","sigma0", "sigmad", "log_lik")
  return(list(para=para, data=bugs.data, inits=inits, n.chains=n.chains,
              theta=theta, cens=c(n, ncens)))
}

