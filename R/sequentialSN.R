## creating a subset of data for a given year
## Grouping sampling events to at least n=n_min
erie_sub <- function(file=eriedata, Yr=2018, n_min = 8,
                     x="Chla", y="part_microcystin"){
    subdata <- file[file$Year==Yr & !is.na(file[,y]),]
    cs <- cumsum(sz <- table(subdata$Rdate))
    stepcs <- 0
    k <- 1
    gr <- rep(1, length(cs))

    for (i in 1:length(cs)){
        if (stepcs >= n_min){
	    k <- k+1
	    stepcs <- 0
        }
        stepcs <- stepcs + sz[i]
        gr[i] <- k
    }
    subdata$gr <- gr[as.numeric(ordered(subdata$Rdate))]
    return(subdata)
}

priorNIG <- function(murv, sdrv, n0=20, setn0=F){
    ## fitrv: an rv object of stan fitted model
    ## setn0: whether to use non-informative prior

    tmp <- summary(murv)
    Ex <- tmp$mean
    Vx <- tmp$sd^2

    tmp <- summary(sdrv^2)
    Ev <- tmp$mean
    Vv <- tmp$sd^2
    if (setn0) alpha <- n0+1
    else alpha <- 2+Ev^2/Vv
    beta <- Ev*(alpha-1)
    lambda <- Ev/Vx
    ## limiting alpha+beta < 100
    while (alpha+beta > 100){
        alpha <- alpha/5
        beta <- beta/5
    }
    return(list(m0=Ex, lmbd=lambda, alp=alpha,  bet=beta))
}

priorIG <- function(varrv, n0=20, setn0=F){
    ## setn0: whether to use non-informative prior
    ##

    tmp <- summary(varrv)
    Ev <- tmp$mean
    Vv <- tmp$sd^2

    if (setn0) alpha <- n0+1
    else alpha <- 2+Ev^2/Vv
    beta <- Ev*(alpha-1)
    while (alpha+beta > 100){
        alpha <- alpha/5
        beta <- beta/5
    }
    return(list(alp=alpha,  bet=beta))
}

priorN <- function(xrv, setn0=F){
    ## setn0: whether to use non-informative prior
    tmp <- summary(xrv)
    return(list(m0=tmp$mean, sd=tmp$sd))
}

## Stan model (reporting $\sigma$)
## 1. Without seasonal effects
### without censored data
stan_model1 <- "
	      data{
	      int N; //the number of observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

	      real m0;
	      real mD;
	      real mP;

	      real lmbd0;
	      real lmbdD;
	      real lmbdP;

	      real al0;
	      real alP;
	      real alD;

	      real bt0;
	      real btP;
	      real btD;

	    }
	    parameters {
	      real beta0; //the regression parameters
	      real<lower=0> delta;
	      real phi; //change point

	      real<lower=0> sigma;

	      real mu0;
	      real muD;
	      real muP;

	      real<lower=0> sigma0sq;
	      real<lower=0> sigmaDsq;
	      real<lower=0> sigmaPsq;
	    }
	    transformed parameters {
	      real<lower=0> sigma0;
	      real<lower=0> sigmad;
	      real<lower=0> sigmap;
	      vector[N] mu;

	      sigma0 = sqrt(sigma0sq);
	      sigmad = sqrt(sigmaDsq);
	      sigmap = sqrt(sigmaPsq);
	      for (i in 1:N)
		mu[i] = beta0 + //beta1 * (x[i]-phi) +
                        delta * theta * log1p(exp((x[i]-phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
	      sigma0sq ~ inv_gamma(al0, bt0);
	      sigmaDsq ~ inv_gamma(alD, btD);
	      sigmaPsq ~ inv_gamma(alP, btP);

	      mu0 ~ normal(m0, sqrt(sigma0sq/lmbd0));
	      muD ~ normal(mD, sqrt(sigmaDsq/lmbdD));
	      muP ~ normal(mP, sqrt(sigmaPsq/lmbdP));

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

	      y ~ normal(mu, sigma);
	    }
            generated quantities{
              vector[N] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
            }
"
### without seasonal effect
### with censored data
stan_model1_cens <- "
	      data{
	      int N; //number of uncensored observations
              int Ncens;
	      vector[N] y; //the response
	      vector[Ncens] ycens; //the response
	      vector[N] x;
	      vector[Ncens] xcens;

	      real theta;
	      //real beta1;

	      real m0;
	      real mD;
	      real mP;

	      real lmbd0;
	      real lmbdD;
	      real lmbdP;

	      real al0;
	      real alP;
	      real alD;

	      real bt0;
	      real btP;
	      real btD;

	    }
	    parameters {
	      real beta0; //the regression parameters
	      real<lower=0> delta;
	      real phi; //change point

	      real<lower=0> sigma;

	      real mu0;
	      real muD;
	      real muP;

	      real<lower=0> sigma0sq;
	      real<lower=0> sigmaDsq;
	      real<lower=0> sigmaPsq;
	    }
	    transformed parameters {
	      real<lower=0> sigma0;
	      real<lower=0> sigmad;
	      real<lower=0> sigmap;
	      vector[N] mu;
	      vector[Ncens] mucens;

	      sigma0 = sqrt(sigma0sq);
	      sigmad = sqrt(sigmaDsq);
	      sigmap = sqrt(sigmaPsq);
	      for (i in 1:N)
		mu[i] = beta0 + //beta1 * (x[i]-phi) +
                        delta * theta * log1p(exp((x[i]-phi)/theta));
       	      for (i in 1:Ncens)
		mucens[i] = beta0 + //beta1 * (xcens[i]-phi) +
                        delta * theta * log1p(exp((xcens[i]-phi)/theta));

	    }
	    model {
	      sigma ~ cauchy(0, 1);
	      sigma0sq ~ inv_gamma(al0, bt0);
	      sigmaDsq ~ inv_gamma(alD, btD);
	      sigmaPsq ~ inv_gamma(alP, btP);

	      mu0 ~ normal(m0, sqrt(sigma0sq/lmbd0));
	      muD ~ normal(mD, sqrt(sigmaDsq/lmbdD));
	      muP ~ normal(mP, sqrt(sigmaPsq/lmbdP));

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

              target += normal_lpdf(y | mu, sigma);
              target += normal_lcdf(ycens | mucens, sigma);
            }
            generated quantities{
              vector[N+Ncens] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
              for (i in 1:Ncens){
                log_lik[i+N] = normal_lcdf(ycens[i]|mucens[i], sigma);
              }
            }
"

### without seasonal effect
### all censored data
stan_model1_5 <- "
	      data{
	      int N; //number of uncensored observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

	      real m0;
	      real mD;
	      real mP;

	      real lmbd0;
	      real lmbdD;
	      real lmbdP;

	      real al0;
	      real alP;
	      real alD;

	      real bt0;
	      real btP;
	      real btD;

	    }
	    parameters {
	      real beta0; //the regression parameters
	      real<lower=0> delta;
	      real phi; //change point

	      real<lower=0> sigma;

	      real mu0;
	      real muD;
	      real muP;

	      real<lower=0> sigma0sq;
	      real<lower=0> sigmaDsq;
	      real<lower=0> sigmaPsq;
	    }
	    transformed parameters {
	      real<lower=0> sigma0;
	      real<lower=0> sigmad;
	      real<lower=0> sigmap;
	      vector[N] mu;

	      sigma0 = sqrt(sigma0sq);
	      sigmad = sqrt(sigmaDsq);
	      sigmap = sqrt(sigmaPsq);
	      for (i in 1:N)
		mu[i] = beta0 + //beta1 * (x[i]-phi) +
                        delta * theta * log1p(exp((x[i]-phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
	      sigma0sq ~ inv_gamma(al0, bt0);
	      sigmaDsq ~ inv_gamma(alD, btD);
	      sigmaPsq ~ inv_gamma(alP, btP);

	      mu0 ~ normal(m0, sqrt(sigma0sq/lmbd0));
	      muD ~ normal(mD, sqrt(sigmaDsq/lmbdD));
	      muP ~ normal(mP, sqrt(sigmaPsq/lmbdP));

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

              target += normal_lcdf(y | mu, sigma);
            }
            generated quantities{
              vector[N] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lcdf(y[i]|mu[i], sigma);
              }
           }
"

## 1. Additive seasonal effects
### without censored data
stan_model3 <- "
	      data{
	      int N; //the number of observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

              real mu_dp;
              real<lower=0> sigma_dp;
              real mu_d0;
              real<lower=0> sigma_d0;
              real mu_dd;
              real<lower=0> sigma_dd;

	      real muP;
              real<lower=0> sigmap;
	      real mu0;
              real<lower=0> sigma0;
	      real<lower=0> muD; 
              real<lower=0> sigmad;
	    }
	    parameters {
	      real beta0; //the regression parameters
              real d0; //seasonal effects
	      real<lower=0> delta;
              real dd; 
	      real phi; //change point
              real dp;
	      real<lower=0> sigma;
	    }
	    transformed parameters {
	      vector[N] mu;
              real Beta0;
              real<lower=0> Delta;
              real Phi;

	      Beta0 = beta0 + d0;
              Delta = delta + dd;
              Phi = phi + dp;
	      for (i in 1:N)
		mu[i] = Beta0 + //beta1 * (x[i]-phi) +
                        Delta * theta * log1p(exp((x[i]-Phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
	      
              dp ~ normal(mu_dp, sigma_dp);
              d0 ~ normal(mu_d0, sigma_d0);
              dd ~ normal(mu_dd, sigma_dd);

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

	      y ~ normal(mu, sigma);
	    }
            generated quantities{
              vector[N] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
            }
"

### with partial censored data
stan_model3_cens <- "
	      data{
	      int N; //number of uncensored observations
              int Ncens;
	      vector[N] y; //the response
	      vector[Ncens] ycens; //the response
	      vector[N] x;
	      vector[Ncens] xcens;

	      real theta;
	      //real beta1;

              real mu_dp;
              real<lower=0> sigma_dp;
              real mu_d0;
              real<lower=0> sigma_d0;
              real mu_dd;
              real<lower=0> sigma_dd;

	      real muP;
              real<lower=0> sigmap;
	      real mu0;
              real<lower=0> sigma0;
	      real<lower=0> muD; 
              real<lower=0> sigmad;
	    }
	    parameters {
	      real beta0; //the regression parameters
              real d0; //seasonal effects
	      real<lower=0> delta;
              real dd; 
	      real phi; //change point
              real dp;
	      real<lower=0> sigma;
	    }
	    transformed parameters {
	      vector[N] mu;
	      vector[Ncens] mucens;

              real Beta0;
              real<lower=0> Delta;
              real Phi;

	      Beta0 = beta0 + d0;
              Delta = delta + dd;
              Phi = phi + dp;

	      for (i in 1:N)
		mu[i] = Beta0 + //beta1 * (x[i]-phi) +
                        Delta * theta * log1p(exp((x[i]-Phi)/theta));
	      for (i in 1:Ncens)
		mucens[i] = Beta0 + //beta1 * (xcens[i]-phi) +
                        Delta * theta * log1p(exp((xcens[i]-Phi)/theta));
	    }
	    model {
              dp ~ normal(mu_dp, sigma_dp);
              d0 ~ normal(mu_d0, sigma_d0);
              dd ~ normal(mu_dd, sigma_dd);

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

	      sigma ~ cauchy(0, 1);

              target += normal_lpdf(y | mu, sigma);
              target += normal_lcdf(ycens | mucens, sigma);
            }
            generated quantities{
              vector[N+Ncens] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
              for (i in 1:Ncens){
                log_lik[i+N] = normal_lcdf(ycens[i]|mucens[i], sigma);
              }
            }
"
## with all censored data
stan_model3_5 <- "
	      data{
	      int N; //the number of observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

              real mu_dp;
              real<lower=0> sigma_dp;
              real mu_d0;
              real<lower=0> sigma_d0;
              real mu_dd;
              real<lower=0> sigma_dd;

	      real muP;
              real<lower=0> sigmap;
	      real mu0;
              real<lower=0> sigma0;
	      real<lower=0> muD; 
              real<lower=0> sigmad;
	    }
	    parameters {
	      real beta0; //the regression parameters
              real d0; //seasonal effects
	      real<lower=0> delta;
              real dd; 
	      real phi; //change point
              real dp;
	      real<lower=0> sigma;
	    }
	    transformed parameters {
	      vector[N] mu;
              real Beta0;
              real<lower=0> Delta;
              real Phi;

	      Beta0 = beta0 + d0;
              Delta = delta + dd;
              Phi = phi + dp;
	      for (i in 1:N)
		mu[i] = Beta0 + //beta1 * (x[i]-phi) +
                        Delta * theta * log1p(exp((x[i]-Phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
	      
              dp ~ normal(mu_dp, sigma_dp);
              d0 ~ normal(mu_d0, sigma_d0);
              dd ~ normal(mu_dd, sigma_dd);

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

              target += normal_lcdf(y | mu, sigma);
	      //y ~ normal(mu, sigma);
	    }
            generated quantities{
              vector[N] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
            }
"

## 2. Nested seasonal effects Seasonal and annual updating
### without censored data
stan_model4 <- "
	      data{
	      int N; //the number of observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

	      real muP;
              real<lower=0> sigmap;
	      real mu0;
              real<lower=0> sigma0;
	      real muD; 
              real<lower=0> sigmad;

	      real<lower=0> al0;
	      real<lower=0> bt0;
	      real<lower=0> alD;
	      real<lower=0> btD;
	      real<lower=0> alP;
 	      real<lower=0> btP;
	    }
	    parameters {
	      real beta0; //the regression parameters
              real d0; //seasonal effects
	      real<lower=0> delta;
              real dd; 
	      real phi; //change point
              real dp;
	      real<lower=0> sigma;
	      real<lower=0> sigma_d0sq;
	      real<lower=0> sigma_ddsq;
	      real<lower=0> sigma_dpsq;
	    }
	    transformed parameters {
	      vector[N] mu;
              real Beta0;
              real<lower=0> Delta;
              real Phi;
              real<lower=0> sigma_d0;
              real<lower=0> sigma_dp;
              real<lower=0> sigma_dd;

	      Beta0 = beta0 + d0;
              Delta = delta + dd;
              Phi = phi + dp;

	      sigma_d0 = sqrt(sigma_d0sq);
	      sigma_dd = sqrt(sigma_ddsq);
	      sigma_dp = sqrt(sigma_dpsq);

	      for (i in 1:N)
		mu[i] = Beta0 + //beta1 * (x[i]-phi) +
                        Delta * theta * log1p(exp((x[i]-Phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
              
	      sigma_d0sq ~ inv_gamma(al0, bt0);
	      sigma_ddsq ~ inv_gamma(alD, btD);
	      sigma_dpsq ~ inv_gamma(alP, btP);
	      
              dp ~ normal(0, sigma_dp);
              d0 ~ normal(0, sigma_d0);
              dd ~ normal(0, sigma_dd);

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

	      y ~ normal(mu, sigma);
	    }
            generated quantities{
              vector[N] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
            }
"

### with censored data
stan_model4_cens <- "
	      data{
	      int N; //number of uncensored observations
              int Ncens;
	      vector[N] y; //the response
	      vector[Ncens] ycens; //the response
	      vector[N] x;
	      vector[Ncens] xcens;

	      real theta;
	      //real beta1;

	      real muP;
              real<lower=0> sigmap;
	      real mu0;
              real<lower=0> sigma0;
	      real muD; 
              real<lower=0> sigmad;

	      real<lower=0> al0;
	      real<lower=0> bt0;
	      real<lower=0> alD;
	      real<lower=0> btD;
	      real<lower=0> alP;
 	      real<lower=0> btP;
	    }
	    parameters {
	      real beta0; //the regression parameters
              real d0; //seasonal effects
	      real<lower=0> delta;
              real dd; 
	      real phi; //change point
              real dp;
	      real<lower=0> sigma;
	      real<lower=0> sigma_d0sq;
	      real<lower=0> sigma_ddsq;
	      real<lower=0> sigma_dpsq;
	    }
	    transformed parameters {
	      vector[N] mu;
	      vector[Ncens] mucens;

              real Beta0;
              real<lower=0> Delta;
              real Phi;

	      real<lower=0> sigma_d0;
	      real<lower=0> sigma_dd;
	      real<lower=0> sigma_dp;

	      Beta0 = beta0 + d0;
              Delta = delta + dd;
              Phi = phi + dp;

	      sigma_d0 = sqrt(sigma_d0sq);
	      sigma_dd = sqrt(sigma_ddsq);
	      sigma_dp = sqrt(sigma_dpsq);

	      for (i in 1:N)
		mu[i] = Beta0 + //beta1 * (x[i]-phi) +
                        Delta * theta * log1p(exp((x[i]-Phi)/theta));
	      for (i in 1:Ncens)
		mucens[i] = Beta0 + //beta1 * (xcens[i]-phi) +
                        Delta * theta * log1p(exp((xcens[i]-Phi)/theta));
	    }
	    model {
	      sigma_d0sq ~ inv_gamma(al0, bt0);
	      sigma_ddsq ~ inv_gamma(alD, btD);
	      sigma_dpsq ~ inv_gamma(alP, btP);
	      
              dp ~ normal(0, sigma_dp);
              d0 ~ normal(0, sigma_d0);
              dd ~ normal(0, sigma_dd);

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

	      sigma ~ cauchy(0, 1);

              target += normal_lpdf(y | mu, sigma);
              target += normal_lcdf(ycens | mucens, sigma);
            }
            generated quantities{
              vector[N+Ncens] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
              for (i in 1:Ncens){
                log_lik[i+N] = normal_lcdf(ycens[i]|mucens[i], sigma);
              }
            }
"

## with all censored data
stan_model4_5 <- "
	      data{
	      int N; //the number of observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

	      real muP;
              real<lower=0> sigmap;
	      real mu0;
              real<lower=0> sigma0;
	      real muD; 
              real<lower=0> sigmad;

	      real<lower=0> al0;
	      real<lower=0> bt0;
	      real<lower=0> alD;
	      real<lower=0> btD;
	      real<lower=0> alP;
 	      real<lower=0> btP;
	    }
	    parameters {
	      real beta0; //the regression parameters
              real d0; //seasonal effects
	      real<lower=0> delta;
              real dd; 
	      real phi; //change point
              real dp;
	      real<lower=0> sigma;
	      real<lower=0> sigma_d0sq;
	      real<lower=0> sigma_ddsq;
	      real<lower=0> sigma_dpsq;
	    }
	    transformed parameters {
	      vector[N] mu;
              real Beta0;
              real<lower=0> Delta;
              real Phi;
              real<lower=0> sigma_d0;
              real<lower=0> sigma_dp;
              real<lower=0> sigma_dd;

	      Beta0 = beta0 + d0;
              Delta = delta + dd;
              Phi = phi + dp;

	      sigma_d0 = sqrt(sigma_d0sq);
	      sigma_dd = sqrt(sigma_ddsq);
	      sigma_dp = sqrt(sigma_dpsq);

	      for (i in 1:N)
		mu[i] = Beta0 + //beta1 * (x[i]-phi) +
                        Delta * theta * log1p(exp((x[i]-Phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
              
	      sigma_d0sq ~ inv_gamma(al0, bt0);
	      sigma_ddsq ~ inv_gamma(alD, btD);
	      sigma_dpsq ~ inv_gamma(alP, btP);
	      
              dp ~ normal(0, sigma_dp);
              d0 ~ normal(0, sigma_d0);
              dd ~ normal(0, sigma_dd);

	      phi ~ normal(muP, sigmap);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmad);

	      target += normal_lcdf(y | mu, sigma);;
	    }
            generated quantities{
              vector[N] log_lik;
              for (i in 1:N){
                log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);
              }
            }
"


stan.in_seq <- function(infile, x="Chla", y="part_microcystin",
                        n.chains=nchains, priorYr_rv, priorSn_rv,
                        additive=T){
	keep <-  (infile[,x] > 0) & (infile[,y] >0)
	infile <- infile[keep & !is.na(keep),]
        cens <- infile$part_microcystin<=0.10
        ncens <- sum(cens)
        logx_bar <- mean(log(infile[,x]), na.rm=T)
        n <- dim(infile)[1]-ncens
        prior_yr0 <- priorN(priorYr_rv$B0)
        prior_yrD <- priorN(priorYr_rv$De)
        prior_yrP <- priorN(priorYr_rv$Ph)
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
        inits <- list()
        if (ncens>0 & n>0){
            if (additive){
                prior_snd0 <- priorN(priorSn_rv$d0)
                prior_sndp <- priorN(priorSn_rv$dp)
                prior_sndd <- priorN(priorSn_rv$dd)
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
                                  muD = prior_yrD$m0, sigmaD = prior_yrD$sd,
                                  muP = prior_yrP$m0, sigmaP = prior_yrP$sd)
            } else {
                prior_snd0 <- priorIG(priorSn_rv$d0)
                prior_sndd <- priorIG(priorSn_rv$dd)
                prior_sndp <- priorIG(priorSn_rv$dp)
                bugs.data <- list(N=n, Ncens=ncens, y=array(y, dim=n),
                                  ycens=array(ycens, dim=ncens), x=array(x, dim=n),
                                  xcens=array(xcens, dim=ncens),
                                  theta=theta, #beta1=0,
                                  al0=prior_snd0$alp, bt0=prior_snd0$bet,
                                  alP=prior_sndp$alp, btP=prior_sndp$bet, 
                                  alD=prior_sndd$alp, btD=prior_sndd$bet,
                                  mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                                  muD = prior_yrD$m0, sigmaD = prior_yrD$sd,
                                  muP = prior_yrP$m0, sigmaP = prior_yrP$sd)
            }
	} else {
            if (additive) {
                prior_snd0 <- priorN(priorSn_rv$d0)
                prior_sndp <- priorN(priorSn_rv$dp)
                prior_sndd <- priorN(priorSn_rv$dd)
                bugs.data <- list(N=ncens+n, y=y, x=x, 
                                  theta=theta, #beta1=0,
                                  mu_dp = prior_sndp$m0,
                                  sigma_dp=prior_sndp$sd,
                                  mu_d0 = prior_snd0$m0,
                                  sigma_d0=prior_snd0$sd,
                                  mu_dd = prior_sndd$m0,
                                  sigma_dd=prior_sndd$sd,
                                  mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                                  muD = prior_yrD$m0, sigmaD = prior_yrD$sd,
                                  muP = prior_yrP$m0, sigmaP = prior_yrP$sd)
            } else {
                prior_snd0 <- priorIG(priorSn_rv$d0)
                prior_sndd <- priorIG(priorSn_rv$dd)
                prior_sndp <- priorIG(priorSn_rv$dp)
                bugs.data <- list(N=ncens+n, y=y, x=x, 
                                  theta=theta, #beta1=0,
                                  al0=prior_snd0$alp, bt0=prior_snd0$bet,
                                  alP=prior_sndp$alp, btP=prior_sndp$bet,
                                  alD=prior_sndd$alp, btD=prior_sndd$bet,
                                  mu0 = prior_yr0$m0, sigma0 = prior_yr0$sd,
                                  muD = prior_yrD$m0, sigmad = prior_yrD$sd,
                                  muP = prior_yrP$m0, sigmap = prior_yrP$sd)
            }
        }
        for (i in 1:n.chains){
            if (additive)
                inits[[i]] <- list(beta0=rnorm(1, prior_yr0$m0, prior_yr0$sd),
                                   delta=abs(rnorm(1,prior_yrD$m0,prior_yrD$sd)),
                                   phi=runif(1, range(x)[1], range(x)[2]),
                                   d0=runif(1), dd=runif(1), dp=runif(1),
                                   sigma=runif(1))
            else
                inits[[i]] <- list(beta0=rnorm(1,  prior_yr0$m0, prior_yr0$sd),
                                   delta=abs(rnorm(1,prior_yrD$m0,prior_yrD$sd)),
                                   phi=runif(1, range(x)[1], range(x)[2]),
                                   d0=runif(1), dd=runif(1), dp=runif(1),
                                   sigma=runif(1), sigma_ddsq=runif(1),
                                   sigma_d0sq=runif(1), sigma_dpsq=runif(1))
        }
	if (additive) para <- c("beta0", "delta", "phi","sigma",
                                "d0", "dd","dp")
        else para <- c("beta0", "delta", "phi","sigma",
                       "d0", "dd","dp",
                       "sigma_dp", "sigma_d0","sigma_dd")
	return(list(para=para, data=bugs.data,
                    inits=inits,n.chains=n.chains, theta=theta,
                    cens=(ncens>0)))
}

