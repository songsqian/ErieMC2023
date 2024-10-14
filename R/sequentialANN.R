## Stan model (reporting $\sigma$)
## 5. Annual updating -- assuming censored values always present

stan_model5_add <- "
          data{
            int N; //number of observations
            int Ncens; //number of censored observations
            vector[N] y; //the response
            vector[N] x;
            vector[Ncens] xcens;
            int S; // number of seasons
            int season[N];
            int sncens[Ncens];

            real theta;
            //real beta1;
            vector[Ncens] ycens; //0.1

            real<lower=0> alR0;
            real<lower=0> btR0;
            real<lower=0> alRd;
            real<lower=0> btRd;
            real<lower=0> alRp;
            real<lower=0> btRp;
            real<lower=0> alS0;
            real<lower=0> btS0;
            real<lower=0> alSd;
            real<lower=0> btSd;
            real<lower=0> alSp;
            real<lower=0> btSp;

            real<lower=0> lmdp;
            real<lower=0> lmdd;
            real<lower=0> lmd0;
            real muP;
            real mu0;
            real muD;
            vector[S] snmu0;
            vector[S] snmuP;
            vector[S] snmuD;
            vector<lower=0>[S] sigmu0;
            vector<lower=0>[S] sigmuP;
            vector<lower=0>[S] sigmuD;
          }
          parameters {
            real beta0; //the regression parameters
            real<lower=0> delta;
            real phi; //change point
            vector[S] sn0;
            vector[S] snP;
            vector[S] snD;
            vector[S] sn0m;
            vector[S] snPm;
            vector[S] snDm;
            real<lower=0> sigma;
            real<lower=0> sigmaR0;
            real<lower=0> sigmaRd;
            real<lower=0> sigmaRp;
            real<lower=0> sigmaS0;
            real<lower=0> sigmaSd;
            real<lower=0> sigmaSp;
          }
          transformed parameters {
            vector[N] mu;
            vector[Ncens] mucens;
            vector[N] Beta0;
            vector[Ncens] B0cens;
            vector[N] Delta;
            vector[Ncens] Decens;
            vector[N] Phi;
            vector[Ncens] Phcens;

            for (i in 1:N){
              Beta0[i] = beta0+sn0[season[i]];
              Phi[i] = phi  +snP[season[i]];
              Delta[i] = delta+snD[season[i]];
              mu[i] = Beta0[i] + //beta1 * (x[i]-Phi[i]) +
                      Delta[i] * theta * log1p(exp((x[i]-Phi[i])/theta));
            }
            for (i in 1:Ncens){
              B0cens[i] =beta0+sn0[sncens[i]];
              Phcens[i] =phi  +snP[sncens[i]];
              Decens[i] =delta+snD[sncens[i]];
              mucens[i] =B0cens[i]+//beta1*(xcens[i]-Phcens[i]) +
                         Decens[i]*theta*log1p(exp((xcens[i]-Phcens[i])/theta));
            }
          }
          model {
            sigma ~ cauchy(0, 1);
            sigmaR0 ~ inv_gamma(alR0, btR0);
            sigmaRd ~ inv_gamma(alRd, btRd);
            sigmaRp ~ inv_gamma(alRp, btRp);
            sigmaS0 ~ inv_gamma(alS0, btS0);
            sigmaSd ~ inv_gamma(alSd, btSd);
            sigmaSp ~ inv_gamma(alSp, btSp);

            phi ~ normal(muP,sigmaRp/sqrt(lmdp));
            beta0 ~ normal(mu0, sigmaR0/sqrt(lmd0));
            delta ~ normal(muD, sigmaRd/sqrt(lmdd));
            sn0m ~ normal(snmu0, sigmaS0);
            snDm ~ normal(snmuD, sigmaSd);
            snPm ~ normal(snmuP, sigmaSp);
            sn0 ~ normal(sn0m, sigmaS0);
            snD ~ normal(snDm, sigmaSd);
            snP ~ normal(snPm, sigmaSp);

            target += normal_lpdf(y | mu, sigma);
            target += normal_lcdf(ycens | mucens, sigma);
          }
          generated quantities {
            real Ph;
            vector[S] delPhS;
            real B0;
            vector[S] delB0S;
            real De;
            vector[S] delDS;
            vector[N+Ncens] log_lik;;
            Ph = phi + mean(snP[]);
            B0 = beta0 + mean(sn0[]);
            De = delta + mean(snD[]);
            for (i in 1:S){
              delPhS[i] = snP[i] - mean(snP[]);
              delB0S[i] = sn0[i] - mean(sn0[]);
              delDS[i]  = snD[i] - mean(snD[]);
            }
            for (i in 1:N){
              log_lik[i] <- normal_lpdf(y[i] | mu[i], sigma);
            }
            for (i in 1:Ncens){
              log_lik[i+N] <- normal_lcdf(ycens[i] | mucens[i], sigma);
            }
          }
"

stan_model5_nstd <- "
          data{
            int N; //number of observations
            int Ncens; //number of censored observations
            vector[N] y; //the response
            vector[N] x;
            vector[Ncens] xcens;
            int S;
            int season[N];
            int sncens[Ncens];

            real theta;
            //real beta1;
            vector[Ncens] ycens; //0.1

            real<lower=0> alR0;
            real<lower=0> btR0;
            real<lower=0> alRd;
            real<lower=0> btRd;
            real<lower=0> alRp;
            real<lower=0> btRp;
            real<lower=0> alS0;
            real<lower=0> btS0;
            real<lower=0> alSd;
            real<lower=0> btSd;
            real<lower=0> alSp;
            real<lower=0> btSp;

            real<lower=0> lmdp;
            real<lower=0> lmdd;
            real<lower=0> lmd0;
            real muP;
            real mu0;
            real muD;
          }
          parameters {
            real beta0; //the regression parameters
            real<lower=0> delta;
            real phi; //change point
            vector[S] sn0;
            vector[S] snP;
            vector[S] snD;
            real<lower=0> sigma;
            real<lower=0> sigmaR0;
            real<lower=0> sigmaRd;
            real<lower=0> sigmaRp;
            real<lower=0> sigmaS0;
            real<lower=0> sigmaSd;
            real<lower=0> sigmaSp;
          }
          transformed parameters {
            vector[N] mu;
            vector[Ncens] mucens;
            vector[N] Beta0;
            vector[Ncens] B0cens;
            vector[N] Delta;
            vector[Ncens] Decens;
            vector[N] Phi;
            vector[Ncens] Phcens;

            for (i in 1:N){
              Beta0[i] = beta0+sn0[season[i]];
              Phi[i] = phi  +snP[season[i]];
              Delta[i] = delta+snD[season[i]];
              mu[i] = Beta0[i] + //beta1 * (x[i]-Phi[i]) +
                      Delta[i] * theta * log1p(exp((x[i]-Phi[i])/theta));
            }
            for (i in 1:Ncens){
              B0cens[i] =beta0+sn0[sncens[i]];
              Phcens[i] =phi  +snP[sncens[i]];
              Decens[i] =delta+snD[sncens[i]];
              mucens[i] =B0cens[i]+//beta1*(xcens[i]-Phcens[i]) +
                         Decens[i]*theta*log1p(exp((xcens[i]-Phcens[i])/theta));
            }
          }
          model {
            sigma ~ cauchy(0, 1);
            sigmaR0 ~ inv_gamma(alR0, btR0);
            sigmaRd ~ inv_gamma(alRd, btRd);
            sigmaRp ~ inv_gamma(alRp, btRp);
            sigmaS0 ~ inv_gamma(alS0, btS0);
            sigmaSd ~ inv_gamma(alSd, btSd);
            sigmaSp ~ inv_gamma(alSp, btSp);

            phi ~ normal(muP,sigmaRp/sqrt(lmdp));
            beta0 ~ normal(mu0, sigmaR0/sqrt(lmd0));
            delta ~ normal(muD, sigmaRd/sqrt(lmdd));
            sn0 ~ normal(0, sigmaS0);
            snD ~ normal(0, sigmaSd);
            snP ~ normal(0, sigmaSp);

            target += normal_lpdf(y | mu, sigma);
            target += normal_lcdf(ycens | mucens, sigma);
          }
          generated quantities {
            real Ph;
            vector[S] delPhS;
            real B0;
            vector[S] delB0S;
            real De;
            vector[S] delDS;
            vector[N+Ncens] log_lik;;
            Ph = phi + mean(snP[]);
            B0 = beta0 + mean(sn0[]);
            De = delta + mean(snD[]);
            for (i in 1:S){
              delPhS[i] = snP[i] - mean(snP[]);
              delB0S[i] = sn0[i] - mean(sn0[]);
              delDS[i]  = snD[i] - mean(snD[]);
            }
            for (i in 1:N){
              log_lik[i] <- normal_lpdf(y[i] | mu[i], sigma);
            }
            for (i in 1:Ncens){
              log_lik[i+N] <- normal_lcdf(ycens[i] | mucens[i], sigma);
            }
          }
"

stan.in_seqA_add <- function(infile, x="Chla", y="part_microcystin",
                              n.chains=nchains, priorYr_rv, priorSn_rv){
  keep <-  (infile[,x] > 0) & (infile[,y] >0)
  infile <- infile[keep & !is.na(keep),]
  cens <- infile[,y]<=0.1
  ncens <- sum(cens)
  n <- dim(infile)[1]-ncens
  seasons <- infile$Mon_ord
  nsns <- max(seasons)
  season <- seasons[!cens]
  sncens <- seasons[cens]
  prior_yr0 <- priorNIG(priorYr_rv$B0, priorYr_rv$sigmaR0)
  prior_yrD <- priorNIG(priorYr_rv$De, priorYr_rv$sigmaRd)
  prior_yrP <- priorNIG(priorYr_rv$Ph, priorYr_rv$sigmaRp)
  prior_sn0 <- priorIG(priorYr_rv$sigmaS0)
  prior_snP <- priorIG(priorYr_rv$sigmaSp)
  prior_snD <- priorIG(priorYr_rv$sigmaSd)
  snDmusig <- NULL
  snPmusig <- NULL
  sn0musig <- NULL
  for (i in 1:nsns){
    sn0musig <- rbind(sn0musig, priorN(priorSn_rv[[i]]$d0))
    snPmusig <- rbind(snPmusig, priorN(priorSn_rv[[i]]$dp))
    snDmusig <- rbind(snDmusig, priorN(priorSn_rv[[i]]$dd))
  }
  xtemp <- infile[,x]
  y <- log(infile[!cens,y])
  ycens <- rep(log(0.01), ncens)
  x <- log(xtemp[!cens])
  xcens <- log(xtemp[cens])
  inits <- list()
  theta <- 0.01*diff(range(x))
  bugs.data <- list(N=n, Ncens=ncens, S=nsns, y=y,
                    ycens=array(ycens, dim=ncens),
                    x=x, xcens=array(xcens,dim=ncens),
                    season=season, sncens=sncens, theta=theta, #beta1=0,
                    alR0=prior_yr0$alp, btR0=prior_yr0$bet,
                    mu0=prior_yr0$m0, lmd0=prior_yr0$lmbd,
                    alRp=prior_yrP$alp, btRp=prior_yrP$bet,
                    muP=prior_yrP$m0, lmdp=prior_yrP$lmbd,
                    alRd=prior_yrD$alp, btRd=prior_yrD$bet,
                    muD=prior_yrD$m0, lmdd=prior_yrD$lmbd,
                    alS0=prior_sn0$alp, btS0=prior_sn0$bet,
                    alSp=prior_snP$alp, btSp=prior_snP$bet,
                    alSd=prior_snD$alp, btSd=prior_snD$bet,   
                    snmu0=unlist(sn0musig[,1]), sigmu0=unlist(sn0musig[,2]), 
                    snmuP=unlist(snPmusig[,1]), sigmuP=unlist(snPmusig[,2]), 
                    snmuD=unlist(snDmusig[,1]), sigmuD=unlist(snDmusig[,2]))
  for (i in 1:n.chains)
    inits[[i]] <- list(beta0=rnorm(1, prior_yr0$m0),
                       delta=abs(rnorm(1,prior_yrD$m0)),
                       phi=runif(1, range(x)[1], range(x)[2]),
                       sn0=rnorm(nsns), snD=rnorm(nsns), snP=rnorm(nsns),
                       sigma=runif(1), sn0m=rnorm(nsns), snPm=rnorm(nsns),
                       snDm=rnorm(nsns), sigmaR0=runif(1), sigmaRd=runif(1), 
                       sigmaRp=runif(1), sigmaS0=runif(1), sigmaSd=runif(1), 
                       sigmaSp=runif(1))
  
  para <- c("Ph","De","B0",
            "delB0S", "delDS", "delPhS",
            "sigma", "sigmaR0", "sigmaRd","sigmaRp",
            "sigmaS0","sigmaSd","sigmaSp")
  return(list(para=para, data=bugs.data,
              inits=inits,n.chains=n.chains, theta=theta,
              cens=(ncens>0)))
}

stan.in_seqA_nstd <- function(infile, x="Chla", y="part_microcystin",
                         n.chains=nchains, priorYr_rv, priorSn_rv,
                         biweek=F){
	keep <-  (infile[,x] > 0) & (infile[,y] >0)
	infile <- infile[keep & !is.na(keep),]
        cens <- infile[,y]<=0.1
        ncens <- sum(cens)
        n <- dim(infile)[1]-ncens
        if (biweek) {seasons <- as.numeric(ordered(infile$sn_nstd2w))
        }else {seasons <- as.numeric(ordered(infile$sn_nstd))}
        nsns <- max(seasons)
        season <- seasons[!cens]
        sncens <- seasons[cens]
        prior_yr0 <- priorNIG(priorYr_rv$B0, priorYr_rv$sigmaR0)
        prior_yrD <- priorNIG(priorYr_rv$De, priorYr_rv$sigmaRd)
        prior_yrP <- priorNIG(priorYr_rv$Ph, priorYr_rv$sigmaRp)
        prior_sn0 <- priorIG(priorSn_rv$sigmaS0)
        prior_snP <- priorIG(priorSn_rv$sigmaSp)
        prior_snD <- priorIG(priorSn_rv$sigmaSd)
        xtemp <- infile[,x]
        y <- log(infile[!cens,y])
        ycens <- rep(log(0.01), ncens)
        x <- log(xtemp[!cens])
        xcens <- log(xtemp[cens])
	inits <- list()
	theta <- 0.01*diff(range(x))
        bugs.data <- list(N=n, Ncens=ncens, S=nsns, y=y,
                          ycens=array(ycens, dim=ncens),
                          x=x, xcens=array(xcens,dim=ncens),
                          season=season, sncens=sncens, theta=theta, #beta1=0,
                          alR0=prior_yr0$alp, btR0=prior_yr0$bet,
                          mu0=prior_yr0$m0, lmd0=prior_yr0$lmbd,
                          alRp=prior_yrP$alp, btRp=prior_yrP$bet,
                          muP=prior_yrP$m0, lmdp=prior_yrP$lmbd,
                          alRd=prior_yrD$alp, btRd=prior_yrD$bet,
                          muD=prior_yrD$m0, lmdd=prior_yrD$lmbd,
                          alS0=prior_sn0$alp, btS0=prior_sn0$bet,
                          alSp=prior_snP$alp, btSp=prior_snP$bet,
                          alSd=prior_snD$alp, btSd=prior_snD$bet)
        for (i in 1:n.chains)
            inits[[i]] <- list(beta0=rnorm(1, prior_yr0$m0),
                               delta=abs(rnorm(1,prior_yrD$m0)),
                               phi=runif(1, range(x)[1], range(x)[2]),
                               sn0=rnorm(nsns), snD=rnorm(nsns), snP=rnorm(nsns),
                               sigma=runif(1))
            para <- c("Ph","De","B0",
                      "delB0S", "delDS", "delPhS",
                      "sigma", "sigmaR0", "sigmaRd","sigmaRp",
                      "sigmaS0","sigmaSd","sigmaSp")
        return(list(para=para, data=bugs.data,
                    inits=inits,n.chains=n.chains, theta=theta,
                    cens=(ncens>0)))
}
