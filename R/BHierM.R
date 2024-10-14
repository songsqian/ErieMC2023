## A hockey stick model
## Censored data version -- with additive/nested seasonal component
stan_model2 <- "
          data{
            int N; //number of observations
            int Ncens; //number of censored observations
            vector[N] y; //the response
            vector[N] x;
            vector[Ncens] xcens;
            int R;
            int S;
            int region[N];
            int regcens[Ncens];
            int season[N];
            int sncens[Ncens];

            real theta;
            //real beta1;
            vector[Ncens] ycens; //0.01
          }
          parameters {
            real beta0; //the regression parameters
            real<lower=0> delta;
            real phi; //change point

            vector[R] zre0;
            vector[R] zreP;
            vector[R] zreD;
            vector[S] zsn0;
            vector[S] zsnP;
            vector[S] zsnD;
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
            vector[R] re0;
            vector[R] reP;
            vector[R] reD;
            vector[S] sn0;
            vector[S] snP;
            vector[S] snD;

            re0 = sigmaR0*zre0;
            reD = sigmaRd*zreD;
            reP = sigmaRp*zreP;
            sn0 = sigmaS0*zsn0;
            snD = sigmaSd*zsnD;
            snP = sigmaSp*zsnP;

            for (i in 1:N){
              Beta0[i] = beta0+re0[region[i]]+sn0[season[i]];
              Phi[i] = phi  +reP[region[i]]+snP[season[i]];
              Delta[i] = delta+reD[region[i]]+snD[season[i]];
              mu[i] = Beta0[i] + //beta1 * (x[i]-Ph[i]) +
                      Delta[i] * theta * log1p(exp((x[i]-Phi[i])/theta));
            }
            for (i in 1:Ncens){
              B0cens[i] =beta0+re0[regcens[i]]+sn0[sncens[i]];
              Phcens[i] =phi  +reP[regcens[i]]+snP[sncens[i]];
              Decens[i] =delta+reD[regcens[i]]+snD[sncens[i]];
              mucens[i] =B0cens[i]+//beta1*(xcens[i]-Phcens[i]) + 
                         Decens[i]*theta*log1p(exp((xcens[i]-Phcens[i])/theta));
            }
          }
          model {
            phi ~ cauchy(0,1);
            beta0 ~ normal(0,5);
            delta ~ normal(0,5);
            sigma ~ cauchy(0, 1);
            sigmaR0 ~ cauchy(0, 1);
            sigmaRd ~ cauchy(0, 1);
            sigmaRp ~ cauchy(0, 1);
            sigmaS0 ~ cauchy(0, 1);
            sigmaSd ~ cauchy(0, 1);
            sigmaSp ~ cauchy(0, 1);

            zre0 ~ std_normal();
            zreD ~ std_normal();
            zreP ~ std_normal();
            zsn0 ~ std_normal();
            zsnD ~ std_normal();
            zsnP ~ std_normal();

            target += normal_lpdf(y | mu, sigma);
            target += normal_lcdf(ycens | mucens, sigma);
          }
          generated quantities {
            vector[N+Ncens] log_lik;
            real Ph;
            vector[R] delPhR;
            vector[S] delPhS;
            real B0;
            vector[R] delB0R;
            vector[S] delB0S;
            real De;
            vector[R] delDR;
            vector[S] delDS;
            Ph = phi + mean(reP[]) + mean(snP[]);
            B0 = beta0 + mean(re0[]) + mean(sn0[]);
            De = delta + mean(reD[]) + mean(snD[]);
            for (i in 1:R){
              delPhR[i] = reP[i] - mean(reP[]);
              delB0R[i] = re0[i] - mean(re0[]);
              delDR[i]  = reD[i] - mean(reD[]);
            }
            for (i in 1:S){
              delPhS[i] = snP[i] - mean(snP[]);
              delB0S[i] = sn0[i] - mean(sn0[]);
              delDS[i]  = snD[i] - mean(snD[]);
            }
            for (i in 1:N){
              log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
            }
            for (i in 1:Ncens){
              log_lik[N+i] = normal_lcdf(ycens[i] | mucens[i], sigma);
            }
          }
"

stan.in <- function(infile=eriedata, x = "POC", grp="Year", seasn="Month",
                    n.chains=nchains, nested=F){
  #infile$part_microcystin[infile$part_microcystin<0.01]<-0.01
    cens <- infile$part_microcystin<=0.1
    xtemp <- infile[,x]
    x <- log(xtemp[!cens])
    xcens <- log(xtemp[cens])
    y <- log(infile$part_microcystin[!cens])
    grtmp <- as.numeric(ordered(infile[, grp]))
    gr <- grtmp[!cens]
    grcens <- grtmp[cens]
    if (nested) sntmp <- as.numeric(ordered(paste(infile[, grp],
                                                  infile[, seasn],
                                                  sep=":")))
    else sntmp <- as.numeric(ordered(infile[, seasn]))
    sn <- sntmp[!cens]
    sncens <- sntmp[cens]
    n <- dim(infile)[1]-sum(cens)
    ncens <- sum(cens)
    ycens <- rep(log(0.1), ncens)
    R <- max(grtmp)
    S <- max(sntmp)
    xlimit <- log(range(xtemp, na.rm=T))

    inits <- list()
    bugs.data <- list(N=n, Ncens=ncens, R=R, S=S, y=y, x=x, region=gr,
                      season=sn, xcens=xcens, ycens=ycens, regcens=grcens,
                      sncens=sncens, theta=0.01*diff(xlimit),
                      phi_low=xlimit[1], phi_up=xlimit[2],
                      beta1=0)
  for (i in 1:n.chains)
    inits[[i]] <- list(beta0=rnorm(1, 0, 0.1), delta=runif(1,0,0.1),
                       phi=runif(1, xlimit[1], xlimit[2]/2),
                       zre0=rep(0,R), zreD=rep(0,R),
                       zreP=rep(0,R), zsn0=rep(0,S),
                       zsnD=rep(0,S), zsnP=rep(0,S),
                       sigma=runif(1),
                       sigmaR0=runif(1),sigmaRd=runif(1),
                       sigmaRp=runif(1),sigmaS0=runif(1),
                       sigmaSd=runif(1),sigmaSp=runif(1))

    para <- c("B0", "De", "Ph", "delB0R", "delDR", "delPhR", "delB0S",
              "delDS", "delPhS","sigma", "sigmaR0","sigmaRd","sigmaRp",
              "sigmaS0","sigmaSd","sigmaSp", "log_lik")
  return(list(para=para, data=bugs.data,
              inits=inits,n.chains=n.chains,
              theta=0.01*diff(xlimit)))
}

running_fun <- function(Data, stanfit=fit, Nst=F, ssn="Month"){
  input.to.bugs <- stan.in(infile=Data, n.chains=nchains,
                           x="Chla", nested=Nst, seasn=ssn)
  fit2keep <- sampling(stanfit, data = input.to.bugs$data,
                       init=input.to.bugs$inits,
                       pars = input.to.bugs$para,
                       iter=niters, thin=nthin,
                       chains=input.to.bugs$n.chains)#,
##                       control = list(adapt_delta = 0.99,
##                                      max_treedepth=15))
  return(list(fit=fit2keep, input=input.to.bugs))
}
