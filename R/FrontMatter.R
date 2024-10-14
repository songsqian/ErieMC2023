packages<-function(x,
repos="https://cloud.r-project.org/", ...){
      x<-as.character(match.call()[[2]])
      if (!require(x,character.only=TRUE)){
          install.packages(pkgs=x, repos=repos, ...)
          require(x,character.only=TRUE)
      }
}

packages(arm)
##packages(dplyr)
packages(lattice)
##devtools::install_github("jsta/rv")
packages(rv)
packages(tikzDevice)
##packages(maptools)
packages(maps)
##packages(mapproj)
##packages(rpart)

base <- getwd()
dataDIR <- paste(base, "Data", sep="/")
plotDIR <- paste(base, "Figures", sep="/")
saveDIR <- paste(base, "Saved", sep="/")
packages(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))
rstan_options(threads_per_chain = 1)

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

rv0 <- function (length = 0)
{
    if (is.numeric(length)) {
        x <- as.list(rep.int(0, length))
    }
    else {
        stop("length must be numeric")
    }
    class(x) <- "rv"
    return(x)
}
to3 <- function(x){
    ifelse (x < 10, paste("00", x, sep=""),
            ifelse(x < 100, paste("0", x, sep=""),
                   as.character(x)))
}

to2 <- function(x)
    return(ifelse(x<10, paste("0", x, sep=""), as.character(x)))

hockey_smooth <- function(x, beta0, beta1=0, delta,
                          phi, theta=NULL){
    if (is.null(theta)) theta=0.01*diff(range(x))
    return(beta0 + beta1 * (x-phi) +
         delta * theta * log1p(exp((x-phi)/theta)))
}
