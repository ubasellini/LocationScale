## ------------------------------------------------------------------- ##
##  R functions accompanying the file "FittingLSmodels.R".
##  Author: Ugofilippo Basellini
##  Last Update: 09/08/2018
## ------------------------------------------------------------------- ##

#################################################################---
## LS, LLS and LS-like parameterization of 12 mortality models  ----
#################################################################---

## Logistic model
LogisticMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## force of mortality
  num <- exp(y)
  den <- 1+exp(y)
  out <- (1/c)*num/den
  return(out)
}

## Normal model
NormalMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## Normal probability density function
  fx <- dnorm(x=x,mean=u,sd=c)
  ## Normal cumulative distribution function
  Fx <- pnorm(q=x,mean=u,sd=c) 
  ## force of mortality
  num <- fx
  den <- 1 - Fx
  out <- num/den
  return(out)
}

## Largest Extreme-Value model
LargestExVaMu_LS <- function(ages,theta) {
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## force of mortality
  num <- exp(-y)
  den <- exp(exp(-y))-1
  out <- (1/c)*num/den
  return(out)
}

## Weibull model
WeibullMu_LLS <- function(ages,theta){
  ## LLS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LLS-standardize 
  x <- ages
  y <- (log(x)-u)/c
  ## force of mortality
  t1 <- 1/(c*x)
  t2 <- exp(y)
  out <- t1*t2
  return(out)
}

## Log-Logistic model
LogLogisticMu_LLS <- function(ages,theta){
  ## LLS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LLS-standardize 
  x <- ages
  y <- (log(x)-u)/c
  ## force of mortality
  t1 <- 1/(c*x)
  num <- exp(y)
  den <- 1+exp(y)
  out <- t1*num/den
  return(out)
}

## Log-Normal model
LogNormalMu_LLS <- function(ages,theta){
  ## LLS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LLS-standardize  
  x <- ages
  y <- (log(x)-u)/c
  ## Normal probability density function
  fx <- dnorm(x=log(x),mean=u,sd=c)
  ## Normal cumulative distribution function
  Fx <- pnorm(q=log(x),mean=u,sd=c) 
  ## force of mortality
  num <- fx
  den <- 1 - Fx
  out <- (1/x)*num/den
  return(out)
}

## Gompertz model (= Smallest Extreme-Value)
GompertzMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## force of mortality
  t1 <- 1/c
  t2 <- exp(y)
  out <- t1*t2
  return(out)
}

## Gamma-Gompertz model
GammaGompertzMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  gam <- theta[3]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## force of mortality
  t1 <- 1/c
  num <- exp(y)
  den <- 1+gam*(exp(y)-exp(-u/c))
  out <- t1*num/den
  return(out)
}

## Kannisto model
KannistoMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## force of mortality
  num <- exp(y)
  den <- 1+exp(y)
  out <- num/den
  return(out)
}

## Minimal Generalized Extreme-Value
MinGEVMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  d <- theta[3]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## force of mortality
  t1 <- (1/c)
  t2 <- (1+d*(-y))^(-1/d-1)
  out <- t1*t2
  return(out)
}

## Maximal Generalized Extreme-Value
MaxGEVMu_LS <- function(ages,theta){
  ## LS parameters
  u <- theta[1]
  c <- theta[2]
  d <- theta[3]
  ## redefine age range and LS-standardize  
  x <- ages
  y <- (x-u)/c
  ## transformation of y
  s <- 1+d*y
  ## force of mortality
  t1 <- (1/c)
  t2 <- exp(-s^(-1/d))
  t3 <- s^(-1/d-1)
  den <- 1-exp(-s^(-1/d))
  out <- t1*t2*t3/den
  return(out)
}


#################################################################---
## Log-likelihood and estimation function for LS parameters  -------
#################################################################---

model.names <- c("Logistic","Normal","LargeEV","Weibull",
                 "LogLogis","LogNormal","Gompertz","GammaGomp","Kannisto",
                 "MinGEV","MaxGEV")

## LS mx according to model choice
LSmx_fun <- function(model,ages,theta){
  mx <- switch(model, Logistic={LogisticMu_LS(ages = ages, theta = theta)},
               Normal={NormalMu_LS(ages = ages, theta = theta)},
               LargeEV={LargestExVaMu_LS(ages = ages, theta = theta)},
               Weibull={WeibullMu_LLS(ages = ages, theta = theta)},
               LogLogis={LogLogisticMu_LLS(ages = ages, theta = theta)},
               LogNormal={LogNormalMu_LLS(ages = ages, theta = theta)},
               Gompertz={GompertzMu_LS(ages = ages, theta = theta)},
               GammaGomp={GammaGompertzMu_LS(ages = ages, theta = theta)},
               Kannisto={KannistoMu_LS(ages = ages, theta = theta)},
               MinGEV={MinGEVMu_LS(ages = ages, theta = theta)},
               MaxGEV={MaxGEVMu_LS(ages = ages, theta = theta)})
  return(mx)
}

## LS Log-likelihood 
Loglike_LS <- function(theta,model,ages,Deaths,Exposures){
  mx <- LSmx_fun(model = model, ages = ages, theta = theta)
  out <- - sum(Deaths * log(mx) - mx * Exposures,na.rm = T)
  return(out)
}

## general function to estimate the LS parameters of the model of interest
LSfit_fun <- function(model,ages,Deaths,Exposures,iter.max=500){
  ## check that model is within the acceptable ones
  if(!model %in% model.names){ 
    par.hat <- cat("Model does not belong to LS family, please choose a model among:",
                   model.names,"")
  }else{
    theta_start <- switch(model, Logistic={c(80,9)},Normal={c(80,10)},
                          LargeEV={c(80,10)},Gompertz={c(80,10)},
                          Kannisto={c(100,10)})
    ## subset of models to be fitted with DEoptim
    DEopt_models <- c("Weibull","LogLogis", "LogNormal","GammaGomp",
                      "MinGEV","MaxGEV")
    ## optimize log likelihood
    if(model %in% DEopt_models){
      ## with DEoptim
      require(DEoptim)
      theta_low <- switch(model,Weibull={c(3,0.08)}, LogLogis={c(3,0.08)}, 
                          LogNormal={c(3,0.08)},GammaGomp={c(50,6,-0.4)},
                          MinGEV={c(60,6,-0.1)},MaxGEV={c(60,6,-0.4)})
      theta_up <- switch(model, Weibull={c(5,0.15)}, LogLogis={c(5,0.15)}, 
                         LogNormal={c(5,0.15)},GammaGomp={c(110,20,0.4)},
                         MinGEV={c(100,30,0.2)},MaxGEV={c(100,30,0.2)})
      opt <- DEoptim(fn=Loglike_LS,lower = theta_low,upper=theta_up,
                     model=model, ages=ages, Deaths=Deaths,Exposures=Exposures,
                     control = DEoptim.control(trace = F,itermax = iter.max))
      par.hat <- opt$optim$bestmem
    }else{
      ## with optim
      opt <- optim(par=theta_start, fn=Loglike_LS, model=model, ages=ages, 
                   Deaths=Deaths, Exposures=Exposures)
      par.hat <- opt$par
    }
  }
  return(par.hat)
}

