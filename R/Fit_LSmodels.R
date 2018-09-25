## ------------------------------------------------------------------- ##
##  R code to fit location-scale mortality models described in: 
##  Basellini U., Canudas-Romo V. and Lenart A. (2018), 
##  "Location-Scale Models in Demography: A Useful Re-parameterization
##  of Mortality Models", European Journal of Population 
##  
##  Author: Ugofilippo Basellini 
##  Last update: 10/08/2018
##  sessionInfo() details:
##  
##  R version 3.4.2 (2017-09-28)
##  Platform: x86_64-apple-darwin15.6.0 (64-bit)
##  Running under: macOS High Sierra 10.13.6
##  
##  locale: en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
##  attached base packages:
##  parallel  stats  graphics  grDevices  utils    
##  datasets  methods  base     
## 
##  other attached packages:
##  DEoptim_2.2-4  demography_1.21  forecast_8.4
##
## ------------------------------------------------------------------- ##

## clean the workspace
rm(list = ls())

## load useful packages
library("demography")
library("DEoptim")

## load Location-Scale functions
source("LSmodels.R")

## -- READING THE DATA -------------

## Read mortality data   
cou <- "DNK"   ## choose any country from HMD
username <- username   ## set your HMD credentials
password <- password   ## set your HMD credentials
FullData <- hmd.mx(country = cou, username, password)

## Specify age range, years and sex of analysis
ages <- x <- 30:110   
m <- length(x)
years <- y <- 1960:2016
n <- length(y)
sex <- "Female"    ## choose "Female", "Male" or "Both"

## Subset data
FittingData <- extract.years(FullData, years=y)
FittingData <- extract.ages(FittingData, ages=x)
if(sex=="Female"){
  ## females
  Exposures <- FittingData$pop$female
  MX.act <- FittingData$rate$female
  Deaths <- Exposures*MX.act
}else if(sex=="Male"){
  ## males
  Exposures <- FittingData$pop$male
  MX.act <- FittingData$rate$male
  Deaths <- Exposures*MX.act
}else if(sex=="Both"){
  ## both sexes
  Exposures <- FittingData$pop$total
  MX.act <- FittingData$rate$total
  Deaths <- Exposures*MX.act
}

## Choose model to fit 
## (type model.names to see possible choices)
model.names
LSmodel <- "MinGEV"

## Data frame to store LS parameters & matrix for fitted mx
LSparsHAT <- data.frame(Year=years,Model=LSmodel,Country=cou,Sex=sex,u = NA, c= NA, d=NA)
LSmxHAT <- matrix(NA,m,n)

## Estimate model's parameters for all years 
iter.max <- 100
for (i in 1:n){
  cat("Fitting year",y[i],"\n")
  ## Fitting
  fitLSmodel <- LSfit_fun(model=LSmodel,ages=x,Deaths=Deaths[,i],Exposures=Exposures[,i])
  ## Save parameters
  LSparsHAT[(LSparsHAT$Year==y[i]),]$u <- fitLSmodel[1]
  LSparsHAT[(LSparsHAT$Year==y[i]),]$c <- fitLSmodel[2]
  if(length(fitLSmodel)==3){
    LSparsHAT[(LSparsHAT$Year==y[i]),]$d <- fitLSmodel[3]
  }else if(length(fitLSmodel)==0){
    ## break loop if name provided incorrect
    break
  }
  ## save mx
  LSmxHAT[,i] <- LSmx_fun(model=LSmodel,ages=x,theta=fitLSmodel)
}

## plotting LS parameters
par(mfrow=c(1,2))
plot(y,LSparsHAT$u,xlab="Year",ylab="",main="Location, u",t="l",lwd=2)
plot(y,LSparsHAT$c,xlab="Year",ylab="",main="Scale, c",t="l",lwd=2)
par(mfrow=c(1,1))

## compute BIC
DXhat <- LSmxHAT*Exposures
DEVt1 <- Deaths * log(ifelse(Deaths==0, 1e-08, Deaths) / 
                        ifelse(DXhat==0, 1e-08, DXhat))
## second term of Deviance
DEVt2 <- Deaths - DXhat
## DEVIANCE
DEV <- 2 * sum(DEVt1 - DEVt2,na.rm = T)
## EFFECTIVE DIMENSION
ED <- length(LSparsHAT$u) + length(LSparsHAT$c) + length(LSparsHAT$d[!is.na(LSparsHAT$d)])
## BIC
BIC <- DEV + log(m*n)*ED
round(BIC/100)

## END