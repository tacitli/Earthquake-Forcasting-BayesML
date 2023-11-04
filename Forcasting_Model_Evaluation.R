rm(list=ls())
# ZIP INGARCH(1,1) time series data
setwd("~/desktop/Earthquake Research")
chain=read.csv("Source Data/chain_e5.csv")[,-1]
chain=chain[-(1:20000),]
a <- read.csv("Source Data/southwest_earthquake_eall.csv")
x <- a$e5

chain0=read.csv("Source Data/iteration60000.csv")
chain=chain0[-(1:15000),]
###
# Convergence diagnostics
library(coda)
geweke.diag(as.mcmc(chain))
geweke.plot(as.mcmc(chain))
z=NULL
for(i in 1:5){
  a=chain[1:10000,i]
  b=chain[10001:dim(chain)[1],i]
  ma=mean(a)
  mb=mean(b)
  sa=var(a)
  sb=var(b)
  na=length(a)
  nb=length(b)
  z=c(z,(ma-mb)/sqrt((1/na)*sa+(1/nb)*sb))
  z
}
z
###

alpha0 <- mean(chain[,1])
alpha1 <- mean(chain[,2])
beta1 <- mean(chain[,3])
rho <- mean(chain[,4])
lambda0 <- mean(chain[,5])
c(alpha0,alpha1,beta1,rho,lambda0)
# Lambda iteration function
nextlambda <- function(alpha0,alpha1,xt,beta1,lambdat){
  alpha0+alpha1*xt+beta1*lambdat
}
# ZIP distribution sampling function
rzip <- function(lambda,rho){
  if(rbinom(1,1,prob=rho)==1){return(0)}
  else{
    return(rpois(1,lambda/(1-rho)))
  }
}

# Mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# ZIP single point probability function
dzip <- function(x,lambda,rho){
  if(x==0){
    return(rho+(1-rho)*dpois(x,lambda/(1-rho)))
  }else{
    return((1-rho)*dpois(x,lambda/(1-rho)))
  }
}
# ZIP cumulative probability function
qzip <- function(x,lambda,rho){
  if(x<=rho){
    return(0)
  }else{
    return(qpois((x-rho)/(1-rho),lambda/(1-rho)))
  }
}

# Function to find the value with the highest probability in ZIP distribution
zipmode <- function(lambda,rho){
  md=NULL
  for (i in 0:100){
    md=c(md,dzip(i,lambda/(1-rho),rho))
  }
  return(which.max(md)-1)
}

# ZIP distribution median function
zipmedian <- function(lambda,rho){
  if(rho>=0.5){return(0)}else{
    return(qpois((0.5-rho)/(1-rho),lambda/(1-rho)))
  }
}

# ZIP distribution mean function
zipmean <- function(lambda,rho){
  return((1-rho)*lambda/(1-rho))
}

# Get lambda series
lambda=lambda0
for (i in 2:564){
  lambda=c(lambda,nextlambda(alpha0,alpha1,x[i-1],beta1,lambda[i-1]))  
}

setwd("/Users/lee/Desktop/")
d2017 <- read.csv("Source Data/level5_earthquake.csv")
e2017 = d2017$n[565:572] # Real number of earthquakes in 2017

# Set the value of lambda for the first month of 2017
lambda_2017=nextlambda(alpha0,alpha1,x[564],beta1,lambda[564])
# Get all lambda values for 2017 based on real values
for (i in 2:8){
  lambda_2017=c(lambda_2017,nextlambda(alpha0,alpha1,e2017[i-1],beta1,lambda_2017[i-1]))
}

#####################
# One-step forecast using the mean, rounded
fore_mean_2017=round(zipmean(lambda_2017[1],rho))
fore_meannoround_2017=zipmean(lambda_2017[1],rho)
for (i in 2:8){
  fore_mean_2017=c(fore_mean_2017,round(zipmean(lambda_2017[i],rho)))
  fore_meannoround_2017=c(fore_meannoround_2017,zipmean(lambda_2017[i],rho))
}
# One-step forecast using the median
fore_median_2017=zipmedian(lambda_2017[1],rho)

for (i in 2:8){
  fore_median_2017=c(fore_median_2017,zipmedian(lambda_2017[i],rho))
}

# One-step forecast using the mode
fore_mode_2017=zipmode(lambda_2017[1],rho)
for (i in 2:8){
  fore_mode_2017=c(fore_mode_2017,zipmode(lambda_2017[i],rho))
}

# Display results
one_step=cbind(e2017,fore_mean_2017,fore_meannoround_2017,fore_median_2017,fore_mode_2017)
write.csv(one_step,"Forecasting Results/one_step_forecast.csv")
# Real lambda values
lambda_2017

### Accurate interval forecast
# Upper 95%
range_2017=qzip(0.95,lambda_2017[1],rho)
for (i in 2:8){
  range_2017=c(range_2017,qzip(0.95,lambda_2017[i],rho))
}
range_2017
# Upper 97.5%
range2_2017=qzip(0.975,lambda_2017[1],rho)
for (i in 2:8){
  range2_2017=c(range2_2017,qzip(0.975,lambda_2017[i],rho))
}
range2_2017

output=cbind(e2017,fore_meannoround_2017,fore_mean_2017,fore_median_2017,fore_mode_2017,lambda_2017,range_2017,range2_2017)
write.csv(output,"Forecasting Results/lambda_star_2017_8_month_5_level_point_and_interval_forecast.csv")
write.csv(fore_meannoround_2017,"Forecasting Results/2017_8_month_5_level_mean_forecast.csv")
