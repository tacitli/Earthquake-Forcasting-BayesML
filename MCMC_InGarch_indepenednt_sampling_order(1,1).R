rm(list=ls())
library("truncnorm")
# ZIP INGARCH(1,1) model
# Time series data
setwd("~/desktop/Different Order Tests") # Setting working directory
a <- read.csv("level5.csv") # Reading the CSV file
x <- a$n[1:564] # Extracting the data from the CSV file

# Likelihood function, given parameters, return the log-likelihood value
likelihood <- function(x,param){
  alpha0 <- param[1]
  alpha1 <- param[2]
  beta1 <- param[3]
  rho <- param[4]
  lambda0 <- param[5]
  alpha2 <- param[6]
  lambda1 <- param[7]
  beta2 <- param[8]
  
  lambda <- c()
  lambda[1] <- lambda0
  lambda[2] <- lambda1
  if(x[1]==0){lnL<-log(rho+(1-rho)*exp(-lambda[1]/(1-rho)))}
  else{lnL<-(1-x[1])*log(1-rho)+x[1]*log(lambda[1])-log(factorial(x[1]))-lambda[1]/(1-rho)}
  if(x[2]==0){lnL<-lnL+log(rho+(1-rho)*exp(-lambda[2]/(1-rho)))}
  else{lnL<-lnL+(1-x[2])*log(1-rho)+x[2]*log(lambda[2])-log(factorial(x[2]))-lambda[2]/(1-rho)}
  for (j in 3:length(x)){
    lambda[j]<-alpha0+alpha1*x[j-1]+alpha2*x[j-2]+beta1*lambda[j-1]+beta2*lambda[j-2]
    if(x[j]==0){lnL<-lnL+log(rho+(1-rho)*exp(-lambda[j]/(1-rho)))}
    else{lnL<-lnL+(1-x[j])*log(1-rho)+x[j]*log(lambda[j])-log(factorial(x[j]))-lambda[j]/(1-rho)}
    if(is.nan(lnL)){print(j)}
  } 
  return(lnL)	
}

# Prior function, given parameters, return log-prior
prior <- function(param){
  alpha0 <- param[1]
  alpha1 <- param[2]
  beta1 <- param[3]
  rho<- param[4]
  lambda0<-param[5]
  alpha2<-param[6]
  lambda1<-param[7]
  beta2<-param[8]
  
  # Non-informative prior for alpha0
  prior_alpha1 <- dunif(alpha1, min=0, max=1, log = T)
  prior_beta1 <- dunif(beta1, min=0, max=1, log = T)#                
  prior_alpha2 <- dunif(alpha2, min=0, max=1, log = T)#
  prior_rho <- dunif(rho, min=0, max=1, log = T)
  prior_lambda0<-dgamma(lambda0,shape=1, rate=1.25,log = T)#Gamma（1,1.25）
  prior_lambda1<-dgamma(lambda1,shape=1, rate=1.25,log = T)#Gamma（1,1.25）
  prior_beta2 <- dunif(beta2, min=0, max=1, log = T)#     
  
  return(prior_alpha1+prior_beta1+prior_rho+prior_lambda0+prior_alpha2+prior_lambda1+prior_beta2)
}


# Posterior function, given parameters, return log-posterior
posterior <- function(x,param){
  return (likelihood(x,param) + prior(param))
}

######## Metropolis algorithm ################
#
# Definition of MH sampling for each parameter


MH_alpha0<-function(x,alpha0_old,alpha1,beta1,rho,lambda0,alpha2,lambda1,beta2){
  # Log form of the proposal density function value, posterior density function value
  #a0<--1
  #while(a0<0){a0<-rnorm(1,mean=0,sd=3)}
  #alpha0_new<-a0
  alpha0_new<-rtruncnorm(1, a=0, b=Inf, mean = alpha0_old, sd = 0.5)
  prop_old<-log(dtruncnorm(alpha0_old, a=0, b=Inf, mean = alpha0_new, sd = 0.5))
  prop_new<-log(dtruncnorm(alpha0_new, a=0, b=Inf, mean = alpha0_old, sd = 0.5))
  post_old<-posterior(x,c(alpha0_old,alpha1,beta1,rho,lambda0,alpha2,lambda1,beta2))
  post_new<-posterior(x,c(alpha0_new,alpha1,beta1,rho,lambda0,alpha2,lambda1,beta2))
  ratio <- post_new - post_old + prop_old - prop_new
  u <- runif(1,0,1)
  if (u < min(1, exp(ratio))){
    r <- alpha0_new
    #k_1<-k_1+1
  }else{r <- alpha0_old}
  r
}

MH_alpha1<-function(x,alpha0,alpha1_old,beta1,rho,lambda0,alpha2,lambda1,beta2){
  r<-alpha1_old
  alpha1_new<-rtruncnorm(1, a=0, b=Inf, mean = alpha1_old, sd = 0.1)
  if(alpha1_new>0){
    prop_old<-log(dtruncnorm(alpha1_old, a=0, b=Inf, mean = alpha1_new, sd = 0.1))
    prop_new<-log(dtruncnorm(alpha1_new, a=0, b=Inf, mean = alpha1_old, sd = 0.1))
    post_old<-posterior(x,c(alpha0,alpha1_old,beta1,rho,lambda0,alpha2,lambda1,beta2))
    post_new<-posterior(x,c(alpha0,alpha1_new,beta1,rho,lambda0,alpha2,lambda1,beta2))
    ratio <- post_new - post_old + prop_old - prop_new
    u <- runif(1,0,1)
    if (u < min(1, exp(ratio))){
      r <- alpha1_new
      #k_2<-k_2+1
    }
  }
  r
}
MH_beta1<-function(x,alpha0,alpha1,beta1_old,rho,lambda0,alpha2,lambda1,beta2){
  r<-beta1_old
  beta1_new<-rtruncnorm(1, a=0, b=1-alpha1, mean = beta1_old, sd = 0.2)
  if(beta1_new>0){
    prop_old<-log(dtruncnorm(beta1_old, a=0, b=1-alpha1, mean = beta1_new, sd = 0.2))
    prop_new<-log(dtruncnorm(beta1_new, a=0, b=1-alpha1, mean = beta1_old, sd = 0.2))
    post_old<-posterior(x,c(alpha0,alpha1,beta1_old,rho,lambda0,alpha2,lambda1,beta2))
    post_new<-posterior(x,c(alpha0,alpha1,beta1_new,rho,lambda0,alpha2,lambda1,beta2))
    ratio <- post_new - post_old + prop_old - prop_new
    u <- runif(1,0,1)
    if (u < min(1, exp(ratio))){
      r <- beta1_new
      #k_3<-k_1+3
    }
  }
  r
}
MH_rho<-function(x,alpha0,alpha1,beta1,rho_old,lambda0,alpha2,lambda1,beta2){
  rho_new<-rtruncnorm(1, a=0, b=Inf, mean = rho_old, sd = 0.1)
  prop_old<-log(dtruncnorm(rho_old, a=0, b=Inf, mean = rho_new, sd = 0.1))
  prop_new<-log(dtruncnorm(rho_new, a=0, b=Inf, mean = rho_old, sd = 0.1))
  post_old<-posterior(x,c(alpha0,alpha1,beta1,rho_old,lambda0,alpha2,lambda1,beta2))
  post_new<-posterior(x,c(alpha0,alpha1,beta1,rho_new,lambda0,alpha2,lambda1,beta2))
  ratio <- post_new - post_old + prop_old - prop_new
  u <- runif(1,0,1)
  if (u < min(1, exp(ratio))){
    r <- rho_new
    #k_4<-k_4+1
  }else{r<-rho_old}
  r
}

MH_lambda0<-function(x,alpha0,alpha1,beta1,rho,lambda0_old,alpha2,lambda1,beta2){
  lambda0_new<-rtruncnorm(1, a=0, b=Inf, mean = lambda0_old, sd = 1)
  prop_old<-log(dtruncnorm(lambda0_old, a=0, b=Inf, mean = lambda0_new, sd = 1))
  prop_new<-log(dtruncnorm(lambda0_new, a=0, b=Inf, mean = lambda0_old, sd = 1))
  post_old<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0_old,alpha2,lambda1,beta2))
  post_new<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0_new,alpha2,lambda1,beta2))
  ratio <- post_new - post_old + prop_old - prop_new
  u <- runif(1,0,1)
  if (u < min(1, exp(ratio))){
    r <- lambda0_new
    #k_5<-k_5+1
  }else{r<-lambda0_old}
  r
}

MH_alpha2<-function(x,alpha0,alpha1,beta1,rho,lambda0,alpha2_old,lambda1,beta2){
  r<-alpha2_old
  alpha2_new<-rtruncnorm(1, a=0, b=Inf, mean = alpha2_old, sd = 0.1)
  if(alpha2_new>0){
    prop_old<-log(dtruncnorm(alpha2_old, a=0, b=1-alpha1-beta1, mean = alpha2_new, sd = 0.1))
    prop_new<-log(dtruncnorm(alpha2_new, a=0, b=1-alpha1-beta1, mean = alpha2_old, sd = 0.1))
    post_old<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0,alpha2_old,lambda1,beta2))
    post_new<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0,alpha2_new,lambda1,beta2))
    ratio <- post_new - post_old + prop_old - prop_new
    u <- runif(1,0,1)
    if (u < min(1, exp(ratio))){
      r <- alpha2_new
      #k_2<-k_2+1
    }
  }
  r
}

MH_lambda1<-function(x,alpha0,alpha1,beta1,rho,lambda0,alpha2,lambda1_old,beta2){
  lambda1_new<-rtruncnorm(1, a=0, b=Inf, mean = lambda1_old, sd = 1)
  prop_old<-log(dtruncnorm(lambda1_old, a=0, b=Inf, mean = lambda1_new, sd = 1))
  prop_new<-log(dtruncnorm(lambda1_new, a=0, b=Inf, mean = lambda1_old, sd = 1))
  post_old<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0,alpha2,lambda1_old,beta2))
  post_new<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0,alpha2,lambda1_new,beta2))
  ratio <- post_new - post_old + prop_old - prop_new
  u <- runif(1,0,1)
  if (u < min(1, exp(ratio))){
    r <- lambda1_new
    #k_5<-k_5+1
  }else{r<-lambda1_old}
  r
}

MH_beta2<-function(x,alpha0,alpha1,beta1,rho,lambda0,alpha2,lambda1,beta2_old){
  r<-beta2_old
  beta2_new<-rtruncnorm(1, a=0, b=1-alpha1-beta1-alpha2, mean = beta2_old, sd = 0.2)
  if(beta2_new>0){
    prop_old<-log(dtruncnorm(beta2_old, a=0, b=1-alpha1-beta1-alpha2, mean = beta2_new, sd = 0.2))
    prop_new<-log(dtruncnorm(beta2_new, a=0, b=1-alpha1-beta1-alpha2, mean = beta2_old, sd = 0.2))
    post_old<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0,alpha2,lambda1,beta2_old))
    post_new<-posterior(x,c(alpha0,alpha1,beta1,rho,lambda0,alpha2,lambda1,beta2_new))
    ratio <- post_new - post_old + prop_old - prop_new
    u <- runif(1,0,1)
    if (u < min(1, exp(ratio))){
      r <- beta2_new
      #k_3<-k_1+3
    }
  }
  r
}


# define block function
block<-function(x,param){
  alpha0_new<-MH_alpha0(x,param[1],param[2],param[3],param[4],param[5],param[6],param[7],param[8])
  alpha1_new<-MH_alpha1(x,alpha0_new,param[2],param[3],param[4],param[5],param[6],param[7],param[8])
  beta1_new<-MH_beta1(x,alpha0_new,alpha1_new,param[3],param[4],param[5],param[6],param[7],param[8])
  rho_new<-MH_rho(x,alpha0_new,alpha1_new,beta1_new,param[4],param[5],param[6],param[7],param[8])
  lambda0_new<-MH_lambda0(x,alpha0_new,alpha1_new,beta1_new,rho_new,param[5],param[6],param[7],param[8])
  alpha2_new<-MH_alpha2(x,alpha0_new,alpha1_new,beta1_new,rho_new,lambda0_new,param[6],param[7],param[8])
  lambda1_new<-MH_lambda1(x,alpha0_new,alpha1_new,beta1_new,rho_new,lambda0_new,alpha2_new,param[7],param[8])
  beta2_new<-MH_beta2(x,alpha0_new,alpha1_new,beta1_new,rho_new,lambda0_new,alpha2_new,lambda1_new,param[8])
  param_new<-c(alpha0_new,alpha1_new,beta1_new,rho_new,lambda0_new,alpha2_new,lambda1_new,beta2_new)
  param_new
  #if(alpha1_new+beta1_new<1){param_new}else{param}
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    chain[i+1,] = block(x,c(chain[i,1],chain[i,2],chain[i,3],chain[i,4],chain[i,5],chain[i,6],chain[i,7],chain[i,8]))#
  }
  return(chain)
}

################a0     a1 b1   p  l0   a2  l1  b2
startvalue = c(0.385, 0.05, 0.5, 0.6, 2.5, 0.2, 0.8, 0.1)
chain = run_metropolis_MCMC(startvalue, 20000) # Run the Metropolis-Hastings MCMC algorithm
write.csv(chain, "Iteration20000_without_annealing.csv")  # Save the results

chain0 = read.csv("Iteration20000_without_annealing.csv")[, -1]

# Annealing
chain = chain0[-(1:6001), ]
colnames(chain) <- c("alpha0", "alpha1", "beta1", "rho", "lambda1", "alpha2", "lambda2", "beta2")
# Convergence diagnostics
library(coda)
geweke.diag(as.mcmc(chain))
geweke.plot(as.mcmc(chain))
cumuplot(as.mcmc(chain))
z = NULL
for (i in 1:7) {
  a = chain[1:3500, i]
  b = chain[3501:dim(chain)[1], i]
  ma = mean(a)
  mb = mean(b)
  sa = var(a)
  sb = var(b)
  na = length(a)
  nb = length(b)
  z = c(z, (ma - mb) / sqrt((1 / na) * sa + (1 / nb) * sb))
}
z

alpha0 <- mean(chain[, 1])
alpha1 <- mean(chain[, 2])
beta1 <- mean(chain[, 3])
rho <- mean(chain[, 4])
lambda0 <- mean(chain[, 5])
alpha2 <- mean(chain[, 6])
lambda1 <- mean(chain[, 7])
beta2 <- mean(chain[, 8])
para = c(alpha0, alpha1, beta1, rho, lambda0, alpha2, lambda1, beta2)
write.csv(para, "parameters_2_2.csv")
# Calculate AIC and BIC
k = 6
AIC = 2 * k - 2 * likelihood(x, para); AIC
BIC = log(564) * k - 2 * likelihood(x, para); BIC

# Function for iteratively computing lambda
nextlambda <- function(alpha0, alpha1, alpha2, xt, xt2, beta1, beta2, lambda1, lambda2) {
  alpha0 + alpha1 * xt + alpha1 * xt2 + beta1 * lambda1 + beta2 * lambda2
}

# Obtain the historical series of lambda
lambda = c(lambda0, lambda1)
for (i in 3:564) {
  lambda = c(lambda, nextlambda(alpha0, alpha1, alpha2, x[i - 1], x[i - 2], beta1, beta2, lambda[i - 1], lambda[i - 2]))
}
d2017 <- read.csv("five_level.csv")
e2017 = d2017$n[565:572] # Real number of earthquakes in 2017

# Set the lambda values for the first and second month of 2017
lambda_2017 = nextlambda(alpha0, alpha1, alpha2, x[564], x[563], beta1, beta2, lambda[564], lambda[563])
lambda_2017 = c(lambda_2017, nextlambda(alpha0, alpha1, alpha2, e2017[1], x[564], beta1, beta2, lambda_2017[1], lambda[564]))
# Obtain all lambda values for 2017 based on real data
for (i in 3:8) {
  lambda_2017 = c(lambda_2017, nextlambda(alpha0, alpha1, alpha2, e2017[i - 1], e2017[i - 2], beta1, beta2, lambda_2017[i - 1], lambda_2017[i - 2]))
}

# ZIP distribution sampling function
rzip <- function(lambda, rho) {
  if (rbinom(1, 1, prob = rho) == 1) {
    return(0)
  } else {
    return(rpois(1, lambda / (1 - rho)))
  }
}

# Mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# ZIP distribution single point probability function
dzip <- function(x, lambda, rho) {
  if (x == 0) {
    return(rho + (1 - rho) * dpois(x, lambda / (1 - rho)))
  } else {
    return((1 - rho) * dpois(x, lambda / (1 - rho)))
  }
}
# ZIP distribution cumulative probability function
qzip <- function(x, lambda, rho) {
  if (x <= rho) {
    return(0)
  } else {
    return(qpois((x - rho) / (1 - rho), lambda / (1 - rho)))
  }
}

# Function for the most probable value of the ZIP distribution
zipmode <- function(lambda, rho) {
  md = NULL
  for (i in 0:100) {
    md = c(md, dzip(i, lambda / (1 - rho), rho))
  }
  return(which.max(md) - 1)
}

# ZIP distribution median function
zipmedian <- function(lambda, rho) {
  if (rho >= 0.5) {
    return(0)
  } else {
    return(qpois((0.5 - rho) / (1 - rho), lambda / (1 - rho)))
  }
}

# ZIP distribution mean function
zipmean <- function(lambda, rho) {
  return((1 - rho) * lambda / (1 - rho))
}

#####################
# One-step forecast using the mean, rounded
fore_mean_2017 = round(zipmean(lambda_2017[1], rho))
fore_meannoround_2017 = zipmean(lambda_2017[1], rho)
for (i in 2:8) {
  fore_mean_2017 = c(fore_mean_2017, round(zipmean(lambda_2017[i], rho)))
  fore_meannoround_2017 = c(fore_meannoround_2017, zipmean(lambda_2017[i], rho))
}
# One-step forecast using the median
fore_median_2017 = zipmedian(lambda_2017[1], rho)

for (i in 2:8) {
  fore_median_2017 = c(fore_median_2017, zipmedian(lambda_2017[i], rho))
}

# One-step forecast using the mode
fore_mode_2017 = zipmode(lambda_2017[1], rho)
for (i in 2:8) {
  fore_mode_2017 = c(fore_mode_2017, zipmode(lambda_2017[i], rho))
}

### Precise interval forecasting
# Upper 95%
range_2017 = qzip(0.95, lambda_2017[1], rho)
for (i in 2:8) {
  range_2017 = c(range_2017, qzip(0.95, lambda_2017[i], rho))
}
range_2017
# Upper 97.5%
range2_2017 = qzip(0.975, lambda_2017[1], rho)
for (i in 2:8) {
  range2_2017 = c(range2_2017, qzip(0.975, lambda_2017[i], rho))
}
range2_2017

yibu = cbind(e2017, fore_meannoround_2017, fore_mean_2017, fore_median_2017, fore_mode_2017, lambda_2017, p95 = range_2017, p975 = range2_2017)
write.csv(yibu, "point_and_interval_forecast_2_2.csv")

mean(abs(fore_meannoround_2017 - e2017))
mean(abs(fore_mean_2017 - e2017))
mean(abs(fore_median_2017 - e2017))
mean(abs(fore_mode_2017 - e2017))

mean(abs(fore_meannoround_2017 - e2017) / (e2017 + 1))
mean(abs(fore_mean_2017 - e2017) / (e2017 + 1))
mean(abs(fore_median_2017 - e2017) / (e2017 + 1))
mean(abs(fore_mode_2017 - e2017) / (e2017 + 1))






