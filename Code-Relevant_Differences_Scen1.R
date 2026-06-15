#####################################################################
#### Code for implementing the Test for relevant differences ###
#####################################################################

### Code for intensive simulations of the level for various thresholds ###
### tuning parameter epsilon can be chosen manually, default is 0.1 ###
### this is the code used for simulating Scenario 1 in the manuscript ###

### Preparation: Code for W in order to obtain the quantile 

m <- 10000
n <- 10000
counter <- rep(0,m)
epsilon <- 0.1 # tuning parameter for W

for(j in 1:m)
{
  set.seed(j)
  #print(j)
  data <- rnorm(n,0,1/sqrt(n)) 
  bbridge <- cumsum(data) 
  numerator <- bbridge[n] 
  vectordenom <- seq(epsilon,1,by=0.1) 
  denominator <- sapply(vectordenom, function(x) bbridge[x*n]/x - bbridge[n]) 
  denominator <- sqrt(mean(denominator^2) * (1-epsilon))
  counter[j] <- numerator/denominator
}
quant <- unname(quantile(counter,0.95)) # 1.889786 for epsilon = 0.1 

#### Simulation: general Settings ###

Nsim <- 1000 # number of simulation runs 
n <- 1000 # sample size per group 
n1 <- n
n2 <- n
sd <- 0.5 # assumed standard deviation 

DeltaVec <- seq(0.00126,0.4,0.0025) #choice for simulating an entire range of thresholds

# true values are
#sd=1: 0.1354556
#sd=2: 0.06598765
#sd= .5: 0.2559087

# I = grid for y
minI <- -1
maxI <- 5
thresholds <- seq(minI,maxI,by=0.1)

# this is for saving some results
test_decision <- matrix(NA,nrow=Nsim,ncol=length(DeltaVec))
res1 <-matrix(NA,nrow=length(thresholds),ncol=3)
res2 <-matrix(NA,nrow=length(thresholds),ncol=3)
res1red <-matrix(NA,nrow=length(thresholds),ncol=3)
res2red <-matrix(NA,nrow=length(thresholds),ncol=3)
list1 <- list()
list2 <- list()
powerVec <- vector()
  
for(j in 1:Nsim){
  set.seed(12345 + j)
  
  # data generation
  # X, U independent distributed
  
  X1 <- rnorm(n1)
  U1 <- rnorm(n1,sd=sd)
  
  X2 <- rnorm(n2)
  U2 <- rnorm(n2,sd=sd)
  
  Y1 <- 1+X1+U1
  Y2 <- 1.7+X2+U2 
  
for(i in 1:length(thresholds)){
dependent1 <- ifelse(Y1 > thresholds[i],0,1)
dependent2 <- ifelse(Y2 > thresholds[i],0,1)

mod1 <- glm(dependent1 ~ X1, family=binomial(link="probit"))
mod2 <- glm(dependent2 ~ X2, family=binomial(link="probit"))

beta1hat <- mod1$coefficients
beta2hat <- mod2$coefficients

res1[i,] <- c(thresholds[i],beta1hat) # coefficients first model
res2[i,] <- c(thresholds[i],beta2hat) # coefficients 2nd model
}
  
# Preparations for V^2 
vectordenom <- seq(epsilon,1,by=0.1) # discretize the integral over [epsilon,1]
for(t in 1:length(vectordenom)){
  
  # "reduced" sample -> taking the first n_l*t observations
  X1red <- X1[1:floor(n1*vectordenom[t])]
  U1red <- U1[1:floor(n1*vectordenom[t])]
  
  X2red <- X2[1:floor(n2*vectordenom[t])]
  U2red <- U2[1:floor(n2*vectordenom[t])]
  
  Y1red <- Y1[1:floor(n1*vectordenom[t])]
  Y2red <- Y2[1:floor(n2*vectordenom[t])]
  
  for(i in 1:length(thresholds)){
    dependent1red <- ifelse(Y1red > thresholds[i],0,1)
    dependent2red <- ifelse(Y2red > thresholds[i],0,1)
    
    mod1red <- glm(dependent1red ~ X1red, family=binomial(link="probit"))
    mod2red <- glm(dependent2red ~ X2red, family=binomial(link="probit"))
    
    beta1hatred <- mod1red$coefficients
    beta2hatred <- mod2red$coefficients
    
    res1red[i,] <- c(thresholds[i],beta1hatred) # coefficients first model
    res2red[i,] <- c(thresholds[i],beta2hatred) # coefficients 2nd model
  }
  list1[[t]] <- res1red
  list2[[t]] <- res2red
}

######################
#########x_2=1########
######################

x <- 1 # we consider the conditional CDFs with x = 1

#true T
DeltaTq <- function(x){(pnorm(x-2, sd=sd)-pnorm(x-2.7, sd=sd))^2}
integrate(DeltaTq,minI,maxI) # this is the true underlying value

#Teststatistic T (only approximation possible, no smooth function)
DeltaThatq <- (pnorm(res1[,2]+res1[,3])-(pnorm(res2[,2]+res2[,3])))^2
df <- data.frame(x=thresholds,y=DeltaThatq)
That <- NA
try(That <- integrate(approxfun(df$x, df$y), min(df$x), max(df$x))$value)
if(is.na(That)){next}

# For obtaining Vhat, we need the estimates for the first n*t observations and calculate DeltaHat for them
# Vsum returns the 10 summands for approximating the integral over [epsilon,1]
Vsum <- vector()
for(t in 1:length(vectordenom)){
  dfred <- data.frame(x=thresholds,y=(pnorm(list1[[t]][,2]+list1[[t]][,3])-(pnorm(list2[[t]][,2]+list2[[t]][,3])))^2)
  try(Vsum[t] <- (integrate(approxfun(dfred$x, dfred$y), min(dfred$x), max(dfred$x))$value-That)^2)
}
if(sum(is.na(Vsum))>1){print("NAs")}
Vhat <- sqrt(mean(Vsum,na.rm=TRUE)) # this approximates the integral

### with That and Vhat we can now get a test decision ###

for (r in 1:length(DeltaVec)){
  
Delta_squared <- DeltaVec[r]
test_decision[j,r] <- ifelse(That > Delta_squared + quant*Vhat,1,0)
} 
}

# show the results
powerVec <- colMeans(test_decision,na.rm=TRUE)
powerVec

