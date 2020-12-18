#This is the r file used in the project.It used a dataset that is included in the zip-file,
#and some c++ functions that also are inlcuded in the zip. 
library(TMB)
library(glmmTMB)
data <- read.delim("data032010.txt") # dataen med elgen

#all functions I made in the project are included below

#some of the data-points are missing. I therefore made a function that 
#only includes the moose were all the covariats of interest are present
#factor is set to be region, year, age, ovulation, weight and day of kill.
#but it can be specified by user to be something else
initialize_data <- function(factors=c(10,17,2,5,6,20),data1){
  df = data1
  keep = rep(TRUE,dim(data1)[1])
  for(i in 1:dim(data1)[1])
    if (sum(is.na(as.numeric(df[i,][factors])))>0||((as.numeric(df[i,][17])==0))){
      keep[i]= FALSE
    }
  df = df[keep,]
  return(df)
}
#this is just a function that calculates AIC given the negative log-likelihood and the number of fixed parameters
my_aic <- function(nll,par){
  return(2*nll+2*par)
}
# mu is a vector of estimated expected ovulation time
# sigma is a vector of estimated standard diviations
# q is a vector estimated ovulation probabilities
# TID is the acctual time of kill for each moose
# the functions take the estimated values for each moose and returnes a value 0 or 1 for their ovulation based on the fitted model
# these values are then meant to be used as responce values when performing bootstrapping
simulation <- function(mu,sigma,q,TID){
  p = q*pnorm(TID,mu,sigma)
  y = rbinom(length(p),1,p)
  return(y)
}
# This model performes B bootsrapps. And returns a matrix where each collumn is the estimated parameter values for each bootstrap run
# mu, sigma, q and TID are as described above
# dataB is the data needed to run the c++ file
# parametersB are the parameters to be estimated
# randomB indicates what parameters are included as random effects
# mapB is a list of what random effeccts are to be ignored
# file is a string giving what c++ file is to be used
bootstrapping <- function(B,mu,sigma,q,TID,dataB,parametersB,randomB,mapB,file){
  boot = cbind()
  for(i in 1:B){
    Y = simulation(mu,sigma,q,TID)
    dataB[[1]] = Y
    objB <- MakeADFun(dataB,parametersB,DLL = file,map=mapB,random = randomB, silent = TRUE)
    resB <- optim(par = objB$par, fn = objB$fn, gr = objB$gr, method = "BFGS")
    repB = sdreport(objB)
    boot = cbind(boot,c(repB$par.fixed,repB$par.random))
  }
  return(boot)
}

# this function is very specific. It is meant to plot the model and observed values for a model with age and region as only effects
# it also plots a 95% confidence interval based on bootstrapping
plotting_regionAndAge <- function(beta,df,alder,region,boot){
  t = c(243:317)
  mu = beta[1]+beta[5+region-4100]
  sigma = exp(beta[2])
  q = plogis(beta[5+8+region-4100]+beta[5+24+alder[1]])
  if(length(alder)==1){
    main = paste("Age=",alder[1], ", Region=",locations[region-4100],", q=",round(q,2))
  }
  else{
    main = paste("Age=",alder[1],":",alder[length(alder)],", Region=",locations[region-4100],", q=",round(q,2))
  }
  plot(t,q*pnorm(t,mu,sigma),
       type = "l",ylim = c(0,1),main=main,
       xlab = "Days into the year",ylab= "probability of having ovulated",lty=1)
  rates = rep(0,length(min(df$daynr):max(df$daynr)))
  size = rep(0,length(min(df$daynr):max(df$daynr)))
  obs =rep(0,length(min(df$daynr):max(df$daynr)))
  for(i in t){
    obs[i-min(df$daynr)+1] = sum(df$ovul[which(df$daynr==i&df$age%in%alder&df$region==region)])
    rates[i-min(df$daynr)+1] = sum(df$ovul[which(df$daynr==i&df$age%in%alder&df$region==region)])/length(df$ovul[which(df$daynr==i&df$age%in%alder&df$region==region)])
    size[i-min(df$daynr)+1] = length(df$ovul[which(df$daynr==i&df$age%in%alder&df$region==region)])
    
  }
  #conf int with bootstrapping 
  
  boot_conf = matrix(0,ncol=length(t),nrow=dim(boot)[2])
  
  for( b in 1:dim(boot_conf)[1]){
    muB = boot[1,b] + boot[5+region-4100,b]
    log_sigmaB = boot[2,b]
    logit_qB =boot[5+8-4100+region,b]+boot[5+24+alder[1],b]
    boot_conf[b,] = plogis(logit_qB)*pnorm(t,muB,exp(log_sigmaB))
  }
  quantiles = matrix(0,nrow=length(t),ncol=2)
  for(tid in 1:(length(t))){
    quantiles[tid,] = quantile(sort(boot_conf[,tid]),c(0.025,0.975))
  }
  #moving average is made here
  points(min(df$daynr):max(df$daynr),rates,cex = 0.25*sqrt(size),pch=19)
  keep =which(!is.nan(rates))
  t = min(df$daynr):max(df$daynr)
  pred = predict(loess(rates[keep]~keep,weights = size[keep]),se=TRUE)
  lines(t[keep],pred$fit+ qt(0.975,pred$df)*pred$se.fit,col="blue")
  lines(t[keep],pred$fit-qt(0.975,pred$df)*pred$se.fit,col="blue")
  lines(t,quantiles[,1],col="red",lty=2)
  lines(t,quantiles[,2],col = "red",lty=2)
  legend = c("Estimated probability","Bootstrapped conf interval","Observed ovulation rate","95% conf.int for ovulationrate")
  if(1 %in% alder){
    legend("topleft",legend,col=c("black","red","black","blue"),lty=c(1,2,NA,1),pch=c(NA,NA,19,NA))
  }
  
  
}
#########################################################################################
###############                                                     #####################
###############         HERE I START MAKING THE MODELS              #####################
###############                                                     #####################
#########################################################################################


#here I start to make the models, the first part makes the model as simple as possible
df = initialize_data(data1= data)#just inizialize the data
df$cweight = as.numeric(df$cweight)
#the first example uses glmTMB 

#first I make the glm with probit link function
glm <- glmmTMB(ovul~1+daynr,data=df,family = binomial(link=probit))
glm$fit$par
#the estimated values for mu and sigma as described in the project-paper
mu = - glm$fit$par[[1]]/glm$fit$par[[2]]
sigma = 1/glm$fit$par[[2]]
#now I also include a constant zero-inflation. z is modeld with a logit link-function
ziglm <- glmmTMB(ovul~1+daynr,data=df,family = binomial(link=probit),ziformula = ~1,
                 start = list(beta=c(-12,0.05),betazi=0.6))
#the estimated biologicaly interesting values
muzi = - ziglm$fit$par[[1]]/ziglm$fit$par[[2]]
sigmazi = 1/ziglm$fit$par[[2]]
z = plogis(ziglm$fit$par[[3]])

#########################################################################################################
#Here I begin using TMB and code written in c++
#########################################################################################################
#The simplest model, but I use TMB and code from c++.
compile("simpleModel.cpp") # this is the c++ file I have written. It calculates the nll for this model
dyn.load(dynlib("simpleModel"))
Y=df$ovul # the responce in the bernulli-function
N = rep(1,length(Y)) 
TID = df$daynr
###the model

data <- list(Y = Y,N = N,TID = TID)
parameters = list(mu = 280, log_sig = log(14))
obj <- MakeADFun(data,parameters,DLL = "simpleModel", silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")#optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
mu = rep$par.fixed[1] # the estimate for mu
sigma = exp(rep$par.fixed[2]) #the estimate for sigma

#a simple model when q does not need to be 1.
compile("simpleModelq.cpp") # this is the c++ file I have written. It calculates the nll for this model
dyn.load(dynlib("simpleModelq"))
Y=df$ovul # the responce in the bernulli-function
N = rep(1,length(Y)) 
TID = df$daynr # the kill date for all the moose

data <- list(Y = Y,N = N,TID = TID)#the data the c++ code depends on
parameters = list(mu = 280, log_sig = log(14),q = 0) # starting values for the parameters
obj <- MakeADFun(data,parameters,DLL = "simpleModelq", silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
repq = sdreport(obj)
muq = repq$par.fixed[1] # the estimate for mu
sigmaq = exp(repq$par.fixed[2]) # the estimate for sigma
q = plogis(repq$par.fixed[3]) # the estimate for q

#The models are fitted and those values will be used to make some plots 
#the next lines make the observed rate of ovulation for each day, and how many observation each day has
rates = rep(0,length(min(df$daynr):max(df$daynr)))
size = rep(0,length(min(df$daynr):max(df$daynr)))
obs =rep(0,length(min(df$daynr):max(df$daynr)))
for(i in min(df$daynr):max(df$daynr) ){
  obs[i-min(df$daynr)+1] = sum(df$ovul[which(df$daynr==i)])
  rates[i-min(df$daynr)+1] = sum(df$ovul[which(df$daynr==i)])/length(df$ovul[which(df$daynr==i)])
  size[i-min(df$daynr)+1] = length(df$ovul[which(df$daynr==i)])
  
}
t = c(min(df$daynr):max(df$daynr))
#plot for q=1
plot(t,pnorm(t,mu,sigma),type="l",xlab="Days into the year",ylab ="Probability of ovulation",main="Simple model without q",ylim=c(0,1.2),col="red")
points(t,rates,cex=0.1*sqrt(size),pch=19)
pred = predict(loess(rates[keep]~keep,weights = size[keep]),se=TRUE)
lines(t[keep],pred$fit+1.96*pred$se.fit,col="blue",lty=2)
lines(t[keep],pred$fit-1.96*pred$se.fit,col="blue",lty=2)
legend = c(paste("The model, mu =",round(mu,1),", sd=",round(sigma,1)), "The observed ovulation rate","95% conf. int. for observed mean")
legend("topright",legend = legend,col=c("red","black","blue"),pch=c(NA,19,NA),lty=c(1,NA,2))

#this is the plot for q<1
plot(t,q*pnorm(t,muq,sigmaq),type="l",xlab="Days into the year",ylab ="Probability of ovulation",main="Simple model with q",ylim=c(0,1.2),col="red")
points(t,rates,cex=0.1*sqrt(size),pch=19)
pred = predict(loess(rates[keep]~keep,weights = size[keep]),se=TRUE)
lines(t[keep],pred$fit+1.96*pred$se.fit,col="blue",lty=2)
lines(t[keep],pred$fit-1.96*pred$se.fit,col="blue",lty=2)
legend = c(paste("The model, mu =",round(mu,1),", sd=",round(sigma,1),", q=",round(q,2)),"The observed ovulation rate","95% conf. int. for observed mean")
legend("topright",legend = legend,col=c("red","black","blue"),pch=c(NA,19,NA),lty=c(1,NA,2))




##########################################################################################################
###################  A model based only on region will be written about now      #########################
##########################################################################################################
compile("simpleModelRegion.cpp")
dyn.load(dynlib("simpleModelRegion"))
#the data that is included in the c++ file

Y=df$ovul # the responce in the bernulli-function
N = rep(1,length(Y)) 
TID = df$daynr # the kill date for all the moose
id = factor(df$region)
XM = cbind(N)
XS = cbind(N)
XQ = cbind(N)
data = list(Y=Y,TID=TID,N=N,XM=XM,XS=XS,XQ=XQ,idm=id,ids=id,idq=id)
parameters = list(BetaM = c(275),BetaS = c(log(10)),BetaQ = c(0),
                  um = rep(0,nlevels(id)),us = rep(0,nlevels(id)),uq = rep(0,nlevels(id)),
                  log_sigmaM = log(0.5), log_sigmaS = log(0.5), log_sigmaQ = log(0.5) )
random = c("um","us","uq")
map = list(log_sigmaS=factor(NA))
obj <- MakeADFun(data,parameters,DLL = "simpleModelRegion",random=random,map=map ,silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
my_aic(res$value,length(rep$par.fixed)) # testing AIC, one needs to run multiple models to compare
#better result without random effect on sigma, so that is the model I choose


#now the model is chosen and the plots will be made
#assumed form of model
# mu ~ 1 + region
# sigma ~ 1
# q ~ 1 + region 
t = c(min(df$daynr):max(df$daynr))
locations = c("Troms","Nordland","Nord-Trondelag","Oppland","Hedmark","Telemark","Aust-Agder","Vest-Agder")
mu = rep$par.fixed[1] + rep$par.random[1:8]#there are 8 locations
sigma = exp(rep$par.fixed[2])
q = plogis(rep$par.fixed[3]+rep$par.random[9:16])

#the plotting will be done here
for(i in 1:8){
  rates = rep(0,length(min(df$daynr):max(df$daynr)))
  size = rep(0,length(min(df$daynr):max(df$daynr)))
  obs =rep(0,length(min(df$daynr):max(df$daynr)))
  for(j in min(df$daynr):max(df$daynr) ){
    obs[j-min(df$daynr)+1] = sum(df$ovul[which(df$daynr==j&df$region==4100+i)])
    rates[j-min(df$daynr)+1] = sum(df$ovul[which(df$daynr==j&df$region==4100+i)])/length(df$ovul[which(df$daynr==j&df$region==4100+i)])
    size[j-min(df$daynr)+1] = length(df$ovul[which(df$daynr==j&df$region==4100+i)])
    
  }
  main = paste("Region = ", locations[i])
  plot(t,q[i]*pnorm(t,mu[i],sigma),type="l",xlab="Days into the year",ylab ="Probability of ovulation",main=main,ylim=c(0,1.2),col="red")
  points(t,rates,cex=0.15*sqrt(size),pch=19)
  keep = which(!is.nan(rates))
  pred = predict(loess(rates[keep]~c(1:75)[keep],weights = size[keep]),se=TRUE)
  lines(t[keep],pred$fit+qt(0.975,pred$df)*pred$se.fit,col="blue",lty=2)
  lines(t[keep],pred$fit-qt(0.975,pred$df)*pred$se.fit,col="blue",lty=2)
  legend = c(paste("The model, mu =",round(mu[i],1),", sd=",round(sigma,1),", q=",round(q[i],2)),"observed ovulation rate", "95% confidence interval for ovulation rate")
  legend("topright",legend = legend,col=c("red","black","blue"),pch=c(NA,19,NA),lty=c(1,NA,2))
  
}

##########################################################################
############            Modell selection                      ############
############                                                  ############
##########################################################################
compile("withoutAge.cpp") # this is file described in the paper
dyn.load(dynlib("withoutAge"))

#The covariats can be changed and this will also chage the model

Y=df$ovul
N = rep(1,length(Y))
TID = df$daynr
###the model
id = factor(df$region)
idy = factor(df$year)
age = factor(df$age-1) # so age goes from 0-19 and can be used as index in the c++ file
XM = cbind(N) 
XS = cbind(N,df$cweight) 
XQ = cbind(df$cweight) 
data <- list(Y = Y,N = N,TID = TID,XM = XM,XS = XS,XQ = XQ,idm = id,idq = id,ids = id,age = age,idym=idy,idys=idy,idyq=idy)
parameters = list(BetaM=c(280),BetaS=c(log(10),0),BetaQ=c(0),
                  um=rep(0,nlevels(id)),uq=rep(0,nlevels(id)),us=rep(0,nlevels(id)),aq = rep(0,nlevels(age)),
                  log_sigmaM=log(0.5),log_sigmaQ=log(0.5),log_sigmaS=log(0.5),log_sigmaAq=log(0.5),
                  ym=rep(0,nlevels(idy)),ys=rep(0,nlevels(idy)),yq=rep(0,nlevels(idy)),
                  log_tauM=log(0.5),log_tauS=log(0.5),log_tauQ=log(0.5))
random = c("us","um","uq","aq","ym","ys","yq")
map <-list()# a list of what random effects should not be included in the model
obj <- MakeADFun(data,parameters,DLL = "withoutAge",map=map,random = random, silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")#optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
rep
my_aic(res$value,length(rep$par.fixed))# based onn this score the best model can be chosen
################################################################
#the best model based on AIC
#mu ~1  + region + year
#log_sig ~1 +year +region+ (weight)
#q ~age + weight + region + year
#aic= 7917.718
################################################################


##############################################################################
#######   The simulation, or bootstrapping, start here    ####################
##############################################################################


############################################################################################################
#this is the part where I define how the parameters are estimated, this will change depending on how they are estimated
#need to run the proper res first, the model needs to be on the form 
#mu~1+region + year
#sigma~1+region + year+weight
#q ~ age + weight + region + year
TID = df$daynr
mu = rep$par.fixed[[1]] + (rep$par.random[df$region-4100]) + rep$par.random[df$year-1980+45]
log_sig = rep$par.fixed[[2]]+(rep$par.random[df$region-4100+8*2]) + rep$par.random[df$year-1980+45+29]
logit_q = rep$par.fixed[[3]]*df$cweight+(rep$par.random[df$region-4100+8]) + rep$par.random[df$year-1980+45+2*29] + rep$par.random[24+df$age]
############################################################################################################

full_model = cbind()
for(i in 1:50){# the parameters, data, map, random and file need to be from the same model
  print(i)
  full_model = cbind(full_model,bootstrapping(1,mu,exp(log_sig),plogis(logit_q),TID,data,parameters,random,map,"withoutAge")) 
}
saveRDS(full_model,"simulations.rds")
full_model = readRDS("simulations.rds")


##############################################################################################
#here I make a model based only on age and region, this makes it easier to plot the effects.
# and since age and weight are heavily correlated the weight effect will be incorperated into age kind of
#a model based only on age and region
locations = c("Troms","Nordland","Nord-Trondelag","Oppland","Hedmark","Telemark","Aust-Agder","Vest-Agder")
compile("regionAndAge.cpp")
dyn.load(dynlib("regionAndAge"))
Y=df$ovul
N = rep(1,length(Y))
TID = df$daynr
###the model
id = factor(df$region)
idy = factor(df$year)
age = factor(df$age-1) # so age goes from 0-19 and can be used as index in the c++ file
XM = cbind(N) 
XS = cbind(N) 
data <- list(Y = Y,N = N,TID = TID,XM = XM,XS = XS,zeros=rep(0,sum(N)),idm = id,idq = id,ids = id,age = age,idym=idy,idys=idy,idyq=idy)
parameters = list(BetaM=c(280),BetaS=c(log(10)),
                  um=rep(0,nlevels(id)),uq=rep(0,nlevels(id)),us=rep(0,nlevels(id)),aq = rep(0,nlevels(age)),
                  log_sigmaM=log(0.5),log_sigmaQ=log(0.5),log_sigmaS=log(0.01),log_sigmaAq=log(0.5),
                  ym=rep(0,nlevels(idy)),ys=rep(0,nlevels(idy)),yq=rep(0,nlevels(idy)),
                  log_tauM=log(0.5),log_tauS=log(0.5),log_tauQ=log(0.5))
random = c("us","um","uq","aq","ym","ys","yq")
map <-list(log_tauS=factor(NA),log_tauQ=factor(NA),log_tauM=factor(NA))#log_sigmaQ=factor(NA),log_sigmaS=factor(NA),log_sigmaM=factor(NA),
obj <- MakeADFun(data,parameters,DLL = "regionAndAge",map=map,random = random, silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")#optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
rep
mu = rep$par.fixed[[1]] + (rep$par.random[df$region-4100]) 
log_sig = rep$par.fixed[[2]]+(rep$par.random[df$region-4100+8*2]) 
logit_q = (rep$par.random[df$region-4100+8]) +  rep$par.random[24+df$age]

# starting some bootstrpping for this model
I_want_to_boot = FALSE
if(I_want_to_boot){
B=1000
for(i in 1:B){
  print(i)
  temp = bootstrapping(1,mu,exp(log_sig),plogis(logit_q),TID,data,parameters,random,map,"regionAndAge")
  regionAndAge = cbind(regionAndAge,temp)
}
}
I_want_to_save_boot = FALSE
if(I_want_to_save_boot){
saveRDS(regionAndAge,"regionsAndAge.rds")
}
regionAndAge = readRDS("regionsAndAge.rds")
beta = rep(0,dim(regionAndAge)[1])
var = rep(0,dim(regionAndAge)[1])
for(i in 1:dim(regionAndAge)[1]){
  beta[i] = mean(regionAndAge[i,])
  var[i] = var(regionAndAge[i,])
}

#this is the part where i plot the results
for(age in 5:5){
  for(region in 4101:4108){
    plotting_regionAndAge(beta,df,c(5:13),region,regionAndAge)
  }
}





#######################################################################################
################################## The  data ##########################################
#######################################################################################


#the weight
hist(as.numeric(df$cweight),freq = F,breaks = 30,main = "Histogram of the weight of moose",xlab = "the weight")
vekt = as.numeric(df$cweight)
lines(50:300,dnorm(50:300,mean(vekt),sqrt(1/(length(vekt)-1)*sum((vekt-mean(vekt))^2))),col="red")
legend("topright",legend = c("Normal distibution","mu=mean(weight)","sigma^2=s^2(weight)" ),col="red",lty=c(1,0,0))

#the weight as a function of age
plot(df$age,df$cweight,xlab = "Age / years", ylab = "Weight / kg",main = "Observed weigth distribution for different ages")
vekt = rep(0,20)
for(i in 1:20){
  vekt[i] = mean(df$cweight[which(df$age==i)])
}
lines(1:20,vekt,col="red")
legend("bottomright",legend="Mean weight for each age",lty=1,col="red")


# Simple Pie Chart
#the the region
slices <- hist(df$region)$counts[c(1:2,4,6,8,10,12,14)]
lbls <- c("Troms","Nordland","Nord-Trondelag","Oppland","Hedemark","Telemark","Aust-Agder","Vest-Agder")
pie(slices, labels = paste(lbls,",",slices), main="Pie Chart of Counties")




######################################################################################
######################################################################################
######################################################################################

#models looking only at q in different ways

#first a model with q as 2RW 
Y=df$ovul
N = rep(1,length(Y))
TID = df$daynr
###the model
id = factor(df$region)
idy = factor(df$year)
age = factor(df$age-1) # so age goes from 0-19 and can be used as index in the c++ file
XM = cbind(N) 
XS = cbind(N) 
data <- list(Y = Y,N = N,TID = TID,XM = XM,XS = XS,zeros=rep(0,sum(N)),idm = id,idq = id,ids = id,age = age,idym=idy,idys=idy,idyq=idy)
parameters = list(BetaM=c(280),BetaS=c(log(10)),
                  um=rep(0,nlevels(id)),uq=rep(0,nlevels(id)),us=rep(0,nlevels(id)),aq = rep(0,nlevels(age)),
                  log_sigmaM=log(0.5),log_sigmaQ=log(0.5),log_sigmaS=log(0.01),log_sigmaAq=log(0.5),
                  ym=rep(0,nlevels(idy)),ys=rep(0,nlevels(idy)),yq=rep(0,nlevels(idy)),
                  log_tauM=log(0.5),log_tauS=log(0.5),log_tauQ=log(0.5))
random = c("us","um","uq","aq","ym","ys","yq")
map <-list(log_sigmaM=factor(NA),log_sigmaS=factor(NA),log_sigmaQ=factor(NA),log_tauS=factor(NA),log_tauQ=factor(NA),log_tauM=factor(NA))#log_sigmaQ=factor(NA),log_sigmaS=factor(NA),log_sigmaM=factor(NA),
obj <- MakeADFun(data,parameters,DLL = "regionAndAge",map=map,random = random, silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")#optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
rep
q_2rw = plogis(rep$par.random[25:44])

# now a model with age as random intercept 
####
Y=df$ovul # the responce in the bernulli-function
N = rep(1,length(Y)) 
TID = df$daynr # the kill date for all the moose
id = factor(df$age-1)
XM = cbind(N)
XS = cbind(N)
XQ = cbind(N)
data = list(Y=Y,TID=TID,N=N,XM=XM,XS=XS,XQ=XQ,idm=id,ids=id,idq=id)
parameters = list(BetaM = c(275),BetaS = c(log(10)),BetaQ = c(0),
                  um = rep(0,nlevels(id)),us = rep(0,nlevels(id)),uq = rep(0,nlevels(id)),
                  log_sigmaM = log(0.5), log_sigmaS = log(0.5), log_sigmaQ = log(0.5) )
random = c("um","us","uq")
map = list(log_sigmaS=factor(NA),log_sigmaM=factor(NA))
obj <- MakeADFun(data,parameters,DLL = "simpleModelRegion",random=random,map=map ,silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
q_ri = plogis(rep$par.fixed[3] + rep$par.random[21:40])


# a model with one intercept for each of the 20 ages
Y=df$ovul # the responce in the bernulli-function
N = rep(1,length(Y)) 
TID = df$daynr # the kill date for all the moose
id = factor(df$age-1)
XM = cbind(N)
XS = cbind(N)
XQ = cbind(as.numeric(df$age ==1))
for(i in 2:20){
  XQ = cbind(XQ,as.numeric(df$age ==i))
}
data = list(Y=Y,TID=TID,N=N,XM=XM,XS=XS,XQ=XQ,idm=id,ids=id,idq=id)
parameters = list(BetaM = c(275),BetaS = c(log(10)),BetaQ = rep(0,20),
                  um = rep(0,nlevels(id)),us = rep(0,nlevels(id)),uq = rep(0,nlevels(id)),
                  log_sigmaM = log(0.5), log_sigmaS = log(0.5), log_sigmaQ = log(0.5) )
random = c("um","us","uq")
map = list(log_sigmaS=factor(NA),log_sigmaM=factor(NA),log_sigmaQ=factor(NA))
obj <- MakeADFun(data,parameters,DLL = "simpleModelRegion",random=random,map=map ,silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
q_20 = plogis(rep$par.fixed[3:22])

#here comes the plot 


plot(q_2rw,xlab = "age",ylab = "ovulation probability/q", main = "Modeling q with different aproches",col="black",type="l")
lines(q_ri,col="blue")
lines(q_20,col="red")
legend("bottom",legend= c("q with age as 2RW. #parameters for q = 1",
                 "q with age as random intercept. #parameters for q = 2",
                 "q with age as factor. #parameters for q = 20"),col =c("black","blue","red"),lty=1)


###################################################################################################
###################################################################################################
#plot ovulation rate for different weights and different ages



#model with q as 2RW and weight as covariat

Y=df$ovul
N = rep(1,length(Y))
TID = df$daynr
###the model
id = factor(df$region)
idy = factor(df$year)
age = factor(df$age-1) # so age goes from 0-19 and can be used as index in the c++ file
XM = cbind(N) 
XS = cbind(N) 
XQ = cbind(df$cweight) 
data <- list(Y = Y,N = N,TID = TID,XM = XM,XS = XS,XQ = XQ,idm = id,idq = id,ids = id,age = age,idym=idy,idys=idy,idyq=idy)
parameters = list(BetaM=c(280),BetaS=c(log(10)),BetaQ=c(0),
                  um=rep(0,nlevels(id)),uq=rep(0,nlevels(id)),us=rep(0,nlevels(id)),aq = rep(0,nlevels(age)),
                  log_sigmaM=log(0.5),log_sigmaQ=log(0.5),log_sigmaS=log(0.5),log_sigmaAq=log(0.5),
                  ym=rep(0,nlevels(idy)),ys=rep(0,nlevels(idy)),yq=rep(0,nlevels(idy)),
                  log_tauM=log(0.5),log_tauS=log(0.5),log_tauQ=log(0.5))
random = c("us","um","uq","aq","ym","ys","yq")
map <-list(log_sigmaM=factor(NA),log_sigmaS=factor(NA),log_sigmaQ=factor(NA),log_tauM=factor(NA),log_tauS=factor(NA),log_tauQ=factor(NA))# a list of what random effects should not be included in the model
obj <- MakeADFun(data,parameters,DLL = "withoutAge",map=map,random = random, silent = TRUE)
res <- optimx::optimr(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")#optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")
rep = sdreport(obj)
rep

plot(1:20,plogis(rep$par.fixed[3]*120+rep$par.random[25:44]),xlab = 
       "age/years",ylab = "q",type="l",ylim = c(0,1),col=120,
     main = "Ovulation probability for different ages and weights")
for(vekt in c(130,140,150,160,170,180,190,200)){
  lines(1:20,plogis(rep$par.fixed[3]*vekt+rep$par.random[25:44]),col=vekt)
}
legend = c(120,130,140,150,160,170,180,190,200)
legend("bottom",legend=paste("weight = ",legend,"kg"), col = legend,lty =1)
















