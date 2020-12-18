// model for moose data
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /*The data to be input
   the data is something the uses includes from the R file, and needs to be on the format 
  you see under here */
  DATA_VECTOR(Y);   // Response vector
  DATA_VECTOR(N);   // a vector of ones of length equal the length of Y
  DATA_VECTOR(TID); // the time of kill
  DATA_MATRIX(XM);  //data  matrix for mu
  DATA_MATRIX(XQ);  // data matrix for q
  DATA_MATRIX(XS);  // data matrix for sigma
  DATA_FACTOR(idm); // factors for random region effect mu
  DATA_FACTOR(idq); //factors for random region effect logit_q
  DATA_FACTOR(ids); //factors for random region effect log_sigma
  DATA_FACTOR(age); // the age of the moose-1 so that it can be used as index
  DATA_FACTOR(idym);//factors for random year effect mu
  DATA_FACTOR(idys);//factors for random year effect log_sigma
  DATA_FACTOR(idyq);//factors for random year effect logit_q

  /*the parameters
  this is the parameters to be estimated, when using the R-code 
  one needs to give starting values to all these parameter
  ideally the starting values should be close to the real solution.
  In the r file you specify what parameters should be included as random,
  those parameters are then integrated out with laplace approximation*/
  PARAMETER_VECTOR(BetaM);//the parameters for mu
  PARAMETER_VECTOR(BetaS);//the parameters for log_sigma
  PARAMETER_VECTOR(BetaQ); // the parameter for logit_q
  PARAMETER_VECTOR(um); //random region effects for mu
  PARAMETER(log_sigmaM); // log sd for random region effect mu
  PARAMETER_VECTOR(uq);//random region effect logit_q
  PARAMETER(log_sigmaQ); // log sd for random region effect logit_q
  PARAMETER_VECTOR(us);//random region effect log_sigma
  PARAMETER(log_sigmaS); // log sd for random region effect log_sigma
  PARAMETER_VECTOR(aq); //random age effect on logit_q
  PARAMETER(log_sigmaAq);//the log of the sd for the 2RW for logit_q
  PARAMETER(log_tauM);  // log sd for random year effect mu
  PARAMETER(log_tauS); // log sd for random year effect log_sigma
  PARAMETER(log_tauQ); // log sd for random year effect logit_q
  PARAMETER_VECTOR(ym);//random year effect mu
  PARAMETER_VECTOR(ys);//random year effect log_sigma
  PARAMETER_VECTOR(yq);//random year effect logit_q
  
  // the creation of the negative log likelihood
  //this is what the function returnes
  Type nll = 0; 
  //adding the random effects to be integrated out
  // all the random effects for region and year are random intercept
  //so they are independent of each other and are therfore just multiplied in the likelihood
  nll -= sum(dnorm(um, Type(0), exp(log_sigmaM), true));
  nll -= sum(dnorm(uq,Type(0),exp(log_sigmaQ),true));
  nll -= sum(dnorm(us,Type(0),exp(log_sigmaS),true));
  nll -= sum(dnorm(ym, Type(0), exp(log_tauM), true));
  nll -= sum(dnorm(yq,Type(0),exp(log_tauQ),true));
  nll -= sum(dnorm(ys,Type(0),exp(log_tauS),true));
  
  // adding the negative log of the pdf of a_t to the nll to be integrated out
  for (int t = 2; t < 20 ;t++){
    //this is the form of the f(a) and it is the negative log likelihood of a 2RW
    nll -= dnorm(aq[t],2*aq[t-1]-aq[t-2],exp(log_sigmaAq),true);
  }
  
  
  //creating the vectors to add the random effects to
  vector<Type> p; //a vector for the probabilities given that they ovulate
  vector<Type> mu = XM*BetaM; // the mu vector
  vector<Type> q = (XQ*BetaQ); // logit_q vector
  vector<Type> lsd = XS*BetaS; // log_sigma vector

  /*In this part the random effects are added if they are to be included.
  Map is an argument given in the r-file, in it one can specify that some random effects should not be added.
  The prior distribution to the random effect is gaussian and will therefore be integrated out exactly with the laplace 
  approximation.*/ 
  for (int i =0; i< N.size();i++){
    // If log_sigma=NA in map argument then don't add random region effect to mu
    if (CppAD::Variable(log_sigmaM) ){ 
      mu[i] += um(idm(i));
    }
    // If log_sigma=NA in map argument then don't add random year effect to mu
    if (CppAD::Variable(log_tauM) ){ 
      mu[i] += ym(idym(i));
    }
    //if log_sigma = factor(NA) in map argument don't add random region effect for logit_q
    if (CppAD::Variable(log_sigmaQ) ) {
      q[i] += uq(idq(i));
    }
    //if log_sigma = factor(NA) in map argument don't add random year effect for logit_q
    if (CppAD::Variable(log_tauQ) ) {
      q[i] += yq(idyq(i));
    }
    //if log_sigma = factor(NA) in map argument don't add random region effect for log_sigma
    if (CppAD::Variable(log_sigmaS) ) {
      lsd[i] += us(ids(i));
    }
    //if log_sigma = factor(NA) in map argument don't add random year effect for log_sigma
    if (CppAD::Variable(log_tauS) ) {
      lsd[i] += ys(idys(i));
    }
  }
   
  
  
  //adding the a_t to the expectation 
  for (int i = 0; i < N.size(); i++){
    q[i] += aq[age[i]];
  }
  //the probabilities for ovulation given that they ovulate 
  p = pnorm(TID,mu,exp(lsd));	
  // creating the nll using the binomial with prob of success p
  for (int i = 0; i < N.size(); i++) {
    // bernulli distribution and p_i given by the model
    nll -= dbinom(Y[i],N[i],1/(1+exp(-q[i]))*p[i],true);
    
  }
  return nll;
}
