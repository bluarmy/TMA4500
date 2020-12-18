// model for moose data
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // The data to be input
  DATA_VECTOR(Y);   // Response vector
  DATA_VECTOR(N);
  DATA_VECTOR(TID); // the time of kill
  DATA_MATRIX(XM);      // matrix for mean
  //DATA_MATRIX(XQ);      //matrix for q
  DATA_MATRIX(XS);      // matrix for sigma
  DATA_FACTOR(idm);     // factors for random effect mu
  DATA_FACTOR(idq);
  DATA_FACTOR(ids);
  DATA_FACTOR(age);
  DATA_FACTOR(idym);
  DATA_FACTOR(idys);
  DATA_FACTOR(idyq);
  DATA_VECTOR(zeros);
  
  //parameterne
  //the mean
  PARAMETER_VECTOR(BetaM);
  PARAMETER_VECTOR(BetaS);
  //PARAMETER_VECTOR(BetaQ); // the parameter for q
  PARAMETER_VECTOR(um); //random effects for mu
  PARAMETER(log_sigmaM); // log sd for random effect mu
  PARAMETER_VECTOR(uq);
  PARAMETER(log_sigmaQ);
  PARAMETER_VECTOR(us);
  PARAMETER(log_sigmaS);
  PARAMETER_VECTOR(aq);
  PARAMETER(log_sigmaAq);
  PARAMETER(log_tauM);
  PARAMETER(log_tauS);
  PARAMETER(log_tauQ);
  PARAMETER_VECTOR(ym);
  PARAMETER_VECTOR(ys);
  PARAMETER_VECTOR(yq);
  
  
  Type nll = 0; 
  //adding the random effects to be integrated out
  nll -= sum(dnorm(um, Type(0), exp(log_sigmaM), true));
  nll -= sum(dnorm(uq,Type(0),exp(log_sigmaQ),true));
  nll -= sum(dnorm(us,Type(0),exp(log_sigmaS),true));
  nll -= sum(dnorm(ym, Type(0), exp(log_tauM), true));
  nll -= sum(dnorm(yq,Type(0),exp(log_tauQ),true));
  nll -= sum(dnorm(ys,Type(0),exp(log_tauS),true));
  
  // adding the log of the pdf of a_t to the nll to be integrated out
  
  for (int t = 2; t < 20 ;t++){
    nll -= dnorm(aq[t],2*aq[t-1]-aq[t-2],exp(log_sigmaAq),true);// a new 2orw for age but for q
  }
  
  
  //creating the vectors to add the random effects to
  vector<Type> p;
  vector<Type> mu = XM*BetaM;
  vector<Type> q = zeros;
  vector<Type> lsd = XS*BetaS;
  //adding the random effects if they are included
  
  for (int i =0; i< N.size();i++){
    if (CppAD::Variable(log_sigmaM) ){ // If log_sigma=NA in map argument then don't add random effect to eta
      mu[i] += um(idm(i));
    }
    if (CppAD::Variable(log_tauM) ){ // If log_sigma=NA in map argument then don't add random effect to eta
      mu[i] += ym(idym(i));
    }
    
    if (CppAD::Variable(log_sigmaQ) ) {
      q[i] += uq(idq(i));
    }
    if (CppAD::Variable(log_tauQ) ) {
      q[i] += yq(idyq(i));
    }
    if (CppAD::Variable(log_sigmaS) ) {
      lsd[i] += us(ids(i));
    }
    if (CppAD::Variable(log_tauS) ) {
      lsd[i] += ys(idys(i));
    }
  }
  //adding the a_t to the expectation   
  
  
  
  for (int i = 0; i < N.size(); i++){
    q[i] += aq[age[i]];
  }
  //the probabilities for ovulation given by the model 
  p = pnorm(TID,mu,exp(lsd));	
  // creating the nll using the binomial with prob of success p
  for (int i = 0; i < N.size(); i++) {
    
    nll -= dbinom(Y[i],N[i],1/(1+exp(-q[i]))*p[i],true);
    
  }
  return nll;
}
