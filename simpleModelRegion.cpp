// model for moose data
/*
 This model used only region as covariates. It will estimate BetaM whitch is a vector of length 1, 
 and will be the mu from the model. The same goes for BeatS and BetaQ. um, us and uq are the random effects from region. 
 Since this model is just used for random regions all the matricies XM, XS, XQ will just be nx1, where n is the number 
 of observed moose.
 */

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // The data to be input
  DATA_VECTOR(Y);   // Response vector
  DATA_VECTOR(N);   // A vector of ones same length as Y
  DATA_VECTOR(TID); // the time of kill
  DATA_MATRIX(XM);      // matrix for mean 
  DATA_MATRIX(XQ);      //matrix for q
  DATA_MATRIX(XS);      // matrix for sigma
  DATA_FACTOR(idm);     // factors for random effect mu
  DATA_FACTOR(idq);     // factors for random effect log(sigma)
  DATA_FACTOR(ids);    // factors for random effect logit(q)
  
  //parameterne
  PARAMETER_VECTOR(BetaM); //mu
  PARAMETER_VECTOR(BetaS);//log(sigma)
  PARAMETER_VECTOR(BetaQ); // logit(q)
  PARAMETER_VECTOR(um); //random effects for mu
  PARAMETER(log_sigmaM); // log sd for random effect mu
  PARAMETER_VECTOR(uq);//random effects for logit(q)
  PARAMETER(log_sigmaQ);// log sd for random effect logit(q)
  PARAMETER_VECTOR(us);//random effects for log(sigma)
  PARAMETER(log_sigmaS);// log sd for random effect log(sigma)

  
  
  Type nll = 0; 
  //adding the random effects to be integrated out
  nll -= sum(dnorm(um, Type(0), exp(log_sigmaM), true));
  nll -= sum(dnorm(uq,Type(0),exp(log_sigmaQ),true));
  nll -= sum(dnorm(us,Type(0),exp(log_sigmaS),true));
  
  
  
  //creating the vectors to add the random effects to
  vector<Type> p;
  vector<Type> mu = XM*BetaM;
  vector<Type> q = (XQ*BetaQ);
  vector<Type> lsd = XS*BetaS;
  //adding the random effects if they are included
  
  for (int i =0; i< N.size();i++){
    if (CppAD::Variable(log_sigmaM) ){ // If log_sigmaM=NA in map argument then don't add random effect to mu
      mu[i] += um(idm(i));
    }
    
    if (CppAD::Variable(log_sigmaQ) ) {// If log_sigmaQ=NA in map argument then don't add random effect to logit(q)
      q[i] += uq(idq(i));
    }
    if (CppAD::Variable(log_sigmaS) ) {// If log_sigmaS=NA in map argument then don't add random effect to log(sigma)
      lsd[i] += us(ids(i));
    }

  }
  
  
  
  //the probabilities for ovulation given by the model 
  p = pnorm(TID,mu,exp(lsd));	
  // creating the nll using the binomial with prob of success p
  for (int i = 0; i < N.size(); i++) {
    
    nll -= dbinom(Y[i],N[i],1/(1+exp(-q[i]))*p[i],true);
    
  }
  return nll;
}