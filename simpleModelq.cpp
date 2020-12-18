// model for moose data
/*
 
 This is the c++ file for the simplest model in the project. It estimated three parameters, mu, log(sigma) and logit(q). 
 This file returns the negative log likelihood for the model. The MLE for mu, logit(q) and log(sigma) will be found using
 optimization functions in R
 
 */
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // The data to be input
  DATA_VECTOR(Y);   // Response vector
  DATA_VECTOR(N);   // a vector of ones of the same length as Y
  DATA_VECTOR(TID); // the time of kill

  PARAMETER(mu); // the mean ovulation time
  PARAMETER(log_sig); // the log standard deviation
  PARAMETER(q); // the logit of the ovulation probability
  
  Type nll = 0;  // the negative log likelihood
  vector<Type> p; // the probability of ovulation at the time of kill given that it would ovulate
  p = pnorm(TID,mu,exp(log_sig));	
  
  // creating the nll using the binomial with prob of success p
  for (int i = 0; i < N.size(); i++) {
    
    nll -= dbinom(Y[i],N[i],1/(1+exp(-q))*p[i],true);
    
  }
 
  return nll;
}