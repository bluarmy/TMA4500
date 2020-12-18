// model for moose data

/*This is the c++ file for the simplest model in the project. It estimated two parameters, mu and log(sigma). 
 This file returns the negative log likelihood for the model. The MLE for mu and log(sigma) will be found using
 optimization functions in R
 */

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // The data to be input
  DATA_VECTOR(Y);   // Response vector
  DATA_VECTOR(N);   //A vector of ones of same length as Y
  DATA_VECTOR(TID); // the time of kill
  
  //THE PARAMETERS TO BE ESTIMATED 
  PARAMETER(mu); // mean ovulation time
  PARAMETER(log_sig); // the logarithm of the standard diviation
  
  Type nll = 0; //the negative log likelihood
  vector<Type> p; // the probability of ovulation for each moose
  p = pnorm(TID,mu,exp(log_sig));	
  
  // creating the nll using the binomial with prob of success p
  for (int i = 0; i < N.size(); i++) {
    
    nll -= dbinom(Y[i],N[i],p[i],true);
    
  }
  
  return nll;
}