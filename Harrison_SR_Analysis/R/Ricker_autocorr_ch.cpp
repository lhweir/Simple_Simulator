#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


template <class Type>
Type minus_one_to_one(Type x)
{
  return Type(2) * invlogit(x) - Type(1);
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_R);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  PARAMETER(a);
  PARAMETER(logbeta);
  PARAMETER(logSigObs);
  PARAMETER(rho);
  //PARAMETER_VECTOR(delta);
  
  int timeSteps=obs_R.size();

  Type rhoo = minus_one_to_one(rho);

  //Type rho = exp(logrho);
  Type beta = exp(logbeta);
  Type SigObs = exp(logSigObs);
  Type Smax  = Type(1.0)/beta;
  Type alpha = exp(a);
  
  Type tau  = Type(1.0)/(SigObs*SigObs);

  Type SigAR  = SigObs*sqrt(1-pow(rhoo,2));
  
  vector<Type> pred_logRS(timeSteps), obs_logRS(timeSteps), residuals(timeSteps), epsilon(timeSteps);
  vector<Type> std_resids(timeSteps);

 

  //priors on precision and variance ratio
  //Type ans= -dbeta(rho,Type(3.0),Type(3.0),true);  
  //Type ans= -dnorm(logSigObs,Type(0.0),Type(5.0),true);   
  Type ans = Type(0);

  obs_logRS(0) = log(obs_R(0)/obs_S(0));
  pred_logRS(0) = a - beta * obs_S(0) ;
  residuals(0) = obs_logRS(0) - pred_logRS(0);
  std_resids(0) = residuals(0)/SigAR;
  epsilon(0) =  residuals(0);
  
  
  ans+= -dnorm(obs_logRS(0),pred_logRS(0),SigObs,true);

  for(int i=1;i<timeSteps;i++){
    if(!isNA(obs_R(i))){
      obs_logRS(i) = log(obs_R(i)/obs_S(i));
      pred_logRS(i) = a - beta * obs_S(i) + epsilon(i-1) * rhoo ;
     // pred_logR(i) = logRS(i) + log(obs_S(i)) + epsilon(i-1) * rhoo ;
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      epsilon(i) = residuals(i);//epsilon(i-1) * rho; //+ delta(i)* sqrt(1-pow(rho,2));
      std_resids(i) = residuals(i)/SigAR;
      
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),SigAR,true);
      
    }
  
  }
  
  Type umsy     = Type(.5) * a - Type(0.07) * (a * a);
  Type Smsy     = a/beta * (Type(0.5) -Type(0.07) * a);  
  


  ADREPORT(alpha);
  ADREPORT(Smax);
  ADREPORT(beta);
  ADREPORT(SigObs);
  ADREPORT(SigAR);
  ADREPORT(rhoo);
  REPORT(pred_logRS);
  REPORT(residuals);
  //REPORT(delta)
  REPORT(umsy)
  REPORT(epsilon)
  REPORT(alpha);
  REPORT(residuals);
  REPORT(std_resids);
  return ans;
}

