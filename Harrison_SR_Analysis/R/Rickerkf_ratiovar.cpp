#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
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

  DATA_VECTOR(obs_R);   // observed recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  DATA_SCALAR(prbeta1);
  DATA_SCALAR(prbeta2);
  //logbeta     -> log of beta from ricker curve
  //a0      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //a       -> Time-varying log(alpha)

  PARAMETER(a0);
  PARAMETER(logbeta);
  PARAMETER(rho);
  PARAMETER(logvarphi);

  PARAMETER_VECTOR(a);
  
  int timeSteps=obs_R.size();

  Type beta = exp(logbeta);
  Type Smax  = Type(1.0)/beta;
  
  //theta       -> total standard deviation
  //sig         -> obs error std
  //tau         -> proc error ( a) std
  
  Type varphi     = exp(logvarphi);
  Type theta     = sqrt(Type(1.0)/varphi);   //varphi is precision, theta is s.d.
  Type sig       = sqrt(rho) * theta;
  Type tau        = sqrt(Type(1.0)-rho) * theta ;
   
  vector<Type> pred_logRS(timeSteps); 
  vector<Type> obs_logRS(timeSteps);
  vector<Type> residuals(timeSteps); 
  vector<Type> std_resids(timeSteps);
  vector<Type> alpha(timeSteps); //umsy(timeSteps), Smsy(timeSteps), ;
  
  alpha = exp(a);

  //priors on precision and variance ratio
  Type ans= -dbeta(rho, prbeta1, prbeta2, true);  
  
  ans+= -dnorm( a(0), a0, tau,true);
  //umsy(0)     = Type(.5) *  a(0) - Type(0.07) * ( a(0) *  a(0));
  //Smsy(0)     =   a(0)/beta * (Type(0.5) -Type(0.07) *  a(0));  

  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm( a(i), a(i-1),tau,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_R(i))){
      obs_logRS(i) = log(obs_R(i)/obs_S(i));
      pred_logRS(i) =  a(i) - beta * obs_S(i) ;
      //pred_logR(i) = logRS(i) + log(obs_S(i));
      //umsy(i) = Type(.5) *  a(i) - Type(0.07) * ( a(i) *  a(i)); 
      //Smsy(i) =   a(i)/beta * (Type(0.5) -Type(0.07) *  a(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      std_resids(i) = residuals(i)/sig;
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),sig,true);
    }
  
  }
  
  // also want avg of last 4 a's
  Type a_avg = (a(timeSteps-1) + a(timeSteps-2) + a(timeSteps-3) + a(timeSteps-4) )/4;
  Type alpha_avg = exp(a_avg);
  
  ADREPORT(alpha);
  ADREPORT(Smax);
  ADREPORT(a_avg);
  ADREPORT(alpha_avg);
  ADREPORT(theta);
  ADREPORT(beta);
  REPORT(alpha);
  REPORT(residuals);
  REPORT(std_resids);

  return ans;
}

