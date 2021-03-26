functions {
  real infected(real t, real I0, real mu, real tau, real eta) {
    if(t<0) return 0;
    if(t<=tau)  return  I0*exp(mu * t);
    else return I0*exp(mu * tau - eta * (t - tau));
  }
   real binomial_lpmf_approx(real x, real n, real p)
   {
     if(x<0 || x>n) return -5000000;
     return lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) +x*log(p)+(n-x)*log(1-p);
   }
}


data {
  int<lower=0> T;            // current time
  int<lower=0> J;
  int<lower=0> time[J];
  int<lower=0> m_p;          // number of positive from the serology test
  int<lower=0> n;            // number of people who did the serology test
  int<lower=0> N;            // population sizes
  int<lower=0> C_cc[J];      // cases counts
  real prior_mu[2];          // prior parameters for mu
  real prior_eta[2]; 	      // prior parameters for eta
  real prior_I0[2];          // prior parameters for I0
  real prior_tau[2];         // prior parameters for tau
  real prior_tau0[2];        // prior parameters for tau0
  real prior_p_cc1[2];       // prior parameters for p_cc1
  real prior_p_cc2[2];       // prior parameters for p_cc2
  real prior_specificity[2]; // prior parameters for p_specificity
  real prior_sensitivity[2]; // prior parameters for p_sensitivity
}


parameters {
  real<lower=0, upper=1> sensitivity; // sensitivity
  real<lower=0, upper=1> specificity; // specificity
  real<lower=1> I0;                   // initial cases
  real<lower=0> mu;                   // growth rate (t<=tau)
  real<lower=1> tau;                  // change pt from growth to decay
  real<lower=0> eta;                  // decay rate (t>tau)
  real<lower=1> tau_0;                // delay parameter
  real<lower=0, upper=1> p_cc1;       // probablity of observing case before April 14
  real<lower=0, upper=1> p_cc2;       // probability of observing case after April 14
}


transformed parameters {
  vector[T] I;
  real<lower=0, upper=1> p_s;
  real<lower=0, upper=1> p_testpositive;
  real<lower=0, upper=1> p_testnegative;
  for(t in 1:T){
    I[t] = infected(t, I0, mu, tau, eta);
  }
  p_s = sum(I)/N;
  p_testpositive = p_s*sensitivity+(1-p_s)*(1-specificity);
  p_testnegative = p_s*(1-sensitivity)+(1-p_s)*specificity;
}


model {
//  target += binomial_lpmf_approx(m_p, n, p_testpositive); 
  target += m_p*log(p_testpositive) + (n-m_p)*log(p_testnegative);
  for(i in 1:J){
    real infected_i = infected(time[i]-tau_0, I0, mu, tau, eta);

    if(time[i]-tau_0 <=0) continue;
    
    if( time[i] < 79 ) target +=  binomial_lpmf_approx(C_cc[i], infected_i, p_cc1);
    else target +=  binomial_lpmf_approx(C_cc[i], infected_i, p_cc2);
  }
  target += beta_lpdf(p_cc1 | prior_p_cc1[1], prior_p_cc1[2]);
  target += beta_lpdf(p_cc2 | prior_p_cc2[1], prior_p_cc2[2]);
  target += beta_lpdf(sensitivity | prior_sensitivity[1], prior_sensitivity[2]);
  target += beta_lpdf(specificity | prior_specificity[1], prior_specificity[2]);
  target += normal_lpdf(I0 | prior_I0[1], prior_I0[2]);
  target += normal_lpdf(tau | prior_tau[1], prior_tau[2]);
  target += normal_lpdf(tau_0 | prior_tau0[1], prior_tau0[2]);
  target += normal_lpdf(mu | prior_mu[1], prior_mu[2]);
  target += normal_lpdf(eta | prior_eta[1], prior_eta[2]);
}


