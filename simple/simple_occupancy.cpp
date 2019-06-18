#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  
  // parameters
  PARAMETER_VECTOR(beta_p);
  PARAMETER_VECTOR(beta_psi);
  
  // data
  DATA_VECTOR(obs);      // observations (not integer??)
  DATA_IVECTOR(site);     // sites corresponding to obs
  DATA_IVECTOR(site_tot); // total obs per site (FIXME: compute this up-front in TMB?)
  DATA_MATRIX(Xpsi);      // occupancy model matrix (Ns x length(beta_psi))
  DATA_MATRIX(Xp);        // detection model matrix (N  x length(beta_p

  
  int N = obs.size();     // total obs

  vector<Type> eta_psi = Xpsi * beta_psi;
  vector<Type> eta_p = Xp * beta_p;

  Type s1, tmp_loglik, nll=0;

  // likelihood (mostly stolen from glmmTMB.cpp)
  for (int i = 0; i < N; i++) {
    // dbinom_robust takes logit_prob rather than prob ...
    tmp_loglik = dbinom_robust(obs(i), Type(1), eta_p(i), true);

    Type logit_pz = eta_psi(site(i)) ;

    // logspace_add computes log(exp(logx)+exp(logy))
    // https://kaskr.github.io/adcomp/group__special__functions.html#ga2d5a07f1202dc628ba18ee9a2af8b9cd

    Type log_pz   = -logspace_add( Type(0) , -logit_pz );
    Type log_1mpz = -logspace_add( Type(0) ,  logit_pz );
    
    if (site_tot(site(i)) == Type(0)) {
      //  equivalent: pz = invlogit(logit_pz);
      //     log( pz(i) + (1.0 - pz(i)) * exp(tmp_loglik) );
      tmp_loglik = logspace_add( log_pz, log_1mpz + tmp_loglik );
    } else {
      // equivalent: log( 1.0 - pz(i))
      tmp_loglik += log_1mpz ;
    }
    nll -= tmp_loglik;
  }
  return nll;
}
