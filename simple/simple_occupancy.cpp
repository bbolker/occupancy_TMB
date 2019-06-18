#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  
  // parameters
  PARAMETER_VECTOR(beta_p);
  PARAMETER_VECTOR(beta_psi);
  
  // data
  DATA_IVECTOR(obs);      // observations
  DATA_IVECTOR(site);     // sites corresponding to obs
  DATA_IVECTOR(site_tot); // total obs per site (FIXME: compute this up-front in TMB?)
  DATA_MATRIX(Xpsi);      // occupancy model matrix
  DATA_MATRIX(Xp);        // detection model matrix
  
  int ns = Xpsi.rows();   // nb of sites
  int N = obs.size();     // total obs
  int npar = b.size();
  int npar_p = Xpsi.rows(); // nb of detection covariates

  // FIXME: do we need to drop dims?
  vector<Type> psi_prob = invlogit(Xpsi * beta_psi);
  vector<Type> p_prob_0 = invlogit(Xp * beta_p);
  
  vector
  // likelihood
  for (int i = 0; i < N; i++) { // loop over sites
    if (obs[
  }

  for (int i = 0; i < nh; i++) {
    vector<int> evennt = ch.row(i);
    ALPHA = PROP * vector<Type>(B.row(e(i))); // element-wise vector product
    REPORT(ALPHA);
    for (int j = 1; j < N; j++) {
      matrix<Type> TEMP = PHI.col(j-1);
      matrix<Type> PHIj(2,2);
      PHIj(0,0) = TEMP(0,0);
      PHIj(0,1) = TEMP(1,0);
      PHIj(1,0) = TEMP(2,0);
      PHIj(1,1) = TEMP(3,0);
      ALPHA = multvecmat(ALPHA,PHIj)* vector<Type>(B.row(evennt(j))); // vector matrix product, then element-wise vector product
    }
    ll += log(sum(ALPHA));
  }
  nll = -ll;
  return nll;
}
