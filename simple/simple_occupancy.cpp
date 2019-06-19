#include <TMB.hpp>

// logspace_add computes log(exp(logx)+exp(logy))
// https://kaskr.github.io/adcomp/group__special__functions.html#ga2d5a07f1202dc628ba18ee9a2af8b9cd
// -log(exp(0) + exp(-log(p/(1-p))) = -log (1 + (1-p)/p) =
// -log (1/p) = log(p)
// convert logit(p) to log(p)
#define invlogit_log(X)    (-logspace_add( Type(0), -X))
// convert logit(p) to log(1-p)
#define invlogit_log1mp(X) (-logspace_add( Type(0),  X))

template<class Type>
Type objective_function<Type>::operator() () {
  
  // parameters
  PARAMETER_VECTOR(beta_p);
  PARAMETER_VECTOR(beta_phi);
  
  // data
  DATA_VECTOR(obs);      // observations (not integer??)
  DATA_IVECTOR(site);     // sites corresponding to obs
  DATA_IVECTOR(site_tot); // total obs per site (FIXME: compute this up-front in TMB?)
  DATA_IVECTOR(site_ind); // site indices (includes 0 and N)
  DATA_MATRIX(Xphi);      // occupancy model matrix (Ns x length(beta_phi))
  DATA_MATRIX(Xp);        // detection model matrix (N  x length(beta_p

  int Ns = site_ind.size()-1; // nb sites

  // linear predictors for occupancy and detection
  vector<Type> eta_phi = Xphi * beta_phi;
  vector<Type> eta_p = Xp * beta_p;

  Type sitesum, nll=0;

  for (int i = 0; i < Ns; i++) { // loop over sites
    // FIXME: more vectorized?
    sitesum=0.0;
    // conditional detection prob for all sites
    // (obs(j) could be all-zero)
    for (int j = site_ind(i); j < site_ind(i+1); j++) {
      sitesum += dbinom_robust(obs(j), Type(1), eta_p(j), true);
    }
    if (site_tot(i) > 0) {
      // non-empty site: log(p_phi(i)) + sum(log(binom detect))
      nll -= (invlogit_log(eta_phi(i)) + sitesum);
    } else {
      // empty site:
      //  unocc = log(1-p_phi(i))
      //  undet = log(p_phi(i)) + sum(log(binom detect)) [all zero]
      nll -= logspace_add(invlogit_log1mp(eta_phi(i)),
			  invlogit_log(eta_phi(i)) + sitesum);
    }
  }
  return nll;
}
