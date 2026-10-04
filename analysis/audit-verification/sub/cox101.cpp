#include <RcppArmadillo.h>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

struct CoxStatsFull {
  arma::vec g_z;   // gradient wrt eta: d - haz_c ∘ h_c  (n×1)
  arma::vec w;     // diag Hessian approx: haz_c ∘ h_c   (n×1)
  double loglik;
};

// eta, d must be ordered by nondecreasing time (risk sets are suffixes)
inline CoxStatsFull cox_suffix_stats_full(const arma::vec& eta,
                                          const arma::vec& d,
                                          const arma::vec& time) {
  const arma::uword n = eta.n_rows;
  arma::vec m(n), rsk_c(n);

  // suffix max
  m[n-1] = eta[n-1];
  for (arma::sword i = (arma::sword)n-2; i >= 0; --i) m[i] = std::max(m[i+1], eta[i]);

  // centered cumulative risk sums: rsk_c[i] = Σ_{j>=i} exp(eta[j]-m[i])
  rsk_c[n-1] = std::exp(eta[n-1] - m[n-1]);
  for (arma::sword i = (arma::sword)n-2; i >= 0; --i) {
    const double scale = std::exp(m[i+1] - m[i]);
    rsk_c[i] = std::exp(eta[i]-m[i]) + rsk_c[i+1]*scale;
  }

  arma::uvec gstart(n);
  gstart[0] = 0;
  for (arma::uword i = 1; i < n; ++i) {
    gstart[i] = (time[i] == time[i-1]) ? gstart[i-1] : i;
  }

  arma::vec logS0_full = arma::log(rsk_c) + m;
  arma::vec logS0 = logS0_full - m[0];
  // cum sum of d / logS0 (centered domain)
  arma::vec h(n);
  // cumulative sum of d / S0^2
  arma::vec B(n);

  const arma::uword g0 = gstart[0];
  h[0] = d[0] * std::exp(-logS0[g0]);
  B[0] = d[0] * std::exp(-2.0 * logS0[g0]);

  for (arma::uword i = 1; i < n; ++i){
    const arma::uword gi = gstart[i];
    const double invS0  = std::exp(-logS0[gi]);
    const double invS0_2 = std::exp(-2.0 * logS0[gi]);
    h[i] = h[i-1] + d[i] * invS0;
    B[i] = B[i-1] + d[i] * invS0_2;
  }


  // haz_c = exp(eta - m0)
  arma::vec eeta = arma::exp(eta - m[0]);

  // gradient and diagonal weight
  arma::vec g_z = d - eeta % h;        // n×1
  arma::vec w = eeta % h - arma::square(eeta) % B;

  // partial log-likelihood
  double ll = 0.0;
  for (arma::uword i = 0; i < n; ++i){
    const arma::uword gi = gstart[i];
    ll += d[i] * (eta[i] - logS0_full[gi]);
  }

  return CoxStatsFull{std::move(g_z), std::move(w), ll};
}
// [[Rcpp::export]]
Rcpp::List cox_stats_101(const arma::vec& eta, const arma::vec& d, const arma::vec& time) {
  auto cs = cox_suffix_stats_full(eta, d, time);
  return Rcpp::List::create(Rcpp::Named("g_z") = Rcpp::NumericVector(cs.g_z.begin(), cs.g_z.end()),
                            Rcpp::Named("w") = Rcpp::NumericVector(cs.w.begin(), cs.w.end()),
                            Rcpp::Named("loglik") = cs.loglik);
}
