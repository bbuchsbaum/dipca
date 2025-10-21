// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// --- finite-difference with zero padding to keep n_rows --- //
static inline arma::mat diff_pad(const arma::mat& X, int d) {
  if (d <= 0) return X;
  arma::mat Y = X;
  const int n = X.n_rows;
  const int p = X.n_cols;
  for (int k = 0; k < d; ++k) {
    arma::mat Z(n, p, arma::fill::zeros);
    for (int j = 0; j < p; ++j) {
      double prev = 0.0;
      for (int i = 0; i < n; ++i) {
        const double curr = Y(i, j);
        Z(i, j) = curr - prev;  // first difference with zero IC
        prev = curr;
      }
    }
    Y = Z;
  }
  return Y;
}

// [[Rcpp::export]]
arma::mat arma_inverse_filter_cpp(const arma::mat& U,
                                  const arma::vec& ar,     // length p (a1..ap)
                                  const arma::vec& ma,     // length q (b1..bq)
                                  const int d) {
  // Implements:  y = H^{-1}(B) x  with  H^{-1}(B) = (1 - a(B)) * (1 - B)^d / (1 + b(B))
  // first difference d times, then ARMA inverse recursion.
  arma::mat X = diff_pad(U, d);
  const int n = X.n_rows, p = X.n_cols;
  const int P = static_cast<int>(ar.n_elem);
  const int Q = static_cast<int>(ma.n_elem);

  arma::vec numer(P + 1); numer(0) = 1.0;
  for (int i = 0; i < P; ++i) numer(i + 1) = -ar(i);  // 1 - sum a_j B^j

  arma::vec denom(Q + 1); denom(0) = 1.0;
  for (int j = 0; j < Q; ++j) denom(j + 1) =  ma(j);  // 1 + sum b_j B^j

  arma::mat Y(n, p, arma::fill::zeros);

  for (int col = 0; col < p; ++col) {
    for (int t = 0; t < n; ++t) {
      // numerator convolution on X
      double accum = 0.0;
      for (int i = 0; i <= P; ++i) {
        const int idx = t - i;
        if (idx >= 0) accum += numer(i) * X(idx, col);
      }
      // denominator recursion on Y
      for (int j = 1; j <= Q; ++j) {
        const int yidx = t - j;
        if (yidx >= 0) accum -= denom(j) * Y(yidx, col);
      }
      Y(t, col) = accum;
    }
  }
  return Y;
}

// [[Rcpp::export]]
arma::vec smallest_eigenvector_crossprod_cpp(const arma::mat& U) {
  arma::mat AtA = U.t() * U; // symmetric PSD
  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, AtA))
    Rcpp::stop("eig_sym failed in smallest_eigenvector_crossprod_cpp");

  arma::uword idx = eigval.index_min();
  arma::vec w = eigvec.col(idx);
  double nrm = arma::norm(w, 2);
  if (nrm > 0) w /= nrm;
  return w;
}

