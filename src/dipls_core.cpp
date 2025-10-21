// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Build X lags aligned with rows s..T-1.
// Returns a list of (T-s) x m matrices: X_s, X_{s-1}, ..., X_0.
// [[Rcpp::export]]
Rcpp::List dipls_build_Xlags_cpp(const arma::mat& X, int s) {
  const arma::uword T = X.n_rows;
  if (s < 0) {
    stop("s must be non-negative.");
  }
  if (T <= static_cast<arma::uword>(s)) {
    stop("Not enough rows for requested lag order.");
  }
  const arma::uword N = T - static_cast<arma::uword>(s);
  Rcpp::List out(s + 1);
  for (int i = 0; i <= s; ++i) {
    const arma::uword r0 = static_cast<arma::uword>(s - i);
    const arma::uword r1 = r0 + N - 1;
    out[i] = X.rows(r0, r1);
  }
  return out;
}

// Build Ts (T-s) x (s+1) with columns [t_s, t_{s-1}, ..., t_0].
// [[Rcpp::export]]
arma::mat dipls_build_Ts_cpp(const arma::vec& t_all, int s) {
  const arma::uword T = t_all.n_rows;
  if (s < 0) {
    stop("s must be non-negative.");
  }
  if (T <= static_cast<arma::uword>(s)) {
    stop("Not enough entries in t for requested lag order.");
  }
  const arma::uword N = T - static_cast<arma::uword>(s);
  arma::mat Ts(N, s + 1, arma::fill::zeros);
  for (int i = 0; i <= s; ++i) {
    const arma::uword r0 = static_cast<arma::uword>(s - i);
    const arma::uword r1 = r0 + N - 1;
    Ts.col(i) = t_all.rows(r0, r1);
  }
  return Ts;
}

// w <- sum_i beta_i * X_{s-i}^T * u_s.
// Xlags: List of (T-s) x m matrices; u_s: length T-s vector; beta: length s+1.
// [[Rcpp::export]]
arma::vec dipls_weight_w_cpp(const Rcpp::List& Xlags,
                             const arma::vec& beta,
                             const arma::vec& u_s) {
  const int s = beta.n_rows - 1;
  if (s < 0) {
    stop("beta must have at least one entry.");
  }
  arma::vec w;
  for (int i = 0; i <= s; ++i) {
    arma::mat Xi = Xlags[i];
    if (Xi.n_rows != u_s.n_rows) {
      stop("Mismatch between X lag rows and length of u_s.");
    }
    arma::vec term = Xi.t() * u_s;
    if (i == 0) {
      w = beta(i) * term;
    } else {
      w += beta(i) * term;
    }
  }
  return w;
}

// Tiny ridge LS solver: coef = (X^T X + lambda I)^{-1} X^T y.
// Returns list(coef=..., fitted=...).
// [[Rcpp::export]]
Rcpp::List dipls_ridge_ls_cpp(const arma::mat& X,
                              const arma::mat& y,
                              double lambda) {
  arma::mat XtX = X.t() * X;
  if (lambda > 0) {
    XtX.diag() += lambda;
  }
  arma::mat Xty = X.t() * y;
  arma::mat coef = arma::solve(XtX, Xty, arma::solve_opts::fast);
  arma::mat fitted = X * coef;
  return Rcpp::List::create(
    Rcpp::Named("coef") = coef,
    Rcpp::Named("fitted") = fitted
  );
}
