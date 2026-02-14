#include <Rcpp.h>
using namespace Rcpp;

// Compute lagged blocks X_i for i=1..s+1; returns a list of NumericMatrix (n x m)
// [[Rcpp::export]]
List form_blocks_cpp(const NumericMatrix& X, const int s) {
  const int T = X.nrow();
  const int m = X.ncol();
  const int n = T - s;
  if (n <= 1) stop("Time series too short for lag order s.");
  List blocks(s + 1);
  for (int i = 0; i < s + 1; ++i) {
    NumericMatrix Xi(n, m);
    // rows: from i .. i + n - 1
    for (int r = 0; r < n; ++r) {
      const int src = i + r;
      for (int c = 0; c < m; ++c) {
        Xi(r, c) = X(src, c);
      }
    }
    blocks[i] = Xi;
  }
  return blocks;
}

// multiply X^T * v, where X is (n x m), v is length n => returns length m
inline NumericVector matTvec(const NumericMatrix& X, const NumericVector& v) {
  const int n = X.nrow();
  const int m = X.ncol();
  NumericVector out(m);
  for (int j = 0; j < m; ++j) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i) {
      acc += X(i, j) * v[i];
    }
    out[j] = acc;
  }
  return out;
}

// multiply X * w, where X is (n x m), w is length m => returns length n
inline NumericVector matvec(const NumericMatrix& X, const NumericVector& w) {
  const int n = X.nrow();
  const int m = X.ncol();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < m; ++j) {
      acc += X(i, j) * w[j];
    }
    out[i] = acc;
  }
  return out;
}

// compute u = X_{s+1} w and vis[,i-1] = X_{s+1-i} w
void compute_u_vis(const List& blocks, const int s,
                   const NumericVector& w,
                   NumericVector& u, NumericMatrix& vis) {
  const NumericMatrix Xsp1 = blocks[s];
  u = matvec(Xsp1, w);
  const int n = Xsp1.nrow();
  for (int i = 1; i <= s; ++i) {
    const NumericMatrix Xi = blocks[s - i];
    NumericVector vi = matvec(Xi, w);
    for (int r = 0; r < n; ++r) vis(r, i - 1) = vi[r];
  }
}

// compute c_i = sum(u * v_i)
NumericVector compute_c(const NumericVector& u, const NumericMatrix& vis) {
  const int s = vis.ncol();
  const int n = vis.nrow();
  NumericVector c(s);
  for (int i = 0; i < s; ++i) {
    double acc = 0.0;
    for (int r = 0; r < n; ++r) acc += u[r] * vis(r, i);
    c[i] = acc;
  }
  return c;
}

// compute d = Y_beta w = 0.5 ( X_{s+1}^T sum beta_i v_i + sum beta_i X_{s+1-i}^T u )
NumericVector compute_d(const List& blocks, const int s,
                        const NumericVector& u, const NumericMatrix& vis,
                        const NumericVector& beta) {
  const NumericMatrix Xsp1 = blocks[s];
  const int n = Xsp1.nrow();
  const int m = Xsp1.ncol();

  // sum_beta_v
  NumericVector sum_beta_v(n);
  for (int i = 0; i < s; ++i) {
    const double b = beta[i];
    for (int r = 0; r < n; ++r) sum_beta_v[r] += b * vis(r, i);
  }
  NumericVector term1 = matTvec(Xsp1, sum_beta_v);

  // term2 = sum beta_i * X_{s+1-i}^T u
  NumericVector term2(m);
  for (int i = 1; i <= s; ++i) {
    const NumericMatrix Xi = blocks[s - i];
    NumericVector tmp = matTvec(Xi, u);
    const double b = beta[i - 1];
    for (int j = 0; j < m; ++j) term2[j] += b * tmp[j];
  }

  NumericVector d(m);
  for (int j = 0; j < m; ++j) d[j] = 0.5 * (term1[j] + term2[j]);
  return d;
}

// matvec for Y_beta: given v, returns Y_beta v without forming matrices
NumericVector ybeta_matvec(const List& blocks, const int s,
                           const NumericVector& beta, const NumericVector& v) {
  const NumericMatrix Xsp1 = blocks[s];
  const int n = Xsp1.nrow();

  // u2, vis2
  NumericVector u2(n);
  NumericMatrix vis2(n, s);
  compute_u_vis(blocks, s, v, u2, vis2);
  // d2
  NumericVector d2 = compute_d(blocks, s, u2, vis2, beta);
  return d2;
}

// One DiPCA component core loop in C++
// [[Rcpp::export]]
List dipca_component_cpp(const NumericMatrix& X, const int s,
                         NumericVector w0, String algorithm,
                         const double tol, const int max_iter,
                         const int inner_power, const double inner_tol,
                         const int verbose) {
  // Build blocks
  List blocks = form_blocks_cpp(X, s);
  const NumericMatrix Xsp1 = blocks[s];
  const int n = Xsp1.nrow();
  const int m = Xsp1.ncol();

  // init
  NumericVector w = clone(w0);
  // normalize w
  double nrm = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
  if (nrm < 1e-15) stop("w0 norm too small.");
  for (int j = 0; j < m; ++j) w[j] /= nrm;

  NumericVector beta(s);
  for (int i = 0; i < s; ++i) beta[i] = 1.0 / std::sqrt((double)s);

  NumericVector u(n);
  NumericMatrix vis(n, s);
  NumericVector c(s);
  NumericVector d(m);

  double lam = 0.0;
  int it = 0;
  std::vector<double> obj_hist;

  for (it = 1; it <= max_iter; ++it) {
    // u, vis
    compute_u_vis(blocks, s, w, u, vis);
    // c and beta update
    c = compute_c(u, vis);
    double lam_new = 0.0;
    for (int i = 0; i < s; ++i) lam_new += c[i] * c[i];
    lam_new = std::sqrt(lam_new);
    if (lam_new <= 1e-15) {
      for (int i = 0; i < s; ++i) beta[i] = 1.0 / std::sqrt((double)s);
    } else {
      for (int i = 0; i < s; ++i) beta[i] = c[i] / lam_new;
    }

    // w update
    if (algorithm == "I") {
      d = compute_d(blocks, s, u, vis, beta);
      // normalize
      double nd = std::sqrt(std::inner_product(d.begin(), d.end(), d.begin(), 0.0));
      if (nd > 1e-15) {
        for (int j = 0; j < m; ++j) w[j] = d[j] / nd;
      } else {
        // keep w unchanged if degenerate
      }
    } else {
      // Algorithm II: power iteration on matvec
      NumericVector v = clone(w);
      for (int k = 0; k < inner_power; ++k) {
        NumericVector Av = ybeta_matvec(blocks, s, beta, v);
        double nv = std::sqrt(std::inner_product(Av.begin(), Av.end(), Av.begin(), 0.0));
        if (nv <= 1e-15) break;
        for (int j = 0; j < m; ++j) Av[j] /= nv;
        // check convergence
        double diff = 0.0;
        for (int j = 0; j < m; ++j) {
          double dtmp = Av[j] - v[j];
          diff += dtmp * dtmp;
        }
        v = Av;
        if (std::sqrt(diff) < inner_tol) break;
      }
      w = v;
      // set d for residual
      d = ybeta_matvec(blocks, s, beta, w);
    }

    // recompute residual and objective with updated w
    compute_u_vis(blocks, s, w, u, vis);
    c = compute_c(u, vis);
    lam = 0.0;
    for (int i = 0; i < s; ++i) lam += c[i] * c[i];
    lam = std::sqrt(lam);
    d = compute_d(blocks, s, u, vis, beta);

    // residual inf-norm
    double res = 0.0;
    for (int j = 0; j < m; ++j) {
      double tmp = d[j] - lam * w[j];
      double a = std::abs(tmp);
      if (a > res) res = a;
    }
    obj_hist.push_back(lam);
    if (verbose && (it % 50 == 0 || res < tol)) {
      Rcpp::Rcout << "[dipca] it=" << it << " obj=" << lam << " res=" << res << std::endl;
    }
    if (res < tol) break;
  }

  // return
  NumericVector obj_hist_R(obj_hist.size());
  for (size_t i = 0; i < obj_hist.size(); ++i) obj_hist_R[i] = obj_hist[i];

  return List::create(
    _["w"] = w,
    _["beta"] = beta,
    _["obj"] = lam,
    _["iters"] = it,
    _["obj_hist"] = obj_hist_R
  );
}

// ------------ DiCCA core (one component) ------------------------------------

// compute crossprod X^T X into a (m x m) matrix
inline NumericMatrix crossprod_matrix(const NumericMatrix& X) {
  const int n = X.nrow();
  const int m = X.ncol();
  NumericMatrix out(m, m);
  for (int j = 0; j < m; ++j) {
    for (int k = j; k < m; ++k) {
      double acc = 0.0;
      for (int i = 0; i < n; ++i) acc += X(i, j) * X(i, k);
      out(j, k) = acc;
      if (k != j) out(k, j) = acc;
    }
  }
  return out;
}

// add matrix B into A in place (A += B)
inline void add_inplace(NumericMatrix& A, const NumericMatrix& B) {
  const int r = A.nrow();
  const int c = A.ncol();
  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < c; ++j) A(i, j) += B(i, j);
  }
}

// matrix-vector product Ts * v where Ts is (n x s), v length s
inline NumericVector matvec_Ts(const NumericMatrix& Ts, const NumericVector& v) {
  const int n = Ts.nrow();
  const int s = Ts.ncol();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < s; ++j) acc += Ts(i, j) * v[j];
    out[i] = acc;
  }
  return out;
}

// t(Ts) * v -> length s
inline NumericVector t_Ts_vec(const NumericMatrix& Ts, const NumericVector& v) {
  const int n = Ts.nrow();
  const int s = Ts.ncol();
  NumericVector out(s);
  for (int j = 0; j < s; ++j) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i) acc += Ts(i, j) * v[i];
    out[j] = acc;
  }
  return out;
}

// build Ts (n x s): columns t_s ... t_1 from blocks and w
inline void build_Ts(const List& blocks, const int s,
                     const NumericVector& w, NumericMatrix& Ts) {
  const NumericMatrix Xsp1 = blocks[s];
  const int n = Xsp1.nrow();
  for (int i = 0; i < s; ++i) {
    const NumericMatrix Xi = blocks[s - 1 - i];
    NumericVector col = matvec(Xi, w);
    for (int r = 0; r < n; ++r) Ts(r, i) = col[r];
  }
}

// Gram G = t(Ts) Ts (s x s)
inline NumericMatrix gram_Ts(const NumericMatrix& Ts) {
  const int n = Ts.nrow();
  const int s = Ts.ncol();
  NumericMatrix G(s, s);
  for (int j = 0; j < s; ++j) {
    for (int k = j; k < s; ++k) {
      double acc = 0.0;
      for (int i = 0; i < n; ++i) acc += Ts(i, j) * Ts(i, k);
      G(j, k) = acc;
      if (k != j) G(k, j) = acc;
    }
  }
  return G;
}

// robust solve(A, b) via base::solve with small ridge fallback
inline NumericVector solve_linear(const NumericMatrix& A_in, const NumericVector& b) {
  using namespace Rcpp;
  Environment base = Environment::namespace_env("base");
  Function solve = base["solve"];
  NumericMatrix A = clone(A_in);
  SEXP res;
  try {
    res = solve(A, b);
  } catch (...) {
    // add small ridge and retry
    const int m = A.nrow();
    for (int i = 0; i < m; ++i) A(i, i) += 1e-10;
    res = solve(A, b);
  }
  return as<NumericVector>(res);
}

// [[Rcpp::export]]
List dicca_component_cpp(const NumericMatrix& X, const int s,
                         NumericVector w0, const double tol,
                         const int max_iter, const int verbose) {
  // lag blocks
  List blocks = form_blocks_cpp(X, s);
  const NumericMatrix Xsp1 = blocks[s];
  const int n = Xsp1.nrow();
  const int m = Xsp1.ncol();

  NumericVector w = clone(w0);
  double nrm = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
  if (nrm < 1e-15) stop("w0 norm too small.");
  for (int j = 0; j < m; ++j) w[j] /= nrm;

  NumericVector beta(s);
  double J = 0.0;
  int it = 0;
  std::vector<double> hist;

  for (it = 1; it <= max_iter; ++it) {
    // ts1 and Ts
    NumericVector ts1 = matvec(Xsp1, w); // (n)
    NumericMatrix Ts(n, s);
    build_Ts(blocks, s, w, Ts);

    // LS for beta: (Ts^T Ts) beta = Ts^T ts1
    NumericMatrix G = gram_Ts(Ts);
    NumericVector rhs = t_Ts_vec(Ts, ts1);
    NumericVector beta_ls = solve_linear(G, rhs);
    for (int i = 0; i < s; ++i) beta[i] = beta_ls[i];

    // normalization beta := beta / sqrt(ts1^T Ts beta)
    NumericVector ts_beta = matvec_Ts(Ts, beta);
    double denom = 0.0; for (int i = 0; i < n; ++i) denom += ts1[i] * ts_beta[i];
    denom = std::sqrt(std::max(1e-15, denom));
    for (int i = 0; i < s; ++i) beta[i] /= denom;

    // Xbeta = sum_i beta_i X_{s - i}
    NumericMatrix Xbeta(n, m);
    for (int i = 0; i < s; ++i) {
      const NumericMatrix Xi = blocks[s - 1 - i];
      const double b = beta[i];
      for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) Xbeta(r, c) += b * Xi(r, c);
      }
    }

    // A = Xsp1^T Xsp1 + Xbeta^T Xbeta
    NumericMatrix A = crossprod_matrix(Xsp1);
    NumericMatrix B = crossprod_matrix(Xbeta);
    add_inplace(A, B);

    // b = Xsp1^T (Ts beta) + Xbeta^T ts1
    NumericVector term1 = matTvec(Xsp1, ts_beta);
    NumericVector term2 = matTvec(Xbeta, ts1);
    NumericVector b(m);
    for (int j2 = 0; j2 < m; ++j2) b[j2] = term1[j2] + term2[j2];

    // solve and normalize w
    NumericVector w_new = solve_linear(A, b);
    double nrmw = std::sqrt(std::inner_product(w_new.begin(), w_new.end(), w_new.begin(), 0.0));
    if (nrmw > 1e-15) for (int j2 = 0; j2 < m; ++j2) w_new[j2] /= nrmw;

    // objective J = ts1_new^T (Ts_new beta)
    NumericVector ts1_new = matvec(Xsp1, w_new);
    NumericMatrix Ts_new(n, s);
    build_Ts(blocks, s, w_new, Ts_new);
    NumericVector tsb_new = matvec_Ts(Ts_new, beta);
    J = 0.0; for (int i2 = 0; i2 < n; ++i2) J += ts1_new[i2] * tsb_new[i2];
    hist.push_back(J);

    // convergence on ||w_new - w||_inf
    double diff = 0.0;
    for (int j2 = 0; j2 < m; ++j2) {
      double d = std::abs(w_new[j2] - w[j2]);
      if (d > diff) diff = d;
    }
    if (verbose && (it % 50 == 0 || diff < tol)) {
      Rcpp::Rcout << "[dicca] it=" << it << "  J=" << J << "  |dw|_inf=" << diff << std::endl;
    }
    w = w_new;
    if (diff < tol) break;
  }

  // R2 using final w, beta
  NumericVector ts1f = matvec(Xsp1, w);
  NumericMatrix Tsf(n, s);
  build_Ts(blocks, s, w, Tsf);
  NumericVector that = matvec_Ts(Tsf, beta);
  // center
  double mean_a = 0.0, mean_c = 0.0;
  for (int i = 0; i < n; ++i) { mean_a += ts1f[i]; mean_c += that[i]; }
  mean_a /= n; mean_c /= n;
  double num = 0.0, da = 0.0, dc = 0.0;
  for (int i = 0; i < n; ++i) {
    double aa = ts1f[i] - mean_a;
    double cc = that[i] - mean_c;
    num += aa * cc; da += aa * aa; dc += cc * cc;
  }
  double R2 = 0.0;
  double denomR = std::sqrt(da * dc);
  if (denomR > 1e-15) { double r = num / denomR; R2 = r * r; }

  NumericVector obj_hist_R(hist.size());
  for (size_t i = 0; i < hist.size(); ++i) obj_hist_R[i] = hist[i];

  return List::create(
    _["w"] = w,
    _["beta"] = beta,
    _["obj"] = J,
    _["iters"] = it,
    _["R2"] = R2,
    _["obj_hist"] = obj_hist_R
  );
}
