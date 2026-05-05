// LRnetSTv2: DAG learning for spatial transcriptomics data.
// Improvements over LRnetST v2.0:
//   - Separate add/delete/reverse score caches with version stamps (no stale entries)
//   - Fixed acyclicUpdate double-if bug in reverse block (if/else if)
//   - mt19937 RNG replaces global srand/rand
//   - indRow precomputed once before HC loop (fixed for each node across all steps)
//   - Loglik-only internal score functions (no unnecessary coefficient/residual output)
//   - Pre-allocated Eigen workspace matrix (no per-call heap allocation in hot path)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"

#include <RcppEigen.h>
#include <RcppNumerical.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Numer;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Ref<const MatrixXd> ConstMatRef;
typedef Eigen::Ref<const VectorXd> ConstVecRef;
typedef std::vector<std::vector<int>> AdjList;

namespace {

const double kInf = std::numeric_limits<double>::infinity();

inline int mat_index(int from, int to, int p) { return from + to * p; }

// ---- Graph utilities ----

// BFS: does adding fromNode->toNode create a cycle?
// Returns true iff toNode is already an ancestor of fromNode.
bool edge_on_loop_cpp(int fromNode, int toNode, const AdjList& parents) {
  const int p = static_cast<int>(parents.size());
  std::vector<char> seen(p, 0);
  std::vector<int> queue;
  queue.reserve(p);
  queue.push_back(fromNode);
  seen[fromNode] = 1;
  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const int cur = queue[pos];
    if (cur == toNode) return true;
    for (int par : parents[cur]) {
      if (!seen[par]) { seen[par] = 1; queue.push_back(par); }
    }
  }
  return false;
}

bool acyclic_add(int from, int to, const AdjList& parents) {
  return !edge_on_loop_cpp(from, to, parents);
}

// Check reverse from->to by temporarily removing from from parents[to].
// Restores parents[to] before returning. In-place to avoid full AdjList copy.
bool acyclic_reverse_inplace(int from, int to, AdjList& parents) {
  std::vector<int>& toParents = parents[to];
  const auto it = std::find(toParents.begin(), toParents.end(), from);
  if (it == toParents.end()) return false;
  const std::ptrdiff_t pos = it - toParents.begin();
  toParents.erase(it);
  const bool result = !edge_on_loop_cpp(to, from, parents);
  toParents.insert(toParents.begin() + pos, from);
  return result;
}

AdjList graph_to_parents(const Rcpp::LogicalMatrix& g, int p) {
  AdjList out(p);
  for (int to = 0; to < p; ++to)
    for (int from = 0; from < p; ++from)
      if (g(from, to)) out[to].push_back(from);
  return out;
}

AdjList graph_to_children(const Rcpp::LogicalMatrix& g, int p) {
  AdjList out(p);
  for (int from = 0; from < p; ++from)
    for (int to = 0; to < p; ++to)
      if (g(from, to)) out[from].push_back(to);
  return out;
}

std::vector<char> compute_ancestors(int i, const AdjList& parents) {
  const int p = static_cast<int>(parents.size());
  std::vector<char> out(p, 0);
  std::vector<int> q;
  q.reserve(p);
  out[i] = 1; q.push_back(i);
  for (std::size_t pos = 0; pos < q.size(); ++pos)
    for (int par : parents[q[pos]])
      if (!out[par]) { out[par] = 1; q.push_back(par); }
  return out;
}

std::vector<char> compute_descendants(int i, const AdjList& children) {
  const int p = static_cast<int>(children.size());
  std::vector<char> out(p, 0);
  std::vector<int> q;
  q.reserve(p);
  out[i] = 1; q.push_back(i);
  for (std::size_t pos = 0; pos < q.size(); ++pos)
    for (int chd : children[q[pos]])
      if (!out[chd]) { out[chd] = 1; q.push_back(chd); }
  return out;
}

void erase_value(std::vector<int>& x, int value) {
  const auto it = std::find(x.begin(), x.end(), value);
  if (it != x.end()) x.erase(it);
}

// ---- Design matrix filling ----

// Fill workspace columns from Y (all n rows).
// Layout: [extraParent if >=0] [existing parents except dropParent] [ones]
// Returns number of columns written.
int fill_design(MatrixXd& ws,
                const Rcpp::NumericMatrix& Y,
                const std::vector<int>& parents,
                int extraParent = -1,
                int dropParent  = -1) {
  const int n = Y.nrow();
  int col = 0;
  if (extraParent >= 0) {
    ws.col(col) = Eigen::Map<const VectorXd>(&Y(0, extraParent), n);
    ++col;
  }
  for (int par : parents) {
    if (par == dropParent) continue;
    ws.col(col) = Eigen::Map<const VectorXd>(&Y(0, par), n);
    ++col;
  }
  ws.col(col).head(n).setOnes();
  ++col;
  return col;
}

// Fill workspace with row-selected columns (continuous nodes: only rows where Y[r,i]!=0).
// Result occupies the first n_sel rows; outer stride is still ws.outerStride() = n.
//
// Constant columns are silently dropped. This handles the SC-specific case where a
// binary indicator parent (whitelisted as indicator_i -> logcount_i) equals 1 on every
// selected row (because selection is logcount_i != 0 ↔ indicator_i == 1), making it
// perfectly collinear with the intercept. Dropping it keeps XtX full-rank.
//
// Returns number of columns written (including the intercept, excluding dropped constants).
int fill_design_sel(MatrixXd& ws,
                    const Rcpp::NumericMatrix& Y,
                    const std::vector<int>& parents,
                    const std::vector<char>& indRow,
                    int n_sel,
                    int extraParent = -1,
                    int dropParent  = -1) {
  const int n = Y.nrow();
  int col = 0;
  // Fill one source column into ws[:,col] using only selected rows.
  // If all selected values are identical (constant column), skip it — constant columns
  // are linearly dependent with the intercept and make XtX singular.
  auto fill_col_if_nonconstant = [&](int src) {
    int r_out = 0;
    double first_val = 0.0;
    bool first_set = false;
    bool is_constant = true;
    for (int r = 0; r < n; ++r) {
      if (!indRow[r]) continue;
      const double v = Y(r, src);
      ws(r_out++, col) = v;
      if (!first_set) { first_val = v; first_set = true; }
      else if (v != first_val) is_constant = false;
    }
    if (!is_constant) ++col;  // keep: advance; otherwise overwrite slot next iteration
  };
  if (extraParent >= 0) fill_col_if_nonconstant(extraParent);
  for (int par : parents) {
    if (par == dropParent) continue;
    fill_col_if_nonconstant(par);
  }
  ws.block(0, col, n_sel, 1).setOnes();
  ++col;
  return col;
}

// ---- Likelihood functions ----

// Gaussian MLE log-likelihood: beta solved via Cholesky of X'X.
// sigma2_hat = ||Y - X*beta||^2 / n  (MLE divisor).
// Returns -Inf on singular X'X, non-finite beta, or sigma2 <= 0.
double continuous_loglik(const ConstMatRef& X, const ConstVecRef& Y) {
  const int n = X.rows();
  const int q = X.cols();
  // X'X via rankUpdate fills only the lower triangle; LLT uses Lower by default.
  MatrixXd XtX(q, q);
  XtX.setZero();
  XtX.selfadjointView<Lower>().rankUpdate(X.adjoint());
  LLT<MatrixXd> llt(XtX);
  if (llt.info() != Eigen::Success) return -kInf;
  const VectorXd beta = llt.solve(X.adjoint() * Y);
  if (!beta.allFinite()) return -kInf;
  const double s2 = (Y - X * beta).squaredNorm() / static_cast<double>(n);
  if (!std::isfinite(s2) || s2 <= 0.0) return -kInf;
  return -static_cast<double>(n) / 2.0 *
         (std::log(2.0 * M_PI) + std::log(s2) + 1.0);
}

// Logistic regression negative log-likelihood via L-BFGS (from RcppNumerical).
class LogisticReg : public MFuncGrad {
 private:
  const ConstMatRef X;
  const ConstVecRef Y;
  const int n;
  VectorXd xbeta, prob;
 public:
  LogisticReg(const ConstMatRef& x_, const ConstVecRef& y_)
      : X(x_), Y(y_), n(X.rows()), xbeta(n), prob(n) {}
  double f_grad(Constvec& beta, Refvec grad) {
    xbeta.noalias() = X * beta;
    const double yxbeta = Y.dot(xbeta);
    for (int i = 0; i < n; ++i) prob[i] = R::log1pexp(xbeta[i]);
    const double f = prob.sum() - yxbeta;
    prob = (xbeta - prob).array().exp();
    grad.noalias() = X.transpose() * (prob - Y);
    return f;
  }
};

double logistic_loglik(const ConstMatRef& X, const ConstVecRef& Y) {
  const int q = X.cols();
  LogisticReg nll(X, Y);
  VectorXd beta = VectorXd::Zero(q);
  double fopt = kInf;
  const int status = optim_lbfgs(nll, beta, fopt, 300, 1e-8, 1e-5);
  if (status < 0 || !std::isfinite(fopt)) return -kInf;
  return -fopt;
}

// BIC = -2*logL + log(n_eff) * (q-1).
// q = total columns (parents + intercept); q-1 = structural parameters.
// For continuous nodes: n_eff = n_sel (nonzero rows); for binary: n_eff = n.
// extraParent / dropParent describe the candidate operation being scored.
//
// topLeftCorner(n_sel, q) has non-unit outer stride (= ws.rows() = n), so
// continuous_loglik receives an explicit MatrixXd copy — cost negligible vs Cholesky.
double bic_score_SC(MatrixXd& ws,
                    VectorXd& y_buf,
                    const Rcpp::NumericMatrix& Y,
                    const std::vector<int>& parents,
                    const std::string& nodeType,
                    int node,
                    const std::vector<char>& indRow,
                    int n_sel,
                    int extraParent = -1,
                    int dropParent  = -1) {
  if (nodeType == "c") {
    const int q = fill_design_sel(ws, Y, parents, indRow, n_sel,
                                  extraParent, dropParent);
    int r_out = 0;
    const int n = Y.nrow();
    for (int r = 0; r < n; ++r)
      if (indRow[r]) y_buf[r_out++] = Y(r, node);
    // Explicit copy: topLeftCorner has outer stride n != n_sel.
    const MatrixXd X(ws.topLeftCorner(n_sel, q));
    const double loglik = continuous_loglik(X, y_buf.head(n_sel));
    if (!std::isfinite(loglik)) return kInf;
    return -2.0 * loglik + std::log(static_cast<double>(n_sel)) * (q - 1);
  } else {
    const int n = Y.nrow();
    const int q = fill_design(ws, Y, parents, extraParent, dropParent);
    Eigen::Map<const VectorXd> y(&Y(0, node), n);
    const double loglik = logistic_loglik(ws.leftCols(q), y);
    if (!std::isfinite(loglik)) return kInf;
    return -2.0 * loglik + std::log(static_cast<double>(n)) * (q - 1);
  }
}

// ---- Score cache ----

// Cache entry for one (from, to) pair. Semantics by operation type:
//   Add from->to:     scoreA = BIC(to after add),   versionA = version[to]
//   Delete from->to:  scoreA = BIC(to after delete), versionA = version[to]
//   Reverse from->to: scoreA = BIC(to after delete), versionA = version[to]
//                     scoreB = BIC(from after add),  versionB = version[from]
// Default versionA/versionB = -1 ensures cache miss on first access
// (version[] initializes to 0).
struct OneCache {
  double value, scoreA, scoreB;
  int versionA, versionB;
  OneCache() : value(kInf), scoreA(kInf), scoreB(kInf), versionA(-1), versionB(-1) {}
};

bool better_delta(double candidate, double currentBest, double tol,
                  bool flip, std::mt19937& rng) {
  if (!std::isfinite(candidate)) return false;
  if (candidate < currentBest - tol) return true;
  if (candidate > currentBest + tol) return false;
  if (flip) {
    std::uniform_int_distribution<int> coin(0, 1);
    return coin(rng) == 1;
  }
  return false;
}

// ---- Acyclicity cache update (fixed) ----

struct LastOpState {
  int from = -1, to = -1, type = -1;
  std::vector<char> fromAn, fromDe, toAn, toDe;
};

// Incrementally update acyStatus for one candidate given the last accepted op.
//
//   acyStatus(from, to) = true  iff adding from->to is acyclic
//   acyStatus(to, from) = true  iff reversing from->to is acyclic
//   NA_LOGICAL           = not yet computed
//
// The lastOper==3 block uses if/else if to avoid both conditions firing when
// the same pair is simultaneously in the "new cycle path" and "freed path" sets.
// (Original LRnetST had two independent if blocks, which could overwrite a
// freshly-set false with a re-check that might return true — incorrect.)
void acyclic_cache_update(const LastOpState& last,
                          int operFrom, int operTo, int operType,
                          AdjList& parents,
                          Rcpp::LogicalMatrix& acyStatus) {
  const int i = operFrom, j = operTo;

  if (last.type == 1) {
    // Add created a new path; some previously-acyclic candidates may now cycle.
    if (operType == 1 && acyStatus(i,j) == true &&
        last.toDe[i] && last.fromAn[j])
      acyStatus(i,j) = false;
    if (operType == 3 && acyStatus(j,i) == true &&
        last.toDe[j] && last.fromAn[i])
      if (i != last.from || j != last.to)
        acyStatus(j,i) = false;

  } else if (last.type == 2) {
    // Delete removed a path; some previously-cyclic candidates may now be acyclic.
    if (operType == 1 && acyStatus(i,j) == false &&
        last.toDe[i] && last.fromAn[j])
      acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
    if (operType == 3 && acyStatus(j,i) == false &&
        last.toDe[j] && last.fromAn[i])
      acyStatus(j,i) = acyclic_reverse_inplace(i, j, parents) ? true : false;

  } else if (last.type == 3) {
    // Reverse both creates a new path and removes an old path; handle with
    // if/else if so that a pair in both condition sets is resolved by BFS once.
    if (operType == 1) {
      if (acyStatus(i,j) == true && last.fromDe[i] && last.toAn[j]) {
        if (last.toDe[i] && last.fromAn[j])
          acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
        else
          acyStatus(i,j) = false;
      } else if (acyStatus(i,j) == false && last.toDe[i] && last.fromAn[j]) {
        acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
      }
    }
    if (operType == 3) {
      const bool isReversedEdge = (i == last.to && j == last.from);
      if (!isReversedEdge && acyStatus(j,i) == true &&
          last.fromDe[j] && last.toAn[i]) {
        if (last.toDe[j] && last.fromAn[i])
          acyStatus(j,i) = acyclic_reverse_inplace(i, j, parents) ? true : false;
        else
          acyStatus(j,i) = false;
      } else if (acyStatus(j,i) == false && last.toDe[j] && last.fromAn[i]) {
        acyStatus(j,i) = acyclic_reverse_inplace(i, j, parents) ? true : false;
      }
    }
  }
}

// ---- Score aggregation helpers (for score_shd_cpp / score_shd_freq_cpp) ----

inline bool has_edge(const std::vector<unsigned char>& g, int from, int to, int p) {
  return g[from + to * p] != 0;
}
inline void set_edge(std::vector<unsigned char>& g, int from, int to, int p, bool v) {
  g[from + to * p] = v ? 1 : 0;
}

std::vector<unsigned char> logical_to_uc(const Rcpp::LogicalMatrix& x, const char* name) {
  const int p = x.nrow();
  if (x.ncol() != p) Rcpp::stop("%s must be square", name);
  std::vector<unsigned char> out(p * p, 0);
  for (int to = 0; to < p; ++to)
    for (int from = 0; from < p; ++from) {
      const int val = x(from, to);
      if (val == NA_LOGICAL) Rcpp::stop("%s must not contain NA values", name);
      out[from + to * p] = val ? 1 : 0;
    }
  return out;
}

AdjList parents_from_uc(const std::vector<unsigned char>& g, int p) {
  AdjList par(p);
  for (int to = 0; to < p; ++to)
    for (int from = 0; from < p; ++from)
      if (has_edge(g, from, to, p)) par[to].push_back(from);
  return par;
}

bool is_dag_uc(const std::vector<unsigned char>& g, int p) {
  std::vector<int> indeg(p, 0);
  for (int to = 0; to < p; ++to)
    for (int from = 0; from < p; ++from)
      if (has_edge(g, from, to, p)) ++indeg[to];
  std::vector<int> q;
  q.reserve(p);
  for (int i = 0; i < p; ++i)
    if (indeg[i] == 0) q.push_back(i);
  int seen = 0;
  for (int qi = 0; qi < (int)q.size(); ++qi) {
    int node = q[qi]; ++seen;
    for (int to = 0; to < p; ++to)
      if (has_edge(g, node, to, p) && --indeg[to] == 0) q.push_back(to);
  }
  return seen == p;
}

Rcpp::IntegerMatrix wrap_int_graph(const std::vector<unsigned char>& g, int p) {
  Rcpp::IntegerMatrix out(p, p);
  for (int to = 0; to < p; ++to)
    for (int from = 0; from < p; ++from)
      out(from, to) = has_edge(g, from, to, p) ? 1 : 0;
  return out;
}

struct EdgeCandidate { int from, to; double sf, gsf; };

Rcpp::IntegerMatrix aggregate_freq(const Rcpp::NumericMatrix& seleFreq,
                                   double alpha, double freqCutoff,
                                   const Rcpp::LogicalMatrix& whiteList,
                                   const Rcpp::LogicalMatrix& blackList,
                                   bool verbose) {
  const int p = seleFreq.nrow();
  if (seleFreq.ncol() != p)  Rcpp::stop("freq must be a square matrix");
  if (p > 46000)             Rcpp::stop("p exceeds safe index limit (46000)");
  if (whiteList.nrow() != p || whiteList.ncol() != p ||
      blackList.nrow() != p || blackList.ncol() != p)
    Rcpp::stop("whiteList and blackList must be p by p matrices");
  if (!std::isfinite(alpha) || alpha <= 0.0)
    Rcpp::stop("alpha must be a positive finite scalar");
  if (!std::isfinite(freqCutoff) || freqCutoff < 0.0 || freqCutoff > 1.0)
    Rcpp::stop("freqCutoff must be a finite scalar in [0, 1]");

  std::vector<unsigned char> white = logical_to_uc(whiteList, "whiteList");
  std::vector<unsigned char> black = logical_to_uc(blackList, "blackList");
  for (int i = 0; i < p; ++i) {
    set_edge(black, i, i, p, true);
    if (has_edge(white, i, i, p))
      Rcpp::stop("whiteList diagonal must be FALSE");
    for (int j = 0; j < p; ++j) {
      if (has_edge(white, i, j, p) && has_edge(black, i, j, p))
        Rcpp::stop("whiteList and blackList conflict");
      if (i != j && has_edge(white, i, j, p) && has_edge(white, j, i, p))
        Rcpp::stop("whiteList cannot contain both directions of an edge");
      const double v = seleFreq(i, j);
      if (!std::isfinite(v) || v < 0.0 || v > 1.0)
        Rcpp::stop("selection frequencies must be finite values in [0, 1]");
    }
  }
  if (!is_dag_uc(white, p)) Rcpp::stop("whiteList must be acyclic");

  std::vector<EdgeCandidate> cands;
  cands.reserve(p * p);
  for (int to = 0; to < p; ++to)
    for (int from = 0; from < p; ++from) {
      if (from == to) continue;
      const double sf  = seleFreq(from, to);
      const double gsf = sf + (1.0 - alpha / 2.0) * seleFreq(to, from);
      if (gsf > freqCutoff && !has_edge(white, from, to, p)) {
        EdgeCandidate c; c.from = from; c.to = to; c.sf = sf; c.gsf = gsf;
        cands.push_back(c);
      }
    }

  std::sort(cands.begin(), cands.end(), [](const EdgeCandidate& a, const EdgeCandidate& b) {
    if (a.gsf != b.gsf) return a.gsf > b.gsf;
    if (a.sf  != b.sf)  return a.sf  > b.sf;
    if (a.from != b.from) return a.from < b.from;
    return a.to < b.to;
  });

  std::vector<unsigned char> result = white;
  AdjList parents = parents_from_uc(result, p);
  for (std::size_t k = 0; k < cands.size(); ++k) {
    const int from = cands[k].from, to = cands[k].to;
    if (verbose)
      Rcpp::Rcout << (k + 1) << "th operation: " << (from + 1) << " -> " << (to + 1) << "\n";
    if (has_edge(black, from, to, p) || has_edge(result, to, from, p)) continue;
    if (!edge_on_loop_cpp(from, to, parents)) {
      set_edge(result, from, to, p, true);
      parents[to].push_back(from);
    } else if (verbose) {
      Rcpp::Rcout << "not pass acyclic check\n";
    }
  }
  return wrap_int_graph(result, p);
}

} // namespace

// ---- Exported: score_shd_freq_cpp ----
// [[Rcpp::export]]
Rcpp::IntegerMatrix score_shd_freq_cpp(const Rcpp::NumericMatrix& freq,
                                        double alpha, double freqCutoff,
                                        const Rcpp::LogicalMatrix& whiteList,
                                        const Rcpp::LogicalMatrix& blackList,
                                        bool verbose = false) {
  Rcpp::NumericMatrix cleanFreq = Rcpp::clone(freq);
  const int p = cleanFreq.nrow();
  if (cleanFreq.ncol() != p) Rcpp::stop("freq must be a square matrix");
  for (int i = 0; i < p; ++i) cleanFreq(i, i) = 0.0;
  return aggregate_freq(cleanFreq, alpha, freqCutoff, whiteList, blackList, verbose);
}

// ---- Exported: score_shd_cpp ----
// [[Rcpp::export]]
Rcpp::IntegerMatrix score_shd_cpp(const Rcpp::NumericVector& bootAdj,
                                   double alpha, double freqCutoff,
                                   const Rcpp::LogicalMatrix& whiteList,
                                   const Rcpp::LogicalMatrix& blackList,
                                   bool verbose = false) {
  Rcpp::IntegerVector dims = bootAdj.attr("dim");
  if (dims.size() != 3) Rcpp::stop("boot.adj must be a p by p by B array");
  const int p = dims[0], p2 = dims[1], nb = dims[2];
  if (p <= 0 || p2 != p || nb <= 0)
    Rcpp::stop("boot.adj must be a p by p by B array with B >= 1");
  Rcpp::NumericMatrix freq(p, p);
  const int sliceSize = p * p;
  for (int b = 0; b < nb; ++b)
    for (int to = 0; to < p; ++to)
      for (int from = 0; from < p; ++from) {
        const double v = bootAdj[from + to * p + b * sliceSize];
        if (!std::isfinite(v) || (v != 0.0 && v != 1.0))
          Rcpp::stop("boot.adj entries must be finite 0/1 values");
        freq(from, to) += v;
      }
  const double nb_d = static_cast<double>(nb);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) freq(i, j) /= nb_d;
    freq(i, i) = 0.0;
  }
  return aggregate_freq(freq, alpha, freqCutoff, whiteList, blackList, verbose);
}

// ---- Exported: edgeOnLoop ----
// Called from the score_shd R function with 0-based node indices and a
// Rcpp::List parent set (each element is an integer vector of 0-based indices).
// [[Rcpp::export]]
bool edgeOnLoop(int fromNode, int toNode, const Rcpp::List& parSet) {
  const int p = parSet.size();
  AdjList parents(p);
  for (int i = 0; i < p; ++i)
    parents[i] = Rcpp::as<std::vector<int>>(parSet[i]);
  return edge_on_loop_cpp(fromNode, toNode, parents);
}

// ---- HC single run ----
// SC adaptation: continuous nodes regress only on rows where Y[r,i] != 0.
// indRow[i] is precomputed once before the loop (Y does not change during HC).
Rcpp::List hcSC1(const Rcpp::NumericMatrix& Y,
                 const Rcpp::CharacterVector& nodeType,
                 const Rcpp::LogicalMatrix& whiteList,
                 const Rcpp::LogicalMatrix& blackList,
                 double tol    = 1e-6,
                 int maxStep   = 500,
                 int seed      = 1,
                 bool flip     = true,
                 bool verbose  = false) {
  const int n = Y.nrow();
  const int p = Y.ncol();

  std::vector<std::string> types(p);
  for (int i = 0; i < p; ++i)
    types[i] = Rcpp::as<std::string>(nodeType[i]);

  // Initialize graph and adjacency lists from whitelist.
  Rcpp::LogicalMatrix curGraph = Rcpp::clone(whiteList);
  AdjList parents  = graph_to_parents(curGraph, p);
  AdjList children = graph_to_children(curGraph, p);

  // Precompute row-selection indicators.
  // indRow[i][r] = 1 iff Y(r,i) != 0; for binary nodes all rows are used.
  // n_sel[i] = number of selected rows (= effective sample size for node i).
  std::vector<std::vector<char>> indRow(p, std::vector<char>(n, 1));
  std::vector<int> n_sel(p, n);
  for (int i = 0; i < p; ++i) {
    if (types[i] == "c") {
      int cnt = 0;
      for (int r = 0; r < n; ++r) {
        indRow[i][r] = (Y(r, i) != 0.0) ? 1 : 0;
        cnt += indRow[i][r];
      }
      n_sel[i] = cnt;
      if (cnt < 2)
        Rcpp::stop("continuous node %d has fewer than 2 nonzero observations", i);
    }
  }

  // Pre-allocated workspace for design matrices (n rows, p+1 cols max).
  MatrixXd ws(n, p + 1);
  VectorXd y_buf(n);

  // Initial BIC scores.
  std::vector<double> curScore(p);
  for (int i = 0; i < p; ++i) {
    curScore[i] = bic_score_SC(ws, y_buf, Y, parents[i], types[i], i,
                               indRow[i], n_sel[i]);
    if (!std::isfinite(curScore[i]))
      Rcpp::stop("initial BIC score for node %d is non-finite", i);
  }

  // Version stamps: incremented whenever a node's parent set changes.
  // Cache entries are valid iff their stored version equals the current stamp.
  std::vector<int> version(p, 0);

  std::vector<OneCache> addCache(p * p);
  std::vector<OneCache> deleteCache(p * p);
  std::vector<OneCache> reverseCache(p * p);

  // Acyclicity cache. NA_LOGICAL = not yet checked.
  //   acyStatus(from, to) = acyclic status of adding from->to
  //   acyStatus(to, from) = acyclic status of reversing from->to
  Rcpp::LogicalMatrix acyStatus(p, p);
  std::fill(acyStatus.begin(), acyStatus.end(), NA_LOGICAL);

  LastOpState lastOp;
  std::mt19937 rng(static_cast<unsigned int>(seed));

  Rcpp::List stepOper;
  Rcpp::NumericVector stepDelta;

  int acceptedSteps = 0;
  while (acceptedSteps < maxStep) {
    double bestDelta  = 0.0;
    double bestScoreA = kInf, bestScoreB = kInf;
    int bestFrom = -1, bestTo = -1, bestType = -1;

    for (int to = 0; to < p; ++to) {
      for (int from = 0; from < p; ++from) {
        if (from == to) continue;
        const int idx = mat_index(from, to, p);

        if (curGraph(from, to)) {
          // ---- Existing edge: consider delete and (if eligible) reverse ----
          if (acyStatus(from, to) != true) acyStatus(from, to) = true;
          if (whiteList(from, to)) continue;

          // Delete from->to (only node `to` score changes).
          OneCache& del = deleteCache[idx];
          if (del.versionA != version[to]) {
            del.scoreA = bic_score_SC(ws, y_buf, Y, parents[to], types[to], to,
                                      indRow[to], n_sel[to], -1, from);
            del.value   = del.scoreA - curScore[to];
            del.versionA = version[to];
          }
          if (verbose)
            Rcpp::Rcout << "delete " << from << "->" << to
                        << ": delta=" << del.value << "\n";
          if (better_delta(del.value, bestDelta, tol, flip, rng)) {
            bestDelta = del.value; bestScoreA = del.scoreA;
            bestFrom = from; bestTo = to; bestType = 2;
          }

          // Reverse from->to into to->from.
          if (!blackList(to, from)) {
            if (acyStatus(to, from) == NA_LOGICAL) {
              acyStatus(to, from) =
                acyclic_reverse_inplace(from, to, parents) ? true : false;
            } else {
              acyclic_cache_update(lastOp, from, to, 3, parents, acyStatus);
              if (acyStatus(to, from) == NA_LOGICAL)
                acyStatus(to, from) =
                  acyclic_reverse_inplace(from, to, parents) ? true : false;
            }
            if (acyStatus(to, from)) {
              OneCache& rev = reverseCache[idx];
              if (rev.versionA != version[to] || rev.versionB != version[from]) {
                // scoreA (node `to` after deletion) reuses del.scoreA — same computation.
                rev.scoreA = del.scoreA;
                rev.scoreB = bic_score_SC(ws, y_buf, Y, parents[from], types[from], from,
                                          indRow[from], n_sel[from], to, -1);
                rev.value   = (rev.scoreA - curScore[to]) + (rev.scoreB - curScore[from]);
                rev.versionA = version[to];
                rev.versionB = version[from];
              }
              if (verbose)
                Rcpp::Rcout << "reverse " << from << "->" << to
                            << ": delta=" << rev.value << "\n";
              if (better_delta(rev.value, bestDelta, tol, flip, rng)) {
                bestDelta = rev.value; bestScoreA = rev.scoreA; bestScoreB = rev.scoreB;
                bestFrom = from; bestTo = to; bestType = 3;
              }
            }
          }

        } else if (!curGraph(to, from) && !blackList(from, to)) {
          // ---- No edge in either direction: consider add from->to ----
          if (acyStatus(from, to) == NA_LOGICAL) {
            acyStatus(from, to) = acyclic_add(from, to, parents) ? true : false;
          } else {
            acyclic_cache_update(lastOp, from, to, 1, parents, acyStatus);
            if (acyStatus(from, to) == NA_LOGICAL)
              acyStatus(from, to) = acyclic_add(from, to, parents) ? true : false;
          }
          if (!acyStatus(from, to)) continue;

          OneCache& add = addCache[idx];
          if (add.versionA != version[to]) {
            add.scoreA = bic_score_SC(ws, y_buf, Y, parents[to], types[to], to,
                                      indRow[to], n_sel[to], from, -1);
            add.value   = add.scoreA - curScore[to];
            add.versionA = version[to];
          }
          if (verbose)
            Rcpp::Rcout << "add " << from << "->" << to
                        << ": delta=" << add.value << "\n";
          if (better_delta(add.value, bestDelta, tol, flip, rng)) {
            bestDelta = add.value; bestScoreA = add.scoreA;
            bestFrom = from; bestTo = to; bestType = 1;
          }
        }
      }
    }

    if (bestDelta >= -tol || bestType < 0) break;

    // Apply accepted operation and bump version stamps to invalidate stale cache entries.
    if (bestType == 1) {            // Add bestFrom -> bestTo
      curGraph(bestFrom, bestTo) = true;
      parents[bestTo].push_back(bestFrom);
      children[bestFrom].push_back(bestTo);
      curScore[bestTo] = bestScoreA;
      ++version[bestTo];
    } else if (bestType == 2) {    // Delete bestFrom -> bestTo
      curGraph(bestFrom, bestTo) = false;
      erase_value(parents[bestTo], bestFrom);
      erase_value(children[bestFrom], bestTo);
      curScore[bestTo] = bestScoreA;
      ++version[bestTo];
    } else {                        // Reverse bestFrom -> bestTo
      curGraph(bestFrom, bestTo) = false;
      erase_value(parents[bestTo], bestFrom);
      erase_value(children[bestFrom], bestTo);
      curScore[bestTo] = bestScoreA;
      ++version[bestTo];

      curGraph(bestTo, bestFrom) = true;
      parents[bestFrom].push_back(bestTo);
      children[bestTo].push_back(bestFrom);
      curScore[bestFrom] = bestScoreB;
      ++version[bestFrom];
    }

    // Update LastOpState. BFS runs on the post-operation graph so that
    // acyclic_cache_update in the next step sees the updated topology.
    lastOp.from = bestFrom; lastOp.to = bestTo; lastOp.type = bestType;
    lastOp.fromAn = compute_ancestors(bestFrom, parents);
    lastOp.fromDe = compute_descendants(bestFrom, children);
    lastOp.toAn   = compute_ancestors(bestTo, parents);
    lastOp.toDe   = compute_descendants(bestTo, children);

    stepOper.push_back(Rcpp::IntegerVector::create(bestFrom, bestTo, bestType));
    stepDelta.push_back(bestDelta);
    ++acceptedSteps;
  }

  Rcpp::IntegerMatrix adjOut(p, p);
  for (int i = 0; i < p * p; ++i) adjOut[i] = curGraph[i] ? 1 : 0;
  return Rcpp::List::create(
    Rcpp::Named("adjacency")  = adjOut,
    Rcpp::Named("score")      = Rcpp::wrap(curScore),
    Rcpp::Named("operations") = stepOper,
    Rcpp::Named("deltaMin")   = stepDelta);
}

// ---- Exported: hcSC_ ----
// Runs hcSC1 `restart` times with distinct tie-breaking seeds (spaced by 101)
// and returns the result with the smallest total BIC score.
// [[Rcpp::export]]
Rcpp::List hcSC_(const Rcpp::NumericMatrix& Y,
                 const Rcpp::CharacterVector& nodeType,
                 const Rcpp::LogicalMatrix& whiteList,
                 const Rcpp::LogicalMatrix& blackList,
                 double tol      = 1e-6,
                 int    maxStep  = 500,
                 int    restart  = 10,
                 int    seed     = 1,
                 bool   verbose  = false) {
  const bool flip = restart > 1;
  Rcpp::List bestRes;
  double bestScore = kInf;

  for (int i = 0; i < restart; ++i) {
    const int iseed = static_cast<int>(
      static_cast<unsigned int>(seed) +
      static_cast<unsigned int>(i) * 101U);
    Rcpp::List curRes = hcSC1(Y, nodeType, whiteList, blackList,
                               tol, maxStep, iseed, flip, verbose);
    Rcpp::NumericVector sv = curRes["score"];
    const double s = Rcpp::sum(sv);
    if (std::isfinite(s) && s < bestScore) {
      bestRes   = curRes;
      bestScore = s;
    }
  }

  if (!std::isfinite(bestScore))
    Rcpp::stop("all HC restarts produced non-finite scores");
  return bestRes;
}

#pragma GCC diagnostic pop
