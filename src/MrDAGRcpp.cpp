// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
bool is_DAG_cpp(const arma::mat& A_in) {
    int n = (int)A_in.n_rows;
    std::vector<int> color(n, 0); 

    for (int start = 0; start < n; ++start) {
        if (color[start] != 0) continue;
        std::vector<std::pair<int, int>> stk;
        stk.push_back({start, 0});
        color[start] = 1;

        while (!stk.empty()) {
            auto& [u, nxt] = stk.back();
            if (nxt == n) {
                color[u] = 2;
                stk.pop_back();
            } else {
                int v = nxt++;
                if (u == v || A_in(u, v) == 0.0) continue; 
                if (color[v] == 1) return false; 
                if (color[v] == 0) {
                    color[v] = 1;
                    stk.push_back({v, 0});
                }
            }
        }
    }
    return true;
}

static inline arma::uvec pa_internal(int j0, const arma::mat& DAG) {
    return arma::find(DAG.col(j0) != 0.0);
}

static std::pair<double, double> schur_and_logdet(int j0, const arma::uvec& pa, const arma::mat& M) {
    arma::mat M_papa = M.submat(pa, pa);
    arma::vec M_paj  = M.submat(pa, arma::uvec{(arma::uword)j0});
    arma::vec sol;
    bool success = arma::solve(sol, arma::symmatu(M_papa), M_paj, arma::solve_opts::likely_sympd);
    if(!success) return {M(j0,j0), 0.0};
    double Sjj = M(j0, j0) - arma::dot(M_paj, sol);
    double ld = arma::log_det_sympd(M_papa); 
    return {Sjj, ld};
}

// [[Rcpp::export]]
double mlnodeDAGWishart_cpp(int node, const arma::mat& DAG, const arma::mat& tXX, int n, double a, const arma::mat& U) {
    int j = node - 1;
    int m = (int)tXX.n_cols;
    arma::mat Upost = U + tXX;
    arma::uvec pa = pa_internal(j, DAG);
    int npa = (int)pa.n_elem;
    double a_star = a + npa - m + 1.0;
    double prior_nc, post_nc;

    if (npa == 0) {
        prior_nc = -std::lgamma(a_star * 0.5) + a_star * 0.5 * std::log(U(j, j) * 0.5);
        post_nc  = -std::lgamma(a_star * 0.5 + n * 0.5) + (a_star * 0.5 + n * 0.5) * std::log(Upost(j, j) * 0.5);
    } else {
        auto [U_jj, ld_U] = schur_and_logdet(j, pa, U);
        auto [Upost_jj, ld_Up] = schur_and_logdet(j, pa, Upost);
        prior_nc = -std::lgamma(a_star * 0.5) + a_star * 0.5 * std::log(U_jj * 0.5) + 0.5 * ld_U;
        post_nc  = -std::lgamma(a_star * 0.5 + n * 0.5) + (a_star * 0.5 + n * 0.5) * std::log(Upost_jj * 0.5) + 0.5 * ld_Up;
    }
    return -n * 0.5 * std::log(2.0 * M_PI) + prior_nc - post_nc;
}

// [[Rcpp::export]]
double logMargLik_cpp(const arma::mat& DAG, const arma::mat& tXX, int n, double a, const arma::mat& U) {
    double lml = 0.0;
    for (int j = 1; j <= (int)DAG.n_cols; ++j)
        lml += mlnodeDAGWishart_cpp(j, DAG, tXX, n, a, U);
    return lml;
}

// [[Rcpp::export]]
Rcpp::List rnodeDAGWishart_cpp(int node, const arma::mat& DAG, double aj, const arma::mat& U) {
    int j = node - 1;
    arma::uvec pa = pa_internal(j, DAG);
    int npa = (int)pa.n_elem;
    double sigmaj;
    arma::vec Lj;

    if (npa == 0) {
        sigmaj = 1.0 / R::rgamma(aj * 0.5, 2.0 / U(j, j));
        Lj = arma::vec(); 
    } else {
        arma::mat U_papa = U.submat(pa, pa);
        arma::vec U_paj  = U.submat(pa, arma::uvec{(arma::uword)j});
        arma::vec sol = arma::solve(arma::symmatu(U_papa), U_paj, arma::solve_opts::likely_sympd);
        double U_jj = U(j, j) - arma::dot(U_paj, sol);
        sigmaj = 1.0 / R::rgamma(aj * 0.5, 2.0 / U_jj);
        arma::mat invU_papa = arma::inv_sympd(U_papa);
        arma::vec mu = -invU_papa * U_paj;
        arma::mat R_mat = arma::chol(sigmaj * invU_papa);
        Lj = mu + R_mat.t() * arma::randn(npa);
    } // <--- This was the missing brace in your error log

    arma::ivec pa_r = arma::conv_to<arma::ivec>::from(pa + 1);
    return Rcpp::List::create(Named("sigmaj") = sigmaj, Named("Lj") = Lj, Named("pa_r") = pa_r);
}

// [[Rcpp::export]]
Rcpp::List rDAGWishart_cpp(const arma::mat& DAG, double a, const arma::mat& U) {
    int m = (int)DAG.n_cols;
    arma::mat L_mat = arma::eye(m, m);
    arma::mat D_mat = arma::zeros(m, m);

    for (int j = 0; j < m; ++j) {
        arma::uvec pa = pa_internal(j, DAG);
        double aj = a + (int)pa.n_elem - m + 1.0;
        Rcpp::List params = rnodeDAGWishart_cpp(j + 1, DAG, aj, U);
        D_mat(j, j) = params["sigmaj"];
        arma::vec Lj = params["Lj"];
        arma::ivec pa_r = params["pa_r"];
        for (int k = 0; k < (int)pa_r.n_elem; ++k) L_mat(pa_r(k) - 1, j) = Lj(k);
    }
    return Rcpp::List::create(Named("D") = D_mat, Named("L") = L_mat);
}

// [[Rcpp::export]]
Rcpp::List proposeDAG_cpp(const arma::mat& DAG, arma::ivec Y_idx, arma::ivec X_idx) {
    bool has_MRCheck = (Y_idx.n_elem > 0 && X_idx.n_elem > 0);
    int m = (int)DAG.n_rows;

    struct Op { int type, x, y; };
    std::vector<Op> O;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == j) continue;
            if (DAG(i, j) == 0.0) {
                O.push_back({1, i + 1, j + 1});
            } else {
                O.push_back({2, i + 1, j + 1});
                O.push_back({3, i + 1, j + 1});
            }
        }
    }
    int nOps = (int)O.size();

    arma::mat ANext;
    int opType_c, opX_c, opY_c;
    int validPropMrDAG = NA_INTEGER;

    while (true) {
        int i = (int)std::floor(R::runif(0.0, (double)nOps));
        i = std::min(i, nOps - 1);
        const Op& op = O[i];
        opType_c = op.type; opX_c = op.x; opY_c = op.y;

        ANext = DAG;
        if (op.type == 1)      ANext(op.x - 1, op.y - 1) = 1.0;
        else if (op.type == 2) ANext(op.x - 1, op.y - 1) = 0.0;
        else {
            ANext(op.x - 1, op.y - 1) = 0.0;
            ANext(op.y - 1, op.x - 1) = 1.0;
        }

        if (!is_DAG_cpp(ANext)) continue;

        if (has_MRCheck) {
            double s = 0.0;
            for (int yi : Y_idx)
                for (int xi : X_idx)
                    s += ANext(yi - 1, xi - 1);
            if (s == 0.0) {
                validPropMrDAG = 1;
            } else {
                ANext = DAG;
                validPropMrDAG = 0;
            }
        }
        break;
    }

    arma::ivec opNode = (opType_c == 3) ?
        arma::ivec{opX_c, opY_c} :
        arma::ivec{opY_c};

    return Rcpp::List::create(
        Named("proposedDAG")    = ANext,
        Named("opType")         = opType_c,
        Named("opNode")         = opNode,
        Named("currentOpCard")  = nOps,
        Named("proposedOpCard") = nOps,
        Named("validPropMrDAG") = validPropMrDAG
    );
}

// [[Rcpp::export]]
Rcpp::List acceptRejectDAG_cpp(const arma::mat& tXX, const arma::mat& currentDAG, const arma::mat& proposedDAG, 
                               arma::ivec node, int opType, int currentOpCard, int proposedOpCard, 
                               int n, double a, const arma::mat& U, double pp, double temp) {
    double diff = arma::accu(arma::abs(currentDAG - proposedDAG));
    if (diff == 0.0) return Rcpp::List::create(Named("isAccepted") = false, Named("logMargLik") = NA_REAL);
    double logpriorRatio = (opType == 1) ? std::log(pp / (1.0 - pp)) : (opType == 2 ? std::log((1.0 - pp) / pp) : 0.0);
    double cur_lml = 0, prop_lml = 0;
    for (int i = 0; i < (int)node.n_elem; ++i) {
        cur_lml  += mlnodeDAGWishart_cpp(node(i), currentDAG, tXX, n, a, U);
        prop_lml += mlnodeDAGWishart_cpp(node(i), proposedDAG, tXX, n, a, U);
    }
    double acceptRatio = std::min(0.0, (1.0 / temp) * (prop_lml - cur_lml + logpriorRatio));
    bool isAccepted = (std::log(R::runif(0.0, 1.0)) < acceptRatio);
    return Rcpp::List::create(
        Named("isAccepted") = isAccepted, 
        Named("logMargLik") = isAccepted ? logMargLik_cpp(proposedDAG, tXX, n, a, U) : NA_REAL
    );
}

// [[Rcpp::export]]
arma::vec genPredEffect_cpp(int response, arma::ivec target, const arma::mat& L, const arma::mat& D) {
    int y = response - 1;
    int nt = (int)target.n_elem;
    int m = (int)L.n_rows;
    arma::mat L_I = L;
    for (int k = 0; k < nt; ++k) {
        int tx = target(k) - 1;
        L_I.col(tx).zeros();
        L_I(tx, tx) = 1.0;
    }
    arma::mat LI_inv = arma::solve(arma::trimatl(L_I), arma::eye(m, m));
    arma::mat Sigma_I = LI_inv.t() * D * LI_inv;
    arma::vec effect(nt);
    for (int k = 0; k < nt; ++k) {
        int tx = target(k) - 1;
        effect(k) = Sigma_I(tx, y) / Sigma_I(tx, tx);
    }
    return effect;
}