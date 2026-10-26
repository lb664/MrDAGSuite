# 135 characters #####################################################################################################################
#' @title invalidIV4MrDAG: Detection of invalid instrumental variables in the Bayesian causal graphical model for joint Mendelian 
#' randomization analysis of multiple exposures and outcomes
#' @description Detection of invalid instrumental variables (IVs) due to unmeasured pleiotropy, in which some genetic variants are 
#' directly associated with the outcomes and not via the exposures. The task is performed by using the conditional predictive 
#' ordinate (CPO). For details, see \insertCite{Zuber2025;textual}{MrDAGSuite}
#'
#' @param data Number of observations (IVs in a summary-level MR design) times the number of traits (both outcomes and exposures).
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param LD,LD_rsids,LD_pop,gencov Must match exactly what was passed to the \code{MrDAG} call that produced \code{output}: an
#' n-by-n LD matrix (or rsIDs to download one) and/or a per-IV sample-overlap \code{gencov}, so that the leave-one-out computation 
#' uses the same corrected sufficient statistic. See \code{\link{MrDAG}}.
#' When either is supplied, the leave-one-out CPO is computed by exact per-IV computation rather than the closed form used when 
#' neither is supplied, since that close form assumes IVs are uncorrelated and gencov-free.
#' @param displayIntro Print information regarding the \code{invIV4MrDAG} model on the screen (default \code{TRUE})

#' @param digits Number of digits (default \code{2}) used to round the frequency of invalid IVs
#'
#' @details
#' For details regarding the model and the algorithm, see \insertCite{Bottolo2026;textual}{MrDAGSuite} and
#' \insertCite{Zuber2025;textual}{MrDAGSuite}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(invalidIV, freqInvalidIV, logMargLik, logMargLik_i, CPO, scaleCPO,
#' timeinvIV4MrDAG)}:
#' \describe{
#'   \item{\code{invalidIV}}{ Detected invalid IVs with a scale CPO < 0.01, see \insertCite{Congdon2005;textual}{MrDAGSuite}. It is a 
#'         list object \code{list(IV, idxIV)} with the names of the invalid IVs and the row indices in \code{data} }
#'   \item{\code{freqInvalidIV}}{ Frequency of invalid IVs defined as outliers or high-leverage and influential IVs. It is defined as 
#'         the frequency of a scale CPO < 0.01 }
#'   \item{\code{logMargLik}}{ Log-marginal likelihood of explored DAGs which belong to the Markov Equivalence Classes whose unique 
#'         representative chain graphs are the EGs learned during MCMC iterations excluding the burn-in and with thinning }
#'   \item{\code{logMargLik_i}}{ Matrix (no. IVs times no. explored DAGs) log-marginal likelihood for the i-th IV of the explored 
#'         DAGs which belong to the Markov Equivalence Classes whose unique representative chain graphs are the EGs learned during 
#'         MCMC iterations excluding burn-in and with thinning }
#'   \item{\code{CPO}}{ Conditional predictive ordinate for each instrumental variable, see 
#'         \insertCite{Ntzoufras2009;textual}{MrDAGSuite} }
#'   \item{\code{scaleCPO}}{ Scaled conditional predictive ordinate for each instrumental variable, see 
#'         \insertCite{Congdon2005;textual}{MrDAGSuite} }
#'   \item{\code{timeinvIV4MrDAG}}{ Time in minutes employed by the \code{invIV4MrDAG} algorithm to analyze the data } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # Detection of invalid instrumental variables (IVs) due to unmeasured pleiotropy, in which some 
#' # genetic variants are directly associated with the outcomes and not via the exposures. To 
#' # estimate genetically predicted effects under intervention on mental health phenotypes that
#' # might impact
#' # lifestyle and behavioural traits, 470 independent IVs were selected to be associated at 
#' # genome-wide significance with seven mental health phenotypes that are considered exposures 
#' # for the risk of six lifestyle and behavioural traits
#'
#' # After loading the dataset:
#'
#' data(MD2LBT_data)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list 
#' # object:
#'
#' MRconstraint <- NULL
#' MRconstraint$Y_idx <- 1 : 6    # Lifestyle and behavioural traits
#' MRconstraint$X_idx <- 7 : 13   # Mental health phenotypes
#'
#' # The MrDAG algorithm is run to generate 2,000 posterior samples of all unknowns with memory 
#' # savings
#'
#' output <- MrDAG(data = MD2LBT_data$data, gencov = MD2LBT_data$gencov, 
#'                 niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                 MRconstraint = MRconstraint, fileName = NULL, seed = 28061971)
#'
#' # The invIV4MrDAG analysis is run to detect invalid IVs. For each saved iteration in MrDAG 
#' # and for each IV, invIV4MrDAG calculates the log-marginal likelihood and derives the IV 
#' # scaled CPO
#'
#' invIV_output <- invIV4MrDAG(data = MD2LBT_data$data, output, gencov = MD2LBT_data$gencov)
#' 
#' print(invIV_output$invalidIV)
#' }


########## invalidIV4MrDAG ##########

invIV4MrDAG <- function(data, output, LD = NULL, LD_rsids = NULL, LD_pop = "EUR", gencov = NULL,
                        displayIntro = TRUE, digits = 2)
{

  if (displayIntro == TRUE)
  {
    cat("\n##\n")
    cat("## Bayesian Mendelian randomization with Directed Acyclic Graphs exploration\n")
    cat("## Version: 0.1.0\n")
    cat("##\n")
    cat("## Copyright (C) 2026 L. Bottolo and V. Zuber\n")
    cat("##\n")
    cat("## Detection of invalid IVs\n")
    if (!is.null(output$hyperPar$MRconstraint))
    {
      cat("##\n")
      cat(paste0("## Outcomes  (traits ", paste(output$hyperPar$MRconstraint$Y_idx, collapse = ", "), "): "),
          paste(colnames(output$hyperPar$U)[output$hyperPar$MRconstraint$Y_idx], collapse = ", "), "\n")
      cat(paste0("## Exposures (traits ", paste(output$hyperPar$MRconstraint$X_idx, collapse = ", "), "): "),
          paste(colnames(output$hyperPar$U)[output$hyperPar$MRconstraint$X_idx], collapse = ", "), "\n")
    }
    cat("##\n")

  } else {
    cat("\n")
  }

  startTime <- Sys.time()

  X <- data
  n <- dim(X)[1]
  m <- dim(X)[2]

  ########## LD validation and optional download (must match the MrDAG() call that produced `output`) ##########

  if (!is.null(LD_rsids) && !is.null(LD))
    stop("Provide either 'LD' or 'LD_rsids', not both.", call. = FALSE)

  if (!is.null(LD_rsids)) {
    if (length(LD_rsids) != n)
      stop("'LD_rsids' must have length equal to nrow(data) (", n, ").", call. = FALSE)
    LD <- .get_LD_from_1000G(LD_rsids, pop = LD_pop)
    keep_idx <- match(rownames(LD), LD_rsids)
    keep_idx <- keep_idx[!is.na(keep_idx)]
    if (length(keep_idx) < n) {
      X <- X[keep_idx, , drop = FALSE]
      n <- nrow(X)
      cat("## Data subsetted to", n, "variants present in the 1000G panel.\n")

      if (!is.null(gencov) && (is.list(gencov) ||
                               (is.array(gencov) && length(dim(gencov)) == 3L))) {
        if (is.list(gencov)) {
          gencov <- gencov[keep_idx]
        } else {
          gencov <- gencov[, , keep_idx, drop = FALSE]
        }
        cat("## Per-IV 'gencov' subsetted to match the retained variants.\n")
      }
    }
  }

  if (!is.null(LD))
    LD <- .validate_LD(LD, n)

  if (displayIntro == TRUE) {
    if (!is.null(LD)) {
      ld_details <- c(
        if (!is.null(LD_rsids)) paste0("population: ", LD_pop),
        paste0("shrunk ", 100 * attr(LD, "ld_shrinkage"), "% towards the identity")
      )
      cat("## LD correction: enabled",
          if (length(ld_details) > 0) paste0(" (", paste(ld_details, collapse = "; "), ")"),
          "\n", sep = "")
    } else {
      cat("## LD correction: none (IVs assumed independent)\n")
    }
  }
  if (!is.null(LD)) attr(LD, "ld_shrinkage") <- NULL

  if (!is.null(gencov))
    gencov <- .validate_gencov(gencov, m, n)

  if (displayIntro == TRUE) {
    if (!is.null(gencov)) {
      cat("## Sample overlap correction (gencov): enabled\n")
    } else {
      cat("## Sample overlap correction (gencov): none (traits assumed estimated from non-overlapping samples)\n")
    }
    cat("##\n\n")
  }

  a   <- m
  U   <- diag(1, m) / n

  ########## Sufficient statistic (global, all n IVs) -- same correction as MrDAG() ##########

  if (!is.null(LD)) {
    R_chol  <- chol(LD)
    LD_inv  <- chol2inv(R_chol)
    weights <- diag(LD_inv)  # see .perIV_sufficientStat()'s comment in MrDAG.R
    L_ld    <- t(R_chol)
    X_w     <- forwardsolve(L_ld, X)
    tXX     <- crossprod(X_w)
  } else {
    weights <- rep(1, n)
    tXX     <- crossprod(X)
  }

  if (!is.null(gencov))
    tXX <- .correct_gencov(tXX, gencov, weights)
  attr(tXX, "alpha") <- NULL
  Upost_full <- U + tXX

  ########## Leave-one-out sufficient statistics ##########
  # Fast closed-form rank-1 shortcut (mlnodeLOO_vec) is valid only when IVs
  # are uncorrelated and gencov-free. Once LD and/or gencov are supplied,
  # each left-out IV's sufficient statistic is instead recomputed exactly:
  # whitened against the remaining (n-1)-IV LD submatrix and gencov-corrected
  # the same way as the full tXX above. This is done once upfront (not per
  # MCMC iteration) since it does not depend on the sampled DAG.

  exact_LOO <- !is.null(LD) || !is.null(gencov)

  if (exact_LOO) {
    tXX_loo   <- vector("list", n)
    Upost_loo <- vector("list", n)
    alphas_reduced <- numeric(0)

    for (i in seq_len(n)) {
      idx_i <- setdiff(seq_len(n), i)
      X_i   <- X[idx_i, , drop = FALSE]

      if (!is.null(LD)) {
        LD_i      <- LD[idx_i, idx_i, drop = FALSE]
        R_chol_i  <- chol(LD_i)
        weights_i <- diag(chol2inv(R_chol_i))
        X_w_i     <- forwardsolve(t(R_chol_i), X_i)
        tXX_i     <- crossprod(X_w_i)
      } else {
        weights_i <- rep(1, n - 1L)
        tXX_i     <- crossprod(X_i)
      }

      if (!is.null(gencov)) {
        # gencov was already standardised to an m x m x n array by .validate_gencov() above
        tXX_i <- .correct_gencov(tXX_i, gencov[, , idx_i, drop = FALSE], weights_i)
        alpha_i <- attr(tXX_i, "alpha"); attr(tXX_i, "alpha") <- NULL
        if (!is.null(alpha_i) && alpha_i < 1) alphas_reduced <- c(alphas_reduced, alpha_i)
      }

      tXX_loo[[i]]   <- tXX_i
      Upost_loo[[i]] <- U + tXX_i
    }

    if (length(alphas_reduced) > 0)
      message("WARNING: in ", length(alphas_reduced), " of ", n, " leave-one-out folds the full sample-overlap correction would make the ",
              "sufficient statistic non-positive-definite. Only a fraction of it was applied (smallest fraction ",
              .display_round(min(alphas_reduced), 3), "). Check the input data and gencov.")
  }

  niter     <- output$samplerPar$niter
  burnin    <- output$samplerPar$burnin
  thin      <- output$samplerPar$thin
  niter_TMP <- (niter - burnin) / thin
  q         <- length(output$hyperPar$MRconstraint$Y_idx)
  p         <- length(output$hyperPar$MRconstraint$X_idx)
  colnames_X <- colnames(output$hyperPar$U)
  
  if (is.null(dim(output$graphs)))
  {
    graphs <- array(dim = c(m, m, niter_TMP),
                    dimnames = list(colnames_X, colnames_X,
                                    paste0("iter", seq(burnin + thin, niter, thin))))
    for (iter in seq_len(niter_TMP))
      graphs[, , iter] <- bd_decode(output$graphs[iter])
  } else {
    graphs <- output$graphs
  }
  
  
  ###########
  # Storing #
  ###########
  
  logmarglik   <- numeric(niter_TMP)
  logmarglik_i <- matrix(0, n, niter_TMP)
  
  PBText <- utils::txtProgressBar(min = 2, max = niter_TMP, style = 3)
  
  for (iter in seq_len(niter_TMP))
  {
    G <- graphs[, , iter]

    for (node in seq_len(m))
      logmarglik[iter] <- logmarglik[iter] +
                          mlnodeDAGWishart(node, G, tXX, n, a, U, Upost_full)

    if (exact_LOO) {

      for (node in seq_len(m))
        for (i in seq_len(n))
          logmarglik_i[i, iter] <- logmarglik_i[i, iter] +
                                   mlnodeDAGWishart(node, G, tXX_loo[[i]], n - 1L, a, U, Upost_loo[[i]])

    } else {

      for (node in seq_len(m))
      {
        pa_j <- which(G[, node] != 0)   # parent set of node in iteration iter
        logmarglik_i[, iter] <- logmarglik_i[, iter] +
                                 mlnodeLOO_vec(node, pa_j, X, a, m)
      }
    }

    utils::setTxtProgressBar(PBText, iter)
  }
  
  close(PBText)
  cat("\n\n")
  
  f_i      <- exp(logmarglik_i)
  ICPO     <- rowMeans(1 / f_i)
  CPO      <- 1 / ICPO
  scaleCPO <- CPO / apply(f_i, 1, max)
  names(ICPO) <- names(CPO) <- names(scaleCPO) <- rownames(X)
  
  outlier      <- which(scaleCPO < 0.01)
  freq_outlier <- round(length(outlier) / length(CPO), digits)
  stability    <- round(sum(log(CPO) < -10, na.rm = TRUE) / length(CPO), digits)
  
  if (length(outlier) == 0)
  {
    outlier      <- "No invalid IVs have been detected"
    freq_outlier <- 0
  }
  
  endTime          <- Sys.time()
  timeinvIV4MrDAG  <- as.numeric(difftime(time1 = endTime, time2 = startTime, units = "min"))
  timeinvIV4MrDAG  <- paste0(round(timeinvIV4MrDAG, digits = 3), "m")
  
  output <- list(invalidIV       = outlier,
                 freqInvalidIV   = freq_outlier,
                 logMargLik      = logmarglik,
                 logMargLik_i    = logmarglik_i,
                 CPO             = CPO,
                 scaleCPO        = scaleCPO,
                 timeinvIV4MrDAG = timeinvIV4MrDAG)

  return(output)
}


########## mlnodeLOO_vec ##########

# Vectorised leave-one-out log-marginal likelihood for a single node j across ALL n IVs simultaneously, exploiting the closed-form 
# simplification that arises when U_i = I_m (identity) and n_i = 1 (one observation)
#
# KEY ANALYTICAL SIMPLIFICATION
# --------------------------------
# With U_i = I_m and tXX_i = x_i x_i^T (rank-1 outer product):
#
#   Upost_i = I_m + x_i x_i^T
#
# For parent set pa (size k) and target node j:
#
# PRIOR (does not depend on i at all, since U_i = I_m always):
#   U[pa, j]   = 0            (off-diagonal of identity)
#   U_jj.pa    = 1            (Schur complement of I_m is 1)
#   det(I_k)   = 1            (log = 0)
#   prior_nc   = -lgamma(a*/2) + (a*/2) * log(1/2)
#
# POSTERIOR (Sherman-Morrison + matrix determinant lemma):
#   Upost[pa, j]   = x_i[pa] * x_i[j]     (rank-1 structure)
#   Upost[pa, pa]  = I_k + x_i[pa] x_i[pa]^T
#
#   Let s_i = ||x_i[pa]||^2 = sum(X[i, pa]^2)
#
#   By Sherman-Morrison:
#     solve(I_k + x_i[pa] x_i[pa]^T) = I_k - x_i[pa] x_i[pa]^T / (1 + s_i)
#
#   Schur complement of Upost at (j, j):
#     Upost_jj.pa = 1 + x_i[j]^2 / (1 + s_i)
#
#   By matrix determinant lemma:
#     det(I_k + x_i[pa] x_i[pa]^T) = 1 + s_i
#
#   post_nc_i = -lgamma(a*/2 + 1/2) + (a*/2 + 1/2) * log(Upost_jj.pa / 2) + 1/2 * log(1 + s_i)
#
# The formula is identical for k=0 (empty parent set): s_i = 0, Upost_jj.pa = 1 + x_i[j]^2, det = 1 (log=0)
#
# ARGUMENTS
#   j          : node index (1-indexed)
#   pa_j       : integer vector of parent indices (possibly empty)
#   X          : n × m data matrix
#   a          : degrees-of-freedom hyperparameter (scalar)
#   m          : number of nodes (= ncol(X))
#
# RETURN
#   Numeric vector of length n: LOO log-ML contribution of node j for each IV, using the current DAG's parent set for j.

mlnodeLOO_vec <- function(j, pa_j, X, a, m)
{
  k      <- length(pa_j)
  a_star <- a + k - m + 1
  prior_nc <- -lgamma(a_star / 2) + (a_star / 2) * log(0.5)
  
  x_j <- X[, j]
  
  if (k == 0L) {
    s_vec <- rep(0, nrow(X))
  } else if (k == 1L) {
    s_vec <- X[, pa_j]^2
  } else {
    s_vec <- rowSums(X[, pa_j, drop = FALSE]^2)
  }
  
  Upost_jj_vec <- 1 + x_j^2 / (1 + s_vec)
  log_det_vec  <- log(1 + s_vec)
  post_nc_vec <- -lgamma(a_star / 2 + 0.5) +
                 (a_star / 2 + 0.5) * log(Upost_jj_vec / 2) +
                 0.5 * log_det_vec

  return(-0.5 * log(2 * pi) + prior_nc - post_nc_vec)
}

