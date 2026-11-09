# 135 characters #####################################################################################################################
#' @title Bootstrap MrDAG: Bootstrap evaluation of the Bayesian causal graphical model for joint Mendelian randomization (MR) analysis 
#' of multiple exposures and outcomes with Bayesian Directed Acyclic Graphs (DAGs) exploration
#' @description Bootstrap assessment of the robustness to weak instrumental variables (IVs) of \code{MrDAG} estimated Posterior 
#' Edge Inclusion Probabilities (PEIPs) and estimated genetically predicted effects by running the \code{MrDAG} algorithm on each bootstrap 
#' summary-level statistic, and the results aggregated by the average vote
#'
#' @param data Number of observations (IVs in a summary-level MR design) times the number of traits (both outcomes and exposures).
#' See also the \code{MRconstraint} argument
#' @param nboot Number of bootstrap samples (default \code{100}) of the summary-level dataset
#' @param CIboot Level (default \code{0.95}) of the confidence interval for the bootstrap frequency of edge inclusion and bootstrap 
#' genetically predicted effects. The frequency and the bootstrap genetically predicted effects are calculated as the average of the estimated PEIPs and the 
#' average of the estimated genetically predicted effects obtained by running the \code{MrDAG} algorithm on each bootstrap summary-level dataset
#' @param niter Number of Markov chain Monte Carlo (MCMC) iterations (including burn-in)
#' @param burnin Number of MCMC iterations to be discarded as burn-in
#' @param thin Parameter that defines how often the MCMC output should be stored, i.e., at every thin-th iteration
#' @param pp Prior edge inclusion probability (default \code{0.01}) between each pair of nodes (vertices) in the graph
#' @param a Degrees of freedom of the Wishart prior distribution on the precision matrix, i.e., the inverse of the covariance matrix 
#' between the traits (default \code{m}, where \code{m} is the total number of outcomes and exposures)
#' @param U (Proportional to the) expected value of the Wishart prior distribution on the precision matrix, i.e., the inverse of the 
#' covariance matrix between the traits. Specifically, it is proportional to an m-dimensional diagonal matrix by default, where 
#' \code{m} is the total number of outcomes and exposures
#' @param tempMax Annealing parameter \code{T} used to facilitate the convergence of the MCMC algorithm to the target distribution 
#' (default \code{10}). Temperature \code{1/T} increases linearly during the burn-in until \code{T=1} at the end of the burn-in 
#' @param MRconstraint List object containing the indices (columns) of the traits in \code{data} that are defined as outcomes and
#' exposures. Outcome indices must precede exposure indices, and this constraint is automatically checked at run time. If \code{NULL}
#' (default), no partition of the traits between outcomes and exposures is applied and the \code{MrDAG} algorithm performs
#' unconstrained structure learning across all traits. Note that this does not correspond to an MR analysis, since instruments valid
#' for exposures are invalid for outcomes and vice versa
#' @param LD,LD_rsids,LD_pop An optional n-by-n LD correlation matrix for the IVs (or rsIDs to download one). As in \code{\link{MrDAG}},
#' LD whitening is applied once, upfront, to the full (original) sample of IVs, and decomposed into one contribution per IV. Each
#' bootstrap replicate resamples IVs with replacement and re-sums the corresponding already-whitened per-IV contributions to build
#' its own sufficient statistic, so a bootstrap draw never needs to re-whiten a resampled (and potentially singular, since a
#' duplicated IV has correlation exactly 1 with itself) LD structure. See \code{\link{MrDAG}}
#' @param gencov An optional per-IV sample-overlap correction. As with \code{LD} above, it is applied once, upfront, to the full
#' sample of IVs and decomposed into one contribution per IV, which each bootstrap replicate resamples and re-sums together with
#' its \code{LD} contribution. See \code{\link{MrDAG}}
#' @param fileName Name of the file for the \code{boot4MrDAG} output object (default \code{"boot4MrDAG_object"}). If \code{NULL},
#' the output is not saved
#' @param filePath Directory where the \code{boot4MrDAG} object is saved. If the directory does not exist, \code{boot4MrDAG} will 
#' create it by using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), 
#' the \code{boot4MrDAG} object is saved in the current working directory
#' @param displayIntro Print information regarding the \code{boot4MrDAG} model on the screen (default \code{TRUE})

#' @param seed Seed to be used in the initialization of the MCMC algorithm (default \code{3112021}) and updates at each bootstrap 
#' (+1) on each summary-level dataset
#'
#' @details
#' For details regarding the model and the algorithm, see \insertCite{Bottolo2026;textual}{MrDAGSuite} and
#' \insertCite{Zuber2025;textual}{MrDAGSuite}. For bootstrap assessment criteria, see \insertCite{Glymour2019;textual}{MrDAGSuite}.
#' \code{LD} whitening and the \code{gencov} correction (if supplied) are applied once, upfront, to the full (original) sample of
#' IVs and decomposed into a per-IV sufficient statistic (see \code{\link{MrDAG}}). Each bootstrap replicate then resamples IVs
#' with replacement directly from this already-corrected per-IV decomposition and re-sums the resampled terms, rather than
#' resampling the raw IVs first and re-deriving the whitening/gencov correction on each (possibly duplicated) bootstrap sample
#'
#' @export
#'
#' @return The value returned is a list object \code{list(bootEdgeFreq, bootEdgeFreq_LL, bootEdgeFreq_UL, bootGenPredEffects,
#' bootGenPredEffects_LL, bootGenPredEffects_UL, bootFreq, bootOutput, timeBoot4MrDAG, hyperPar, samplerPar, correctionPar, opt)}:
#' \describe{
#'   \item{\code{bootEdgeFreq}}{ Frequency of edge inclusion calculated as the average of the estimated PEIPs obtained by running 
#'         the \code{MrDAG} algorithm on each bootstrap summary-level dataset and the results aggregated by the average vote }
#'   \item{\code{bootEdgeFreq_LL}}{ Lower limit of the \code{CIboot} for the bootstrap frequency of edge inclusion. It is calculated 
#'         as the (1-\code{CIboot})/2 % quantile of the estimated PEIPs obtained by running the \code{MrDAG} algorithm on each 
#'         bootstrap summary-level dataset }
#'   \item{\code{bootEdgeFreq_UL}}{ Upper limit of the \code{CIboot} for the bootstrap frequency of edge inclusion. It is calculated 
#'         as the 1-\[(1-\code{CIboot})/2\] % quantile of the estimated PEIPs obtained by running the \code{MrDAG} algorithm on each 
#'         bootstrap summary-level dataset }
#'   \item{\code{bootGenPredEffects}}{ Bootstrap genetically predicted effects calculated as the average of the estimated genetically predicted effects obtained by 
#'         running the \code{MrDAG} algorithm on each bootstrap summary-level dataset and the results aggregated by the average vote }
#'   \item{\code{bootGenPredEffects_LL}}{ Lower limit of the \code{CIboot} for the bootstrap genetically predicted effects. It is calculated as the 
#'         (1-\code{CIboot})/2 % quantile of the estimated genetically predicted effects obtained by running the \code{MrDAG} algorithm on each 
#'         bootstrap summary-level dataset }
#'   \item{\code{bootGenPredEffects_UL}}{ Upper limit of the \code{CIboot} for the bootstrap genetically predicted effects. It is calculated as the 
#'         1-\[(1-\code{CIboot})/2\] % quantile of the estimated genetically predicted effects obtained by running the \code{MrDAG} algorithm on 
#'         each bootstrap summary-level dataset }
#'   \item{\code{bootFreq}}{ The frequency of unique IVs sampled in the bootstrap analysis where IVs are randomly sampled with 
#'         replacement }
#'   \item{\code{bootOutput}}{ For each bootstrap sample, the outcome produced by the \code{MrDAG} algorithm. See \code{MrDAG} for 
#'         details }
#'   \item{\code{timeBoot4MrDAG}}{ Time in minutes employed by the \code{boot4MrDAG} algorithm to analyze the data }
#'   \item{\code{hyperPar}}{ List of the hyper-parameters \code{list(pp, a, U, tempMax)} and, if specified, the \code{MRconstraint} 
#'         list }
#'   \item{\code{samplerPar}}{ List of parameters used in the bootstrap analysis and MCMC algorithm \code{list(nboot, CIboot,
#'         niter, burnin, thin)} }
#'   \item{\code{correctionPar}}{ List recording whether LD and/or gencov corrections were applied, \code{list(LD_used, LD_pop,
#'         gencov_used)} }
#'   \item{\code{opt}}{ List of options used in the bootstrap analysis \code{list(seed)} } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # Assessment of the robustness to weak instrumental variables (IVs) of the estimated Posterior
#' # Edge Inclusion Probabilities (PEIPs) and the estimated genetically predicted effects under
#' # intervention
#' # on mental health phenotypes that might impact lifestyle and behavioural traits. 470
#' # independent IVs were selected to be associated at genome-wide significance with seven mental
#' # health phenotypes that are considered exposures for the risk of six lifestyle and behavioural
#' # traits
#'
#' # After loading the dataset
#'
#' data(MD2LBT_data)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list object
#'
#' MRconstraint <- NULL
#' MRconstraint$Y_idx <- 1 : 6    # Lifestyle and behavioural traits
#' MRconstraint$X_idx <- 7 : 13   # Mental health phenotypes
#'
#' # The boot4MrDAG algorithm is run to generate the bootstrap genetically predicted effects
#' # calculated as the
#' # average of the estimated genetically predicted effects obtained by running the
#' # \code{MrDAG} algorithm on 5
#' # bootstrap summary-level datasets and the results aggregated by the average vote. Each
#' # instance of the MrDAG algorithm is run to generate 2,000 posterior samples of all unknowns,
#' # with an overall 5,000 posterior samples
#'
#' output <- boot4MrDAG(data = MD2LBT_data$data, gencov = MD2LBT_data$gencov, nboot = 50, CIboot = 0.9,
#'                      niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                      MRconstraint = MRconstraint, fileName = NULL)
#' }


########## boot4MrDAG ##########

boot4MrDAG <- function(data, nboot = 100, CIboot = 0.95,
                       niter, burnin, thin,
                       pp = 0.05, a = NULL, U = NULL, tempMax = 10,
                       MRconstraint = NULL,
                       LD = NULL, LD_rsids = NULL, LD_pop = "EUR",
                       gencov = NULL,
                       fileName = "boot4MrDAG_object", filePath = NULL, displayIntro = TRUE,
                       seed = 31122021)
{
  if (!is.null(MRconstraint))
  {
    if (max(MRconstraint$Y_idx) >= min(MRconstraint$X_idx))
      stop("Outcome indices (Y_idx) must precede exposure indices (X_idx)")
  }

  if (displayIntro == TRUE)
  {
    cat("\n##\n")
    cat("## Bayesian Mendelian randomization with Directed Acyclic Graphs exploration\n")
    cat("## Version: 0.1.0\n")
    cat("##\n")
    cat("## Copyright (C) 2026 L. Bottolo and V. Zuber\n")
    cat("##\n")
    cat("## The model selected is:\n")
    cat("##\n")
    
    if (is.null(MRconstraint))
    {
      cat("## Bootstrap for Directed Acyclic Graph exploration\n")
    } else {
      cat("## Bootstrap for Mendelian randomization with Directed Acyclic Graphs exploration\n")
      cat("##\n")
      cat(paste0("## Outcomes  (traits ", paste(MRconstraint$Y_idx, collapse = ", "), "): "),
          paste(colnames(data)[MRconstraint$Y_idx], collapse = ", "), "\n")
      cat(paste0("## Exposures (traits ", paste(MRconstraint$X_idx, collapse = ", "), "): "),
          paste(colnames(data)[MRconstraint$X_idx], collapse = ", "), "\n")
    }
    
    cat("##\n")
  } else {
    cat("\n")
  }

  startTime <- Sys.time()

  ########## LD validation and optional download (once, before the bootstrap loop) ##########

  if (!is.null(LD_rsids) && !is.null(LD))
    stop("Provide either 'LD' or 'LD_rsids', not both.", call. = FALSE)

  if (!is.null(LD_rsids)) {
    if (length(LD_rsids) != nrow(data))
      stop("'LD_rsids' must have length equal to nrow(data) (", nrow(data), ").", call. = FALSE)
    LD <- .get_LD_from_1000G(LD_rsids, pop = LD_pop)
    keep_idx <- match(rownames(LD), LD_rsids)
    keep_idx <- keep_idx[!is.na(keep_idx)]
    if (length(keep_idx) < nrow(data)) {
      data <- data[keep_idx, , drop = FALSE]
      cat("## Data subsetted to", nrow(data), "variants present in the 1000G panel.\n")

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
    LD <- .validate_LD(LD, nrow(data))

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
    if (!is.null(gencov)) {
      cat("## Sample overlap correction (gencov): enabled\n")
    } else {
      cat("## Sample overlap correction (gencov): none (traits assumed estimated from non-overlapping samples)\n")
    }
    cat("##\n\n")
  }
  if (!is.null(LD)) attr(LD, "ld_shrinkage") <- NULL

  X <- data
  n <- dim(X)[1]
  m <- dim(X)[2]

  if (!is.null(gencov))
    gencov <- .validate_gencov(gencov, m, n)   # standardised once to an m x m x n array

  ########## Per-IV sufficient statistic (once, on the original data) ##########
  # LD whitening and the gencov correction (if supplied) are applied once, upfront, to the full sample of n IVs, and
  # decomposed into one contribution per IV -- see .perIV_sufficientStat() in MrDAG.R. Each bootstrap replicate then resamples
  # IVs with replacement and re-sums the corresponding already-corrected per-IV terms (.sumTXX()) to build its own sufficient
  # statistic, so IV_weights are unnecessary: a duplicated IV's contribution to tXX is simply added again

  tXXi <- .perIV_sufficientStat(X, LD, gencov)

  if (is.null(colnames(X)))
  {
    
    if (is.null(MRconstraint))
    {
      colnames(X) <- paste0("EO", seq(1, m))
    } else {
      colnames(X) <- c(paste0("O", seq(1, length(MRconstraint$Y_idx))), 
                       paste0("E", seq(1, length(MRconstraint$X_idx))))
    }
    
    colnames_YX <- toupper(substr(colnames(X), 1, 4))
  
  } else {
    colnames_YX <- toupper(substr(colnames(X), 1, 4))
    
    if (length(colnames_YX) != length(unique(colnames_YX)))
    {
      idx <- which(duplicated(colnames_YX))
      idx <- unique(c(idx, match(colnames_YX[idx], colnames_YX)))
      colnames_YX[idx] <- toupper(colnames(X)[idx])
    }
    
    colnames(X) <- colnames_YX
  }
  
  if (is.null(rownames(X)))
  {
    rownames(X) <- paste0("IV", seq(1, n))
    rowX_names <- rownames(X)
  }
  
  if (is.null(MRconstraint))
  {
    q <- m
    p <- m
  } else {
    q <- length(MRconstraint$Y_idx)
    p <- length(MRconstraint$X_idx)
  }
  
  Y_names <- colnames_YX[1 : q]
  X_names <- colnames_YX[q + (1 : p)]

  hyperPar <- NULL
  hyperPar$pp <- pp
  hyperPar$a <- a
  hyperPar$U <- U
  hyperPar$tempMax <- tempMax
  hyperPar$MRconstraint <- MRconstraint

  if (nboot <= 2)
  {
    warning("nboot must be at least 3; setting nboot = 3")
    nboot <- 3
  }

  samplerPar <- NULL
  samplerPar$nboot <- nboot
  samplerPar$CIboot <- CIboot
  samplerPar$niter <- niter
  samplerPar$burnin <- burnin
  samplerPar$thin <- thin

  ########## Per-replicate hyperPar/samplerPar/correctionPar ##########
  # Identical for every bootstrap replicate (standard bootstrap always resamples n IVs, so the effective sample size used in
  # the Wishart update is n regardless of duplicates), so they are resolved and built once here rather than inside the loop.
  # Mirrors what a bare MrDAG() call would have stored in its own output for the same (a, U) defaults

  a_boot <- if (is.null(a)) m else a
  U_boot <- if (is.null(U)) diag(1, m) / n else U

  hyperPar_boot             <- list(pp = pp, a = a_boot, U = U_boot, tempMax = tempMax, MRconstraint = MRconstraint)
  rownames(hyperPar_boot$U) <- colnames_YX
  colnames(hyperPar_boot$U) <- colnames_YX

  samplerPar_boot <- list(niter = niter, burnin = burnin, thin = thin)

  correctionPar_boot             <- list()
  correctionPar_boot$LD_used     <- !is.null(LD)
  correctionPar_boot$LD_pop      <- if (!is.null(LD)) LD_pop else NULL
  correctionPar_boot$gencov_used <- !is.null(gencov)

  ########## Storing ##########
  
  bootOutput <- NULL
  bootFreq <- rep(NA, nboot)
  names(bootFreq) <- paste0("boot", seq(1, nboot, 1))
  bootEdgeFreq <- array(0, dim = c(m, m, nboot), dimnames = list(colnames_YX, colnames_YX, 
                                                                 paste0("boot", seq(1, nboot, 1))))
  edgeFreq <- matrix(0, nrow = m, ncol = m, dimnames = list(colnames_YX, colnames_YX))
  edgeFreq_LL <- matrix(0, nrow = m, ncol = m, dimnames = list(colnames_YX, colnames_YX))
  edgeFreq_UL <- matrix(0, nrow = m, ncol = m, dimnames = list(colnames_YX, colnames_YX))
  
  if (!is.null(MRconstraint))
  {
    bootGenPredEffects <- array(0, dim = c(p, q, nboot), dimnames = list(X_names, Y_names, 
                                                                        paste0("boot", seq(1, nboot, 1))))
    genPredEffects <- matrix(0, nrow = p, ncol = q, dimnames = list(X_names, Y_names))
    genPredEffects_LL <- matrix(0, nrow = p, ncol = q, dimnames = list(X_names, Y_names))
    genPredEffects_UL <- matrix(0, nrow = p, ncol = q, dimnames = list(X_names, Y_names))
  } else {
    bootGenPredEffects <- array(0, dim = c(p, q, nboot), dimnames = list(colnames_YX, colnames_YX, 
                                                                        paste0("boot", seq(1, nboot, 1))))
    genPredEffects <- matrix(0, nrow = p, ncol = q, dimnames = list(colnames_YX, colnames_YX))
    genPredEffects_LL <- matrix(0, nrow = p, ncol = q, dimnames = list(colnames_YX, colnames_YX))
    genPredEffects_UL <- matrix(0, nrow = p, ncol = q, dimnames = list(colnames_YX, colnames_YX))
  }
  
  PBText <- utils::txtProgressBar(min = 2, max = nboot, style = 3)
  alphas_reduced <- numeric(0)

  for (b in 1 : nboot)
  {
    set.seed(seed + b)

    # Resample and re-sum the already whitened/gencov-corrected per-IV terms directly -- a duplicated IV's contribution is
    # simply added as many times as it was drawn, so a repeated IV never needs a literal (and potentially singular) LD.
    # The corrected sufficient statistic of the resample is safeguarded: if the full correction would make it non-positive-definite,
    # only the fraction that keeps it positive-definite is subtracted (see .safeguardPD() in MrDAG.R).
    idx <- sample(seq_len(n), replace = TRUE)
    tXX_b <- .sumTXX(tXXi, idx)
    alpha_b <- attr(tXX_b, "alpha"); attr(tXX_b, "alpha") <- NULL
    if (!is.null(alpha_b) && alpha_b < 1) alphas_reduced <- c(alphas_reduced, alpha_b)
    bootFreq[b] <- length(unique(idx)) / n

    utils::capture.output(
    mcmcOut_b <- .runMrDAG_MCMC(tXX = tXX_b, n = n, a = a_boot, U = U_boot, tempMax = tempMax, pp = pp,
                                MRconstraint = MRconstraint, niter = niter, burnin = burnin, thin = thin,
                                colnames_YX = colnames_YX, seed = seed + b))

    bootOutput[[b]] <- list(graphs         = mcmcOut_b$graphs,
                            L              = mcmcOut_b$L,
                            D              = mcmcOut_b$D,
                            logMargLik     = mcmcOut_b$logMargLik,
                            validPropMrDAG = mcmcOut_b$validPropMrDAG,
                            acceptPropDAG  = mcmcOut_b$acceptPropDAG,
                            hyperPar       = hyperPar_boot,
                            samplerPar     = samplerPar_boot,
                            correctionPar  = correctionPar_boot,
                            opt            = list(seed = seed + b))

    bootEdgeFreq[, , b] <- get_postEdgeProb(bootOutput[[b]])$postEdgeProb
    bootGenPredEffects[, , b] <- get_genPredEffects(bootOutput[[b]], progress = NULL)$genPredEffects
    utils::setTxtProgressBar(PBText, b)
  }
  close(PBText)
  cat("\n\n")

  if (length(alphas_reduced) > 0)
    message("WARNING: in ", length(alphas_reduced), " of ", nboot, " bootstrap resamples the full sample-overlap correction would make the ",
            "sufficient statistic non-positive-definite. Only a fraction of it was applied (smallest fraction ",
            .display_round(min(alphas_reduced), 3), "). Check the input data and gencov.")

  edgeFreq <- apply(bootEdgeFreq, c(1, 2), mean, na.rm = TRUE)
  edgeFreq_LL <- apply(bootEdgeFreq, c(1, 2), stats::quantile, (1 - CIboot) / 2, na.rm = TRUE)
  edgeFreq_UL <- apply(bootEdgeFreq, c(1, 2), stats::quantile, 1 - ((1 - CIboot) / 2), na.rm = TRUE)
  genPredEffects <- apply(bootGenPredEffects, c(1, 2), mean, na.rm = TRUE)
  genPredEffects_LL <- apply(bootGenPredEffects, c(1, 2), stats::quantile, (1 - CIboot) / 2, na.rm = TRUE)
  genPredEffects_UL <- apply(bootGenPredEffects, c(1, 2), stats::quantile, 1 - ((1 - CIboot) / 2), na.rm = TRUE)
  
  opt <- list(seed = seed)
  endTime <- Sys.time()
  timeBoot4MrDAG <- as.numeric(difftime(time1 = endTime, time2 = startTime, units = "min"))
  timeBoot4MrDAG <- paste0(round(timeBoot4MrDAG, digits = 3), "m")

  output <- list(bootEdgeFreq = edgeFreq, bootEdgeFreq_LL = edgeFreq_LL, bootEdgeFreq_UL = edgeFreq_UL,
                 bootGenPredEffects = genPredEffects, bootGenPredEffects_LL = genPredEffects_LL, bootGenPredEffects_UL = genPredEffects_UL,
                 bootFreq = bootFreq,
                 bootOutput = bootOutput,
                 timeBoot4MrDAG = timeBoot4MrDAG,
                 hyperPar = hyperPar,
                 samplerPar = samplerPar,
                 correctionPar = correctionPar_boot,
                 opt = opt)
  
  if (!is.null(fileName))
  {
    saveas <- paste(fileName, "RData", sep = ".")
    
    if (is.null(filePath))
    {
      filePath <- getwd()
    } else {
      
      if (!dir.exists(filePath))
      {
        dir.create(filePath, recursive = TRUE)
        cat("\nThe directory", paste0("\"", filePath, "\""), "does not exist. It is now created.\n")
      }
    
    }
    
    if (.Platform$OS.type == "windows")
    {
      filePath <- paste0(normalizePath(filePath), "\\")
    } else {
      filePath <- paste0(normalizePath(filePath), "/")
    }
    
    save(output, file = paste0(filePath, saveas))
  }
  
  return(output)
}
