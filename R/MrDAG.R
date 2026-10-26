utils::globalVariables(c("prior_nc", "value_TMP", "value_LL_UL", "Color"))
# 135 characters #####################################################################################################################
#' @title MrDAG: Bayesian causal graphical model for joint Mendelian randomization analysis of multiple exposures and outcomes
#' @description Markov chain Monte Carlo (MCMC) implementation of a Bayesian multivariable, multi-response Mendelian randomization 
#' (MR) model for summary-level data with structure learning and genetically predicted effect estimation. The directionality of the genetically predicted effects 
#' between the exposures and the outcomes is assumed known, i.e., the exposures can only be potential causes of the outcomes, and no 
#' reverse causation is allowed
#'
#' @param data Number of observations (IVs in a summary-level MR design) times the number of traits (both outcomes and exposures).
#' See also the \code{MRconstraint} argument
#' @param niter Number of MCMC iterations (including burn-in)
#' @param burnin Number of MCMC iterations to be discarded as burn-in
#' @param thin Parameter that defines how often the MCMC output should be stored, i.e., at every \code{thin}-th iteration
#' @param pp Prior edge inclusion probability (default \code{0.01}) between each pair of nodes (vertices) in the graph
#' @param a Degrees of freedom of the Wishart prior distribution on the precision matrix, i.e., the inverse of the covariance matrix 
#' between the traits (default \code{m}, where \code{m} is the total number of outcomes and exposures)
#' @param U (Proportional to the) expected value of the Wishart prior distribution on the precision matrix, i.e., the inverse of the 
#' covariance matrix between the traits. Specifically, it is proportional to an \code{m}-dimensional diagonal matrix by default, where 
#' \code{m} is the total number of outcomes and exposures
#' @param tempMax Annealing parameter \code{T} used to facilitate the convergence of the MCMC algorithm to the target distribution 
#' (default \code{10}). Temperature \code{1/T} increases linearly during the burn-in until \code{T=1} at the end of the burn-in
#' @param MRconstraint List object containing the indices (columns) of the traits in \code{data} that are defined as outcomes and
#' exposures. Outcome indices must precede exposure indices, and this constraint is automatically checked at run time. If \code{NULL}
#' (default), no partition of the traits between outcomes and exposures is applied and the \code{MrDAG} algorithm performs
#' unconstrained structure learning across all traits. Note that this does not correspond to an MR analysis, since instruments valid
#' for exposures are invalid for outcomes and vice versa
#' @param LD An optional n-by-n LD correlation matrix for the IVs, where n is the number of IVs. When
#' supplied, the sufficient statistic is corrected for linkage disequilibrium via GLS whitening, replacing \code{t(X) \%*\% X} with
#' \code{t(X) \%*\% solve(LD) \%*\% X}. Must be symmetric; it is shrunk 1\% toward the identity automatically (this makes it
#' positive-definite and bounds the per-IV GLS-whitening weights; it is needed when the reference panel has fewer individuals than
#' there are IVs). If the shrunk matrix is still not positive-definite the algorithm stops.
#' @param LD_rsids Optional character vector of rsIDs used to download the LD matrix automatically
#' from the 1000 Genomes Project via \code{ieugwasr::ld_matrix}. Requires the \code{ieugwasr} package and an OpenGWAS API token
#' stored as \code{OPENGWAS_JWT} in \code{.Renviron} (see \url{https://api.opengwas.io/}). Rows corresponding to variants absent
#' from the reference panel are silently dropped
#' @param LD_pop 1000 Genomes super-population code used when \code{LD_rsids} is supplied (\code{"EUR"} default). Common options
#' are \code{"EUR"}, \code{"EAS"}, \code{"AFR"}, \code{"AMR"}, \code{"SAS"}
#' @param gencov An optional error covariance term capturing correlation between summary statistic estimation errors induced by
#' overlapping participants across the GWAS/regressions used to estimate the traits in \code{data} (outcomes AND exposures).
#' Accepts two equivalent forms: (i) a list of n m-by-m matrices, one per IV (e.g. the direct output of
#' \code{MVMR::phenocov_mvmr()}); or (ii) an m-by-m-by-n array equivalent to (i). A single matrix broadcast to every IV
#' (homoscedastic assumption) is deliberately not accepted — see \code{.validate_gencov()} for why. Only the off-diagonal
#' entries (sample-overlap covariance between traits) are used: the main diagonal (each trait's squared standard error) is not
#' corrected and is ignored. If the full correction would make the corrected sufficient statistic non-positive-definite, only the fraction of it that keeps the statistic positive-definite is subtracted and a warning reports the fraction
#' @param fileName Name of the file for the \code{MrDAG} output object (default \code{"MrDAG_object"}). If \code{NULL}, the output is
#' not saved
#' @param filePath Directory where the \code{MrDAG} object is saved. If the directory does not exist, the \code{MrDAG} algorithm will 
#' create it by using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), 
#' the \code{MrDAG} object is saved in the current working directory
#' @param displayIntro Print information regarding the \code{MrDAG} model on the screen (default \code{TRUE})

#' @param seed Seed to be used in the initialization of the MCMC algorithm (default \code{3112021})
#'
#' @details
#' For details regarding the model and the algorithm, see \insertCite{Bottolo2026;textual}{MrDAGSuite} and
#' \insertCite{Zuber2025;textual}{MrDAGSuite}. When \code{LD} and/or \code{gencov} are supplied, they are applied once, upfront, to
#' the full sample of IVs to build the (whitened, gencov-corrected) sufficient statistic, which is then decomposed into one
#' contribution per IV. \code{\link{boot4MrDAG}} reuses this same per-IV decomposition to build each bootstrap replicate's
#' sufficient statistic directly, by resampling and re-summing the already-corrected per-IV terms, rather than re-deriving the
#' correction on a resampled sample of raw IVs
#'
#' @export
#'
#' @return The value returned is a list object \code{list(graphs, L, D, logMargLik, validPropMrDAG, acceptPropDAG, timeMrDAG,
#' hyperPar, samplerPar, correctionPar, opt)}:
#' \describe{
#'   \item{\code{graphs}}{ Explored DAGs that belong to the learned Completed Partially Directed Acyclic Graphs (CPDAGs) or Essential 
#'         Graphs (EGs). The number of explored DAGs corresponds 
#'         to the number of thinned (\code{thin}) MCMC iterations (excluding burn-in) }
#'   \item{\code{L}}{ Posterior samples of the lower triangular matrix of the modified Cholesky decomposition (\code{L,D}) }
#'   \item{\code{D}}{ Posterior samples of the diagonal matrix of the modified Cholesky decomposition (\code{L,D}) }
#'   \item{\code{logMargLik}}{ Log-marginal likelihood of explored DAGs which belong to the Markov Equivalence Classes whose unique 
#'         representative chain graphs are the EGs learned during MCMC iterations (including burn-in) without thinning }
#'   \item{\code{validPropMrDAG}}{ If \code{MRconstraint} is different from \code{NULL}, the proportion of proposed DAGs that comply 
#'         with the partial ordering implied by the partition of the traits between exposures and outcomes, with genetically predicted effects only 
#'         from the former to the latter }
#'   \item{\code{acceptPropDAG}}{ Proportion of proposed DAGs that are accepted by the Metropolis-Hastings ratio after 
#'         burn-in. If \code{MRconstraint} is different from \code{NULL}, \code{acceptPropDAG} corresponds to the proposed DAGs that 
#'         comply with the partial ordering }
#'   \item{\code{timeMrDAG}}{ Time in minutes employed by the \code{MrDAG} algorithm to analyze the data }
#'   \item{\code{hyperPar}}{ List of the hyper-parameters \code{list(pp, a, U, tempMax)} and, if specified, the 
#'         \code{MRconstraint} list }
#'   \item{\code{samplerPar}}{ List of parameters used in the MCMC algorithm \code{list(niter, burnin, thin)} }
#'   \item{\code{correctionPar}}{ List recording whether LD and/or gencov corrections were applied,
#'         \code{list(LD_used, LD_pop, gencov_used)} }
#'   \item{\code{opt}}{ List of options used in the \code{MrDAG} algorithm \code{list(seed)} } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # Analysis of lifestyle and behavioural exposures that might impact mental health phenotypes.
#' # 708 independent Instrumental Variables (IVs) were selected to be associated at genome-wide
#' # significance with six lifestyle and behavioural traits after pruning or clumping, which are
#' # considered exposures for the risk of seven mental health outcome phenotypes
#'
#' # After loading the dataset:
#'
#' data(LBT2MD_data)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list object:
#'
#' MRconstraint <- NULL
#' MRconstraint$Y_idx <- 1 : 7    # Mental health phenotypes
#' MRconstraint$X_idx <- 8 : 13   # Lifestyle and behavioural traits
#'
#' # The MrDAG algorithm is run to generate 2,000 posterior samples of all unknowns
#'
#' output <- MrDAG(data = LBT2MD_data$data, gencov = LBT2MD_data$gencov, 
#'                 niter = 150000, burnin = 50000, thin = 50, pp = 0.10, tempMax = 20, 
#'                 MRconstraint = MRconstraint, fileName = NULL, 
#'                 seed = 28061971)
#' }
#'
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @useDynLib MrDAGSuite, .registration = TRUE

MrDAG <- function(data,
                  niter, burnin, thin,
                  pp = 0.05, a = NULL, U = NULL, tempMax = 10,
                  MRconstraint = NULL,
                  LD = NULL, LD_rsids = NULL, LD_pop = "EUR",
                  gencov = NULL,
                  fileName = "MrDAG_object", filePath = NULL, displayIntro = TRUE,
                  seed = 31122021)
{

  if (!is.null(MRconstraint))
  {
    if (max(MRconstraint$Y_idx) >= min(MRconstraint$X_idx))
      stop("Outcome indices (Y_idx) must precede exposure indices (X_idx)")
  }

if (displayIntro == TRUE)
  {
    cat("\n")
    cat("##\n")
    cat("## Bayesian Mendelian randomization with Directed Acyclic Graphs exploration\n")
    cat("## Version: 0.1.0\n")
    cat("##\n")
    cat("## Copyright (C) 2026 L. Bottolo and V. Zuber\n")
    cat("##\n")
    cat("## The model selected is:\n")
    cat("##\n")
    
    if (is.null(MRconstraint))
    {
      cat("## Directed Acyclic Graphs exploration\n")
      cat("##\n")
      cat("## Dataset analysed:\n")
      cat("## All traits\n")
    } else {
      cat("## Mendelian randomization with Directed Acyclic Graphs exploration\n")
      
      cat("##\n")
      cat("## Dataset analysed:\n")
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
  set.seed(seed)
  
  X <- data
  n <- dim(X)[1]
  m <- dim(X)[2]

  ########## LD validation and optional download ##########

  if (!is.null(LD_rsids) && !is.null(LD))
    stop("Provide either 'LD' or 'LD_rsids', not both.", call. = FALSE)

  if (!is.null(LD_rsids)) {
    if (length(LD_rsids) != n)
      stop("'LD_rsids' must have length equal to nrow(data) (", n, ").", call. = FALSE)
    LD <- .get_LD_from_1000G(LD_rsids, pop = LD_pop)
    # Subset rows of X (and, if per-IV, gencov) to variants present in the reference panel
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

  ld_shr <- if (!is.null(LD)) attr(LD, "ld_shrinkage") else NULL

  if (displayIntro == TRUE) {
    if (!is.null(LD)) {
      ld_details <- c(
        if (!is.null(LD_rsids)) paste0("population: ", LD_pop),
        paste0("shrunk ", 100 * ld_shr, "% towards the identity")
      )
      cat("## LD correction: enabled (", paste(ld_details, collapse = "; "), ")\n", sep = "")
    } else {
      cat("## LD correction: none (IVs assumed independent)\n")
    }
  }
  if (!is.null(LD)) attr(LD, "ld_shrinkage") <- NULL

  ########## gencov validation ##########

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

  if (is.null(colnames(X)))
  {
    if (is.null(MRconstraint))
      colnames(X) <- paste0("EO", seq(1, m))
    else
      colnames(X) <- c(paste0("O", seq(1, length(MRconstraint$Y_idx))), 
                       paste0("E", seq(1, length(MRconstraint$X_idx))))
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
  
  colnames_YX <- colnames(X)
  
  if (is.null(rownames(X)))
  {
    rownames(X) <- paste0("IV", seq(1, n))
    rownames_X  <- rownames(X)
  }
  
  if (is.null(a)) a <- m
  if (is.null(U)) U <- diag(1, m) / n

  hyperPar             <- list()
  hyperPar$pp          <- pp
  hyperPar$a           <- a
  hyperPar$U           <- U
  rownames(hyperPar$U) <- colnames_YX
  colnames(hyperPar$U) <- colnames_YX
  hyperPar$tempMax     <- tempMax
  hyperPar$MRconstraint  <- MRconstraint

  samplerPar        <- list()
  samplerPar$niter  <- niter
  samplerPar$burnin <- burnin
  samplerPar$thin   <- thin

  correctionPar               <- list()
  correctionPar$LD_used       <- !is.null(LD)
  correctionPar$LD_pop        <- if (!is.null(LD)) LD_pop else NULL
  correctionPar$gencov_used   <- !is.null(gencov)

  ########## Sufficient statistic ##########
  # LD whitening and the gencov correction are both applied once, upfront, to
  # the full set of n IVs, and decomposed into one contribution per IV:
  # tXXi[, , i] is IV i's whitened, gencov-corrected contribution, and
  # tXX = sum_i tXXi[, , i] recovers the usual sufficient statistic exactly
  # (see .perIV_sufficientStat() and .sumTXX() below). Decomposing tXX this
  # way is what lets boot4MrDAG() build each bootstrap replicate's sufficient
  # statistic directly, by resampling and re-summing these already-corrected
  # per-IV terms, instead of re-deriving the whitening/gencov correction from
  # a resampled (and potentially singular) LD structure on every draw

  tXXi <- .perIV_sufficientStat(X, LD, gencov)
  tXX  <- .sumTXX(tXXi)
  gencov_alpha <- attr(tXX, "alpha"); attr(tXX, "alpha") <- NULL
  correctionPar$gencov_alpha <- gencov_alpha
  if (!is.null(gencov_alpha) && gencov_alpha < 1)
    message(strrep("!", 70), "\nWARNING: the full sample-overlap correction would make the sufficient statistic non-positive-definite. ",
            "Only ", .display_round(100 * gencov_alpha, 1), "% of it was applied. Check the input data and gencov.\n", strrep("!", 70))

  mcmcOut <- .runMrDAG_MCMC(tXX = tXX, n = n, a = a, U = U, tempMax = tempMax, pp = pp,
                            MRconstraint = MRconstraint, niter = niter, burnin = burnin,
                            thin = thin, colnames_YX = colnames_YX, seed = seed)

  opt       <- list(seed = seed)
  endTime   <- Sys.time()
  timeMrDAG <- as.numeric(difftime(time1 = endTime, time2 = startTime, units = "min"))
  timeMrDAG <- paste0(round(timeMrDAG, digits = 3), "m")

  output <- list(graphs          = mcmcOut$graphs,
                 L               = mcmcOut$L,
                 D               = mcmcOut$D,
                 logMargLik      = mcmcOut$logMargLik,
                 validPropMrDAG  = mcmcOut$validPropMrDAG,
                 acceptPropDAG   = mcmcOut$acceptPropDAG,
                 timeMrDAG       = timeMrDAG,
                 hyperPar        = hyperPar,
                 samplerPar      = samplerPar,
                 correctionPar   = correctionPar,
                 opt             = opt)
  
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
        cat("\nThe directory", paste0("\"", filePath, "\""),
            "does not exist. It is now created.\n")
      }
    }
    
    if (.Platform$OS.type == "windows")
      filePath <- paste0(normalizePath(filePath), "\\")
    else
      filePath <- paste0(normalizePath(filePath), "/")
    
    save(output, file = paste0(filePath, saveas))
  }
  
  return(output)
}


########## .runMrDAG_MCMC ##########
# The MCMC engine shared by MrDAG() and, per bootstrap replicate, boot4MrDAG():
# given an already-built sufficient statistic tXX (and its effective number of
# IVs, n), runs the DAG structure-learning chain and returns the raw sampler
# output. Called by MrDAG() on the full-sample tXX, and by boot4MrDAG() on
# each bootstrap replicate's resampled tXX -- so the two are guaranteed to
# analyse a sufficient statistic in exactly the same way

.runMrDAG_MCMC <- function(tXX, n, a, U, tempMax, pp, MRconstraint,
                           niter, burnin, thin, colnames_YX, seed)
{
  set.seed(seed)
  m     <- ncol(tXX)
  Upost <- U + tXX

  if (!is.null(MRconstraint))
  {
    Y_idx_cpp <- as.integer(MRconstraint$Y_idx)
    X_idx_cpp <- as.integer(MRconstraint$X_idx)
  } else {
    Y_idx_cpp <- integer(0)
    X_idx_cpp <- integer(0)
  }

  iterCounter    <- 1L
  validPropMrDAG <- 0
  acceptPropDAG  <- 0

  ########## Storing ##########

  graphs <- vector("double", (niter - burnin) / thin)
  L      <- vector("double", (niter - burnin) / thin)
  D      <- vector("double", (niter - burnin) / thin)
  names(graphs) <- paste0("iter", seq(burnin + thin, niter, thin))
  names(L)      <- paste0("iter", seq(burnin + thin, niter, thin))
  names(D)      <- paste0("iter", seq(burnin + thin, niter, thin))

  temp       <- c(seq(tempMax, 1, length.out = burnin), rep(1, niter))
  logMargLik <- vector("double", niter)
  names(logMargLik) <- paste0("iter", seq(1, niter, 1))

  currentDAG <- matrix(0, ncol = m, nrow = m, dimnames = list(colnames_YX, colnames_YX))
  PBText     <- utils::txtProgressBar(min = 2, max = niter, style = 3)

  ########## MCMC iterations ##########

  for (iter in 1L : niter)
  {
    proposed <- proposeDAG_cpp(currentDAG, Y_idx_cpp, X_idx_cpp)

    if (iter > burnin & !is.null(MRconstraint))
      validPropMrDAG <- validPropMrDAG + proposed$validPropMrDAG

    acceptReject <- acceptRejectDAG_cpp(
                      tXX, currentDAG,
                      proposed$proposedDAG, proposed$opNode,
                      proposed$opType, proposed$currentOpCard, proposed$proposedOpCard,
                      n, a, U, pp, temp[iter])

    if (acceptReject$isAccepted == TRUE)
    {
      currentDAG     <- proposed$proposedDAG
      logMargLik[iter] <- acceptReject$logMargLik

      if (iter > burnin)
        acceptPropDAG <- acceptPropDAG + 1

    } else {

      if (iter == 1L)
      {
        # Replicate original behaviour: initialise full vector with
        # the log-marginal likelihood of the empty starting DAG.
        lml_init   <- logMargLik_cpp(currentDAG, tXX, n, a, U)
        logMargLik <- logMargLik + lml_init
      } else {
        logMargLik[iter] <- logMargLik[iter - 1L]
      }
    }

    postPar <- rDAGWishart_cpp(currentDAG, a + n, Upost)

    if (iter > burnin & iter %% thin == 0L)
    {
      graphs[iterCounter] <- bd_encode(currentDAG)
      L[iterCounter]      <- bd_encode(postPar$L)
      D[iterCounter]      <- bd_encode(postPar$D)
      iterCounter <- iterCounter + 1L
    }

    utils::setTxtProgressBar(PBText, iter)
  }

  close(PBText)

  cat("\n\n")

  if (!is.null(MRconstraint))
    validPropMrDAG <- validPropMrDAG / niter
  else
    validPropMrDAG <- NULL

  acceptPropDAG <- acceptPropDAG / niter

  return(list(graphs = graphs, L = L, D = D, logMargLik = logMargLik,
             validPropMrDAG = validPropMrDAG, acceptPropDAG = acceptPropDAG))
}


########## acceptRejectDAG ##########

acceptRejectDAG <- function(tXX, currentDAG, proposedDAG,
                            node, opType, currentOpCard, proposedOpCard,
                            n, a, U, pp, temp,
                            currentLogML = NULL, Upost = NULL)
{
  if (sum(abs((currentDAG - proposedDAG))) > 0)
  {
    logpriorRatios   <- c(log(pp / (1 - pp)), log((1 - pp) / pp), log(1))
    logpriorRatio    <- logpriorRatios[opType]
    logproposalRatio <- log(currentOpCard) - log(proposedOpCard)
    
    if (opType != 3L)
    {
      currentlogMargLik  <- mlnodeDAGWishart(node, currentDAG,  tXX, n, a, U, Upost)
      proposedlogMargLik <- mlnodeDAGWishart(node, proposedDAG, tXX, n, a, U, Upost)
    } else {
      currentlogMargLik  <- mlnodeDAGWishart(node[1], currentDAG,  tXX, n, a, U, Upost) +
                            mlnodeDAGWishart(node[2], currentDAG,  tXX, n, a, U, Upost)
      proposedlogMargLik <- mlnodeDAGWishart(node[1], proposedDAG, tXX, n, a, U, Upost) +
                            mlnodeDAGWishart(node[2], proposedDAG, tXX, n, a, U, Upost)
    }
    
    acceptRatio <- min(0, 1 / temp * (proposedlogMargLik - currentlogMargLik +
                                      logpriorRatio + logproposalRatio))
    isAccepted  <- log(stats::runif(1)) < acceptRatio
    
    if (isAccepted)
    {
      logMargLik <- if (!is.null(currentLogML))
        currentLogML - currentlogMargLik + proposedlogMargLik
      else {
        lml <- 0
        for (nd in seq_len(nrow(U)))
          lml <- lml + mlnodeDAGWishart(nd, proposedDAG, tXX, n, a, U, Upost)
        lml
      }
    } else {
      logMargLik <- NULL
    }
  } else {
    isAccepted <- FALSE
    logMargLik <- NULL
  }
  
  return(list(isAccepted = isAccepted, logMargLik = logMargLik))
}

########## bd_decode ##########

bd_decode <- function(string, separator = ";") 
{
  vec4mat <- as.numeric(strsplit(string, separator)[[1]])
  q <- length(vec4mat)
  matrix(vec4mat, ncol = q)
  return(vec4mat)
}

########## bd_encode ##########

bd_encode <- function(matrix, separator = ";")
{
  matrix <- paste(matrix, collapse = separator)
  return(matrix)
}

########## genPredEffect_R ##########
# Pure-R version retained; use genPredEffect_cpp() directly for speed.

genPredEffect_R <- function(response, target, L, D)
{
  y    <- response
  L_I  <- L
  L_I[, target] <- 0
  diag(L_I)     <- 1
  Sigma_I <- solve(t(L_I)) %*% D %*% solve(L_I)
  effect  <- sapply(target, function(x) Sigma_I[x, y] / Sigma_I[x, x])
  return(effect)
}

########## fa ##########

fa <- function(node, DAG)
{
  pa <- which(DAG[, node] != 0)
  fa <- c(node, pa)
  return(fa)
}

########## mlnodeDAGWishart ##########

mlnodeDAGWishart <- function(node, DAG, tXX, n, a, U, Upost = NULL)
{
  j      <- node
  pa     <- pa(j, DAG)
  m      <- ncol(tXX)
  a.star <- (a + length(pa) - m + 1)
  if (is.null(Upost)) Upost <- U + tXX
  
  if (length(pa) == 0)
  {
    U_jj     <- U[j, j]
    Upost_jj <- Upost[j, j]
    prior.normcost <- -lgamma(a.star / 2) + a.star / 2 * log(U_jj / 2)
    post.normcost  <- -lgamma(a.star / 2 + n / 2) + (a.star / 2 + n / 2) * log(Upost_jj / 2)
  } else {
    U_paj.j    <- U[pa, j]
    U_jj       <- U[j, j] - t(U_paj.j) %*% solve(U[pa, pa]) %*% U_paj.j
    Upost_paj.j <- Upost[pa, j]
    Upost_jj   <- Upost[j, j] - t(Upost_paj.j) %*% solve(Upost[pa, pa]) %*% Upost_paj.j
    prior.normcost <- -lgamma(a.star / 2) + a.star / 2 * log(U_jj / 2) +
                      0.5 * log(det(as.matrix(U[pa, pa])))
    post.normcost  <- -lgamma(a.star / 2 + n / 2) + (a.star / 2 + n / 2) * log(Upost_jj / 2) +
                      0.5 * log(det(as.matrix(Upost[pa, pa])))
  }
  
  nodelml <- -n / 2 * log(2 * pi) + prior.normcost - post.normcost
  return(nodelml)
}

########## pa ##########

pa <- function(node, DAG)
{
  pa <- which(DAG[, node] != 0)
  return(pa)
}

########## rDAG ##########

rDAG <- function(m, w)
{
  DAG <- matrix(0, nrow = m, ncol = m)
  DAG[lower.tri(DAG)] <- stats::rbinom(n = m * (m - 1) / 2, size = 1, prob = w)
  return(DAG)
}

########## rDAGWishart ##########

rDAGWishart <- function(n, DAG, a, U)
{
  m    <- ncol(DAG)
  ajs  <- vapply(seq_len(m), function(j) { a + sum(DAG[, j] == 1) - m + 1 }, numeric(1))
  L.array <- array(0, dim = c(m, m, n))
  D.array <- array(0, dim = c(m, m, n))
  
  for (i in seq_len(n))
  {
    params <- lapply(seq_len(m), function(j) { rnodeDAGWishart(j, DAG, ajs[j], U) })
    sigmas <- vapply(seq_len(m), function(x) { params[[x]]$sigmaj }, numeric(1))
    L      <- lapply(seq_len(m), function(x) { params[[x]]$Lj })
    D.array[, , i] <- diag(sigmas)
    
    for (j in seq_len(m))
    {
      whc <- which(DAG[, j] == 1)
      L.array[whc, j, i] <- as.numeric(L[[j]])
    }
    diag(L.array[, , i]) <- 1
  }
  
  if (n == 1) { D.array <- D.array[, , 1]; L.array <- L.array[, , 1] }
  return(list(D = D.array, L = L.array))
}

########## rnodeDAGWishart ##########

rnodeDAGWishart <- function(node, DAG, aj, U)
{
  j   <- node
  pa  <- pa(j, DAG)
  out <- list(sigmaj = 0, Lj = 0)
  
  if (length(pa) == 0)
  {
    U_jj      <- U[j, j]
    out$sigmaj <- stats::rgamma(1, shape = aj / 2, rate = U_jj / 2) ^(-1)
  } else {
    U_paj.j   <- U[pa, j]
    invU_papa  <- solve(U[pa, pa])
    U_jj       <- U[j, j] - t(U_paj.j) %*% invU_papa %*% U_paj.j
    out$sigmaj  <- stats::rgamma(1, shape = aj / 2, rate = U_jj / 2) ^(-1)
    out$Lj      <- MASS::mvrnorm(1, -invU_papa %*% U_paj.j, out$sigmaj * invU_papa)
  }
  return(out)
}

########## .get_LD_from_1000G ##########

.get_LD_from_1000G <- function(rsids, pop = "EUR") {
  if (!requireNamespace("ieugwasr", quietly = TRUE))
    stop(
      "Package 'ieugwasr' is required for automatic LD download. ",
      "Install it with: install.packages('ieugwasr')",
      call. = FALSE
    )
  jwt <- ieugwasr::get_opengwas_jwt()
  if (identical(jwt, "") || is.null(jwt))
    stop(
      "No OpenGWAS API token found (required since 1 May 2024 for LD downloads via 'LD_rsids'). ",
      "To fix this:\n",
      "  1. Register/login at https://api.opengwas.io/ and generate a token from your profile page.\n",
      "  2. Add the line OPENGWAS_JWT=<your token> to your .Renviron file ",
      "(find its location with Sys.getenv(\"R_ENVIRON_USER\")).\n",
      "  3. Restart your R session.\n",
      "  4. Confirm it is recognised with ieugwasr::get_opengwas_jwt() (should return your token, not \"\").\n",
      "Alternatively, supply a pre-computed LD matrix directly via the 'LD' argument instead of 'LD_rsids'.",
      call. = FALSE
    )
  cat("## Downloading LD matrix for", length(rsids),
      "variants from 1000 Genomes (population:", pop, ")...\n")
  R <- tryCatch(
    suppressMessages(ieugwasr::ld_matrix(rsids, pop = pop, with_alleles = FALSE)),
    error = function(e)
      stop("LD download failed: ", conditionMessage(e),
           "\nCheck your rsIDs and internet connection.", call. = FALSE)
  )
  missing_ids <- setdiff(rsids, rownames(R))
  if (length(missing_ids) > 0L)
    warning(
      length(missing_ids), " variant(s) not found in the 1000G panel and were dropped: ",
      paste(missing_ids[seq_len(min(5L, length(missing_ids)))], collapse = ", "),
      if (length(missing_ids) > 5L) " ...",
      call. = FALSE
    )
  R
}

########## .display_round ##########

#' @title Round for on-screen display, falling back to signif() only where round() hits zero
#' @description Rounds to \code{digits} decimal places, which is more readable than significant-figure
#' rounding for typical-magnitude numbers. Wherever that rounds an element to exactly 0 (and the value
#' itself isn't 0), \code{signif(x, digits)} is used instead for that element only, so a small-but-real
#' value like 0.000106 is never displayed as a misleading 0
#' @param x A numeric vector, matrix, or data frame column
#' @param digits Decimal places for \code{round()} (also the significant figures for the signif() fallback)
#' @return \code{x} with the same shape, rounded for display

.display_round <- function(x, digits = 3) {
  r <- round(x, digits)
  zeroed <- !is.na(x) & x != 0 & r == 0
  if (any(zeroed)) r[zeroed] <- signif(x[zeroed], digits)
  r
}

#' @title Format a vector/matrix of .display_round() values for printing, without forced scientific notation
#' @description \code{print()} picks ONE shared format for an entire vector/matrix/data-frame column: if any single
#' element needs many decimal places (e.g. a .display_round() signif() fallback like 0.00000012), R silently switches
#' EVERY element in that same print call to scientific notation -- so normal values like 0.295 would suddenly show as
#' 2.95e-01 too. This formats every element independently (each keeps its own natural width) so the rest of the table
#' stays in clean, ordinary decimal form regardless of what any one outlier needs
#' @param x A numeric vector, matrix, or data frame column
#' @param digits Passed to \code{.display_round()}
#' @return A character vector/matrix, same shape as \code{x} (print matrices with \code{quote = FALSE})

.print_num <- function(x, digits = 3) {
  r <- .display_round(x, digits)
  s <- vapply(r, function(v) if (is.na(v)) NA_character_ else format(v, scientific = FALSE, trim = TRUE), character(1))
  if (!is.null(dim(r))) {
    dim(s) <- dim(r)
    dimnames(s) <- dimnames(r)
  } else {
    names(s) <- names(r)
  }
  s
}

########## .validate_LD ##########
## The LD matrix must be symmetric, with one row per IV. It is shrunk 1% towards the identity, R <- 0.99 R + 0.01 I,
## which makes it positive-definite (it is a positive semi-definite correlation matrix computed from genotypes) and
## bounds the weights (Omega^-1)_ii used below. The shrinkage is needed in practice: a reference panel with fewer
## individuals than IVs gives a singular matrix. If the matrix is still not positive-definite, the algorithm stops.

.validate_LD <- function(R, n_data) {
  p <- nrow(R)
  if (!isSymmetric(R, tol = 1e-8))
    stop("LD matrix is not symmetric.", call. = FALSE)
  if (p != n_data)
    stop("LD matrix dimension (", p, ") does not match the number of IVs (", n_data, ").", call. = FALSE)
  diag(R) <- 1

  ld_shrinkage <- 0.01
  R <- (1 - ld_shrinkage) * R + ld_shrinkage * diag(p)

  if (inherits(try(chol(R), silent = TRUE), "try-error"))
    stop("The LD matrix is not positive-definite, even after shrinkage. Check the LD panel and the list of SNPs ",
         "(the rsIDs and their alleles).", call. = FALSE)

  attr(R, "ld_shrinkage") <- ld_shrinkage
  R
}

########## .validate_gencov ##########
## Accepts gencov as:
##   (a) a list of n m x m matrices -> one per IV, e.g. the direct output of
##       MVMR::phenocov_mvmr() (heteroscedastic, matches per-SNP standard errors)
##   (b) an m x m x n array         -> equivalent to (a) in a single object
## A single matrix broadcast to every IV is deliberately not accepted: it is numerically
## equivalent to (a)/(b) for the main-analysis sufficient statistic under uniform per-IV
## weights, but invIV4MrDAG()'s leave-one-out folds and boot4MrDAG()'s bootstrap resamples
## touch individual IVs' own gencov terms separately, not just the aggregate sum, so a
## broadcast single matrix silently loses per-IV heterogeneity there.
## Only the OFF-DIAGONAL entries are used (the cross-trait sample-overlap covariance, from LD Score
## Regression). The main diagonal (each trait's squared standard error) is NOT corrected and is set to
## zero here, whatever was supplied: subtracting it distorts the results (it overstates the information in
## the data for traits with weak signal) without making the corrected matrix non-positive-definite.
## Returns a standardised m x m x n array so downstream code
## (.perIV_sufficientStat/.correct_gencov) has one shape to handle.

.validate_gencov <- function(gencov, m_data, n_data) {

  as_array_from_list <- function(lst, m_data, n_data) {
    if (length(lst) != n_data)
      stop(
        "'gencov' list/array must have one m x m matrix per IV (expected ", n_data,
        ", got ", length(lst), ").", call. = FALSE
      )
    arr <- array(0, dim = c(m_data, m_data, n_data))
    for (i in seq_len(n_data)) {
      Gi <- lst[[i]]
      if (!is.matrix(Gi) || nrow(Gi) != m_data || ncol(Gi) != m_data)
        stop("'gencov[[", i, "]]' is not an ", m_data, " x ", m_data, " matrix.", call. = FALSE)
      if (!isSymmetric(Gi, tol = 1e-8))
        stop("'gencov[[", i, "]]' is not symmetric.", call. = FALSE)
      arr[, , i] <- Gi
    }
    arr
  }

  if (is.list(gencov)) {

    gencov <- as_array_from_list(gencov, m_data, n_data)

  } else if (is.array(gencov) && length(dim(gencov)) == 3L) {

    if (!all(dim(gencov) == c(m_data, m_data, n_data)))
      stop(
        "'gencov' array must have dimension c(", m_data, ", ", m_data, ", ", n_data, ").",
        call. = FALSE
      )
    for (i in seq_len(n_data))
      if (!isSymmetric(gencov[, , i], tol = 1e-8))
        stop("'gencov[,,", i, "]' is not symmetric.", call. = FALSE)

  } else {
    stop(
      "'gencov' must be a list of n m x m matrices, or an m x m x n array ",
      "(one matrix per IV; a single matrix broadcast to every IV is not accepted).",
      call. = FALSE
    )
  }

  for (j in seq_len(m_data)) gencov[j, j, ] <- 0

  gencov
}

########## .correct_gencov ##########

.correct_gencov <- function(tXX, gencov_arr, weights) {
  n_data <- dim(gencov_arr)[3L]
  weighted_sum <- matrix(0, dim(gencov_arr)[1L], dim(gencov_arr)[2L])
  for (i in seq_len(n_data))
    weighted_sum <- weighted_sum + weights[i] * gencov_arr[, , i]
  .safeguardPD(tXX, weighted_sum)
}


########## .perIV_sufficientStat ##########
## Applies LD whitening and the gencov correction once, upfront, to the full
## set of n IVs, decomposed into one contribution per IV rather than already
## summed. Summing over all n slices (.sumTXX()) recovers the usual
## sufficient statistic tXX exactly; summing over a resampled (with
## replacement) subset of slices gives the sufficient statistic for that
## resample -- which is how boot4MrDAG() builds each bootstrap replicate
## directly from the corrected per-IV terms, instead of re-whitening and
## re-applying the gencov correction on a resampled (and potentially
## singular) LD structure for every draw
##
## With LD:    IV i's contribution is tcrossprod(X_w[i, ]), where
##             X_w = L_ld^{-1} %*% X is the GLS-whitened data (L_ld the lower
##             Cholesky factor of LD)
## Without LD: IV i's contribution is tcrossprod(X[i, ])
## gencov (if supplied) is subtracted per IV, weighted by (LD^{-1})_ii (or 1
## when LD is not used) -- it acts on the trait (column) axis, not the IV
## (row) axis:
##   tXXi[, , i] = tcrossprod(X_w[i, ]) - weights_i * gencov[, , i]

.perIV_sufficientStat <- function(X, LD = NULL, gencov = NULL) {
  n <- nrow(X)
  m <- ncol(X)

  if (!is.null(LD)) {
    R_chol  <- chol(LD)          # upper Cholesky: t(R_chol) %*% R_chol = LD
    LD_inv  <- chol2inv(R_chol)  # solve(LD), reusing the same factorisation
    weights <- diag(LD_inv)      # (LD^{-1})_ii, one weight per IV -- see
                                  # gencov_derivation.pdf: B ~ MN(0, LD, Sigma*)
                                  # (true signal correlates via LD) but
                                  # E ~ MN(0, I_n, Xi) (estimation noise
                                  # independent across IVs), so
                                  # E[t(E) LD^-1 E] = tr(LD^-1) * Xi in the
                                  # homoscedastic case, sum_i (LD^-1)_ii * Xi_i
                                  # per-IV -- NOT weight 1.
    L_ld    <- t(R_chol)         # lower Cholesky: L_ld %*% t(L_ld) = LD
    X_w     <- forwardsolve(L_ld, X)
  } else {
    weights <- rep(1, n)
    X_w     <- X
  }

  tXXi <- array(0, dim = c(m, m, n))
  for (i in seq_len(n))
    tXXi[, , i] <- tcrossprod(X_w[i, ])

  if (!is.null(gencov)) {
    corr <- array(0, dim = c(m, m, n))
    for (i in seq_len(n)) {
      corr[, , i] <- weights[i] * gencov[, , i]
      tXXi[, , i] <- tXXi[, , i] - corr[, , i]
    }
    attr(tXXi, "correction") <- corr   # the per-IV correction, kept so that a resample can be safeguarded (see .safeguardPD())
  }

  tXXi
}


########## .sumTXX ##########
## Sums the per-IV sufficient statistic contributions (see .perIV_sufficientStat()) over a set of IVs, possibly with repeats as in a
## bootstrap draw (idx), and safeguards the result (see .safeguardPD())

.sumTXX <- function(tXXi, idx = NULL) {
  if (is.null(idx)) idx <- seq_len(dim(tXXi)[3L])
  tXX  <- rowSums(tXXi[, , idx, drop = FALSE], dims = 2)
  corr <- attr(tXXi, "correction")
  if (is.null(corr)) return(.safeguardPD(tXX, NULL))
  C <- rowSums(corr[, , idx, drop = FALSE], dims = 2)
  .safeguardPD(tXX + C, C)
}


########## .safeguardPD ##########
## The corrected sufficient statistic is T = T_raw - C, where T_raw is the cross-product of the data and C the sample-overlap
## correction (off-diagonal terms only). T_raw is positive-definite, and so is T in every case we analysed. As a safeguard, if
## T were not positive-definite, only the fraction alpha of the correction that keeps it positive-definite is subtracted:
## T_raw - alpha C is positive-definite if and only if alpha * lambda_max(R^{-T} C R^{-1}) < 1, where T_raw = R'R (Cholesky),
## so alpha = 0.99 / lambda_max, one eigenvalue. alpha = 1 (the full correction) whenever T is positive-definite. The
## fraction applied is returned in attr "alpha", for the caller to report.

.safeguardPD <- function(Traw, C = NULL) {
  Traw <- (Traw + t(Traw)) / 2
  is_pd <- function(M) !inherits(try(chol(M), silent = TRUE), "try-error")
  if (!is_pd(Traw))
    stop("The uncorrected cross-product of the data is not positive-definite. Check the input data.", call. = FALSE)
  if (is.null(C) || all(C == 0)) { attr(Traw, "alpha") <- 1; return(Traw) }
  C  <- (C + t(C)) / 2
  Tc <- Traw - C
  if (is_pd(Tc)) { attr(Tc, "alpha") <- 1; return(Tc) }
  R <- chol(Traw)
  A <- forwardsolve(t(R), C)
  M <- t(forwardsolve(t(R), t(A)))
  lam   <- max(eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  alpha <- min(1, 0.99 / lam)
  Tc    <- Traw - alpha * C
  while (!is_pd(Tc) && alpha > 1e-8) { alpha <- alpha / 2; Tc <- Traw - alpha * C }
  attr(Tc, "alpha") <- alpha
  Tc
}
