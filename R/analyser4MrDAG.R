# 135 characters #####################################################################################################################
#' @title Analyzer MrDAG: Detection of invalid instrumental variables and bootstrap evaluation of the Bayesian causal graphical model 
#' for joint Mendelian randomization (MR) analysis of multiple exposures and outcomes
#' @description Detection of invalid instrumental variables (IVs) due to unmeasured pleiotropy, in which some genetic variants are 
#' directly associated with the outcomes and not via the exposures. After removing invalid IVs, a bootstrap assessment of the 
#' robustness to weak IVs of \code{MrDAG} estimated Posterior Edge Inclusion Probabilities (PEIPs) and estimated genetically predicted effects is 
#' performed
#'
#' @param data Number of observations (IVs in a summary-level MR design) times the number of traits (both outcomes and exposures).
#' See also the \code{MRconstraint} argument
#' @param ninvIV Number of invalid IV iterations. At each iteration, the detection of invalid IVs is performed by using the 
#' conditional predictive ordinate (CPO) conditionally on the current set of IVs. If, after an iteration, some IVs are removed, in 
#' the next iteration, invalid IVs are detected conditionally on the set of remaining IVs. Default (\code{ninvIV=1}) is when only 
#' one iteration for the invalid IV detection procedure is performed, although other invalid IVs might be present. If 
#' \code{ninvIV=NULL}, the procedure terminates when no more invalid IVs are detected
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
#' @param LD,LD_rsids,LD_pop An optional n-by-n LD correlation matrix for the IVs (or rsIDs to download one, resolved once
#' upfront). See \code{\link{MrDAG}}
#' @param gencov An optional per-IV sample-overlap correction. See \code{\link{MrDAG}}
#' @param fileName Name of the file for the \code{analyzer4MrDAG} output object (default \code{"analyzer4MrDAG_object"}). If
#' \code{NULL}, the output is not saved
#' @param filePath Directory where the \code{analyzer4MrDAG} object is saved. If the directory does not exist, 
#' \code{analyzer4MrDAG} will create it by using the current working directory as the root directory. If the directory is not 
#' specified (default \code{NULL}), the \code{analyzer4MrDAG} object is saved in the current working directory
#' @param displayIntro Print information regarding the \code{analyzer4MrDAG} model on the screen (default \code{TRUE})

#' @param seed Seed to be used in the initialization of the MCMC algorithm (default \code{3112021}) and updates at each bootstrap 
#' (+1) on each summary-level dataset
#'
#' @details
#' For details regarding the model and the algorithm, see \insertCite{Bottolo2026;textual}{MrDAGSuite} and 
#' \insertCite{Zuber2025;textual}{MrDAGSuite}. For the detection of outliers or high-leverage and influential IVs, see 
#' \insertCite{Congdon2005;textual}{MrDAGSuite}. For bootstrap assessment criteria, see \insertCite{Glymour2019;textual}{MrDAGSuite} 
#'
#' @export
#'
#' @return The value returned is a list object \code{list(outputMrDAG, outputInvIV, outputBoot, invalidIV, freqInvalidIV,
#' timeAnalyzer4MrDAG, hyperPar, samplerPar, correctionPar, opt)}:
#' \describe{
#'   \item{\code{outputMrDAG}}{ List object \code{list(graphs, L, D, logMargLik, validPropMrDAG, acceptPropDAG)} obtained after 
#'         removing invalid IVs, for a specified number of iterations defined in \code{ninvIV}. See the \code{MrDAG} function for 
#'         details}
#'   \item{\code{outputInvIV}}{ List object \code{list(invalidIV, freqInvalidIV, logMargLik, logMargLik_i, CPO, scaleCPO, opt)} of 
#'         the last iteration of the invalid IV detection procedure, for a specified number of iterations defined in \code{ninvIV}. 
#'         See the \code{invalidIV4MrDAG} function for details}
#'   \item{\code{outputBoot}}{ List object \code{list(bootEdgeFreq, bootEdgeFreq_LL, bootEdgeFreq_UL, bootGenPredEffects, 
#'         bootGenPredEffects_LL, bootGenPredEffects_UL, bootFreq)} obtained after removing invalid IVs, as specified in \code{ninvIV}. 
#'         See the \code{boot4MrDAG} function for details}
#'   \item{\code{invalidIV}}{ Detected invalid IVs with a scale CPO < 0.01 (see \insertCite{Congdon2005;textual}{MrDAGSuite}) in 
#'         each iteration of the invalid IV detection procedure }
#'   \item{\code{freqInvalidIV}}{ Frequency of invalid IV iterations defined as outliers or high-leverage and influential IVs 
#'         across all iterations of the invalid IV detection procedure }
#'   \item{\code{timeAnalyzer4MrDAG}}{ Time in minutes employed by the \code{analyzer4MrDAG} algorithm to analyze the data }
#'   \item{\code{hyperPar}}{ List of the hyper-parameters \code{list(pp, a, U, tempMax)} and, if specified, the
#'         \code{MRconstraint} list }
#'   \item{\code{samplerPar}}{ List of parameters used in the invalid IV detection, bootstrap analysis and MCMC algorithm
#'         \code{list(ninvIV, nboot, CIboot, niter, burnin, thin)} }
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
#' # Detection of invalid instrumental variables (IVs) due to unmeasured pleiotropy, in which some
#' # genetic variants are directly associated with the outcomes and not via the exposures.
#' # Subsequent assessment of the robustness to weak IVs of the estimated Posterior Probabilities
#' # of Edge Inclusion (PEIPs) and the estimated genetically predicted effects under
#' # intervention on mental
#' # health phenotypes that might impact lifestyle and behavioural traits. 470 independent IVs
#' # were selected to be associated at genome-wide significance with seven mental health
#' # phenotypes that are considered exposures for the risk of six lifestyle and behavioural traits
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
#' # analyzer4MrDAG is run to detect invalid IVs with two iterations of the procedure. After
#' # removing invalid IVs, it generates the bootstrap frequencies of edge inclusion and the
#' # bootstrap genetically predicted effects as the average of the estimated PEIPs and the
#' # estimated genetically predicted
#' # effects, respectively, by running the \code{MrDAG} algorithm on five bootstrap summary-
#' # level datasets. Results are aggregated by the average vote. Each instance of the
#' # MrDAG algorithm generates 2,000 posterior samples of all unknowns, with an overall 5,000
#' # posterior samples
#'
#' output <- analyzer4MrDAG(data = MD2LBT_data$data, gencov = MD2LBT_data$gencov,
#'                          ninvIV = 5, nboot = 50,
#'                          niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20,
#'                          MRconstraint = MRconstraint, fileName = NULL)
#'
#' # Finally, PEIPs are shown with an order of appearance that is changed within the outcomes and
#' # within the exposures compared to the order of the traits in the analyzer4MrDAG output
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' plot_postEdgeProb(output$outputMrDAG, ord = ord,
#'               fileName = "Example_alternative_plot_postEdgeProb")
#'
#' # as well as bootstrap frequencies of edge inclusion
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' plot_bootEdgeFreq(output$outputBoot, ord = ord, 
#'                   fileName = "Example_plot_bootEdgeFreq")
#' }


########## analyzer4MrDAG ##########

analyzer4MrDAG <- function(data, ninvIV = 1, nboot = 100, CIboot = 0.95,
                           niter, burnin, thin,
                           pp = 0.05, a = NULL, U = NULL, tempMax = 10,
                           MRconstraint = NULL,
                           LD = NULL, LD_rsids = NULL, LD_pop = "EUR",
                           gencov = NULL,
                           fileName = "analyzer4MrDAG_object", filePath = NULL, displayIntro = TRUE,
                           seed = 31122021)
{
  if (!is.null(MRconstraint))
  {
    if (max(MRconstraint$Y_idx) >= min(MRconstraint$X_idx))
      stop("Outcome indices (Y_idx) must precede exposure indices (X_idx)")
  }

  if (displayIntro)
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
      cat("## Analyzer for Directed Acyclic Graph exploration\n")
    } else {
      cat("## Analyzer for Mendelian randomization with Directed Acyclic Graphs exploration\n")
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

  displayIntro <- FALSE
  startTime <- Sys.time()

  ########## LD validation and optional download (once, before the pipeline) ##########

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

  if (!is.null(LD)) {
    ld_shr   <- attr(LD, "ld_shrinkage")
    attr(LD, "ld_shrinkage") <- NULL
    ld_details <- c(
      if (!is.null(LD_rsids)) paste0("population: ", LD_pop),
      paste0("shrunk ", 100 * ld_shr, "% towards the identity")
    )
    cat("## LD correction: enabled",
        if (length(ld_details) > 0) paste0(" (", paste(ld_details, collapse = "; "), ")"),
        "\n", sep = "")
  } else {
    cat("## LD correction: none (IVs assumed independent)\n")
  }

  X <- data
  n <- dim(X)[1]
  m <- dim(X)[2]

  if (!is.null(gencov))
    gencov <- .validate_gencov(gencov, m, n)   # standardised once to an m x m x n array

  if (!is.null(gencov)) {
    cat("## Sample overlap correction (gencov): enabled\n")
  } else {
    cat("## Sample overlap correction (gencov): none (traits assumed estimated from non-overlapping samples)\n")
  }
  cat("##\n")

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
  
  samplerPar <- NULL
  samplerPar$ninvIV <- ninvIV
  samplerPar$nboot <- nboot
  samplerPar$CIboot <- CIboot
  samplerPar$niter <- niter
  samplerPar$burnin <- burnin
  samplerPar$thin <- thin
  
  X_TMP <- X
  LD_TMP <- LD
  gencov_TMP <- gencov
  orig_idx_TMP <- seq_len(n)   # which original rows of X are still present in X_TMP
  invIV_round <- 1
  invalidIV <- NULL
  freqInvalidIV <- NULL
  freqInvalidIV_TMP <- 1

  if (is.null(ninvIV))
  {
    ninvIV <- Inf
  }

  startTime <- Sys.time()

  cat("\n############################\n")
  cat("## Step 1: MrDAG analysis ##\n")
  cat("############################\n")

  outputMrDAG <- MrDAG(data = X_TMP,
                       niter = niter, burnin = burnin, thin = thin,
                       pp = pp, a = a, U = U, tempMax = tempMax,
                       MRconstraint = MRconstraint,
                       LD = LD_TMP, gencov = gencov_TMP,
                       fileName = NULL, filePath = NULL, displayIntro = displayIntro,
                       seed = seed)

  while (freqInvalidIV_TMP > 0 & invIV_round <= ninvIV)
  {


    cat("#################################################\n")
    cat("## Step 2.1: Detection of invalid IVs, round", invIV_round, "##\n")
    cat("#################################################\n")

    outputInvIV <- invIV4MrDAG(data = X_TMP, outputMrDAG, LD = LD_TMP, gencov = gencov_TMP,
                               displayIntro = displayIntro)
    freqInvalidIV_TMP <- outputInvIV$freqInvalidIV

    if (freqInvalidIV_TMP > 0)
    {
      # outputInvIV$invalidIV are indices into this round's X_TMP, not the
      # original X -- translate through orig_idx_TMP before accumulating, or
      # rounds after the first would remove the wrong rows from X.
      invalidIV <- c(invalidIV, orig_idx_TMP[outputInvIV$invalidIV])
      freqInvalidIV <- c(freqInvalidIV, outputInvIV$freqInvalidIV)
      orig_idx_TMP <- orig_idx_TMP[-outputInvIV$invalidIV]
      X_TMP <- X[orig_idx_TMP, ]
      # LD/gencov re-subsetted the same way, from the same original LD/gencov,
      # so they stay row-aligned with X_TMP
      LD_TMP <- if (!is.null(LD)) LD[orig_idx_TMP, orig_idx_TMP, drop = FALSE] else NULL
      gencov_TMP <- if (!is.null(gencov)) gencov[, , orig_idx_TMP, drop = FALSE] else NULL

      cat("##################################################################\n")
      cat("## Step 2.2: MrDAG analysis after removing invalid IVs, round", invIV_round, "##\n")
      cat("##################################################################\n")

      outputMrDAG <- MrDAG(data = X_TMP,
                           niter = niter, burnin = burnin, thin = thin,
                           pp = pp, a = a, U = U, tempMax = tempMax,
                           MRconstraint = MRconstraint,
                           LD = LD_TMP, gencov = gencov_TMP,
                           fileName = NULL, filePath = NULL, displayIntro = displayIntro,
                           seed = seed)
    }

  invIV_round <- invIV_round + 1
  }
  
  if (is.null(invalidIV))
  {
    invalidIV <- "No invalid IVs have been detected"
    freqInvalidIV <- 0
  }
  
  cat("########################################################################\n")
  cat("## Step 3: Bootstrap analysis for weak IVs after removing invalid IVs ##\n")
  cat("########################################################################\n")
  
  outputBoot <- boot4MrDAG(data = X_TMP, nboot = nboot, CIboot = CIboot,
                           niter, burnin, thin,
                           pp = pp, a = NULL, U = NULL, tempMax = tempMax,
                           MRconstraint = MRconstraint,
                           LD = LD_TMP, gencov = gencov_TMP,
                           fileName = NULL, filePath = NULL, displayIntro = displayIntro,
                           seed = seed)
  
  opt <- list(seed = seed)
  endTime <- Sys.time()
  timeAnalyzer4MrDAG <- as.numeric(difftime(time1 = endTime, time2 = startTime, units = "min"))
  timeAnalyzer4MrDAG <- paste0(round(timeAnalyzer4MrDAG, digits = 3), "m")

  correctionPar             <- list()
  correctionPar$LD_used     <- !is.null(LD)
  correctionPar$LD_pop      <- if (!is.null(LD)) LD_pop else NULL
  correctionPar$gencov_used <- !is.null(gencov)

  output <- list(outputMrDAG = outputMrDAG, outputInvIV = outputInvIV, outputBoot = outputBoot,
                 invalidIV = invalidIV, freqInvalidIV = freqInvalidIV,
                 timeAnalyzer4MrDAG = timeAnalyzer4MrDAG,
                 hyperPar = hyperPar,
                 samplerPar = samplerPar,
                 correctionPar = correctionPar,
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
