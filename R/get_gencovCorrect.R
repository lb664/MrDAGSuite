# 135 characters #####################################################################################################################
#' @title Magnitude of the sample-overlap (gencov) correction
#' @description Quantifies how much of the raw sufficient statistic (crossprod of the data, and if \code{LD} is supplied,
#' GLS-whitened) is attributable to the \code{gencov} sample-overlap correction actually applied inside \code{MrDAG}: overall
#' (Frobenius norm) and the mean per-IV \code{gencov} matrix expressed as a
#' correlation (which trait pairs actually share samples). Built directly on the package's own \code{.validate_LD()},
#' \code{.validate_gencov()}, and \code{.perIV_sufficientStat()} internals, so the reported correction is exactly the one \code{MrDAG}
#' itself would apply, not an approximation of it
#'
#' @param data Number of observations (IVs) times the number of traits (both outcomes and exposures), as passed to \code{MrDAG}
#' @param LD Optional n-by-n LD correlation matrix for the IVs, as passed to \code{MrDAG}. \code{NULL} means IVs are independent
#' @param gencov Sample-overlap correction, as passed to \code{MrDAG}: a single m-by-m matrix, a list of n m-by-m matrices, or an
#' m-by-m-by-n array. \code{NULL} means no correction (a message is printed and \code{NULL} is returned)
#' @param digits Number of decimal places used to round every numeric value in the returned summary. If \code{NULL} (default),
#' the calling script's own \code{digits} variable is used if one exists in its environment; otherwise the summary is returned
#' at full, unrounded precision
#'
#' @return An object of class \code{gencovCorrect}, a list \code{list(overall, mean_gencov_correlation)}, or \code{NULL} if \code{gencov} is \code{NULL}. Printing the object shows the
#' overall correction and the mean per-IV sample-overlap correlations between the traits:
#' \describe{
#'   \item{\code{overall}}{ \code{list(frobenius_no_gencov, frobenius_with_gencov, frobenius_correction, pct_of_signal)}: Frobenius
#'         norm of the sufficient statistic without and with the correction applied, the norm of the correction itself, and the
#'         correction as a percentage of the uncorrected sufficient statistic's magnitude }
#'   \item{\code{mean_gencov_correlation}}{ The mean per-IV \code{gencov} matrix, expressed as a correlation matrix -- shows which
#'         trait pairs actually share samples, independent of any LD reweighting }
#' }
#'
#' @param x Object returned by \code{get_gencovCorrect}
#' @param label Optional label (for example the name of the dataset) shown in the heading when the object is printed
#' @param ... Not used
#'
#' @examples
#' \dontrun{
#' data(CardioNHS_data)
#' gencovCorrect <- get_gencovCorrect(CardioNHS_data$data, gencov = CardioNHS_data$gencov)
#' gencovCorrect
#' print(gencovCorrect, label = "CardioNHS")
#' }
#'
#' @export
#' @rdname get_gencovCorrect
get_gencovCorrect <- function(data, LD = NULL, gencov, digits = NULL) {

  if (is.null(gencov)) {
    message("No 'gencov' supplied -- no sample-overlap correction applied.")
    return(invisible(NULL))
  }

  # NB: can't just do exists("digits", inherits = TRUE) here -- `digits` is
  # itself a formal argument of this function, so that would always find
  # this very parameter (currently NULL) before ever reaching the caller's
  # own `digits`. Look in the CALLER's frame specifically instead. If the
  # caller has no `digits` of its own, leave it NULL -- no rounding, full
  # precision is returned.
  if (is.null(digits)) {
    caller_env <- parent.frame()
    if (exists("digits", envir = caller_env, inherits = FALSE))
      digits <- get("digits", envir = caller_env)
  }

  X <- as.matrix(data)
  n <- nrow(X)
  m <- ncol(X)

  ## Exactly mirrors MrDAG.R's own sufficient-statistic computation (same LD validation, same
  ## gencov validation, same .perIV_sufficientStat() whitening/correction), so the "correction"
  ## reported here is the one actually applied, not a proxy.
  gencov_raw <- gencov
  if (!is.null(LD)) LD <- .validate_LD(LD, n)
  gencov <- .validate_gencov(gencov, m, n)

  tXXi_no_gencov   <- .perIV_sufficientStat(X, LD = LD, gencov = NULL)
  tXXi_with_gencov <- .perIV_sufficientStat(X, LD = LD, gencov = gencov)

  tXX_no_gencov   <- rowSums(tXXi_no_gencov,   dims = 2)
  tXX_with_gencov <- rowSums(tXXi_with_gencov, dims = 2)
  correction      <- tXX_no_gencov - tXX_with_gencov

  overall <- list(
    frobenius_no_gencov   = norm(tXX_no_gencov, "F"),
    frobenius_with_gencov = norm(tXX_with_gencov, "F"),
    frobenius_correction  = norm(correction, "F"),
    pct_of_signal         = 100 * norm(correction, "F") / norm(tXX_no_gencov, "F")
  )

  ## Mean per-IV gencov, as a correlation -- unweighted by LD, since this is
  ## about which traits share samples, not about the reweighted correction.
  ## se_diag comes from the raw gencov's diagonal (each trait's squared standard error).
  raw_diag_mean <- if (is.list(gencov_raw)) {
    vapply(seq_len(m), function(j)
      mean(vapply(gencov_raw, function(G) G[j, j], numeric(1))), numeric(1))
  } else {
    vapply(seq_len(m), function(j) mean(gencov_raw[j, j, ]), numeric(1))
  }
  mean_gencov <- apply(gencov, c(1, 2), mean)
  dimnames(mean_gencov) <- list(colnames(X), colnames(X))
  se_diag     <- sqrt(raw_diag_mean)
  mean_gencov_correlation <- mean_gencov / outer(se_diag, se_diag)
  dimnames(mean_gencov_correlation) <- dimnames(mean_gencov)

  if (!is.null(digits)) {
    overall   <- lapply(overall, round, digits = digits)
    mean_gencov_correlation <- round(mean_gencov_correlation, digits = digits)
  }

  structure(list(overall = overall, mean_gencov_correlation = mean_gencov_correlation), class = "gencovCorrect")
}


#' @rdname get_gencovCorrect
#' @method print gencovCorrect
#' @export
print.gencovCorrect <- function(x, label = NULL, ...) {
  cat("\n----  ", if (!is.null(label)) paste0(label, ": "), "sample-overlap (gencov) correction ----\n", sep = "")

  cat("Correction is", .display_round(x$overall$pct_of_signal, 3),
      "% of the uncorrected sufficient statistic's magnitude (off-diagonal sample overlap only; Frobenius norm",
      .display_round(x$overall$frobenius_correction, 3), "of", .display_round(x$overall$frobenius_no_gencov, 3), ")\n")

  cat("\nMean per-IV gencov, as a correlation (which trait pairs share samples):\n")
  print(.print_num(x$mean_gencov_correlation, 3), quote = FALSE)

  invisible(x)
}
