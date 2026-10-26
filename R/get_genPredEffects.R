# 135 characters #####################################################################################################################
#' @title Estimation of the genetically predicted effects under intervention on the exposures
#' @description Estimation of the genetically predicted effects under intervention on the exposures (targets) and measured on the outcomes 
#' (responses) based on the Directed Acyclic Graphs (DAGs) explored by the \code{MrDAG} algorithm
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param BMA If \code{TRUE} (default), Bayesian Model Averaging of the estimated genetically predicted effects across all explored DAGs in the 
#' visited Completed Partially Directed Acyclic Graph (CPDAG) or Essential Graph (EG) is performed
#' @param CI Level (default \code{0.95}) of the credible intervals (CIs) for the genetically predicted effects. If \code{BMA = TRUE}, they are 
#' calculated based on suitable quantiles of the estimated genetically predicted effects across all explored DAGs
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when printing the genetically predicted effects. If \code{ord = NULL}, the ordering follows \code{MRconstraint} in the \code{MrDAG} output.
#' @param progress Logical value set to \code{FALSE} by default to print on the screen the progress of the genetically predicted effects estimation
#'
#' @export
#'
#' @return The value returned is a list object \code{list(genPredEffects, genPredEffects_LL, genPredEffects_UL, group, BMA, CI, ord)}:
#' \describe{
#'   \item{\code{genPredEffects}}{ Estimate of the genetically predicted effects under intervention on the exposures based on the DAGs explored by 
#'         the \code{MrDAG} algorithm. A feedback loop genetically predicted effect of the same outcome is not allowed (NA) }
#'   \item{\code{genPredEffects_LL}}{ Lower limit of the \code{CI} for the genetically predicted effects. It is calculated as the 
#'         (1-\code{CIboot})/2 % quantile of the genetically predicted effects across all explored DAGs }
#'   \item{\code{genPredEffects_UL}}{ Upper limit of the \code{CI} for the genetically predicted effects. It is calculated as the 
#'         1-\[(1-\code{CIboot})/2\] % quantile of the genetically predicted effects across all explored DAGs }
#'   \item{\code{group}}{ If \code{MRconstraint} in the \code{MrDAG} algorithm contains the indices of the traits in \code{data} that 
#'         define the outcomes and the exposures, \code{group} coincides with this list }
#'   \item{\code{BMA}}{ Logical option }
#'   \item{\code{CI}}{ Level of the credible interval option }
#'   \item{\code{ord}}{ Indices of the traits in \code{data} that specify the outcomes and the exposures. It might differ from 
#'         \code{group} if \code{ord} has been specified and it is different from \code{MRconstraint} } }
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # Estimation of the genetically predicted effects under intervention on lifestyle and
#' # behavioural exposures,
#' # and measured on mental health phenotypes. 708 independent instrumental variables (IVs) were 
#' # selected to be associated at genome-wide significance with six lifestyle and behavioural 
#' # traits after pruning or clumping, which are considered exposures for the risk of seven mental 
#' # health outcome phenotypes
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
#'                 niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                 MRconstraint = MRconstraint, fileName = NULL, 
#')
#'
#' # Finally, the posterior genetically predicted effects and 90% credible intervals are
#' # estimated with an order
#' # of appearance that is changed within the outcomes and within the exposures compared to the 
#' # order of the traits in the MrDAG output
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' genPredEffects <- get_genPredEffects(output, CI = 0.90, ord = ord, progress = TRUE)
#' }


########## get_genPredEffects ##########

get_genPredEffects <- function(output, BMA = TRUE, CI = 0.95, ord = NULL, progress = FALSE)
{
  niter <- output$samplerPar$niter
  burnin <- output$samplerPar$burnin
  thin <- output$samplerPar$thin
  niter_TMP <- (niter - burnin) / thin
  q <- length(output$hyperPar$MRconstraint$Y_idx)
  p <- length(output$hyperPar$MRconstraint$X_idx)
  colnames_YX <- colnames(output$hyperPar$U)
  
  if (is.null(dim(output$graphs)))
  {
    m <- sqrt(length(bd_decode(output$graphs[1])))
  } else {
    m <- dim(output$graphs)[1]
  }
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    group <- NULL
    group$Outcomes <- output$hyperPar$MRconstraint$Y_idx
    group$Exposures <- output$hyperPar$MRconstraint$X_idx
  } else {
    group <- NULL
  }
  
  if (!is.null(ord))
  {
    
    if (!is.null(output$hyperPar$MRconstraint))
    {
      Y_ord <- ord[ord <= q]
      X_ord <- ord[ord > q]
    } else {
      Y_ord <- ord
      X_ord <- NULL
    }
  
  } else {
    
    if (!is.null(output$hyperPar$MRconstraint))
    {
      Y_ord <- output$hyperPar$MRconstraint$Y_idx
      X_ord <- output$hyperPar$MRconstraint$X_idx
    } else {
      Y_ord <- seq(1, m, 1)
      X_ord <- NULL
    }
  
  }
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    dimnames <- list(colnames_YX, colnames_YX[output$hyperPar$MRconstraint$Y_idx])
    
    if (BMA == TRUE)
    {
      genPredEffects <- matrix(NA, nrow = m, ncol = q, dimnames = dimnames)
      genPredEffects_LL <- matrix(NA, nrow = m, ncol = q, dimnames = dimnames)
      genPredEffects_UL <- matrix(NA, nrow = m, ncol = q, dimnames = dimnames)
    } else {
      genPredEffects <- array(NA, dim = c((niter - burnin) / thin, m, q), dimnames = list(paste0("iter", seq(burnin + thin, niter, thin)), 
                                                                                         colnames_YX, colnames_YX[output$hyperPar$MRconstraint$Y_idx])) 
    }
    
    for (k in 1 : q)
    {
      response_idx <- k
      target_idx <- (1 : m)[-k]
      genPredEffect_TMP <- get_genPredEffect(output, response = response_idx, target = target_idx, BMA = BMA, CI = CI)
      
      if (BMA == TRUE)
      {
        genPredEffects[, k] <- append(genPredEffect_TMP$genPredEffect, NA, after = (k - 1))
        genPredEffects_LL[, k] <- append(genPredEffect_TMP$genPredEffect_LL, NA, after = (k - 1))
        genPredEffects_UL[, k] <- append(genPredEffect_TMP$genPredEffect_UL, NA, after = (k - 1))
      } else {
        genPredEffects[, target_idx, k] <- genPredEffect_TMP$genPredEffect
      }
      
      if (!is.null(progress))
      {
      
        if (progress == TRUE)
        {
          
          if (k == 1)
          {
            cat("\nEstimation (average) genetically predicted effect under intervention on exposure X and measured on outcome", colnames_YX[k], "\n")
          } else if (k > 1 & k < q) {
            cat("Estimation (average) genetically predicted effect under intervention on exposure X and measured on outcome", colnames_YX[k], "\n")
          } else if (k == q) {
            cat("Estimation (average) genetically predicted effect under intervention on exposure X and measured on outcome", colnames_YX[k], "\n\n")
          }
        
        } else {
          
          if (k == 1)
          {
            cat("\n")
            cat("Estimation genetically predicted effects under intervention on the exposure\n\n")
          }
        
        }
      
      }
    
    }
    
    if (BMA == TRUE)
    {
      genPredEffects <- genPredEffects[X_ord, Y_ord]
      genPredEffects_LL <- genPredEffects_LL[X_ord, Y_ord]
      genPredEffects_UL <- genPredEffects_UL[X_ord, Y_ord]
    } else {
      genPredEffects <- genPredEffects[, X_ord, Y_ord]
    }
  
  } else {
    dimnames <- list(colnames_YX, colnames_YX)
    
    if (BMA == TRUE)
    {
      genPredEffects <- matrix(NA, nrow = m, ncol = m, dimnames = dimnames)
      genPredEffects_LL <- matrix(NA, nrow = m, ncol = m, dimnames = dimnames)
      genPredEffects_UL <- matrix(NA, nrow = m, ncol = m, dimnames = dimnames)
    } else {
      genPredEffects <- array(NA, dim = c((niter - burnin) / thin, m, m), dimnames = list(paste0("iter", seq(burnin + thin, niter, thin)), 
                                                                                         colnames_YX, colnames_YX)) 
    }
    
    for (k in 1 : m)
    {
      response_idx <- k
      target_idx <- (1 : m)[-k]
      genPredEffect_TMP <- get_genPredEffect(output, response = response_idx, target = target_idx, BMA = BMA, CI = CI)
      
      if (BMA == TRUE)
      {
        genPredEffects[, k] <- append(genPredEffect_TMP$genPredEffect, NA, after = (k - 1))
        genPredEffects_LL[, k] <- append(genPredEffect_TMP$genPredEffect_LL, NA, after = (k - 1))
        genPredEffects_UL[, k] <- append(genPredEffect_TMP$genPredEffect_UL, NA, after = (k - 1))
      } else {
        genPredEffects[, target_idx, k] <- genPredEffect_TMP$genPredEffect
      }
      
      if (!is.null(progress))
      {
        
        if (progress == TRUE)
        {
        
          if (k == 1)
          {
            cat("\nEstimation genetically predicted effects under intervention on a trait and measured on trait", colnames_YX[k], "\n")
          } else if (k > 1 & k < m) {
            cat("Estimation genetically predicted effects under intervention on a trait and measured on trait", colnames_YX[k], "\n")
          } else if (k == m) {
            cat("Estimation genetically predicted effects under intervention on a trait and measured on trait", colnames_YX[k], "\n\n")
          }
        
        } else {
           
          if (k == 1)
          {
            cat("\n")
            cat("Estimation genetically predicted effects under intervention on the traits\n\n")
          }
         
        }
      
      }
    
    }
    
    if (BMA == TRUE)
    {
      genPredEffects <- genPredEffects[Y_ord, Y_ord]
      genPredEffects_LL <- genPredEffects_LL[Y_ord, Y_ord]
      genPredEffects_UL <- genPredEffects_UL[Y_ord, Y_ord] 
    } else {
      genPredEffects <- genPredEffects[, Y_ord, Y_ord]
    }
  
  }
  
  if (BMA == TRUE)
  {
    return(list(genPredEffects = genPredEffects, genPredEffects_LL = genPredEffects_LL, genPredEffects_UL = genPredEffects_UL, group = group, BMA = BMA, CI = CI, ord = ord))
  } else {
    return(list(genPredEffects = genPredEffects, group = group, BMA = BMA, CI = CI, ord = ord))
  }

}
