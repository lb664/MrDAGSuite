# 135 characters #####################################################################################################################
#' @title Estimation of the marginal genetically predicted effect under intervention on an exposure
#' @description Estimation of the marginal genetically predicted effect under intervention on a trait (target) and measured on another one (outcome) 
#' based on the Directed Acyclic Graphs (DAGs) explored by the \code{MrDAG} algorithm
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param response Trait (outcome) where the effect of the intervention on another trait (target) is measured. The index refers to 
#' the position of the trait (outcome) in the output produced by the \code{MrDAG} algorithm
#' @param target Trait (target) under intervention. The index refers to the position of the trait (target) in the output produced by 
#' the \code{MrDAG} algorithm
#' @param BMA If (default) \code{TRUE}, Bayesian Model Averaging of the estimated (average) genetically predicted effect across all explored DAGs in
#' the visited Completed Partially Directed Acyclic Graph (CPDAG) or Essential Graph (EG) is performed
#' @param CI Level (default \code{0.95}) of the credible interval (CI) for the (average) genetically predicted effect. If \code{BMA = TRUE}, it is
#' calculated based on suitable quantiles of the estimated (average) genetically predicted effect across all explored DAGs
#'
#' @export
#'
#' @return The value returned is a list object \code{list(genPredEffect, genPredEffect_LL, genPredEffect_UL, group, BMA, CI)}:
#' \describe{
#'   \item{\code{genPredEffect}}{ Estimate of the (average) genetically predicted effect under intervention on a trait (target) on another one 
#'         (outcome) based on the DAGs explored by the \code{MrDAG} algorithm. A feedback loop genetically predicted effect of the same trait is not 
#'         allowed (NA). Note that the estimated genetically predicted effect corresponds to the marginal genetically predicted effect of a target on an outcome 
#'         *without* considering the effects of the other traits }
#'   \item{\code{genPredEffect_LL}}{ Lower limit of the \code{CI} for the (average) genetically predicted effect. It is calculated as the 
#'         (1-\code{CIboot})/2 % quantile of the (average) genetically predicted effect across all explored DAGs }
#'   \item{\code{genPredEffect_UL}}{ Upper limit of the \code{CI} for the (average) genetically predicted effect. It is calculated as the 
#'         1-\[(1-\code{CIboot})/2\] % quantile of the (average) genetically predicted effect across all explored DAGs }
#'   \item{\code{group}}{ If \code{MRconstraint} in the \code{MrDAG} algorithm contains the indices of the traits in \code{data} that 
#'         define the outcomes and the exposures, \code{group} coincides with this list }
#'   \item{\code{BMA}}{ Logical option }
#'   \item{\code{CI}}{ Level of the credible interval option } }
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # Estimation of the (average) genetically predicted effect under intervention on lifestyle
#' # and behavioural
#' # exposures, and measured on mental health phenotypes. 708 independent instrumental variables 
#' # (IVs) were selected to be associated at genome-wide significance with six lifestyle and 
#' # behavioural traits after pruning or clumping, which are considered exposures for the risk of 
#' # seven mental health outcome phenotypes
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
#' # Finally, the posterior (average) genetically predicted effect and credible interval (CI)
#' # of the intervention
#' # on ALC and measured on SCZ are estimated
#'
#' genPredEffect <- get_genPredEffect(output, 7, 8)   # Intervention on ALC(8) and measured on SCZ(7)
#'
#' genPredEffect <- get_genPredEffect(output, 7, c(4, 8))   # Intervention on ALC(8) and measured on 
#'                                                        # SCZ(7) after conditioning on BD(4)
#' }


########## get_genPredEffect ##########

get_genPredEffect <- function(output, response, target, BMA = TRUE, CI = 0.95)
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
  
  if (is.null(dim(output$graphs)))
  {
    graphs <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                              paste0("iter", seq(burnin + thin, niter, thin))))
    L <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                         paste0("iter", seq(burnin + thin, niter, thin))))
    D <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                         paste0("iter", seq(burnin + thin, niter, thin))))
    
    for (iter in 1 : niter_TMP)
    {
      graphs[, , iter] <- bd_decode(output$graphs[iter])
      L[, , iter] <- bd_decode(output$L[iter])
      D[, , iter] <- bd_decode(output$D[iter])
    }
  
  } else {
    graphs <- output$graphs
    L <- output$L
    D <- output$D
  }
  
  genPredEffect <- matrix(0, ncol = length(target), nrow = niter_TMP)
  colnames(genPredEffect) <- paste("Genetically predicted effect under intervention on", colnames_YX[target], "and measured on", colnames_YX[response], sep = " ")
  rownames(genPredEffect) <- paste0("iter", seq(burnin + thin, niter, thin))
  
  for (iter in 1 : niter_TMP)
  {
    genPredEffect[iter, ] <- genPredEffect_R(response, target, L[, , iter], D[, , iter])
  }
  
  if (BMA == TRUE)
  {
    genPredEffect_LL <- apply(genPredEffect, 2, stats::quantile, (1 - CI) / 2)
    genPredEffect_UL <- apply(genPredEffect, 2, stats::quantile, 1 - ((1 - CI) / 2))
    genPredEffect <- colMeans(genPredEffect)
    
    return(list(genPredEffect = genPredEffect, genPredEffect_LL = genPredEffect_LL, genPredEffect_UL = genPredEffect_UL, group = group, BMA = BMA, CI = CI))
  } else {
    return(list(genPredEffect = genPredEffect, group = group, BMA = BMA, CI = CI))
  }
  
}
