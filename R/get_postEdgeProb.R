# 135 characters #####################################################################################################################
#' @title Estimation of the Posterior Edge Inclusion Probabilities
#' @description Estimation of the Posterior Edge Inclusion Probabilities (PEIPs) based on the Directed Acyclic Graphs (DAGs) 
#' explored by the \code{MrDAG} algorithm
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when printing the estimated PEIPs. If (default) \code{ord = NULL}, the ordering follows \code{MRconstraint} in the
#' \code{MrDAG} output. Note that despite \code{ord} ordering, the outcomes always precede the exposures
#'
#' @export
#'
#' @return The value returned is a list object \code{list(postEdgeProb, group, ord)}:
#' \describe{
#'   \item{\code{postEdgeProb}}{ Estimate of the PEIPs between each pair of traits. Feedback loop of the same trait is not allowed (NA). 
#'          If a partition of the traits between outcomes and exposures is specified, the estimated PEIPs between outcomes and 
#'          exposures are not considered (NA) }
#'   \item{\code{group}}{ If \code{MRconstraint} in the \code{MrDAG} algorithm contains the indices of the traits in \code{data} that 
#'         define the outcomes and the exposures, \code{group} coincides with this list }
#'   \item{\code{ord}}{ Indices of the traits in \code{data} that specify the outcomes and the exposures. It might differ from 
#'         \code{group} if \code{ord} has been specified and it is different from \code{MRconstraint} } }
#'
#' @examples
#' \dontrun{
#' # Example
#' # Estimation of the Posterior Edge Inclusion Probabilities (PEIPs) of lifestyle and 
#' # behavioural traits, and mental health phenotypes. 708 independent instrumental variables 
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
#'                 MRconstraint = MRconstraint, fileName = NULL)
#'
#' # Finally, PEIPs are estimated and presented with an order of appearance that is changed within
#' # the outcomes and within the exposures compared to the order of the traits in the MrDAG output
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' postEdgeProb <- get_postEdgeProb(output, ord = ord)
#' }


########## get_postEdgeProb ##########

get_postEdgeProb <- function(output, ord = NULL)
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
  
  if (is.null(dim(output$graphs)))
  {
    graphs <- array(dim = c(m, m, niter_TMP), dimnames = list(colnames_YX, colnames_YX, 
                                                              paste0("iter", seq(burnin + thin, niter, thin))))
    
    for (iter in 1 : niter_TMP)
    {
      graphs[, , iter] <- bd_decode(output$graphs[iter])
    }
  
  } else {
    graphs <- output$graphs
  }
  
  postEdgeProb <- apply(graphs, c(1, 2), mean)
  diag(postEdgeProb) <- NA
  postEdgeProb[1 : q, q + (1 : p)] <- NA
  colnames(postEdgeProb) <- colnames_YX
  rownames(postEdgeProb) <- colnames_YX
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    postEdgeProb <- postEdgeProb[c(Y_ord, X_ord), c(Y_ord, X_ord)]
  } else {
    postEdgeProb <- postEdgeProb[Y_ord, Y_ord]
  }
  
  return(list(postEdgeProb = postEdgeProb, group = group, ord = ord))
}
