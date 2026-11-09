# 135 characters #####################################################################################################################
#' @title Maximum a Posteriori model
#' @description List of unique visited models by the \code{MrDAG} algorithm ordered by the frequency of visits or the log-marginal 
#' likelihood (\insertCite{Castelletti2018;textual}{MrDAGSuite}). The Maximum a Posteriori (MAP) model is therefore identified as 
#' the top-ordered model. The models can be either Directed Acyclic Graphs (DAGs), Completed Partially DAGs (CPDAGs), or Essential 
#' Graphs (EGs)
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when printing the Posterior Edge Inclusion Probabilities (PEIPs). If \code{ord = NULL}, the ordering follows
#' \code{MRconstraint} in the \code{MrDAG} output.
#' @param model Either (default) "DAG" or "CPDAG" / "EG"
#' @param sortby Unique models sorted in descending order of importance according to (default) \code{"freq"} or \code{"logMargLik"}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(map_XXX, freq_map_XXX, logMargLik_map_XXX, ord)}, where XXX can be either 
#' DAG or CPDAG/EG depending on \code{model}:
#' \describe{
#'   \item{\code{map_XXX}}{ Adjacency matrices of the list of unique visited models ordered by the frequency of visits (in 
#'         descending order) }
#'   \item{\code{freq_map_XXX}}{ Frequency of visits (in descending order) }
#'   \item{\code{logMargLik_map_XXX}}{ Log-marginal likelihood for each unique visited model ordered by the frequency of visits (in 
#'         descending order) }
#'   \item{\code{ord}}{ Indices of the traits in \code{data} that specify the outcomes and the exposures } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # Estimation of the Maximum a Posteriori model of lifestyle and behavioural traits, and mental
#' # health phenotypes. 708 independent instrumental variables (IVs) were selected to be
#' # associated at genome-wide significance with six lifestyle and behavioural traits after
#' # pruning or clumping, which are considered exposures for the risk of seven mental health
#' # outcome phenotypes
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
#' # Finally, the Maximum a Posteriori CPDAG model is estimated with an order of appearance that
#' # is changed within the outcomes and within the exposures compared to the order of the traits
#' # in the MrDAG output
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' MAP <- get_MAP(output, ord = ord, model = "CPDAG", sortby = "logMargLik")
#' }


########## get_MAP ##########

get_MAP <- function(output, ord = NULL, model = "DAG", sortby = "freq")
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
    dag_code <- as.vector(sapply(output$graphs, gsub, pattern = ";", replacement = ""))
    m <- sqrt(nchar(dag_code[[1]]))
  } else {
    graphs <- output$graphs
    dag_code <- apply(graphs, 3, bd_encode, separator = "")
    m <- dim(graphs)[1]
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
  
  if (model == "DAG")
  {
    logMargLik <- output$logMargLik[seq(burnin + thin, burnin + niter, thin)]
    uniq_dag <- unique(dag_code)
    map_dagcode <- uniq_dag[which.max(tabulate(match(dag_code, uniq_dag)))]
    map_dagcode <- bd_decode(map_dagcode, separator = "")
    s_tab <- sort(tabulate(match(dag_code, uniq_dag)), decreasing = TRUE, index.return = TRUE)
    map_dag <- NULL
    freq_map_dag <- NULL
    logMargLik_map_dag <- NULL
    
    for (d in 1 : length(uniq_dag))
    {
      map_dagcode <- uniq_dag[s_tab$ix[d]]
      logMargLik_map_dag <- c(logMargLik_map_dag, logMargLik[which(match(dag_code, map_dagcode) == 1)][1])
      map_dagcode <- bd_decode(map_dagcode, separator = "")
      map_dag_TMP <- matrix(map_dagcode, ncol = m)
      diag(map_dag_TMP) <- NA
      map_dag_TMP[1 : q, q + (1 : p)] <- NA
      colnames(map_dag_TMP) <- colnames_YX
      rownames(map_dag_TMP) <- colnames_YX
      
      if (!is.null(output$hyperPar$MRconstraint))
      {
        map_dag_TMP <- map_dag_TMP[c(Y_ord, X_ord), c(Y_ord, X_ord)]
      } else {
        map_dag_TMP <- map_dag_TMP[Y_ord, Y_ord]
      }
      
      map_dag[[d]] <- map_dag_TMP
      freq_map_dag <- c(freq_map_dag, s_tab$x[d] / sum(s_tab$x))
    }
    
    if (tolower(sortby) == "logmarglik")
    {
      s_logMargLik_map_dag <- sort(logMargLik_map_dag, decreasing = TRUE, index.return = TRUE)
      map_dag <- map_dag[s_logMargLik_map_dag$ix]
      freq_map_dag <- freq_map_dag[s_logMargLik_map_dag$ix]
      logMargLik_map_dag <- logMargLik_map_dag[s_logMargLik_map_dag$ix]
    }
  
  } else {
    cpdag_code <- dag_code
    
    for (d in 1 : length(dag_code))
    {
      map_dagcode <- dag_code[d]
      map_dagcode <- bd_decode(map_dagcode, separator = "")
      cpdag_code[d] <- paste(c(pcalg::dag2cpdag(matrix(map_dagcode, ncol = m)) * 1), collapse = "")
    }
    
    logMargLik <- output$logMargLik[seq(burnin + thin, burnin + niter, thin)]
    uniq_cpdag <- unique(cpdag_code)
    map_cpdagcode <- uniq_cpdag[which.max(tabulate(match(cpdag_code, uniq_cpdag)))]
    map_cpdagcode <- bd_decode(map_cpdagcode, separator = "")
    s_tab <- sort(tabulate(match(cpdag_code, uniq_cpdag)), decreasing = TRUE, index.return = TRUE)
    map_cpdag <- NULL
    freq_map_cpdag <- NULL
    logMargLik_map_cpdag <- NULL
    
    for (d in 1 : length(uniq_cpdag))
    {
      map_cpdagcode <- uniq_cpdag[s_tab$ix[d]]
      logMargLik_map_cpdag <- c(logMargLik_map_cpdag, logMargLik[which(match(cpdag_code, map_cpdagcode) == 1)][1])
      map_cpdagcode <- bd_decode(map_cpdagcode, separator = "")
      map_cpdag_TMP <- matrix(map_cpdagcode, ncol = m)
      diag(map_cpdag_TMP) <- NA
      map_cpdag_TMP[1 : q, q + (1 : p)] <- NA
      colnames(map_cpdag_TMP) <- colnames_YX
      rownames(map_cpdag_TMP) <- colnames_YX
      
      if (!is.null(output$hyperPar$MRconstraint))
      {
        map_cpdag_TMP <- map_cpdag_TMP[c(Y_ord, X_ord), c(Y_ord, X_ord)]
      } else {
        map_cpdag_TMP <- map_cpdag_TMP[Y_ord, Y_ord]
      }
      
      map_cpdag[[d]] <- map_cpdag_TMP
      freq_map_cpdag <- c(freq_map_cpdag, s_tab$x[d] / sum(s_tab$x))
    }
    
    if (tolower(sortby) == "logmarglik")
    {
      s_logMargLik_map_cpdag <- sort(logMargLik_map_cpdag, decreasing = TRUE, index.return = TRUE)
      map_cpdag <- map_cpdag[s_logMargLik_map_cpdag$ix]
      freq_map_cpdag <- freq_map_cpdag[s_logMargLik_map_cpdag$ix]
      logMargLik_map_cpdag <- logMargLik_map_cpdag[s_logMargLik_map_cpdag$ix]
    }
  
  }
  
  if (model == "DAG")
  {
    return(list(map_DAG = map_dag, freq_map_DAG = freq_map_dag, logMargLik_map_DAG = logMargLik_map_dag, ord = ord))
  } else {
    return(list(map_CPDAG = map_cpdag, freq_map_CPDAG = freq_map_cpdag, logMargLik_map_CPDAG = logMargLik_map_cpdag, ord = ord))
  }
 
}
