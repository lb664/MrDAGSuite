# 135 characters #####################################################################################################################
#' @title Partially DAG (PDAG) of the estimated genetically predicted effects under intervention on the exposures
#' @description qgraph of the Partially DAG (PDAG) of the estimated (average) genetically predicted effects under intervention on the exposures 
#' (targets) and measured on the outcomes (responses) based on the Directed Acyclic Graphs (DAGs) explored by the \code{MrDAG} 
#' algorithm. Feedback loop genetically predicted effects of the same trait are not allowed
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param donotplot Logical value (default \code{FALSE}) indicating if the qgraph is run but not plotted. Useful for saving the 
#' output (i.e., layout) without plotting
#' @param progress Logical value set to \code{FALSE} by default to print on the screen the progress of the genetically predicted effects estimation
#' @param layout Layout of the graph (default \code{"spring"}). For details, see \insertCite{Epskamp2012;textual}{MrDAGSuite} and 
#' \insertCite{Epskamp2023;textual}{MrDAGSuite}
#' @param threshold Threshold (default \code{0}) on the estimated genetically predicted effects. Genetically predicted effects in absolute value below 
#' \code{threshold} are not plotted
#' @param labels Labels of the traits when plotting the PDAG of the genetically predicted effects. If \code{labels = NULL}, the labels are 
#' retrieved from the \code{MrDAG} output
#' @param legend Logical value indicating if a legend should be plotted (default \code{FALSE}) when the traits are partitioned 
#' into outcomes and exposures
#' @param palette Palette used for colouring nodes when the traits are partitioned into outcomes and exposures. Default is 
#' \code{"colorblind"}. For alternative \code{palette}s, see \insertCite{Epskamp2012;textual}{MrDAGSuite} and 
#' \insertCite{Epskamp2023;textual}{MrDAGSuite}. If the traits are not partitioned, the colour for each node is white
#' @param fileType Type of plot file. Default is \code{pdf}. Other file types include \code{jpeg}, \code{tiff}, and \code{png}
#' @param fileName Name of the plot (default \code{"qgraph_genPredEffects"})
#' @param filePath Directory where the plot is saved. If the directory does not exist, \code{qgraph_genPredEffects} will create 
#' it by using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), 
#' the plot is saved in the current working directory
#' @param qgraphPar List of parameters that can be changed to modify the appearance of the graph. Default values are \code{list(
#' label.cex = 1.25, node.width.height = 0.9, edge.label.position = 0.675, edge.label.cex = 1.125, edge.label.color = "dimgrey", 
#' asize = 4, legend.cex = 1)}. For details, see \insertCite{Epskamp2012;textual}{MrDAGSuite} and 
#' \insertCite{Epskamp2023;textual}{MrDAGSuite}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example
#' # Partially DAG (PDAG) of the estimated genetically predicted effects under intervention on
#' # lifestyle and
#' # behavioural exposures, and measured on mental health phenotypes. 708 independent instrumental 
#' # variables (IVs) were selected to be associated at genome-wide significance with six lifestyle 
#' # and behavioural traits after pruning or clumping, which are considered exposures for the risk 
#' # of seven mental health outcome phenotypes
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
#' # Finally, the PDAG of the estimated genetically predicted effects is plotted with default values
#'
#' qgraph_genPredEffects(output, threshold = 0.01, 
#'                      fileName = "Example_qgraph_genPredEffects")
#'
#' # and the figure is saved as a .pdf file in the directory
#'
#' print(getwd())
#'
#' # An alternative "circle" representation can be selected
#'
#' qgraph_genPredEffects(output, progress = TRUE, 
#'                      layout = "circle", threshold = 0.01, legend = TRUE, palette = "pastel", 
#'                      fileName = "Example_alternative_qgraph_genPredEffects", 
#'                      qgraphPar = list(edge.label.position = 0.5, edge.label.cex = 1))
#' }

########## qgraph_genPredEffects ##########

qgraph_genPredEffects <- function(output, donotplot = FALSE, progress = FALSE, 
                                 layout = NULL, threshold = 0, labels = NULL, legend = FALSE, palette = NULL, 
                                 fileType = "pdf", fileName = "qgraph_genPredEffects", filePath = NULL, 
                                 qgraphPar = list(label.cex = 1.25, node.width.height = 0.9, edge.label.position = 0.675, edge.label.cex = 1.125, edge.label.color = "dimgrey", asize = 4, legend.cex = 1))
{
  
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
    
    if (is.null(palette))
    {
      palette <- "colorblind"
    }
  
  } else {
    group <- NULL
    legend <- FALSE
  }
  
  if (is.null(layout))
  {
    layout <- "spring"
  }
  
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
  
  qgraphParDefault <- list(label.cex = 1.25, node.width.height = 0.9, edge.label.position = 0.675, edge.label.cex = 1.125, edge.label.color = "dimgrey", asize = 4, legend.cex = 1)
  
  if (length(match(names(qgraphPar), names(qgraphParDefault))) > 0)
  {
    
    for (p in 1 : length(names(qgraphPar)))
    {
      qgraphParDefault[[match(names(qgraphPar)[p], names(qgraphParDefault))]] <- qgraphPar[[p]]
    }
  
  }
  
  genPredEffects <- get_genPredEffects(output, progress = progress)$genPredEffects
  genPredEffects[is.na(genPredEffects)] <- 0
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    q <- length(output$hyperPar$MRconstraint$Y_idx)
    p <- length(output$hyperPar$MRconstraint$X_idx)
    colnames_YX <- colnames(output$hyperPar$U)
    genPredEffects <- rbind(matrix(0, nrow = q, ncol = m), cbind(genPredEffects, matrix(0, nrow = p, ncol = p)))
    colnames(genPredEffects) <- colnames_YX
    rownames(genPredEffects) <- colnames_YX
  }
  
  if (is.null(labels))
  {
    suppressMessages(
    output <- qgraph::qgraph(genPredEffects, 
                             layout = layout, 
                             DoNotPlot = donotplot, 
                             threshold = threshold, 
                             groups = group, 
                             palette = palette, 
                             posCol = "#D55E00", 
                             negCol = "#009E73", 
                             label.cex = qgraphParDefault$label.cex, 
                             node.height = qgraphParDefault$node.width.height, 
                             node.width = qgraphParDefault$node.width.height,
                             edge.labels = TRUE, 
                             edge.label.position = qgraphParDefault$edge.label.position, 
                             edge.label.color = qgraphParDefault$edge.label.color, 
                             edge.label.cex = qgraphParDefault$edge.label.cex, 
                             asize = qgraphParDefault$asize, 
                             legend = legend, 
                             legend.cex = qgraphParDefault$legend.cex, 
                             filetype = fileType, 
                             filename = paste0(filePath, fileName)))
  } else {
    suppressMessages(
    output <- qgraph::qgraph(genPredEffects, 
                             layout = layout, 
                             DoNotPlot = donotplot, 
                             threshold = threshold, 
                             groups = group, 
                             palette = palette, 
                             posCol = "#D55E00", 
                             negCol = "#009E73", 
                             labels = labels, 
                             label.cex = qgraphParDefault$label.cex, 
                             node.height = qgraphParDefault$node.width.height, 
                             node.width = qgraphParDefault$node.width.height,
                             edge.labels = TRUE, 
                             edge.label.position = qgraphParDefault$edge.label.position, 
                             edge.label.color = qgraphParDefault$edge.label.color, 
                             edge.label.cex = qgraphParDefault$edge.label.cex, 
                             asize = qgraphParDefault$asize, 
                             legend = legend, 
                             legend.cex = qgraphParDefault$legend.cex, 
                             filetype = fileType, 
                             filename = paste0(filePath, fileName)))
  }
  
  if (donotplot == TRUE)
  {
    return(output)
  }

}
