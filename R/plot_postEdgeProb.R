# 135 characters #####################################################################################################################
#' @title Plot of the estimated Posterior Edge Inclusion Probabilities
#' @description ggplot of the estimated Posterior Edge Inclusion Probabilities (PEIPs) based on the Directed Acyclic Graphs (DAGs) 
#' explored by the \code{MrDAG} algorithm. Feedback loop genetically predicted effects of the same trait are not allowed and they are represented by 
#' black-and-white stripes. If a partition of the traits between outcomes and exposures is specified, the PEIPs between outcomes and 
#' exposures are not considered and replaced by black-and-white stripes
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param donotplot Logical value (default \code{FALSE}) indicating if the plot is run but not plotted
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when plotting the estimated PEIPs. If \code{ord = NULL}, the ordering follows \code{MRconstraint} in the \code{MrDAG} output.
#' @param threshold Threshold (default \code{0}) on the estimated PEIPs. PEIPs below \code{threshold} are not 
#' plotted
#' @param labels Labels of the traits when plotting the estimated PEIPs. If \code{labels = NULL}, the labels are retrieved from 
#' the \code{MrDAG} output
#' @param digits Number of digits (default \code{2}) in the plot for the estimated PEIPs
#' @param fileType Type of plot file. Default is \code{pdf}. Other file types include \code{jpeg}, \code{tiff}, and \code{png}
#' @param fileName Name of the plot (default \code{"plot_postEdgeProb"})
#' @param filePath Directory where the plot is saved. If the directory does not exist, \code{plot_postEdgeProb} will create it by 
#' using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), the plot is 
#' saved in the current working directory
#' @param plotPar List of parameters that can be changed to modify the appearance of the plot, with elements
#' \code{cex.axis} (size of the axis labels), \code{cex.lab} (size of the axis titles) and \code{geom_text.size} (size
#' of the numbers in the tiles). Elements that are \code{NULL} (default) are set automatically from the number of traits
#' and the number of characters of the labels, so that the numbers fit in the tiles and the labels remain readable. A
#' number overrides the automatic choice
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example
#' # Plotting the Posterior Edge Inclusion Probabilities (PEIPs) of intervention on lifestyle 
#' # and behavioural exposures, and measured on mental health phenotypes. 708 independent 
#' # instrumental variables (IVs) were selected to be associated at genome-wide significance 
#' # with six lifestyle and behavioural traits after pruning or clumping, which are considered 
#' # exposures for the risk of seven mental health outcome phenotypes
#'
#' # After loading the dataset:
#'
#' data(LBT2MD_data)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list 
#' # object:
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
#' # Finally, PEIPs are plotted
#'
#' plot_postEdgeProb(output, 
#'               fileName = "Example_plot_postEdgeProb")
#'
#' # and the figure is saved as a .pdf file in the directory
#'
#' print(getwd())
#'
#' # Alternatively, the estimated PEIPs are displayed with an order of appearance that is 
#' # changed within the outcomes and within the exposures compared to the order of the traits 
#' # in the MrDAG output
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' plot_postEdgeProb(output, ord = ord, 
#'               fileName = "Example_alternative_plot_postEdgeProb")
#' }


########## plot_postEdgeProb ##########

plot_postEdgeProb <- function(output, donotplot = FALSE, ord = NULL, threshold = 0, labels = NULL, digits = 2, 
                          fileType = "pdf", fileName = "plot_postEdgeProb", filePath = NULL, 
                          plotPar = list(cex.axis = NULL, cex.lab = NULL, geom_text.size = NULL))
{
  NotAllowed <- NULL
  value <- NULL
  variable <- NULL
  
  q <- length(output$hyperPar$MRconstraint$Y_idx)
  p <- length(output$hyperPar$MRconstraint$X_idx)
  
  if (!is.null(labels))
  {
    colnames_YX <- labels
  } else {
    colnames_YX <- colnames(output$hyperPar$U)
  }
  
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
    Y_intercept <- p
    X_intercept <- q
  } else {
    group <- NULL
    Y_intercept <- -m
    X_intercept <- -m
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
  
  Y_names <- colnames_YX[Y_ord]
  X_names <- colnames_YX[X_ord]
  
  postEdgeProb <- get_postEdgeProb(output)$postEdgeProb
  postEdgeProb[abs(postEdgeProb) < threshold] <- 0
  
  if (!is.null(labels))
  {
    rownames(postEdgeProb) <- labels
    colnames(postEdgeProb) <- labels
  }
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    postEdgeProb <- postEdgeProb[c(Y_ord, X_ord), c(Y_ord, X_ord)]
  } else {
    postEdgeProb <- postEdgeProb[Y_ord, Y_ord]
  }
  
  tt <- data.frame(round(postEdgeProb, digits))
  tt$row <- c(Y_names, X_names)
  ttmelt <- reshape2::melt(tt, measure.vars = c(Y_names, X_names))
  ttmelt$row <- factor(ttmelt$row, levels = unique(rev(c(Y_names, X_names))))
  ttmelt$variable <- factor(ttmelt$variable, levels = unique(c(Y_names, X_names)))
  
  saveas <- paste(fileName, fileType, sep = ".")
  
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
  
  if (donotplot == FALSE)
  {
  
    if (fileType == "pdf")
    {
      grDevices::pdf(file = paste0(filePath, saveas))
    } else if (fileType == "jpeg") {
      grDevices::jpeg(file = paste0(filePath, saveas))
    } else if (fileType == "tiff") {
      grDevices::tiff(file = paste0(filePath, saveas))
    } else if (fileType == "png") {
      grDevices::png(file = paste0(filePath, saveas))
    } else {
      stop("\nfileType ", paste0(".", fileType), " not currently supported\n")
    }
  
  }
  
  plotParDefault <- list(cex.axis = NULL, cex.lab = NULL, geom_text.size = NULL)
  
  for (nm in intersect(names(plotPar), names(plotParDefault))) plotParDefault[nm] <- list(plotPar[[nm]])
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    xlab <- "Outcomes and Exposures"
    ylab <- "Outcomes and Exposures"
  } else {
    xlab <- "Traits"
    ylab <- "Traits"
  }
  
  ttmelt <- cbind(ttmelt, rep("No", nrow(ttmelt)))
  ttmelt[which(is.na(ttmelt$value)), 4] <- "Yes"
  colnames(ttmelt)[4] <- "NotAllowed"
  ttmelt[ ttmelt[, 3] == 0 & !is.na(ttmelt[, 3]), 3] <- NA
  ttmelt <- cbind(ttmelt, rep("white", nrow(ttmelt)))
  ttmelt[which(ttmelt$value <= 0.25), 5] <- "black"
  colnames(ttmelt)[5] <- "Color"
  
  if (length(unique(ttmelt$Color)) == 1)
  {
    values = c("white", "white")
  } else {
    values = c("black", "white")
  }
  
  plotParDefault <- .autoPlotPar(plotParDefault, nRows = length(unique(ttmelt$row)), nCols = length(unique(ttmelt$variable)), labelChars = (digits + 2))

  fig <- ggplot2::ggplot(data = ttmelt, ggplot2::aes(x = variable, y = row, pattern = NotAllowed, fill = value)) + 
         ggplot2::xlab(xlab) + 
         ggplot2::ylab(ylab) + 
         ggpattern::geom_tile_pattern(pattern_color = NA, 
                                      pattern_fill = "black", 
                                      pattern_angle = 45, 
                                      pattern_density = 0.4, 
                                      pattern_spacing = 0.02, 
                                      pattern_key_scale_factor = 1) + 
         ggplot2::geom_hline(yintercept = Y_intercept + 0.5, linetype = "dashed", color = "black", linewidth = 1) + 
         ggplot2::geom_vline(xintercept = X_intercept + 0.5, linetype = "dashed", color = "black", linewidth = 1) + 
         ggpattern::scale_pattern_manual(values = c(Yes = "stripe", No = "none")) + 
         ggplot2::scale_fill_gradient(low = "white", high = "black", na.value = "white", limits = c(0, 1), breaks = seq(0, 1, length.out = 5)) + 
         ggplot2::geom_text(ggplot2::aes(label = value, color = Color), size = plotParDefault$geom_text.size, show.legend = FALSE) + 
         ggplot2::scale_colour_manual(values = values) + 
         ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis)), 
                        axis.title = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.lab)), 
                        axis.text.y = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), color = 'black'), 
                        axis.text.x = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), angle = 45, hjust = 1, color = 'black'), 
                        legend.title = ggplot2::element_blank(), 
                        legend.text = ggplot2::element_text(size = ggplot2::rel(1)), 
                        legend.key.size = ggplot2::unit(1, 'cm')) + 
         ggplot2::guides(pattern = "none")
  
  if (donotplot == TRUE)
  {
    return(fig)
  } else {
    suppressWarnings(print(fig))
    invisible(invisible(grDevices::dev.off()))
  }

}
