# 135 characters #####################################################################################################################
#' @title Plot of the bootstrap frequencies of edge inclusion
#' @description ggplot of the bootstrap frequencies of edge inclusion calculated as the average of the estimated Posterior 
#' Edge Inclusion Probabilities (PEIPs) obtained by running the \code{MrDAG} algorithm on each bootstrap summary-level dataset 
#' and the results aggregated by the average vote. Feedback loop genetically predicted effects of the same trait are not allowed and they are 
#' represented by black-and-white stripes. If a partition of the traits between outcomes and exposures is specified, the bootstrap 
#' frequencies of edge inclusion between outcomes and exposures are not considered and replaced by black-and-white stripes
#'
#' @param output Output produced by the \code{boot4MrDAG} algorithm
#' @param donotplot Logical value (default \code{FALSE}) indicating if the plot is run but not plotted
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when plotting the bootstrap frequencies of edge inclusion. If \code{ord = NULL}, the ordering follows \code{MRconstraint}
#' in the \code{MrDAG} output.
#' @param threshold Threshold (default \code{0}) on the bootstrap frequencies of edge inclusion. Frequencies below \code{threshold} 
#' are not plotted
#' @param labels Labels of the traits when plotting the bootstrap frequencies of edge inclusion. If \code{labels = NULL}, the labels 
#' are retrieved from the \code{MrDAG} output
#' @param digits Number of digits (default \code{2}) used in the plot for the bootstrap frequencies of edge inclusion 
#' @param fileType Type of plot file. Default is \code{pdf}. Other file types include \code{jpeg}, \code{tiff}, and \code{png}
#' @param fileName Name of the plot (default \code{"plot_bootEdgeFreq"})
#' @param filePath Directory where the plot is saved. If the directory does not exist, \code{plot_bootEdgeFreq} will create it by 
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
#' # Plotting bootstrap assessment of the robustness to weak instrumental variables (IVs) of the 
#' # estimated Posterior Edge Inclusion Probabilities (PEIPs) of mental health phenotypes that 
#' # might impact lifestyle and behavioural traits. 470 independent IVs were selected to be 
#' # associated at genome-wide significance with seven mental health phenotypes that are 
#' # considered exposures for the risk of six lifestyle and behavioural traits
#'
#' # After loading the dataset:
#'
#' data(MD2LBT_data)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list 
#' # object:
#'
#' MRconstraint <- NULL
#' MRconstraint$Y_idx <- 1 : 6    # Lifestyle and behavioural traits
#' MRconstraint$X_idx <- 7 : 13   # Mental health phenotypes
#'
#' # The boot4MrDAG algorithm is run to generate the bootstrap frequencies of edge inclusion 
#' # calculated as the average of the estimated PEIPs obtained by running the \code{MrDAG} 
#' # algorithm on 5 bootstrap summary-level datasets and the results aggregated by the average 
#' # vote. Each instance of the MrDAG algorithm is run to generate 2,000 posterior samples of 
#' # all unknowns, with an overall 5,000 posterior samples
#'
#' output <- boot4MrDAG(data = MD2LBT_data$data, gencov = MD2LBT_data$gencov,
#'                      niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                      MRconstraint = MRconstraint, fileName = NULL)
#'
#' # Bootstrap frequencies of edge inclusion are plotted with an order of appearance that is 
#' # changed within the outcomes and within the exposures compared to the order of the traits 
#' # in the boot4MrDAG output
#'
#' ord <- c(c(2, 4, 6, 1, 5, 3), c(12, 8, 7, 9, 10, 13, 11))
#' plot_bootEdgeFreq(output, ord = ord, 
#'                    fileName = "Example_plot_bootEdgeFreq")
#'
#' # and the figure is saved as a .pdf file in the directory
#'
#' print(getwd())
#' }


########## plot_bootEdgeFreq ##########

plot_bootEdgeFreq <- function(output, donotplot = FALSE, 
                              ord = NULL, threshold = 0, labels = NULL, digits = 2, 
                              fileType = "pdf", fileName = "plot_bootEdgeFreq", filePath = NULL, 
                              plotPar = list(cex.axis = NULL, cex.lab = NULL, geom_text.size = NULL))
{
  NotAllowed <- NULL
  value <- NULL
  variable <- NULL
  
  edgeFreq <- output$bootEdgeFreq
  output <- output$bootOutput[[1]]
  
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
  
  edgeFreq[abs(edgeFreq) < threshold] <- 0
  
  if (!is.null(labels))
  {
    rownames(edgeFreq) <- labels
    colnames(edgeFreq) <- labels
  }
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    edgeFreq <- edgeFreq[c(Y_ord, X_ord), c(Y_ord, X_ord)]
  } else {
    edgeFreq <- edgeFreq[Y_ord, Y_ord]
  }
  
  tt <- data.frame(round(edgeFreq, digits))
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
