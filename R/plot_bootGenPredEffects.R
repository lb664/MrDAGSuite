# 135 characters #####################################################################################################################
#' @title Plot of the bootstrap genetically predicted effects
#' @description ggplot of the bootstrap genetically predicted effects calculated as the average of the estimated genetically predicted effects obtained by running 
#' the \code{MrDAG} algorithm on each bootstrap summary-level dataset and the results aggregated by the average vote. Feedback loop 
#' genetically predicted effects of the same trait are not allowed and they are represented by black-and-white stripes
#'
#' @param output Output produced by the \code{boot4MrDAG} algorithm
#' @param donotplot Logical value (default \code{FALSE}) indicating if the plot is run but not plotted
#' @param CI If logical \code{CI == TRUE}, credible intervals (CIs) are plotted at the same level chosen in \code{boot4MrDAG}. 
#' By default, they are not plotted (\code{CI = FALSE}), unless otherwise specified
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when plotting the bootstrap genetically predicted effects. If \code{ord = NULL}, the ordering follows \code{MRconstraint} in the
#' \code{MrDAG} output.
#' @param threshold Threshold (default \code{0}) on the bootstrap genetically predicted effects. Genetically predicted effects in absolute value below 
#' \code{threshold} are not plotted
#' @param labels Labels of the traits when plotting the bootstrap genetically predicted effects. If \code{labels = NULL}, the labels 
#' are retrieved from the \code{MrDAG} output
#' @param digits Number of digits (default \code{2}) used in the plot for the bootstrap genetically predicted effects
#' @param fileType Type of plot file. Default is \code{pdf}. Other file types include \code{jpeg}, \code{tiff}, and \code{png}
#' @param fileName Name of the plot (default \code{"plot_bootGenPredEffects"})
#' @param filePath Directory where the plot is saved. If the directory does not exist, \code{plot_bootGenPredEffects} will 
#' create it by using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), 
#' the plot is saved in the current working directory
#' @param plotPar List of parameters that can be changed to modify the appearance of the plot, with elements
#' \code{cex.axis} (size of the axis labels), \code{cex.lab} (size of the axis titles) and \code{geom_text.size} (size
#' of the numbers in the tiles). Elements that are \code{NULL} (default) are set automatically from the number of traits
#' and the number of characters of the labels, so that the numbers fit in the tiles and the labels remain readable. A
#' number overrides the automatic choice. In addition, \code{limit} (default \code{NULL}, automatic) is the range of the
#' plotted genetically predicted effects
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example
#' # Plotting bootstrap assessment of the robustness to weak instrumental variables (IVs) of the 
#' # estimated genetically predicted effects of mental health phenotypes that might impact
#' # lifestyle and
#' # behavioural traits. 470 independent IVs were selected to be associated at genome-wide 
#' # significance with seven mental health phenotypes that are considered exposures for the 
#' # risk of six lifestyle and behavioural traits
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
#' # The boot4MrDAG algorithm is run to generate the bootstrap genetically predicted effects
#' # calculated as the
#' # average of the estimated genetically predicted effects obtained by running the
#' # \code{MrDAG} algorithm on 5
#' # bootstrap summary-level datasets and the results aggregated by the average vote. Each 
#' # instance of the MrDAG algorithm is run to generate 2,000 posterior samples of all 
#' # unknowns, with an overall 5,000 posterior samples
#'
#' output <- boot4MrDAG(data = MD2LBT_data$data, gencov = MD2LBT_data$gencov, CIboot = 0.9, 
#'                      niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                      MRconstraint = MRconstraint, fileName = NULL)
#'
#' # Finally, the bootstrap genetically predicted effects are plotted with an order of
#' # appearance that is
#' # changed within the outcomes and within the exposures compared to the order of the traits 
#' # in the boot4MrDAG output
#'
#' ord <- c(c(2, 4, 6, 1, 5, 3), c(12, 8, 7, 9, 10, 13, 11))
#' plot_bootGenPredEffects(output, ord = ord, 
#'                        fileName = "Example_plot_bootGenPredEffects")
#'
#' # and the figure is saved as a .pdf file in the directory
#'
#' print(getwd())
#'
#' # Alternatively, bootstrap genetically predicted effects are plotted along with the
#' # confidence intervals of
#' # the bootstrap genetically predicted effects selected in the boot4MrDAG algorithm
#'
#' plot_bootGenPredEffects(output, CI = TRUE, ord = ord, 
#'                        fileName = "Example_plot_bootGenPredEffects_w_CI")
#' }


########## plot_bootGenPredEffects ##########

plot_bootGenPredEffects <- function(output, donotplot = FALSE, CI = FALSE, ord = NULL, 
                                   threshold = 0, labels = NULL, digits = 2, 
                                   fileType = "pdf", fileName = "plot_bootGenPredEffects", filePath = NULL, 
                                   plotPar = list(limit = NULL, cex.axis = NULL, cex.lab = NULL, geom_text.size = NULL))
{
  NotAllowed <- NULL
  value <- NULL
  variable <- NULL
  
  genPredEffects <- output$bootGenPredEffects
  genPredEffects_LL <- output$bootGenPredEffects_LL
  genPredEffects_UL <- output$bootGenPredEffects_UL
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
    Y_intercept <- -m
    X_intercept <- -m
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
      X_ord <- ord
    }
  
  } else {
    
    if (!is.null(output$hyperPar$MRconstraint))
    {
      Y_ord <- output$hyperPar$MRconstraint$Y_idx
      X_ord <- output$hyperPar$MRconstraint$X_idx
    } else {
      Y_ord <- seq(1, m, 1)
      X_ord <- seq(1, m, 1)
    }
  
  }
  
  Y_names <- colnames_YX[Y_ord]
  X_names <- colnames_YX[X_ord]
  
  genPredEffects[abs(genPredEffects) < threshold] <- 0
  genPredEffects_LL[abs(genPredEffects) < threshold] <- 0
  genPredEffects_UL[abs(genPredEffects) < threshold] <- 0
  
  if (!is.null(labels))
  {
    rownames(genPredEffects) <- rownames(genPredEffects_LL) <- rownames(genPredEffects_UL) <- labels
    colnames(genPredEffects) <- colnames(genPredEffects_LL) <- colnames(genPredEffects_UL) <- labels
  }
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    genPredEffects <- genPredEffects[X_ord - q, Y_ord]
    genPredEffects_LL <- genPredEffects_LL[X_ord - q, Y_ord]
    genPredEffects_UL <- genPredEffects_UL[X_ord - q, Y_ord]
  } else {
    genPredEffects <- genPredEffects[Y_ord, Y_ord]
    genPredEffects_LL <- genPredEffects_LL[Y_ord, Y_ord]
    genPredEffects_UL <- genPredEffects_UL[Y_ord, Y_ord]
  }

  if (is.null(dim(genPredEffects)))
  {
    genPredEffects <- as.matrix(genPredEffects)
    genPredEffects_LL <- as.matrix(genPredEffects_LL)
    genPredEffects_UL <- as.matrix(genPredEffects_UL)
    colnames(genPredEffects)[1] <- Y_names
    colnames(genPredEffects_LL)[1] <- Y_names
    colnames(genPredEffects_UL)[1] <- Y_names
  }
  
  tt <- data.frame(round(genPredEffects, digits))
  tt$row <- X_names
  ttmelt <- reshape2::melt(tt, measure.vars = Y_names)
  ttmelt$row <- factor(ttmelt$row, levels = unique(rev(c(Y_names, X_names))))
  ttmelt$variable <- factor(ttmelt$variable, levels = unique(c(Y_names, X_names)))
  
  tt_LL <- data.frame(round(genPredEffects_LL, digits))
  tt_LL$row <- X_names
  ttmelt_LL <- reshape2::melt(tt_LL, measure.vars = Y_names)
  ttmelt_LL$row <- factor(ttmelt_LL$row, levels = unique(rev(c(Y_names, X_names))))
  ttmelt_LL$variable <- factor(ttmelt_LL$variable, levels = unique(c(Y_names, X_names)))
  
  tt_UL <- data.frame(round(genPredEffects_UL, digits))
  tt_UL$row <- X_names
  ttmelt_UL <- reshape2::melt(tt_UL, measure.vars = Y_names)
  ttmelt_UL$row <- factor(ttmelt_UL$row, levels = unique(rev(c(Y_names, X_names))))
  ttmelt_UL$variable <- factor(ttmelt_UL$variable, levels = unique(c(Y_names, X_names)))
  
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
  
  plotParDefault <- list(limit = NULL, cex.axis = NULL, cex.lab = NULL, geom_text.size = NULL)
  
  for (nm in intersect(names(plotPar), names(plotParDefault))) plotParDefault[nm] <- list(plotPar[[nm]])
  
  if (!is.null(output$hyperPar$MRconstraint))
  {
    xlab <- "Outcomes"
    ylab <- "Intervention on Exposures"
  } else {
    xlab <- "Traits"
    ylab <- "Intervention on Traits"
  }
  
  if (is.null(plotParDefault$limit))
  {
    ## Rounded up to the next 0.5 (not just any 2-sig-fig value) so the color scale always
    ## lands on a "clean" bound -- e.g. an actual max of 1.193 or 1.273 both round up to
    ## 1.5, not to two different values like 1.2/1.3 that would happen to differ only
    ## because of where the raw estimates landed.
    plotParDefault$limit <- max(1, ceiling(max(abs(genPredEffects), na.rm = TRUE) / 0.5) * 0.5)
  }
  
  ttmelt <- cbind(ttmelt, rep("No", nrow(ttmelt)))
  ttmelt[which(is.na(ttmelt$value)), 4] <- "Yes"
  colnames(ttmelt)[4] <- "NotAllowed"
  ttmelt[ttmelt[, 3] == 0 & !is.na(ttmelt[, 3]), 3] <- NA
  ttmelt_LL[ttmelt_LL[, 3] == 0 & !is.na(ttmelt_LL[, 3]), 3] <- 0
  ttmelt_UL[ttmelt_UL[, 3] == 0 & !is.na(ttmelt_UL[, 3]), 3] <- 0
  ttmelt <- cbind(ttmelt, rep(NA, nrow(ttmelt)))
  colnames(ttmelt)[5] <- "value_TMP"
  ttmelt <- cbind(ttmelt, rep(NA, nrow(ttmelt)))
  colnames(ttmelt)[6] <- "value_LL_UL"
  
  if (CI != FALSE)
  {
    
    for (i in 1 : nrow(ttmelt))
    {
      
      if (!is.na(ttmelt[i, 3]))
      {
        ttmelt[i, 5] <- paste0(ttmelt[i, 3], "\n")
        ttmelt[i, 6] <- paste0("\n(", 
                               ttmelt_LL[i, 3], 
                               ";", 
                               ttmelt_UL[i, 3], 
                               ")")
      }
    
    }
  
  } else {
    ttmelt[, 5] <- ttmelt[, 3]
  }
  
  autoTextSize <- is.null(plotParDefault$geom_text.size)
  plotParDefault <- .autoPlotPar(plotParDefault, nRows = length(unique(ttmelt$row)), nCols = length(unique(ttmelt$variable)), labelChars = if (is.null(CI)) digits + 3 else 0.7 * 8)

  fig <- ggplot2::ggplot(data = ttmelt, ggplot2::aes(x = variable, y = row, pattern = NotAllowed, fill = value)) + 
         ggplot2::xlab(xlab) + 
         ggplot2::ylab(ylab) + 
         ggpattern::geom_tile_pattern(pattern_color = NA, 
                                      pattern_fill = "black", 
                                      pattern_angle = 45, 
                                      pattern_density = 0.4, 
                                      pattern_spacing = 0.02, 
                                      pattern_key_scale_factor = 1) + 
         ggpattern::scale_pattern_manual(values = c(Yes = "stripe", No = "none")) + 
         ggplot2::scale_fill_gradient2(low = "#009E73", mid = "white", high = "#D55E00", na.value = "white", limits = c(-plotParDefault$limit, plotParDefault$limit), breaks = as.numeric(seq(-plotParDefault$limit, plotParDefault$limit, length.out = 5))) + 
         ggplot2::geom_text(ggplot2::aes(label = value_TMP), color = 'black', size = plotParDefault$geom_text.size) + 
         ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis)), 
                        axis.title = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.lab)), 
                        axis.text.y = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), color = 'black'), 
                        axis.text.x = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), angle = 45, hjust = 1, color = 'black'), 
                        legend.title = ggplot2::element_blank(), 
                        legend.text = ggplot2::element_text(size = ggplot2::rel(1)), 
                        legend.key.size = ggplot2::unit(1, 'cm')) + 
         ggplot2::guides(pattern = "none")
  
  if (CI != FALSE)
  {
    ## Confidence-interval text: 0.7 of the size of the estimates. With automatic sizes, intervals with both limits
    ## negative (two minus signs, the longest ones) are drawn slightly smaller (factor 0.85) so that they stay within the tile
    ciSize <- plotParDefault$geom_text.size * 0.7
    ciData <- ttmelt
    ciData$ciSize <- ciSize
    if (autoTextSize)
    {
      nMinus <- lengths(regmatches(as.character(ttmelt$value_LL_UL), gregexpr("[-\u2212]", as.character(ttmelt$value_LL_UL))))
      ciData$ciSize <- ifelse(nMinus >= 2, 0.85 * ciSize, ciSize)
    }
    fig <- fig + ggplot2::geom_text(data = ciData, ggplot2::aes(label = value_LL_UL, size = ciSize), color = 'black') + 
                 ggplot2::scale_size_identity(guide = "none")
  }
  
  if (donotplot == TRUE)
  {
    return(fig)
  } else {
    suppressWarnings(print(fig))
    invisible(invisible(grDevices::dev.off()))
  }
  
}
