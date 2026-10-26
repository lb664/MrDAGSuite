# 135 characters #####################################################################################################################
#' @title Plot of the estimated genetically predicted effects under intervention on the exposures
#' @description ggplot of the estimated genetically predicted effects under intervention on the exposures (targets) and measured on the responses 
#' (outcomes) based on the Directed Acyclic Graphs (DAGs) explored by the \code{MrDAG} algorithm. If a partition of the traits 
#' between outcomes and exposures is not specified, feedback loop genetically predicted effects of the same trait are not allowed and they are 
#' represented by black-and-white stripes
#
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param donotplot Logical value (default \code{FALSE}) indicating if the plot is run but not plotted
#' @param CI Level of the credible intervals (CIs) of the genetically predicted effects. Default is \code{NULL}, and no CIs are plotted
#' @param ord Indices of the traits in \code{data} that define the outcomes and exposures. They determine the order in which traits
#' appear when plotting the estimated genetically predicted effects. If \code{ord = NULL}, the ordering follows \code{MRconstraint} in the
#' \code{MrDAG} output.
#' @param progress Logical value set to \code{FALSE} by default to print on the screen the progress of the genetically predicted effects estimation
#' @param threshold Threshold (default \code{0}) on the estimated genetically predicted effects. Genetically predicted effects in absolute value below 
#' \code{threshold} are not plotted
#' @param labels Labels of the traits when plotting the estimated genetically predicted effects. If \code{labels = NULL}, the labels are retrieved 
#' from the \code{MrDAG} output
#' @param digits Number of digits (default \code{2}) used in the plot for the estimated genetically predicted effects
#' @param fileType Type of plot file. Default is \code{pdf}. Other file types include \code{jpeg}, \code{tiff}, and \code{png}
#' @param fileName Name of the plot (default \code{"plot_genPredEffects"})
#' @param filePath Directory where the plot is saved. If the directory does not exist, \code{plot_genPredEffects} will create it by 
#' using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), the plot is 
#' saved in the current working directory
#' @param plotPar List of parameters that can be changed to modify the appearance of the plot, with elements
#' \code{cex.axis} (size of the axis labels), \code{cex.lab} (size of the axis titles) and \code{geom_text.size} (size
#' of the numbers in the tiles). Elements that are \code{NULL} (default) are set automatically from the number of traits
#' and the number of characters of the labels, so that the numbers fit in the tiles and the labels remain readable. A
#' number overrides the automatic choice. In addition, \code{limit} (default \code{NULL}, automatic) is the range of the
#' plotted genetically predicted effects
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example
#' # Plotting the estimated genetically predicted effects under intervention on lifestyle and
#' # behavioural
#' # exposures, and measured on mental health phenotypes. 708 independent instrumental variables 
#' # (IVs) were selected to be associated at genome-wide significance with six lifestyle and 
#' # behavioural traits after pruning or clumping, which are considered exposures for the risk 
#' # of seven mental health outcome phenotypes
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
#' # Finally, the estimated genetically predicted effects are plotted
#'
#' plot_genPredEffects(output, progress = TRUE, 
#'                    fileName = "Example_plot_genPredEffects")
#'
#' # and the figure is saved as a .pdf file in the directory
#'
#' print(getwd())
#'
#' # Alternatively, the estimated genetically predicted effects are plotted with an order of
#' # appearance that
#' # is changed within the outcomes and within the exposures compared to the order of the 
#' # traits in the MrDAG output
#'
#' ord <- c(c(6, 2, 1, 3, 4, 7, 5), c(9, 11, 13, 8, 12, 10))
#' plot_genPredEffects(output, ord = ord, 
#'                    fileName = "Example_alternative_plot_genPredEffects")
#'
#' # or along with 95% CIs
#'
#' plot_genPredEffects(output, CI = 0.95, ord = ord, 
#'                    fileName = "Example_alternative_plot_genPredEffects_w_CIs", 
#'                    plotPar = list(geom_text.size = 4))
#' }


########## plot_genPredEffects ##########

plot_genPredEffects <- function(output, donotplot = FALSE, CI = NULL, ord = NULL, progress = FALSE, 
                               threshold = 0, labels = NULL, digits = 2, 
                               fileType = "pdf", fileName = "plot_genPredEffects", filePath = NULL, 
                               plotPar = list(limit = NULL, cex.axis = NULL, cex.lab = NULL, geom_text.size = NULL))
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
  
  if (is.null(CI) == FALSE)
  {
    genPredEffects_TMP <- get_genPredEffects(output, ord = ord, progress = progress)
  } else {
    genPredEffects_TMP <- get_genPredEffects(output, CI = CI, ord = ord, progress = progress)
  }
  
  genPredEffects <- genPredEffects_TMP$genPredEffects
  genPredEffects_LL <- genPredEffects_TMP$genPredEffects_LL
  genPredEffects_UL <- genPredEffects_TMP$genPredEffects_UL
  genPredEffects[abs(genPredEffects) < threshold] <- 0
  genPredEffects_LL[abs(genPredEffects) < threshold] <- 0
  genPredEffects_UL[abs(genPredEffects) < threshold] <- 0
  rm(genPredEffects_TMP)
  
  if (!is.null(labels))
  {
    rownames(genPredEffects) <- rownames(genPredEffects_LL) <- rownames(genPredEffects_UL) <- labels
    colnames(genPredEffects) <- colnames(genPredEffects_LL) <- colnames(genPredEffects_UL) <- labels
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
  
  if (is.null(CI) == FALSE)
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
  
  if (is.null(CI) == FALSE)
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
