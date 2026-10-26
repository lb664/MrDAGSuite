# 135 characters #####################################################################################################################
#' @title Convergence diagnostic plots of the \code{MrDAG} algorithm
#' @description Convergence diagnostic plots of the trace of the log-marginal likelihood evaluated during Markov chain Monte Carlo 
#' (MCMC) iterations including the burn-in (grey dotted line) and the trace of the posterior model size, i.e., the number of directed 
#' edges selected during MCMC iterations (after burn-in), both thinned at every thin-th iteration. For details, see 
#' \insertCite{Zuber2025;textual}{MrDAGSuite}
#'
#' @param output Output produced by the \code{MrDAG} algorithm
#' @param donotplot Logical value (default \code{FALSE}) indicating if the plot is run but not plotted
#' @param fileType Type of plot file. Default is \code{pdf}. Other file types include \code{jpeg}, \code{tiff}, and \code{png}
#' @param fileName Name of the plot (default \code{"plot_convDiagnostic"})
#' @param filePath Directory where the plot is saved. If the directory does not exist, \code{plot_convdiagnostic} will create it by  
#' using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), the plot is 
#' saved in the current working directory
#' @param plotPar List of parameters that can be changed to modify the appearance of the plot. Default values are \code{list(cex.axis 
#' = 1.2, cex.lab = 1.5)}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example
#' # Plotting convergence diagnostics of the MrDAG algorithm on the effect of lifestyle and 
#' # behavioural exposures, and measured on mental health phenotypes. 708 independent 
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
#' # Finally, convergence diagnostics are plotted
#'
#' plot_convdiagnostic(output, 
#'                     fileName = "Example_plot_convdiagnostic")
#'
#' # and the figure is saved as a .pdf file in the directory
#'
#' print(getwd())
#' }


########## plot_convdiagnostic ##########

plot_convdiagnostic <- function(output, donotplot = FALSE, 
                                fileType = "pdf", fileName = "plot_convdiagnostic", filePath = NULL, 
                                plotPar = list(cex.axis = 1.2, cex.lab = 1.5))
{
  NotAllowed <- NULL
  value <- NULL
  variable <- NULL
  
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
                                                              paste0("iter", seq(1, niter - burnin, thin))))
    for (iter in 1 : niter_TMP)
    {
      graphs[, , iter] <- bd_decode(output$graphs[iter])
    }
  
  } else {
    graphs <- output$graphs
  }
  
  logMargLik <- output$logMargLik[seq(1, niter, thin)]
  logMargLik[1 : 10] <- rep(logMargLik[11], 10)
  modelSize <- apply(apply(graphs, c(1, 3), sum), 2, sum)
  
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
  
  plotParDefault <- list(cex.axis = 1.2, cex.lab = 1.5)
  
  if (length(match(names(plotPar), names(plotParDefault))) > 0)
  {
    
    for (p in 1 : length(names(plotPar)))
    {
      plotParDefault[[match(names(plotPar)[p], names(plotParDefault))]] <- plotPar[[p]]
    }
  
  }
  
  iter <- seq(1, length(logMargLik), 1)
  data <- data.frame(iter, logMargLik)
  
  xlab <- "Thinned MCMC Iters"
  ylab <- "Marginal log-Likelihood"
  X_intercept <- burnin / thin
  
  fig1 <- ggplot2::ggplot(data = data, ggplot2::aes(x = iter, y = logMargLik)) +
          ggplot2::scale_x_continuous(xlab, 
                                      breaks = seq(0, niter / thin, length.out = 6),
                                      label = paste0(seq(0, niter / (thin * 1000), length.out = 6), "K")) + 
          ggplot2::ylab(ylab) + 
          ggplot2::geom_line() + 
          ggplot2::geom_vline(xintercept = X_intercept + 0.5, linetype = "dashed", color = "grey", linewidth = 1) + 
          ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis)), 
                         axis.title = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.lab)), 
                         axis.text.y = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), color = 'black'), 
                         axis.text.x = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), angle = 45, hjust = 1, color = 'black'), 
                         legend.title = ggplot2::element_blank())
  
  iter <- seq(1, length(modelSize), 1)
  data <- data.frame(iter, modelSize)
  
  xlab <- "Thinned MCMC Iters"
  ylab <- "Model Size"
  
  fig2 <- ggplot2::ggplot(data = data, ggplot2::aes(x = iter, y = modelSize)) +
          ggplot2::scale_x_continuous(xlab, 
                                      breaks = seq(0, niter_TMP, length.out = 5),
                                      label = paste0(seq(burnin / (thin * 1000), niter / (thin * 1000), length.out = 5), "K")) + 
          ggplot2::ylab(ylab) + 
          ggplot2::geom_line() + 
          ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis)), 
                         axis.title = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.lab)), 
                         axis.text.y = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), color = 'black'), 
                         axis.text.x = ggplot2::element_text(size = ggplot2::rel(plotParDefault$cex.axis), angle = 45, hjust = 1, color = 'black'),           
                         legend.title = ggplot2::element_blank())
  
  fig <- ggpubr::ggarrange(fig1, fig2,
                           ncol = 2, nrow = 1)
  
  if (donotplot == TRUE)
  {
    return(fig)
  } else {
    suppressWarnings(print(fig))
    invisible(invisible(grDevices::dev.off()))
  }

}
