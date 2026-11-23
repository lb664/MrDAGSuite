# Internal helper: automatic appearance parameters for the heatmap plots ------------------------------------------------
#
# Any element of plotPar that is NULL (the default) is set from the size of the plot, so that the numbers fit inside the
# tiles and the axis labels stay readable. Elements supplied as numbers are left untouched.
#   cex.axis      relative size of the axis labels, decreasing with the number of traits (at most 1.2, about 1.1 for 15 traits, at least 0.6)
#   cex.lab       relative size of the axis titles, 1.3 times cex.axis
#   geom_text.size size of the numbers in the tiles, chosen so that a label of labelChars characters fits the tile
#                 width of a plotting area of about 12.5 cm (the usable part of a 7 inch square device)
#' @noRd
.autoPlotPar <- function(plotPar, nRows, nCols, labelChars, panelMM = 125)
{
  nMax <- max(nRows, nCols)
  tileMM <- panelMM / nMax
  auto <- list(cex.axis = min(1.2, max(0.6, 4.3 / sqrt(nMax))))
  auto$cex.lab <- min(1.5, max(0.8, 1.3 * auto$cex.axis))
  auto$geom_text.size <- min(5, max(1.5, 0.85 * (panelMM / nCols) / (labelChars * 0.502)))
  for (nm in names(auto))
  {
    if (nm %in% names(plotPar) && is.null(plotPar[[nm]])) plotPar[[nm]] <- auto[[nm]]
  }
  plotPar
}
