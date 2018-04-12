#' GG plot of conditioned LOD scores against genomic positions
#'
#' Plot LOD statistics calculated by \code{mediation.scan} against genomic positions
#' using ggplot2.
#'
#' @param x mediation object
#' @param col color of points (default "firebrick4")
#' @param cex character expansion (default 2)
#' @param ylab Y axis label (default "Conditioned LOD")
#' @param ... additional arguments to \code{\link[qtl2ggplot]{ggplot_scan1}}
#' 
#' @seealso \code{\link{identify.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' # Find and remove Tmem68 from mediators because it is target.
#' m <- match("Tmem68", Tmem68$annotation$symbol)
#' Tmem68$annotation[m,]
#' med_scan <- mediation_scan(target = Tmem68$target,
#'                       mediator = Tmem68$mediator[,-m],
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = Tmem68$annotation[-m,],
#'                       covar = Tmem68$covar,
#'                       method = "double-lod-diff")
#' ggplot2::autoplot(med_scan)
#' 
#' @export
#' @importFrom ggplot2 aes autoplot element_blank element_rect facet_grid geom_hline geom_point ggplot theme
#' @importFrom grid unit

ggplot_mediation_scan <- function(x, 
                           col="firebrick4",
                           cex = 1,
                           ylab = "Conditioned LOD",
                           col_target = "blue",
                           gap = 25,
                           ...) {
  if(!is.factor(x$chr))
    x$chr <- factor(x$chr, unique(x$chr))
  x <- dplyr::arrange(x, chr, pos)

  p <- ggplot2::ggplot(x) +
    ggplot2::aes(pos, lod, symbol = symbol) +
    ggplot2::geom_point(col = col, alpha = 0.5) +
    ggplot2::facet_grid( ~ chr, scales = "free_x", space = "free")

  # gap between chromosomes
  p <- p +
    ggplot2::theme(panel.spacing = grid::unit(gap / 10000, "npc")) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black",
                                           fill=NA))
  
  # X axis
  if(length(unique(x$chr)) > 1) {
    p <- p + 
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(), 
        axis.ticks.x = ggplot2::element_blank())
  }
  
  targetFit <- attr(x, "targetFit")
  if(!is.null(targetFit))
    p <- p + ggplot2::geom_hline(yintercept = targetFit, col = col_target)
  p
}
#' @export
#' @rdname ggplot_mediation_scan
#'
autoplot.mediation_scan <- function(x, ...) {
  ggplot_mediation_scan(x, ...)
}
#' @export
#' @rdname ggplot_mediation_scan
#'
plot.mediation_scan <- function(x, ...) {
  ggplot_mediation_scan(x, ...)
}
