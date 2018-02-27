#' GG plot of conditioned LOD scores against genomic positions
#'
#' Plot LOD statistics calculated by \code{mediation.scan} against genomic positions
#' using ggplot2.
#'
#' @param med mediation object
#' @param col color of points (default "firebrick4")
#' @param cex character expansion (default 2)
#' @param ylab Y axis label (default "Conditioned LOD")
#' @param ... additional arguments to \code{\link[qtl2ggplot]{ggplot_scan1}}
#' 
#' @seealso \code{\link{identify.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' med <- mediation.scan(target=Tmem68$target,
#'                       mediator=Tmem68$mediator,
#'                       annotation=Tmem68$annotation,
#'                       covar=Tmem68$covar,
#'                       qtl.geno=Tmem68$qtl.geno)
#' ggplot2::autoplot(med)
#' 
#' @export
#' @importFrom qtl2ggplot ggplot_scan1

ggplot_mediation <- function(med, 
                           col="firebrick4",
                           cex = 2,
                           ylab = "Conditioned LOD",
                           ...){
  map <- split(med$pos, factor(med$chr, unique(med$chr)))
  qtl2ggplot::ggplot_scan1(as.matrix(med[,"LOD", drop = FALSE]), 
                           map, 
                           lines = FALSE,
                           col = col, cex = cex, 
                           ylab = ylab,
                           ...)
}
#' @export
#' @export autoplot.mediation
#' @method autoplot mediation
#' @rdname ggplot_mediation
#'
#' @importFrom ggplot2 autoplot
#'
autoplot.mediation <- function(x, ...) {
  ggplot_mediation(x, ...)
}
