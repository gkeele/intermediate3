#' Plot of conditioned LR scores against index
#'
#' Plot LR statistics calculated by [mediation_scan()] against index
#' using ggplot2.
#'
#' @param x mediation object
#' @param col color of points (default "firebrick4")
#' @param size character size (default 2)
#' @param xlab,ylab X and Y axis label (default `index_name` and "Conditioned LR")
#' @param col_target color for target LR line
#' @param gap gap between facets (default `25`)
#' 
#' @export
#' @importFrom ggplot2 aes autoplot element_blank element_rect facet_grid geom_hline geom_point ggplot theme
#' @importFrom grid unit
#' @importFrom dplyr arrange_at
#' @rdname mediation_scan

ggplot_mediation_scan <- function(x, 
                           col="firebrick4",
                           size = 1,
                           xlab = index_name,
                           ylab = "Conditioned LR",
                           col_target = "blue",
                           gap = 25) {
  
  facet_name <- attr(x, "facet_name")
  index_name <- attr(x, "index_name")
  if(index_name != "index" & "index" %in% names(x)) {
    # Make sure we don't clash with column named index.
    x$index <- NULL
  }

  if(!is.factor(x[[facet_name]]))
    x[[facet_name]] <- factor(x[[facet_name]], unique(x[[facet_name]]))
  
  x <- dplyr::rename(x, index = index_name)
  
  x <- dplyr::arrange_at(x, c(facet_name, "index"))

  p <- ggplot2::ggplot(x) +
    ggplot2::aes(.data$index, .data$LR, symbol = .data$symbol) +
    ggplot2::facet_grid(stats::formula(paste("~", facet_name)),
                        scales = "free_x", space = "free") +
    ggplot2::xlab(xlab)

  if(!is.null(x$col)) {
    p <- p +
      ggplot2::aes(col = col) +
      ggplot2::geom_point(alpha = 0.5, size = size)
    
  } else {
    p <- p +
      ggplot2::geom_point(col = col, alpha = 0.5, size = size)
  }
  # gap between facets
  p <- p +
    ggplot2::theme(panel.spacing = grid::unit(gap / 10000, "npc")) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black",
                                           fill=NA))
  
  # X axis
  if(length(unique(x[[facet_name]])) > 1) {
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
#' @rdname mediation_scan
#'
autoplot.mediation_scan <- function(x, ...) {
  ggplot_mediation_scan(x, ...)
}
#' @export
#' @rdname mediation_scan
#'
plot.mediation_scan <- function(x, ...) {
  ggplot_mediation_scan(x, ...)
}
