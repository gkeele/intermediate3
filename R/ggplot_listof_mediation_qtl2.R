#' @export
#' @rdname mediation_qtl2
plot.listof_mediation_qtl2 <- function(x, ...)
  ggplot_listof_mediation_qtl2(x, ...)
#' @export
#' @rdname mediation_qtl2
autoplot.listof_mediation_qtl2 <- function(x, ...)
  ggplot_listof_mediation_qtl2(x, ...)
#' @export
#' @rdname mediation_qtl2
ggplot_listof_mediation_qtl2 <- function(x, 
                                         plot_type = c("all","causal","reactive","independent","undecided"),
                                         minpvalue = 0.05,
                                         id_name = "mediator_id", ...) {
  # Remake listof as mediation_qtl2 object.
  out <- bind_mediation_index(x, id_name)
  out$best <-
    dplyr::filter(
      out$best,
      .data$pvalue <= minpvalue)
  
  plot_type <- match.arg(plot_type)
  
  switch(
    plot_type,
    all = {
      ggplot_mediation_qtl2(out)
    },
    causal =,
    reactive =,
    independent =,
    undecided = {
      out$best <-
        dplyr::mutate(
          dplyr::filter(
            out$best,
            .data$triad == plot_type),
          symbol = reorder(.data$symbol, -.data$pvalue))
      if(!nrow(out$best))
        return(NULL)
      ggplot_mediation_qtl2(out, pattern_name = "symbol")
    })
}