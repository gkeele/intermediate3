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
                                         minpvalue = 0.05, ...) {
  # Remake listof as mediation_qtl2 object.
  out <- bind_mediation_qtl2(x)
  out$best <-
    dplyr::filter(
      out$best,
      pvalue <= minpvalue)
  
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
          dplyr::rename(
            dplyr::select(
              dplyr::filter(
                out$best,
                triad == plot_type),
              -pattern),
            pattern = "symbol"),
          pattern = reorder(pattern, -pvalue))
      if(!nrow(out$best))
        return(NULL)
      ggplot_mediation_qtl2(out) +
        ggplot2::guides(color = ggplot2::guide_legend("Mediator"))
    })
}