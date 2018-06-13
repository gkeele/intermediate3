#' @export
#' @rdname mediation_index
plot.listof_mediation_index <- function(x, ...)
  ggplot_listof_mediation_index(x, ...)
#' @export
#' @rdname mediation_index
autoplot.listof_mediation_index <- function(x, ...)
  ggplot_listof_mediation_index(x, ...)
#' @export
#' @rdname mediation_index
ggplot_listof_mediation_index <- function(x, 
                                          plot_type = c("all","causal","reactive","independent","undecided"),
                                          minpvalue = 0.05, ...) {
  # Remake listof as mediation_index object.
  out <- bind_mediation_index(x)

  plot_type <- match.arg(plot_type)

  switch(
    plot_type,
    all = {
      ggplot_mediation_index(out)
    },
    causal =,
    reactive =,
    independent =,
    undecided = {
      out$best <-
        dplyr::rename(
          dplyr::select(
            dplyr::filter(
              out$best,
              triad == plot_type),
            -pattern),
          pattern = "symbol")
      if(!nrow(out$best))
        return(NULL)
      ggplot_mediation_index(out) +
        ggplot2::guides(color = ggplot2::guide_legend("Mediator"))
    })
}