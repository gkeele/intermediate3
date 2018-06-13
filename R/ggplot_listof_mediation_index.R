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
  # Remake one element as mediation_index object.
  out <- x[[1]]
  out$best <-
    dplyr::filter(
      dplyr::bind_rows(
        purrr::map(
          x,
          function(x) {
            minbest <- min(x$best$pvalue, na.rm = TRUE)
            if(!is.na(minbest) && minbest <= minpvalue)
              x$best
            else
              NULL
          }),
        .id = "mediator_id"),
      !is.na(pvalue))

  plot_type <- match.arg(plot_type)
  target_index <- out$params$target_index

  switch(
    plot_type,
    all = {
      ggplot_mediation_index(out) +
        ggplot2::geom_vline(xintercept = target_index)
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
        ggplot2::geom_vline(xintercept = target_index) +
        ggplot2::guides(color = ggplot2::guide_legend("Mediator"))
    })
}