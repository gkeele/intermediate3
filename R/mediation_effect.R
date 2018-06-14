#' Driver effect on target and mediator
#' 
#' Compare driver effect on mediator, target and target adjusted by mediator.
#' 
#' @param object data frame with effects
#' @param driver_levels levels of driver
#' 
#' @export
#' 
#' @examples
#' example("mediation_test")
#' out <- mediation_effect(med_test, "symbol")
#' ggplot_mediation_effect(out)
#' 
mediation_effect <- function(object,
                          id_name = "^id$",
                          driver_levels = unique(unlist(object$driver))) {
  coefs <- 
    tidyr::gather(
      dplyr::rename(
        dplyr::inner_join(
          dplyr::select(
            object$best,
            id, triad, pvalue, dplyr::matches(id_name)),
          dplyr::mutate(
            dplyr::select(
              object$fit,
              id, response, dplyr::one_of(driver_levels)),
            response = ifelse(response == "mediation", "adjusted", response),
            response = factor(response, c("mediator","target","adjusted"))),
          by = "id"),
        mediator = "id"),
      level, effect, -mediator, -dplyr::matches(id_name), -triad, -pvalue, -response)
  
  if(length(m <- grep(id_name, names(coefs)))) {
    coefs$mediator <- coefs[[m]]
    coefs[[m]] <- NULL
  }

  class(coefs) <- c("mediation_effect", class(coefs))
  coefs
}
#' @export
ggplot_mediation_effect <- function(object,
                              colors = qtl2::CCcolors,
                              max_facet = 12) {
  udriver <- unique(object$level)
  if(length(colors) != length(udriver)) {
    colors <- seq_along(udriver)
    names(colors) <- udriver
  }
  driver_names <- names(colors)
  names(driver_names) <- udriver
  
  object <- 
    dplyr::mutate(
      dplyr::filter(
        dplyr::mutate(
          dplyr::arrange(
            object,
            pvalue),
          mediator = paste0(mediator, " (", signif(pvalue, 2), " ", triad, ")")),
        mediator %in% unique(mediator)[seq_len(max_facet)]),
      level = factor(driver_names[level], names(colors)),
      mediator = factor(mediator, unique(mediator)))
  
  tmpfn <- function(object) {
    if(!nrow(object))
      return(NULL)
    sum_object <- 
      purrr::map(
        split(object, object$response),
        function(x) c(mean = mean(x$effect), sd = sd(x$effect)))
    sum_object$adjusted <- sum_object$target
    object <- split(object, object$response)
    for(i in names(object)) {
      object[[i]] <-
        dplyr::mutate(
          object[[i]],
          effect = (effect - sum_object[[i]]["mean"]) / sum_object[[i]]["sd"])
    }
    dplyr::bind_rows(object)
  }
  object <- 
    dplyr::bind_rows(
      purrr::map(
        split(object, object$mediator),
        function(x) {
          dplyr::bind_rows(
            purrr::map(
              split(object, object$triad),
              tmpfn))
        }))
  
  ggplot2::ggplot(object) +
    ggplot2::aes(response, effect, color = level, group = level) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size=2) +
    ggplot2::facet_wrap(~ mediator) +
    ggplot2::scale_color_manual(name = "driver",
                                values = colors)
  
}
#' @export
autoplot.mediation_effect <- function(...) ggplot_mediation_effect(...)
