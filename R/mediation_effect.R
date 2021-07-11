#' Driver effect on target and mediator
#' 
#' Compare driver effect on mediator, target and target adjusted by mediator.
#' 
#' @param test_object object of class \code{mediation_test}
#' @param id_name name of identifier column
#' @param driver_levels levels of driver
#' 
#' @importFrom tidyr pivot_longer
#' 
#' @export
#' 
#' @examples
#' example("mediation_test")
#' out <- mediation_effect(med_test, "symbol")
#' ggplot_mediation_effect(out)
#' 
mediation_effect <- function(test_object,
                          id_name = "id",
                          driver_levels = unique(unlist(test_object$driver_levels))) {
  
  if(id_name != "id" & id_name %in% names(test_object$best)) {
    m <- match(test_object$fit$id, test_object$best$id)
    test_object$fit$id <- test_object$best[[id_name]][m]
    test_object$best$id <- test_object$best[[id_name]]
    test_object$best[[id_name]] <- NULL
  }

  coefs <- 
    tidyr::pivot_longer(
      dplyr::rename(
        dplyr::inner_join(
          dplyr::select(
            test_object$best,
            .data$id, .data$triad, .data$pvalue),
          dplyr::mutate(
            dplyr::select(
              test_object$fit,
              .data$id, .data$response, dplyr::one_of(driver_levels)),
            response = ifelse(.data$response == "mediation",
                              "adjusted",
                              .data$response),
            response = factor(.data$response,
                              c("mediator","target","adjusted"))),
          by = "id"),
        mediator = "id"),
      -(1:4), names_to = "level", values_to = "effect")
  
  if(length(m <- grep(id_name, names(coefs)))) {
    coefs$mediator <- coefs[[m]]
    coefs[[m]] <- NULL
  }

  class(coefs) <- c("mediation_effect", class(coefs))
  coefs
}
#' @param object object of class \code{mediation_effect}
#' @param colors colors for plotting
#' @param max_facet maximum number of facets
#' @param size_geom geometry size for points and lines
#' @param ... additional parameters
#' 
#' @export
#' @rdname mediation_effect
ggplot_mediation_effect <- function(object,
                              colors = colors_default,
                              max_facet = 12,
                              size_geom = 2, ...) {
  # Levels are used for color, but user might change level names by specifying
  # names to the colors. If names for colors not provided, then use unique levels.
  udriver <- sort(unique(object$level))
  
  # Default colors if not provided.
  colors_default = seq_len(length(udriver))

  # Make sure colors have same length as udriver.
  # Repeat if necessary, but then provide names as udriver.
  if(length(colors) != length(udriver)) {
    colors <- rep(colors, length = length(udriver))
    names(colors) <- udriver
  }
  if(is.null(names(colors))) {
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
            .data$pvalue),
          mediator = paste0(.data$mediator, " (", signif(.data$pvalue, 2), " ", .data$triad, ")")),
        .data$mediator %in% unique(.data$mediator)[seq_len(max_facet)]),
      level = factor(driver_names[.data$level], driver_names),
      mediator = factor(.data$mediator, unique(.data$mediator)))
  
  tmpfn <- function(object) {
    if(!nrow(object))
      return(NULL)
    sum_object <- 
      purrr::map(
        split(object, object$response),
        function(x) c(mean = mean(x$effect), sd = stats::sd(x$effect)))
    # Use SD of target for adjusted
    sum_object$adjusted["sd"] <- sum_object$target["sd"]
    object <- split(object, object$response)
    for(i in names(object)) {
      object[[i]] <-
        dplyr::mutate(
          object[[i]],
          effect = (.data$effect - sum_object[[i]]["mean"]) / sum_object[[i]]["sd"])
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
    ggplot2::aes(.data$response, .data$effect, color = .data$level, group = .data$level) +
    ggplot2::geom_line(size = size_geom) +
    ggplot2::geom_point(size = size_geom) +
    ggplot2::facet_wrap(~ .data$mediator) +
    ggplot2::scale_color_manual(name = "driver",
                                values = colors)
  
}
#' @rdname mediation_effect
#' @export
autoplot.mediation_effect <- function(object, ...)
  ggplot_mediation_effect(object, ...)
#' @param x object of class \code{mediation_effect}
#' @param ... additional parameters
#' 
#' @export
#' @rdname mediation_effect
plot.mediation_effect <- function(x, ...)
  ggplot_mediation_effect(x, ...)
