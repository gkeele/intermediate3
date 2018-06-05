# Mediation tests across index
#
#' Test mediation across set of indexed drivers
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix of mediators
#' @param driver vector or matrix with driver values
#' @param annotation optional annotation data frame for mediators
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param driver_med driver array for mediators
#' @param driver_index index to driver array
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom qtl2 decomp_kinship fit1 get_common_ids
#' @importFrom dplyr arrange bind_rows desc filter group_by left_join mutate one_of rename ungroup
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes autoplot element_blank facet_grid facet_wrap 
#' geom_hline geom_point geom_vline ggplot 
#' ggtitle scale_color_manual scale_shape_manual theme xlab ylab
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#'
mediation_index <- function(target, mediator, driver = NULL,
                            annotation = NULL, covar_tar = NULL, covar_med = NULL, kinship = NULL,
                            driver_med = NULL, driver_index = NULL,
                            index_name = "pos", ...) {
  # Mediation test over interval
  
  nmed <- dim(driver_med)[3]
  
  # Propagate mediator over third dimension of driver_med.
  mediator <- mediator[, rep(1, nmed), drop = FALSE]
  colnames(mediator) <- dimnames(driver_med)[[3]]
  
  # Propagate annotation over 
  annotation <- annotation[rep(1, nmed),, drop = FALSE]
  annotation$id <- colnames(mediator)
  annotation$driver <- colnames(mediator)
  stopifnot(length(driver_index) == nmed)
  annotation[[index_name]] <- driver_index
  
  #   run mediation test and find the best models (using BIC among the four models)
  out <- intermediate::mediation_test(
    target = target,
    mediator = mediators,
    annotation = annotation,
    covar_tar = covar,
    covar_med = covar,
    kinship = kinship,
    driver = driver,
    driver_med = driver_med,
    index_name = index_name, ...)
  
  class(out) <- c("mediation_index", class(out))
  out
}
