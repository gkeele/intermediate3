#' Mediator Likelihood Ratio
#'
#' For a given QTL haplotype probabilities \code{driver} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LR statistic. A low LR value indicates \code{driver} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{driver} to \code{target}.
#'
#'
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param driver A matrix, haplotype probabilities at QTL we try to mediate
#' @param annotation A data frame with mediators' annotation with columns for `facet_name` and `index_name`
#' @param covar_med A matrix with additive covariates
#' @param intcovar A matrix of covariate interacting with driver
#' @param facet_name name of facet column (default `chr`)
#' @param index_name name of index column (default `pos`)
#' @param verbose If TRUE display information about the progress
#' @param fitFunction function to fit models with driver, target and mediator
#' @param ... additional parameters
#' 
#' @examples
#' data(Tmem68)
#' 
#' med_LR <- mediator_LR(mediator = Tmem68$mediator,
#'                       driver = Tmem68$driver,
#'                       annotation = Tmem68$annotation,
#'                       covar_med = NULL)
#' ggplot2::autoplot(med_LR) +
#'   ggplot2::geom_hline(yintercept = 5, col = "blue")
#' 
#' @export

mediator_LR <- function(mediator, 
                         driver, 
                         annotation, 
                         covar_med=NULL, 
                         intcovar = NULL,
                         facet_name = "chr",
                         index_name = "pos",
                         verbose=TRUE, 
                         fitFunction = fitDefault,
                         ...) {

  # Make sure covariates are numeric
  covar_med <- covar_df_mx(covar_med)
  intcovar <- covar_df_mx(intcovar)
  
  # Get common data.
  commons <- common_data(NULL, mediator, driver, NULL, covar_med, intcovar = intcovar)
  if(is.null(commons))
    return(NULL)
  
  mediator <- commons$mediator
  driver <- commons$driver
  covar_med <- commons$covar_med
  common <- commons$common
  intcovar <- commons$intcovar
  rm(commons)

  # check input
  stopifnot(NROW(driver) == NROW(mediator))
  stopifnot(is.null(covar_med) | NROW(mediator) == NROW(covar_med))
  if(!is.null(covar_med)) {
    stopifnot(!any(is.na(covar_med)))
  }
  if(!is.null(driver)) {
    stopifnot(!any(is.na(driver)))
  }
  stopifnot(c(facet_name, index_name) %in% tolower(names(annotation)))

  # Match up annotation with mediators
  stopifnot(all(!is.na(m <- match(colnames(mediator), annotation$id))))
  annotation <- annotation[m,]
  
  # Transpose mediator and annotation.
  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(cbind(mediator)),
           annotation = split(annotation, rownames(annotation))))
  
  # Fit likelihood on subset with no missing mediator data for a given mediator.
  mapfn <- function(x, driver, covar_med, intcovar, med_rownames, ...) {
    mediator <- x$mediator
    names(mediator) <- med_rownames
    fitFunction(driver, mediator, addcovar = covar_med, intcovar = intcovar, ...)$LR
  }
  
  output <- annotation
  output$LR <- unlist(
    purrr::map(med_pur, mapfn, driver, covar_med, intcovar,
               rownames(mediator), ...))
  attr(output, "targetFit") <- min(output$LR)
  attr(output, "facet_name") <- facet_name
  attr(output, "index_name") <- index_name
  class(output) <- c("mediation_scan", "data.frame")
  return(output)
}
