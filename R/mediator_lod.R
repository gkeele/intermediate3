#' Mediator LOD
#'
#' For a given QTL haplotype probabilities \code{driver} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates \code{driver} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{driver} to \code{target}.
#'
#'
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param driver A matrix, haplotype probabilities at QTL we try to mediate
#' @param annotation A data frame with mediators' annotation, must include columns "chr" and "pos"
#' @param covar_med A matrix with additive covariates
#' @param intcovar A matrix of covariate interacting with driver
#' @param kinship kinship object
#' @param facet_name name of facet column (default `chr`)
#' @param verbose If TRUE display information about the progress
#' 
#' @seealso \code{\link{plot.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' med_lod <- mediator_lod(mediator = Tmem68$mediator,
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = Tmem68$annotation,
#'                       covar_med = NULL)
#' ggplot2::autoplot(med_lod) +
#'   ggplot2::geom_hline(yintercept = 5, col = "blue")
#' 
#' @export

mediator_lod <- function(mediator, 
                         driver, 
                         annotation, 
                         covar_med=NULL, 
                         intcovar = NULL,
                         kinship = NULL,
                         facet_name = "chr",
                         verbose=TRUE) {

  # Make sure covariates are numeric
  covar_med <- covar_df_mx(covar_med)
  
  # Get common data.
  commons <- common_data(, mediator, driver, intcovar, covar_med)
  if(is.null(commons))
    return(NULL)
  
  mediator <- commons$mediator
  driver <- commons$driver
  intcovar <- commons$intcovar
  covar_med <- commons$covar_med
  common <- commons$common
  rm(commons)

  # check input
  stopifnot(NROW(driver) == NROW(mediator))
  stopifnot(is.null(covar_med) | NROW(mediator) == NROW(covar_med))
  if(!is.null(covar_med)) {
    stopifnot(all(is.numeric(covar_med)))
    stopifnot(!any(is.na(covar_med)))
  }
  stopifnot(!any(is.na(driver)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(driver)))
  stopifnot(c(facet_name, "pos") %in% tolower(names(annotation)))

  # Match up annotation with mediators
  stopifnot(all(!is.na(m <- match(colnames(mediator), annotation$id))))
  annotation <- annotation[m,]
  
  # Transpose mediator and annotation.
  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(cbind(mediator)),
           annotation = split(annotation, rownames(annotation))))
  
  # Fit likelihood on subset with no missing mediator data for a given mediator.
  mapfn <- function(x, driver, covar_med, intcovar, kinship) {
    no.na <- !is.na(x$mediator)
    fitDefault(driver[no.na,],
               x$mediator[no.na],
               kinship = kinship[no.na, no.na],
               addcovar = covar_med[no.na,],
               intcovar = intcovar[no.na,])$LR
  }
  
  output <- annotation
  output$lod <- unlist(purrr::map(med_pur, mapfn, driver, covar_med, intcovar, kinship)) / log(10)
  attr(output, "targetFit") <- min(output$lod)
  attr(output, "facet_name") <- facet_name
  class(output) <- c("mediation_scan", "data.frame")
  return(output)
}
