#' Mediation Scan
#'
#' For a given QTL haplotype probabilities \code{driver} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates \code{driver} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{driver} to \code{target}.
#'
#' @param target A numeric vector with gene/protein expression
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param driver A matrix, haplotype probabilities at QTL we try to mediate
#' @param annotation A data frame with mediators' annotation, must include columns `chr`` and `pos``
#' @param covar A matrix with additive covariates
#' @param intcovar A matrix of covariate interacting with driver
#' @param kinship kinship object
#' @param method A method to handle missing cases
#' @param fitFunction function to fit models
#' @param verbose If TRUE display information about the progress
#' 
#' @seealso [ggplot_mediation()]
#' 
#' @examples
#' data(Tmem68)
#' # Find and remove Tmem68 from mediators because it is target.
#' m <- match("Tmem68", Tmem68$annotation$symbol)
#' Tmem68$annotation[m,]
#' med_scan <- mediation_scan(target = Tmem68$target,
#'                       mediator = Tmem68$mediator[,-m],
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = Tmem68$annotation[-m,],
#'                       covar = Tmem68$covar,
#'                       method = "double-lod-diff")
#' ggplot2::autoplot(med_scan)
#' ggplot2::autoplot(subset(med_scan, "4")) +
#'   ggplot2::geom_vline(xintercept = Tmem68$annotation[m,"pos"], linetype = "dashed")
#' 
#' @export

mediation_scan <- function(target, 
                           mediator, 
                           driver, 
                           annotation,
                           covar=NULL,
                           intcovar=NULL,
                           kinship = NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff"), 
                           fitFunction = fitDefault,
                           facet_name = "chr",
                           verbose=TRUE) {
  
  # *** Something is wrong! ***

  # Make sure covariates are numeric
  covar <- covar_df_mx(covar)
  
  # Get common data.
  commons <- common_data(target, mediator, driver, covar, intcovar)
  if(is.null(commons))
    return(NULL)
  
  target <- commons$target
  mediator <- commons$mediator
  driver <- commons$driver
  covar <- commons$covar_tar
  intcovar <- commons$covar_med
  common <- commons$common
  rm(commons)

  # check input
  stopifnot(all(is.numeric(target)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(driver)))
  stopifnot(all(is.numeric(covar)))
  stopifnot(all(is.numeric(intcovar)))
  stopifnot(c(facet_name, "pos") %in% tolower(names(annotation)))
  method = match.arg(method)

  # data preparation
  mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
  N <- nrow(mediator) # number of individuals
  if (is.null(covar)) covar <- cbind(rep(1, N)) # if no covariates, use just intercept
  
  # Match up annotation with mediators
  stopifnot(all(!is.na(m <- match(colnames(mediator), annotation$id))))
  annotation <- annotation[m,]
  
  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(mediator),
           annotation = split(annotation, rownames(annotation))))
  lodfn <- 
    switch(
      method,
      ignore            = function(loglik, loglik0) loglik[1],
      "lod-diff"        = function(loglik, loglik0) loglik[2] - loglik[1],
      "double-lod-diff" = function(loglik, loglik0) loglik0 - (loglik[2] - loglik[1]))
  
  no.na <- !is.na(target)
  loglik0 <- fitFunction(driver[no.na,], target[no.na,], kinship[no.na, no.na],
                        covar[no.na,], intcovar[no.na,])$LR

  mapfn <- function(x, target, covar, driver, loglik0) {
    no.na <- !is.na(target) & !is.na(x$mediator)
    loglik <- c(fitFunction(driver[no.na,], target[no.na,], kinship[no.na, no.na],
                           cbind(covar[no.na,], x$mediator[no.na]),
                           intcovar[no.na,])$LR,
                fitFunction(driver[no.na,], target[no.na,], kinship[no.na, no.na],
                           covar[no.na,], intcovar[no.na,])$LR)
    lodfn(loglik, loglik0)
  }
  output <- annotation
  # Compute LOD (fitFunction provides LL, so divide by log(10))
  output$lod <- unlist(purrr::map(med_pur, mapfn, target, covar, driver, loglik0)) /
    log(10)
  attr(output, "targetFit") <- loglik0 / log(10)
  attr(output, "facet_name") <- facet_name
  class(output) <- c("mediation_scan", "data.frame")
  return(output)
}
#'
#' @export
subset.mediation_scan <- function(object, chrs=NULL, ...) {
  facet_name <- attr(object, "facet_name")
  if(is.null(chrs))
    return(object)
  
  new_object <- object[object[[facet_name]] %in% chrs,]
  modify_object(object, new_object)
}
