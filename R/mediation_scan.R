#' Mediation Scan
#'
#' For a given QTL haplotype probabilities \code{qtl.geno} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates \code{qtl.geno} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{qtl.geno} to \code{target}.
#'
#'
#' @param target A numeric vector with gene/protein expression
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param annotation A data frame with mediators' annotation, must include columns "chr" and "pos"
#' @param qtl.geno A matrix, haplotype probabilities at QTL we try to mediate
#' @param covar A matrix with additive covariates
#' @param method A method to handle missing cases
#' @param verbose If TRUE display information about the progress
#' 
#' @seealso \code{\link{plot.mediation}}, \code{\link{kplot}}
#' 
#' @examples
#' data(Tmem68)
#' med <- mediation_scan(target = Tmem68$target,
#'                       mediator = Tmem68$mediator,
#'                       annotation = Tmem68$annotation,
#'                       covar = Tmem68$covar,
#'                       qtl.geno = Tmem68$qtl.geno,
#'                       method = "double-lod-diff")
#' plot(med)
#' @export

mediation_scan <- function(target, 
                           mediator, 
                           annotation, 
                           qtl.geno, 
                           covar=NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff", 
                                    "lod-ratio"), 
                           verbose=TRUE) {

  # calculates log10-Likelihood of linear model y ~ 1 + X
  LL <- function(y, X) {
    -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
  }

  # check input
  stopifnot(NROW(target) == NROW(mediator))
  stopifnot(NROW(annotation) == NCOL(mediator))
  stopifnot(NROW(qtl.geno) == NROW(target))
  stopifnot(is.null(covar) | NROW(target) == NROW(covar))
  stopifnot(!any(is.na(covar)))
  stopifnot(!any(is.na(qtl.geno)))
  stopifnot(all(is.numeric(target)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(qtl.geno)))
  stopifnot(all(is.numeric(covar)))
  stopifnot(c("CHR", "POS") %in% toupper(names(annotation)))
  method = match.arg(method)

  # data preparation
  mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
  N <- ncol(mediator) # number of points to scan
  if (is.null(covar)) covar <- cbind(rep(1, N)) # if no covariates, use just intercept

  if (method == "double-lod-diff") {
    no.na <- !is.na(target)
    LOD0 <- LL(target[no.na], cbind(covar, qtl.geno)[no.na,]) - LL(target[no.na], covar[no.na,])
  }

  med_pur <- 
    purrr::transpose(
      list(mediator = as.data.frame(mediator),
           annotation = split(annotation, rownames(annotation))))
  lodfn <- 
    switch(
      method,
      ignore = function(loglik, LOD0) {
        loglik[2] - loglik[1]
      },
      "lod-diff" = function(loglik, LOD0) {
        loglik[4] - loglik[3] - (loglik[2] - loglik[1])
      },
      "double-lod-diff" = function(loglik, LOD0) {
        LOD0 - (loglik[4] - loglik[3] - (loglik[2] - loglik[1]))
      },
      "lod-ratio" = function(loglik, LOD0) {
        (10^loglik[2]-10^loglik[1]) / (10^loglik[4] - 10^loglik[3])
      })
  
  mapfn <- function(x, target, covar, qtl.geno, LOD0) {
    no.na <- !is.na(target) & !is.na(x$mediator)
    loglik <- c(LL(target[no.na], cbind(covar[no.na,], x$mediator[no.na])),
                LL(target[no.na], cbind(covar[no.na,], x$mediator[no.na], 
                                        qtl.geno[no.na,])),
                LL(target[no.na], covar[no.na,]),
                LL(target[no.na], cbind(covar[no.na,], qtl.geno[no.na,])))
    lodfn(loglik, LOD0)
  }
  output <- annotation
  output$LOD <- unlist(purrr::map(med_pur, mapfn, target, covar, qtl.geno, LOD0))
  class(output) <- c("mediation", "data.frame")
  return(output)
}
