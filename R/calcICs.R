#' Calculate Z scores
#' 
#' @param models Object with model information
#' @param Shat Matrix of variance and covariance estimates
#' @param ICs Calculated information criteria
#' @param flavor Flavor of penalty from \code{c("B","A","N")}
#' @param ... additional parameters
#' 
#' @export
#'
#' @importFrom purrr simplify_all transpose
#' @importFrom stringr str_split
#' @importFrom dplyr mutate
#'
calcZ <- function(models,
                  Shat = calcShat(models),
                  ICs = calcICs(models, flavor),
                  flavor = "B",
                  ...) {
  
  n_ind <- nrow(models$indLR)
  
  LRt <- outer(ICs, ICs, function(x,y) y-x)
  LRt <- LRt[lower.tri(LRt)]
  LRt <- -LRt / (2 * sqrt(n_ind))
  LRt / sqrt(diag(Shat))
}
#' Calculate information criteria
#' 
#' @param models Object with model information
#' @param flavor Flavor of penalty from \code{c("B","A","N")}
#' 
#' @export
#'
calcICs <- function(models, flavor = "B") {
  n_ind <- nrow(models$indLR)
  -2 * models$LR + models$df * penalty(n_ind, flavor)
}
penalty <- function(n_ind, flavor = c("B","A","N")) {
  flavor <- match.arg(flavor)
  switch(flavor, A = 2, B = log(n_ind), N = 0)
}
