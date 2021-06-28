#' Joint Normal IUCMST
#' 
#' @param models Object with model information
#' @param Zscores Calculated Z scores from \code{calcZ}
#' @param Shat Matrix of variance and covariance estimates
#' @param ... additional parameters
#' 
#' @export
#'
normJointIUCMST <- function(models,
                                Zscores = calcZ(models, ...),
                                Shat = calcShat(models),
                                ...) {

  if(length(models$LR) == 1)
    return(data.frame(ref = names(models$LR),
                     alt = "",
                     pv = 1))

  # Expand to data frame with ref, alt, Z.
  Zscores <- left_right(Zscores)

  # Compare pvalues. this is different from comp_pv().
  # Need to verify at some time that it works properly!
  dplyr::mutate(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::mutate(Zscores,
                        ref = factor(.data$ref, unique(.data$ref))),
          .data$ref),
        pv = as.vector(1 - mnormt::pmnorm(
          rep(min(.data$Z),
              length(.data$Z)),
          varcov = corHat(.data$alt, .data$pair, .data$Shat))),
        alt = .data$alt[which.min(.data$Z)][1])),
    ref = as.character(.data$ref))
}
# Compare reference model with all others and get max pvalue.
corHat <- function(alt, pair, Shat) {
  # Goofy way to get sign right on Shat.
  ss <- 2 * (sapply(
    stringr::str_split(pair, ":"),
    function(x) x[[2]]) == alt) - 1
  Shat <- (ss %o% ss) * Shat[pair,pair] 
  if (!corpcor::is.positive.definite(Shat)) {
    Shat <- corpcor::make.positive.definite(Shat)
  }
  stats::cov2cor(Shat)
}
