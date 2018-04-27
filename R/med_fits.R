med_fits <- function(driver, target, mediator, fitFunction,
                     kinship=NULL, cov_tar=NULL, cov_med=NULL,
                     driver_med = NULL, intcovar = NULL,
                     common = FALSE, ...) {
  
  # Probably want to move this to fitDefault as fitQtl2 does this.
  # Need to look where common_data is found to adjust.
  if(!common) {
    commons <- common_data(target, mediator, driver,
                           cov_tar, cov_med, kinship, driver_med, intcovar)
    driver <- commons$driver
    target <- commons$target
    mediator <- commons$mediator
    kinship <- commons$kinship
    cov_tar <- commons$cov_tar
    cov_med <- commons$cov_med
    driver_med <- commons$driver_med
    intcovar <- commons$intcovar
  }
  if(is.null(driver_med)) {
    driver_med <- driver
    perp_tar <- perp_med <- NULL
  } else {
    # Find perpendicular matrices.
    perp_tar <- driver
    for(i in seq_len(ncol(driver)))
      perp_tar[,i] <- fitFunction(driver, perp_tar[,i, drop = FALSE], kinship, cov_tar, intcovar)
    perp_med <- driver_med
    for(i in seq_len(ncol(driver_med)))
      perp_tar[,i] <- fitFunction(driver_med, perp_med[,i, drop = FALSE], kinship, cov_med, intcovar)
  }
  
  # Fit mediation models.
  # Transpose list of model fits
  fits <- purrr::transpose(list(
    t.d_t    = fitFunction(driver, target, kinship, cov_tar, intcovar),
    m.d_m    = fitFunction(driver_med, mediator, kinship, cov_med, intcovar),
    t.m_t    = fitFunction(cbind(1, mediator, perp_tar), target, kinship, cov_tar, intcovar),
    m.t_m    = fitFunction(cbind(1, target, perp_med), mediator, kinship, cov_med, intcovar),
    t.md_t.m = fitFunction(driver, target, kinship, cbind(cov_tar, mediator, perp_tar), intcovar)))
  fits$LR <- unlist(fits$LR)
  fits$indLR <- as.matrix(as.data.frame(fits$indLR))
  fits$df <- unlist(fits$df)
  
  fits
}
