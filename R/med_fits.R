med_fits <- function(driver, target, mediator, fitFunction,
                     covar_tar=NULL, covar_med=NULL,
                     driver_med = NULL, intcovar = NULL,
                     frobenius = 0.01,
                     fit_list = c("t.d_t", "m.d_m", "t.m_t", "m.t_m", "t.md_t.m","t.md_t"),
                     verbose = FALSE,
                     fitFunction_med = fitFunction, 
                     common = FALSE,
                     ...) {
  
  # Need to look where common_data is found to adjust.
  if(!common) {
    commons <- common_data(target, mediator, driver,
                           covar_tar, covar_med, driver_med, intcovar)
    driver <- commons$driver
    target <- commons$target
    mediator <- commons$mediator
    covar_tar <- commons$covar_tar
    covar_med <- commons$covar_med
    driver_med <- commons$driver_med
    intcovar <- commons$intcovar
  }
  
  # Make sure covariates are numeric
  covar_med <- covar_df_mx(covar_med)
  covar_tar <- covar_df_mx(covar_tar)
  intcovar <- covar_df_mx(intcovar)
  
  if(is.null(driver_med)) {
    driver_med <- driver
    normF <- c(target = NA, mediator = NA)
    perp_tar <- perp_med <- NULL
  } else {
    # Find perpendicular matrices.
    if(verbose) cat("target", file = stderr())
    perp_tar <- perp_frob(driver, driver_med, frobenius, verbose)
    if(verbose) cat("mediator", file = stderr())
    perp_med <- perp_frob(driver_med, driver, frobenius, verbose)
    if(verbose) cat("\ n", file = stderr())
    normF <- c(target = perp_tar$normF, mediator = perp_med$normF)
    perp_tar <- perp_tar$driver
    perp_med <- perp_med$driver
  }
  
  # Fit mediation models.
  fits <- list()
  for(i in fit_list) {
    fits[[i]] <- 
      switch(
        i,
        t.d_t    = fitFunction(driver, target,
                               covar_tar, intcovar, ...),
        m.d_m    = fitFunction_med(driver_med, mediator,
                                   covar_med, intcovar, ...),
        t.m_t    = fitFunction_med(bind_stuff(1, mediator, perp_tar), target,
                                   covar_tar, intcovar, ...),
        m.t_m    = fitFunction(bind_stuff(1, target, perp_med), mediator,
                               covar_med, intcovar, ...),
        t.md_t.m = fitFunction(driver, target, 
                               bind_stuff(covar_tar, mediator, perp_tar), intcovar, ...),
        t.md_t   = fitFunction(bind_stuff(driver, mediator, perp_tar), target, 
                               covar_tar, intcovar, ...))
  }

  # Transpose list of model fits
  fits <- purrr::transpose(fits)
  fits$LR <- unlist(fits$LR)
  fits$indLR <- as.matrix(as.data.frame(fits$indLR))
  fits$df <- unlist(fits$df)
  
  fits$normF <- normF
  
  fits
}
perp_frob <- function(driver, driver_med, frobenius = 0.01, verbose = FALSE) {
  # Find part of driver perpendicular to other driver (driver_med)
  m <- intersect(rownames(driver_med), rownames(driver))
  qrX <- qr(driver_med[m,])
  driver <- driver[m,]
  for(i in seq_len(ncol(driver))) {
    driver[,i] <- qr.resid(qrX, driver[,i])
  }

  # Check Frobenius norm. If too small, then nix this perpendicular piece.
  m <- apply(driver, 1, function(x)!any(is.na(x)))
  if(sum(m) < 3)
    return(NULL)
  driver <- driver[m,]
  normF <- NA
  if(frobenius > 0) {
    normF <- norm(driver, "F") / sqrt(length(driver))
    if(verbose)
      cat(normF, file = stderr())
    if(normF < frobenius)
      return(NULL)
  }
  list(driver = driver, normF = normF)
}
bind_stuff <- function(...) {
  # There has to be a better way to do this.
  stuff <- list(...)
  nr <- unlist(sapply(stuff, function(x) ifelse(is.null(x), 0, nrow(as.matrix(x)))))
  if(all(nr == 1 | nr == max(nr)))
    return(cbind(...))
  
  # Take care of stuff with different number of rows
  nms <- lapply(stuff, function(x) {
    if(is.null(x))
      NULL
    else
      rownames(as.matrix(x))
    })
  nms <- nms[sapply(nms, length) > 0]
  if(length(nms) > 1)
    nmsi <- intersect(nms[[1]],nms[[2]])
  else
    nmsi <- nms[[1]]
  if(length(nms) > 2) for(i in 3:length(nms)) nmsi <- intersect(nmsi, nms[i])
  
  out <- NULL
  for(i in seq_along(stuff)) {
    if(!is.null(stuff[[i]])) {
      if(nr[i] == 1) {
        if(length(out) == 0)
          out <- stuff[[1]]
        else
          out <- cbind(out, stuff[[i]])
      }
      else
        out <- cbind(out, as.matrix(stuff[[i]])[nmsi,, drop = FALSE])
    }
  }
  out
}
