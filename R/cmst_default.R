cmst_default <- function(object, driver, target, 
                         kinship, covar_tar, covar_med,
                         driver_med, intcovar,
                         fitFunction, testFunction,
                         common = TRUE, 
                         flavor = "B",
                         fitRelate = TRUE,
                         undecided_penalty = FALSE,
                         ...) {
  
  # Make sure we have driver or driver_med.
  driver_med <- get_driver_med(driver_med, object)
  if(is.null(driver)) {
    if(!is.null(driver_med))
      driver <- driver_med
    else {
      stop("must supply driver or driver_med")
    }
  }

  # Force x (= mediator column) to be matrix.
  mediator <- as.matrix(object[[1]])
  colnames(mediator) <- "mediator"
  rownames(mediator) <- rownames(driver)
  
  # Make sure covariates are numeric
  covar_med <- covar_df_mx(covar_med)
  covar_tar <- covar_df_mx(covar_tar)
  intcovar <- covar_df_mx(intcovar)
  
  # Fit models
  fits <- med_fits(driver, target, mediator, fitFunction,
                   kinship, covar_tar, covar_med, driver_med,
                   intcovar, common = common, ...)
  
  combos <- combo_models(fitRelate)
  models <-
    purrr::transpose(
      purrr::map(combos[,1:4],
                 combine_models, fits[c("LR", "indLR", "df")]))
  
  models$LR <- unlist(models$LR)
  models$indLR <- as.data.frame(models$indLR)
  models$df <- unlist(models$df)
  
  n_ind <- length(models$indLR[[1]])
  
  # P-values
  if(undecided_penalty | flavor == "N") {
    pvalues <- testFunction(models, flavor = flavor)
  } else {
    # Models without "undecided"
    modelu <- 
      purrr::transpose(
        purrr::transpose(models)[c("causal", "reactive", "independent")])
    modelu$LR <- unlist(modelu$LR)
    modelu$indLR <- as.data.frame(modelu$indLR)
    modelu$df <- unlist(modelu$df)
    
    pvalues <- testFunction(models, flavor = "N")
    pvalue <- testFunction(modelu, flavor = flavor)
    m <- match(pvalue$ref, pvalues$ref)
    pvalues[m,] <- pvalue
  }
  
  # causal model tests
  test <- 
    dplyr::select(
      dplyr::mutate(
        dplyr::rename(
          pvalues,
          model = "ref",
          pvalue = "pv"),
        alt = factor(.data$alt, .data$model),
        model = factor(.data$model, .data$model),
        LR = models$LR,
        df = models$df,
        IC = .data$LR - .data$df * penalty(n_ind, flavor) / 2),
      .data$model, .data$LR, .data$df, .data$IC, .data$pvalue, .data$alt)

  # target and mediator fits
  fit <- apply(fits$LR * combos[,5:7], 2, sum)
  coefs <- fits$coef[c(1,2,5)]
  names(coefs) <- names(fit)
  coef_names <- names(coefs[[1]])
  coef_names <- coef_names[coef_names != ""]
  coefs <- sapply(coefs, function(x, coef_names) {
    x[match(coef_names, names(x))]
  }, coef_names)
  coefs <- t(coefs)
  coefs <- 
    dplyr::select(
      dplyr::mutate(
        as.data.frame(coefs, stringsAsFactors = FALSE),
        response = rownames(coefs),
        LR = fit),
      .data$response, .data$LR, dplyr::everything())
  
  list(test = test, fit = coefs, fitsLR = fits$LR, normF = fits$normF)
}
combo_models <- function(fitRelate) {
  combos <- 
    matrix(
      0, 6, 7,
      dimnames = list(
        c("t.d_t", "m.d_m", "t.m_t", "m.t_m", "t.md_t.m","t.md_t"),
        c("causal", "reactive", "independent", "undecided",
          "target", "mediator", "mediation")))
  combos[  c(2,3), 1] <- 1 # causal: m.d_t.m
  combos[  c(1,4), 2] <- 1 # reactive: t.d_m.t
  combos[  c(1,2), 3] <- 1 # independent: t.d_m.d
  if(fitRelate)
    combos[  c(2,6), 4] <- 1 # undecided: t.md_m.d
  else
    combos[c(2,3,5), 4] <- 1 # undecided: t.md_m.d
  combos[       1, 5] <- 1 # target contrast: t.d_t
  combos[       2, 6] <- 1 # mediator contrast: m.d_m
  combos[       5, 7] <- 1 # mediation contrast: t.md_t.m
  as.data.frame(combos)
}
combine_models <- function(combos, fits) {
  list(LR = sum(fits$LR * combos),
       indLR = fits$indLR %*% combos,
       df = sum(fits$df * combos),
       coef = fits$coef)
}
get_driver_med <- function(driver_med, object) {
  if(!is.null(driver_med)) {
    driver_names <- object$driver_names
    if(driver_names == "")
      stop("must supply driver_names in annotation if including driver_med")
    
    if(is.array(driver_med))
      driver_med <- driver_med[,, driver_names]
    else # must be list
      driver_med <- driver_med[[driver_names]]
  }
  driver_med
}
