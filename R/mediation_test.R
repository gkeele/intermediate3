# Mediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix of mediators
#' @param driver vector or matrix with driver values
#' @param annotation optional annotation data frame for mediators
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param driver_med optional driver matrix for mediators
#' @param intcovar optional interactive covariates (assumed same for `mediator`` and `target``)
#' @param test Type of CMST test.
#' @param pos Position of driver.
#' @param fitFunction function to fit models with driver, target and mediator
#' @param data_type Type of mediator data.
#' @param verbose verbose messages if `TRUE`
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
#' 
#' @examples
#' data(Tmem68)
#' m <- match("Tmem68", Tmem68$annotation$symbol) 
#' # Add noise to target, which is exactly Tmem68Rmediator[,"Tmem68"]
#' target <- Tmem68$target
#' target <- target + rnorm(length(target), sd = 0.5)
#' # Find mediators with significant effect
#' med_lod <- mediator_lod(mediator = Tmem68$mediator,
#'                         driver = Tmem68$qtl.geno,
#'                         annotation = Tmem68$annotation,
#'                         covar_med = NULL)
#' med_signif <- med_lod$id[med_lod$lod >= 5]
#' # Add info column.
#' med_lod$info <- paste("chr =", med_lod$chr)
#' 
#' med_test <- mediation_test(target = target,
#'                       mediator = Tmem68$mediator[, med_signif, drop = FALSE],
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = med_lod,
#'                       covar_tar = Tmem68$covar,
#'                       method = "double-lod-diff")
#' summary(med_test)
#' ggplot2::autoplot(med_test)
#'
#' @export
#'
mediation_test <- function(target, mediator, driver, annotation,
                          covar_tar=NULL, covar_med=NULL, kinship=NULL,
                          driver_med = NULL, intcovar = NULL,
                          test = c("wilcoxon","binomial","joint","normal"),
                          pos = NULL,
                          fitFunction = fitQtl2,
                          data_type = c("phenotype","expression"),
                          verbose = FALSE,
                          ...) {
  
  ## Need to enable different covariates for different mediators.

  if(is.null(mediator))
    return(NULL)
  
  data_type = match.arg(data_type)
  # Need following in annotation for each data_type
  #   id = identifier of mediator
  #   biotype = type of biological measurement
  cmstfn <- switch(data_type,
                   expression = cmst_default,
                   phenotype = cmst_pheno)
  
  test <- match.arg(test)
  testfn <- switch(test,
                   wilcoxon = wilcIUCMST,
                   binomial = binomIUCMST,
                   joint    = normJointIUCMST,
                   normal   = normIUCMST)

  pos_tar <- pos
  
  # Make sure covariates are numeric
  covar_tar <- covar_df_mx(covar_tar)
  intcovar <- covar_df_mx(intcovar)
  
  scan_max <- fitFunction(driver, target, kinship, covar_tar)
  LR_tar <- scan_max$LR
  
  use_1_driver <- is.null(annotation$driver) | is.null(driver_med)
  if(use_1_driver & !is.null(driver_med))
    driver_med <- NULL
  
  # Get common data.
  commons <- common_data(target, mediator, driver,
                         covar_tar, NULL, kinship, intcovar = intcovar,
                         common = use_1_driver)
  if(is.null(commons))
    return(NULL)
  
  target <- commons$target
  mediator <- commons$mediator
  driver <- commons$driver
  kinship <- commons$kinship
  covar_tar <- commons$covar_tar
  intcovar <- commons$intcovar
  common <- commons$common
  rm(commons)

  # Two reasons not to put covar_med in common_data call:
  # 1: different mediators may have different covariates
  # 2: covar_med is data frame, so need to be careful.
  # Fix up covar_med to match the rest
  m <- match(rownames(driver), rownames(covar_med), nomatch = 0)
  covar_med <- covar_med[m,, drop = FALSE]
  
  # Reorganize annotation and mediator data.
  # Need to make sure elements of mediator have same ids.
  annotation <- dplyr::filter(
    annotation,
    id %in% colnames(mediator))
  # Make sure annotation is in same order as mediator.
  m <- match(colnames(mediator), annotation$id)
  if(any(is.na(m))) {
    cat("mediator and annotation do not match\n", file = stderr())
    return(NULL)
  }
  annotation <- annotation[m,]

  # Workhorse: CMST on each mediator.
  mediator <- as.data.frame(mediator)
  best <- purrr::map(
    purrr::transpose(list(
      mediator = mediator,
      annotation = 
        purrr::transpose(
          # Make sure order is maintained to match mediator.
          dplyr::mutate(
            annotation,
            id = factor(id, id))))),
    cmstfn, driver, target, 
    kinship, covar_tar, covar_med,
    driver_med, intcovar,
    fitFunction, testfn, common, verbose)

  best <- dplyr::rename(
    dplyr::bind_rows(best, .id = "id"),
    triad = ref)
  
  # Kludge until I figure out why last level.
#  relabel <- c("causal", "reactive", "independent", "correlated")
#  names(relabel) <- c("m.d_t.m", "t.d_m.t", "t.d_m.d", "t.md_m.d")
#  best$triad <- factor(relabel[best$triad], relabel)
#  best$alt <- factor(relabel[best$alt], relabel)
  
  result <- list(
    best = dplyr::arrange(
      dplyr::left_join(best, annotation, by = "id"),
      pvalue),
    params = list(pos = pos_tar,
                  LR = LR_tar,
                  target = colnames(target),
                  data_type = data_type),
    targetFit = scan_max)
    
  class(result) <- c("mediation_test", class(result))
  result
}
#' @export
subset.mediation_test <- function(object, not_type, ...) {
  attrc <- class(object)
  object$best <- dplyr::filter(object$best, 
                               biotype != not_type)
  class(object) <- attrc
  
  object
}
#' @export
summary.mediation_test <- function(object, ..., lod = FALSE) {
  out <- dplyr::select(
    dplyr::mutate(
      dplyr::arrange(
        object$best,
        pvalue),
      mediation = mediation / log(10),
      pvalue = signif(pvalue, 3),
      pos = round(pos, 2),
      mediation = signif(mediation, 3),
      LRmed = signif(LRmed, 3)),
    id, symbol, chr, pos, mediation, triad, pvalue, LRmed,
    dplyr::everything())
  
  if(lod) {
    out <- 
      dplyr::rename(
        dplyr::mutate(
          out,
          mediation = mediation / log(10),
          LRmed = LRmed / log(10)),
        lodMediator = "LRmed")
  }
  
  out
}
