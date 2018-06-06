# Mediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix of mediators
#' @param driver vector or matrix with driver values
#' @param annotation A data frame with mediators' annotation with columns for `facet_name` and `index_name`
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param driver_med optional driver matrix for mediators
#' @param intcovar optional interactive covariates (assumed same for `mediator` and `target`)
#' @param test Type of CMST test.
#' @param fitFunction function to fit models with driver, target and mediator
#' @param facet_name name of facet column (default `chr`)
#' @param index_name name of index column (default `pos`)
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
#' @return List with elements:
#' - best best fit table
#' - test causal test results in table
#' - driver list of driver names for target and mediator(s)
#' - normF Frobenius norm if using both target and mediator drivers
#' - params list of parameter settings for use by summary and plot methods
#' 
#' @examples
#' data(Tmem68)
#'  
#' # Add noise to target, which is exactly Tmem68$mediator[,"Tmem68"]
#' target <- Tmem68$target
#' target <- target + rnorm(length(target), sd = 0.25)
#' 
#' # Reconstruct 8-allele genotype probabilities.
#' driver <- cbind(A = 1 - apply(Tmem68$qtl.geno, 1, sum), Tmem68$qtl.geno)
#' rownames(driver) <- rownames(Tmem68$qtl.geno)
#' 
#' # Find mediators with significant effect
#' med_lod <- mediator_lod(mediator = Tmem68$mediator,
#'                         driver = driver,
#'                         annotation = Tmem68$annotation,
#'                         covar_med = Tmem68$covar)
#' med_signif <- med_lod$id[med_lod$lod >= 5]
#' # Add info column.
#' med_lod$info <- paste("chr =", med_lod$chr)
#' 
#' med_test <- mediation_test(target = target,
#'                       mediator = Tmem68$mediator[, med_signif, drop = FALSE],
#'                       driver = driver,
#'                       annotation = med_lod,
#'                       covar_tar = Tmem68$covar,
#'                       covar_med = Tmem68$covar)
#' summary(med_test)
#' ggplot2::autoplot(med_test)
#'
#' @export
#'
mediation_test <- function(target, mediator, driver, annotation,
                          covar_tar=NULL, covar_med=NULL, kinship=NULL,
                          driver_med = NULL, intcovar = NULL,
                          test = c("wilcoxon","binomial","joint","normal"),
                          fitFunction = fitQtl2,
                          facet_name = "chr",
                          index_name = "pos",
                          ...) {
  
  ## Need to enable different covariates for different mediators.

  if(is.null(mediator))
    return(NULL)
  
  if(!is.null(driver)) {
    scan_max <- fitFunction(driver, target, kinship, covar_df_mx(covar_tar))
  } else {
    scan_max <- NULL
  }
  target_LR <- scan_max$LR
  
  test <- match.arg(test)
  testFunction <- switch(test,
                   wilcoxon = wilcIUCMST,
                   binomial = binomIUCMST,
                   joint    = normJointIUCMST,
                   normal   = normIUCMST)
  
  result <- mediation_test_internal(target, mediator, driver, annotation,
                                    covar_tar, covar_med, kinship,
                                    driver_med, intcovar,
                                    fitFunction, testFunction,
                                    cmst_default,
                                    ...)
  if(is.null(result))
    return(NULL)
    
  result <-
    purrr::map(
      purrr::transpose(result),
      function(x) {
        if(is.data.frame(x[[1]])) {
          dplyr::bind_rows(x, .id = "id")
        } else {
          isnt <- !sapply(x, is.null)
          if(any(isnt))
            as.data.frame(t(as.data.frame(x[isnt])))
          else
            NULL
        }
      })
  
  if(!is.null(result$normF) && all(is.na(result$normF)))
    result$normF <- NULL
  
  result$driver <- 
    list(
      target = colnames(driver),
      mediator = {
        if(is.null(driver_med)) colnames(driver)
        else {
          if(is.array(driver_med)) colnames(driver_med)
          else names(driver_med)
        }
      })
  
  result$best <-
    dplyr::arrange(
      dplyr::rename(
        dplyr::mutate(
          dplyr::left_join(
            dplyr::ungroup(
              dplyr::filter(
                dplyr::group_by(
                  result$test,
                  id),
                pvalue == min(pvalue))),
            annotation, by = "id"),
          mediation = dplyr::filter(result$fit, response == "mediation")$LR,
          LRmed = dplyr::filter(result$fit, response == "mediator")$LR),
        triad = model),
      pvalue)
  
  result$params <-
    list(target_LR = target_LR,
         target_name = colnames(target),
         facet_name = facet_name,
         index_name = index_name,
         test = test)
  result$targetFit <- scan_max
  
  class(result) <- c("mediation_test", class(result))
  result
}
mediation_test_internal <- function(target, mediator, driver, annotation,
                                    covar_tar, covar_med, kinship,
                                    driver_med, intcovar,
                                    fitFunction,
                                    testFunction,
                                    cmstfn = cmst_default,
                                    ...) {
                                    
  # Make sure covariates are numeric
  covar_tar <- covar_df_mx(covar_tar)
  intcovar <- covar_df_mx(intcovar)
  
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
  if(!is.null(covar_med)) {
    m <- match(rownames(mediator), rownames(covar_med), nomatch = 0)
    covar_med <- covar_med[m,, drop = FALSE]
  }
  
  # If we have driver_med, reduce to same subset as others
  if(!is.null(driver_med)) {
    if(is.array(driver_med)) {
      m <- match(rownames(mediator), rownames(driver_med), nomatch = 0)
      driver_med <- driver_med[m,,, drop = FALSE]
    } else {
      if(is.list(driver_med)) {
        m <- match(rownames(mediator), rownames(driver_med[[1]]), nomatch = 0)
        driver_med <- lapply(driver_med, function(x) x[m,, drop = FALSE])
      } else {
        stop("driver_med is neither an array nor a list")
      }
    }
  }
  
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
  purrr::map(
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
    fitFunction, testFunction, common, ...)
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
  # Change out facet and index names for now.
  facet_name <- object$params$facet_name
  index_name <- object$params$index_name
  best <-
    dplyr::rename(
      object$best,
      facet = facet_name,
      index = index_name)
  
  out <- dplyr::select(
    dplyr::mutate(
      dplyr::arrange(
        best,
        pvalue),
      mediation = mediation / log(10),
      pvalue = signif(pvalue, 3),
      mediation = signif(mediation, 3),
      LRmed = signif(LRmed, 3)),
    id, symbol, facet, index, mediation, triad, pvalue, LRmed,
    dplyr::everything())
  
  # Change facet and index names back in now.
  names(out)[3:4] <- c(facet_name, index_name)

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
