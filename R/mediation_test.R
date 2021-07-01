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
#' @importFrom dplyr any_of arrange bind_rows desc filter group_by left_join
#' mutate one_of rename ungroup
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 aes autoplot element_blank facet_grid facet_wrap 
#' geom_hline geom_point geom_vline ggplot 
#' ggtitle scale_color_manual scale_shape_manual theme xlab ylab
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
#' 
#' @return List with elements:
#' - best best fit table
#' - test causal test results in table
#' - driver list of driver names for target and mediator(s)
#' - normF Frobenius norm if using both target and mediator drivers
#' - params list of parameter settings for use by summary and plot methods
#' 
#' @examples
#' data(Tmem68, package = "Tmem68")
#' # Focus on chromosome 13
#' Tmem68 <- Tmem_helper(Tmem68, "13")
#'  
#' target <- Tmem68$target
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
#' 
#' ggplot_mediation_test(med_test)
#'
#' @export
#' @importFrom dplyr arrange bind_rows everything filter left_join matches
#'             mutate rename select
#' @importFrom purrr map transpose
#'
mediation_test <- function(target, mediator, driver, annotation = NULL,
                          covar_tar=NULL, covar_med=NULL, kinship=NULL,
                          driver_med = NULL, intcovar = NULL,
                          test = c("wilcoxon","binomial","joint","normal"),
                          fitFunction = fitDefault,
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
  
  if(!is.matrix(mediator)) {
    mediator <- as.matrix(mediator)
    if(is.null(colnames(mediator)))
      colnames(mediator) <- "mediator"
  }
  
  # Convert any blank driver names to V1, V2, ...
  driver <- driver_blank_names(driver)
  driver_med <- driver_blank_names(driver_med)

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
  
  result$driver_names <- {
    if(is.null(driver_med)) NULL
    else {
      if(is.array(driver_med)) dimnames(driver_med)[[3]]
      else names(driver_med)
    }
  }
  
  result$driver_levels <- {
    if(!is.null(driver))
      colnames(driver)
    else {
      if(is.null(driver_med)) NULL
      else {
        if(is.array(driver_med)) colnames(driver_med)
        else colnames(driver_med[[1]])
      }
    }
  }
  
  if(is.null(annotation))
    annotation <- data.frame(id = colnames(mediator),
                             stringsAsFactors = FALSE)
  
  decided <- dplyr::filter(
    result$test,
    .data$model_test != "undecided")
  
  result$best <-
    dplyr::select(
      dplyr::arrange(
        dplyr::rename(
          dplyr::left_join(
            dplyr::left_join(
              # join undecided with other tests
              dplyr::rename(
                dplyr::select(
                  dplyr::filter(
                    result$test,
                    .data$model_test == "undecided"),
                  .data$id, .data$pvalue),
                undecided = "pvalue"),
              dplyr::left_join(
                # Join best test with annotation.
                dplyr::bind_rows(
                  purrr::map(
                    split(
                      decided,
                      decided$id),
                      function(x) x[which.min(x$pvalue)[1],, drop = FALSE])),
                annotation,
               by = "id"),
              by = "id"),
            # Join with mediation and mediator LR.
            dplyr::rename(
              tidyr::pivot_wider(
                dplyr::filter(
                  dplyr::select(
                    result$fit,
                    .data$id, .data$response, .data$LR),
                  .data$response %in% c("mediation", "mediator")),
                names_from = "response", values_from = "LR"),
              LRmed = "mediator"),
            by = "id"),
          triad = "model_test"),
        .data$pvalue),
      .data$id, 
      dplyr::matches("symbol"),
      dplyr::matches(facet_name),
      dplyr::matches(index_name),
      .data$triad, .data$pvalue, .data$alt,
      .data$undecided, .data$mediation, .data$LRmed, 
      dplyr::everything())
  
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
  
  use_1_driver <- is.null(annotation$driver_names) | is.null(driver_med)
  if(use_1_driver & !is.null(driver_med))
    driver_med <- NULL
  
  # Get common data.
  commons <- common_data(target, mediator, driver,
                         covar_tar, NULL, kinship, intcovar = intcovar,
                         common = use_1_driver, ...)
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
  if(!is.null(annotation)) {
    annotation <- dplyr::filter(
      annotation,
      .data$id %in% colnames(mediator))
    # Make sure annotation is in same order as mediator.
    m <- match(colnames(mediator), annotation$id)
    if(any(is.na(m))) {
      cat("mediator and annotation do not match\n", file = stderr())
      return(NULL)
    }
    driver_names <- annotation$driver_names
    if(!is.null(driver_names))
      driver_names <- driver_names[m]
  } else {
    driver_names <- NULL
  }
  if(is.null(driver_names))
    driver_names <- rep("", ncol(mediator))
  names(driver_names) <- colnames(mediator)

  # Workhorse: CMST on each mediator.
  mediator <- as.data.frame(mediator, make.names = FALSE)

  purrr::map(
    purrr::transpose(
      list(
        mediator = mediator,
        driver_names = driver_names)),
    cmstfn, driver, target, 
    kinship, covar_tar, covar_med,
    driver_med, intcovar,
    fitFunction, testFunction, common, ...)
}

#' @export
#' @param x object of class \code{mediation_test}
#' @param not_type biotypes to not include
#' @rdname mediation_test
#' 
subset.mediation_test <- function(x, not_type, ...) {
  attrc <- class(x)
  x$best <- dplyr::filter(x$best, 
                          .data$biotype != not_type)
  class(x) <- attrc
  
  x
}
#' @param object object of class \code{mediation_test}
#' @export
#' @rdname mediation_test
#' 
summary.mediation_test <- function(object, ..., lod = FALSE) {
  out <-
    dplyr::mutate(
      dplyr::arrange(
        object$best,
        pmin(.data$pvalue, .data$undecided), .data$id),
      mediation = .data$mediation / log(10),
      pvalue = signif(.data$pvalue, 3),
      mediation = signif(.data$mediation, 3),
      LRmed = signif(.data$LRmed, 3))

  if(lod) {
    out <- 
      dplyr::rename(
        dplyr::mutate(
          out,
          mediation = .data$mediation / log(10),
          LRmed = .data$LRmed / log(10)),
        lodMediator = "LRmed")
  }
  
  out
}
driver_blank_names <- function(driver) {
  if(is.null(driver)) 
    return(NULL)
  
  if(is.array(driver)) {
    is_blank <- ("" == colnames(driver))
    if(any(is_blank)) {
      colnames(driver)[is_blank] <- paste0("V", seq_len(sum(is_blank)))
    }
  } else {
    for(i in seq_along(driver)) {
      is_blank <- ("" == colnames(driver[[i]]))
      if(any(is_blank)) {
        colnames(driver[[i]])[is_blank] <- paste0("V", seq_len(sum(is_blank)))
      }
    }
  }
  driver
}
