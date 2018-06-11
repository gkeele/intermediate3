# Joint likelihood ratio for target and mediator
#
#' Fit joint likelihood for `target` and `mediator` given `driver`.
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
#' @return Data frame with `id` and `LR` as well as `annotation` columns.
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
#' med_joint <- mediation_joint(target = target,
#'                       mediator = Tmem68$mediator[, med_signif, drop = FALSE],
#'                       driver = driver,
#'                       annotation = med_lod,
#'                       covar_tar = Tmem68$covar,
#'                       covar_med = Tmem68$covar)
#' plot(med_joint)
#' @export
#'
mediation_joint <- function(target, mediator, driver, annotation,
                          covar_tar=NULL, covar_med=NULL, kinship=NULL,
                          driver_med = NULL, intcovar = NULL,
                          fitFunction = fitQtl2,
                          index_name = "pos",
                          ...) {
  
  if(is.null(mediator))
    return(NULL)
  
  # If only one mediator, replicate it to match annotation.
  mediator <- as.matrix(mediator)
  if(ncol(mediator) == 1) {
    mediator <- mediator[, rep(1, nrow(annotation)), drop = FALSE]
    colnames(mediator) <- annotation$id
  } else {
    stopifnot(ncol(mediator) == nrow(annotation))
  }
  
  result <- mediation_test_internal(target, mediator, driver, annotation,
                                    covar_tar, covar_med, kinship,
                                    driver_med, intcovar,
                                    fitFunction, NULL,
                                    fit_joint,
                                    ...)
    
  out <- dplyr::left_join(
    dplyr::bind_rows(
      result,
      .id = "id"),
    annotation,
    by = "id")
  attr(out, "index_name") <- index_name
  class(out) <- c("mediation_joint", class(out))
  
  out
}

fit_joint <- function(object, driver, target, 
                      kinship, covar_tar, covar_med,
                      driver_med, intcovar,
                      fitFunction, testFunction,
                      common = TRUE, 
                      ...) {
  
  # Make sure we have driver or driver_med.
  if(!is.null(driver_med)) {
    if(is.array(driver_med))
      driver_med <- driver_med[,, object[[2]]$driver]
    else # must be list
      driver_med <- driver_med[[object[[2]]$driver]]
  }
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
                   intcovar, common = common,
                   fit_list = c("m.d_m","t.md_t"), ...)
  data.frame(LR = sum(fits$LR))
}

#' @export
#' @rdname mediation_joint
plot.mediation_joint <- function(x, ...)
  ggplot_mediation_joint(x, ...)
#' @export
#' @rdname mediation_joint
autoplot.mediation_joint <- function(x, ...)
  ggplot_mediation_joint(x, ...)
#' @export
#' @rdname mediation_joint
ggplot_mediation_joint <- function(x, lod = TRUE,
                                   xlab = index_name, ylab = ylab_name, ...) {
  index_name <- attr(x, "index_name")
  if(index_name != "index" & "index" %in% names(x)) {
    # Make sure we don't clash with column named index.
    x$index <- NULL
  }
  x <- dplyr::rename(x, index = index_name)
  if(lod) {
    x <- dplyr::mutate(x, LR = LR / log(10))
    ylab_name <- "LOD"
  } else {
    ylab_name <- "LR"
  }
  p <- ggplot2::ggplot(x) +
    ggplot2::aes(index, LR)
  if("pattern" %in% names(x))
    p <- p + ggplot2::aes(col = pattern)
  p +
    ggplot2::geom_point() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
}
