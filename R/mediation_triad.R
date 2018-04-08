#' Triad scatter plot for mediator and target
#' 
#' Triad plot. Currently relies on \code{sdp} to provide lines, but want to use
#' coefficients from model fit with \code{\link{mediation_test}} to get lines for
#' each column of driver. Note that the plot uses column \code{info} to provide
#' additional information, which here is the \code{chr} of mediator. The plot uses
#' the mediator position on its home chromosome, which is not really what is wanted.
#' See package \code{qtl2shiny} for a more elegant use.
#' 
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param driver vector or matrix with driver values
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param fitFunction function to fit models with driver, target and mediator
#' @param sdp SNP distribution pattern for plot colors
#' @param allele Driver has alleles if \code{TRUE}, otherwise allele pairs.
#' 
#' @examples
#' data(Tmem68)
#' # Pick Abhd17a as strongest mediator.
#' m <- match("Abhd17a", Tmem68$annotation$symbol)
#' Tmem68$annotation[m,]
#' 
#' med_triad <- mediation_triad(target = Tmem68$target,
#'                       mediator = Tmem68$mediator[,m, drop = FALSE],
#'                       driver = Tmem68$qtl.geno,
#'                       covar_tar = Tmem68$covar,
#'                       sdp = 1)
#' ggplot2::autoplot(med_triad)
#' 
#' @export
#' 
#' @importFrom stringr str_split
#' @importFrom ggplot2 aes autoplot facet_wrap geom_hline geom_smooth 
#' geom_text ggplot ggtitle scale_color_discrete xlab ylab
#' 
mediation_triad <- function(target, mediator, driver,
                        covar_tar = NULL, covar_med = NULL,
                        kinship = NULL, 
                        fitFunction = fitQtl2,
                        sdp = NULL,
                        allele = TRUE,
                        label_fn = pattern_label,
                        group_fn = pattern_sdp) {
  
  # Make sure covariates are numeric
  covar_tar <- covar_df_mx(covar_tar)
  covar_med <- covar_df_mx(covar_med)

  commons <- common_data(target, mediator, driver, 
                         covar_tar, covar_med, kinship)
  
  if(!is.null(covar_med)) {
    cov_names <- colnames(covar_med)[!(colnames(covar_med) %in% colnames(covar_tar))]
    commons$covar_med <- commons$covar_med[,cov_names, drop = FALSE]
  } else {
    commons$covar_med <- matrix(NA, length(commons$target), 0)
  }
  
  for(i in c("target","mediator"))
    colnames(commons[[i]]) <- i

  label <- label_fn(commons$driver, allele)
  group <- as.character(group_fn(label, sdp, colnames(commons$driver)))
  dat <- data.frame(commons$driver, commons$target, commons$mediator,
                    commons$covar_tar, commons$covar_med,
                    label = label, group = group)
  
  # Would like to have option to have line per haplo.
  # But that requires some regression style approach, such as dividing up data
  # or fitting allele model. Guess is this would involve fitting allele model,
  # getting estimates of slopes for each ellele interacted with mediator,
  # creating data frame, and adding this to ggplot object.
  # Probably signal this with sdp = NULL option?

  if(!is.null(dat$sex))
    dat$Sex <- c("Female", "Male")[1 + dat$sex]
  
  # Fit target and target|mediator models
  fit <- med_fits(driver, target, mediator,
                  fitFunction, kinship, covar_tar, covar_med)
  for(i in names(fit$coef)[1:2]) {
    tmp <- fit$coef[[i]][seq_len(ncol(driver))]
    dat[[i]] <- c(as.matrix(dat[names(tmp)]) %*% tmp)
  }
  
  class(dat) <- c("mediation_triad", class(dat))
  
  dat
}
#' @param x object of class \code{mediation_triad}
#' @param \dots additional parameters for plotting
#' @param type type of plot: one of \code{("by_mediator", "by_target", "driver_offset", "driver")}
#' @param tname target name (default \code{"target"})
#' @param mname mediator name (default \code{"mediator"})
#' @param dname driver name (default \code{"driver"})
#' @param centerline horizontal line at value (default = \code{0}); set to \code{NA} for no line or \code{NULL} for mean
#' @param main main title (defautl \code{tname})
#' 
#' @rdname mediation_triad
#' @export
ggplot_mediation_triad <- function(x, ..., 
                             type = c("by_mediator", "by_target", "driver_offset", "driver"),
                             tname = "target", mname = "mediator", dname = "driver",
                             centerline = 0,
                             main = tname) {
  
  type <- match.arg(type)
  
  p <- ggplot2::ggplot(x) +
    ggplot2::aes(col = group, label = label) +
    ggplot2::scale_color_discrete(name = dname)

  switch(type,
         by_mediator = {
           p <- p + 
             ggplot2::aes(mediator, target) +
             ggplot2::xlab(mname) +
             ggplot2::ylab(tname)
         },
         by_target = {
           p <- p + 
             ggplot2::aes(target, mediator) +
             ggplot2::xlab(tname) +
             ggplot2::ylab(mname)
         },
         driver_offset = {
           p <- p + 
             ggplot2::aes(mediator, t.d_t - t.md_t.m) +
             ggplot2::xlab(mname) +
             ggplot2::ylab(paste(dname, "effect offset"))
         },
         driver = {
           p <- p + 
             ggplot2::aes(t.d_t, t.d_t - t.md_t.m) +
             ggplot2::xlab(paste(dname, "effect")) +
             ggplot2::ylab(paste(dname, "effect offset"))
         })
  if(is.null(centerline)) {
    switch(type,
           by_mediator = {
             centerline <- mean(x$target, na.rm = TRUE)
           },
           by_target = {
             centerline <- mean(x$mediator, na.rm = TRUE)
           },
           driver_offset, driver = {
             centerline <- mean(x$t.d_t - x$t.md_t.m, na.rm = TRUE)
           })
  }
  
  if(!is.na(centerline)) {
    p <- p +
      ggplot2::geom_hline(yintercept = centerline)
  }
    
  p + 
    ggplot2::geom_smooth(method = "lm", se=FALSE) +
    ggplot2::geom_text(size=3) +
    ggplot2::facet_wrap(~ Sex) +
    ggplot2::ggtitle(main)
}
#' @export
autoplot.mediation_triad <- function(x, ...) {
  ggplot_mediation_triad(x, ...)
}

# from qtl2pattern
sdp_to_logical <- function(sdp, haplos = LETTERS[1:8]) {
  sapply(sdp, function(x, haplos) {
    as.logical(intToBits(x)[seq_along(haplos)])
  }, haplos)
}

