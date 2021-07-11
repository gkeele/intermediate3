#' Triad scatter plot for mediator and target
#' 
#' Triad plot. Currently relies on \code{sdp} to provide lines, but want to use
#' coefficients from model fit with \code{\link{mediation_test}} to get lines for
#' each column of driver. Note that the plot uses column \code{info} to provide
#' additional information, which here is the \code{chr} of mediator. The plot uses
#' the mediator position on its home chromosome, which is not really what is wanted.
#' 
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param driver vector or matrix with driver values
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param fitFunction function to fit models with driver, target and mediator
#' @param ... additional arguments
#' 
#' @examples
#' data(Tmem68)
#' 
#' target <- cbind(Tmem68 = Tmem68$target)
#' 
#' # Pick strongest mediator that is not target.
#' m <- match("Nnt", Tmem68$annotation$symbol)
#' mediator <- Tmem68$mediator[, m, drop = FALSE]
#' colnames(mediator) <- "Nnt"
#' 
#' med_triad <- mediation_triad(target = target,
#'                       mediator = mediator,
#'                       driver = Tmem68$driver,
#'                       covar_tar = Tmem68$covar,
#'                       sdp = 2)
#' 
#' summary(med_triad)
#' 
#' ggplot_mediation_triad(med_triad, tname = "Tmem68", mname = "Nnt")
#' 
#' @export
#' 
#' @importFrom stringr str_split
#' @importFrom ggplot2 aes autoplot facet_wrap geom_hline geom_smooth 
#' geom_text ggplot ggtitle scale_color_discrete xlab ylab
#' @importFrom broom tidy
#' @importFrom stats lm
#' 
mediation_triad <- function(target, mediator, driver,
                        covar_tar = NULL, covar_med = NULL,
                        fitFunction = fitDefault,
                        ...) {
  
  # Make sure covariates are numeric
  covar_tar <- covar_df_mx(covar_tar)
  covar_med <- covar_df_mx(covar_med)
  
  # Convert any blank driver names to V1, V2, ...
  driver <- driver_blank_names(driver)

  # Would like to have option to have line per haplo.
  # But that requires some regression style approach, such as dividing up data
  # or fitting allele model. Guess is this would involve fitting allele model,
  # getting estimates of slopes for each ellele interacted with mediator,
  # creating data frame, and adding this to ggplot object.
  # Probably signal this with sdp = NULL option?
  
  # Fit target and target|mediator models
  fit <- med_fits(driver, target, mediator,
                  fitFunction, covar_tar, covar_med)
  
  dat <- triad_data(target, mediator, driver, 
                    covar_tar, covar_med, ...)
  
  for(i in c("t.d_t","t.md_t.m")) {
    tmp <- fit$coef[[i]][seq_len(ncol(driver))]
    dat[[i]] <- c(as.matrix(dat[names(tmp)]) %*% tmp)
  }
  
  # Need to account for covariates and sex.
  out <- list(data = dat,
              coef = fit$coef[["t.d_t"]], coef_med = fit$coef[["t.md_t.m"]],
              drivers = colnames(driver), med_name = colnames(mediator))
  
  class(out) <- c("mediation_triad", class(out))
  
  out
}
triad_data <- function(target, mediator, driver, 
                       covar_tar, covar_med,
                       sdp = NULL,
                       allele = TRUE,
                       label_fn = pattern_label,
                       group_fn = pattern_sdp,
                       ...) {
  
  # Find common data.
  commons <- common_data(target, mediator, driver, 
                         covar_tar, covar_med)
  
  # Get covariations from covar_med that are not in covar_tar
  if(!is.null(covar_med)) {
    cov_names <- colnames(covar_med)[!(colnames(covar_med) %in% colnames(covar_tar))]
    commons$covar_med <- commons$covar_med[,cov_names, drop = FALSE]
  } else {
    commons$covar_med <- matrix(NA, length(commons$target), 0)
  }
  
  # Set names for target and mediator (can change in plot)
  for(i in c("target","mediator"))
    colnames(commons[[i]]) <- i
  
  # Set up point labels and groups.
  if(is.null(label_fn))
    label_fn <- function(driver, allele)
      toupper(substr(colnames(driver), 1, 1))[apply(driver, 1, function(x) which.max(x)[1])]
  label <- label_fn(commons$driver, allele)
  if(is.null(group_fn))
    group_fn = function(label, a, b) label
  group <- as.character(group_fn(label, sdp, colnames(commons$driver)))
  
  dat <- data.frame(commons$driver, commons$target, commons$mediator)
  if(length(commons$covar_tar))
    dat <- data.frame(dat, commons$covar_tar)
  if(length(commons$covar_med))
    dat <- data.frame(dat, commons$covar_med)
  dat <- data.frame(dat, label = label, group = group)
  
  if(!is.null(dat$sex))
    dat$Sex <- c("Female", "Male")[1 + dat$sex]

  dat
}
#' @param object object of class \code{mediation_triad}
#' 
#' @rdname mediation_triad
#' @export
#' 
summary.mediation_triad <- function(object, ...) {
  lm_tidy <- function(object, driver) {
    form <- formula(paste("target ~ 0 + mediator",
      ifelse(match("Sex", colnames(object$data), nomatch = 0),
             "* Sex +",
             "+"),
      driver))
    broom::tidy(stats::lm(form, object$data))
  }
  dplyr::bind_rows(
    driver = lm_tidy(object, "group"),
    allele = lm_tidy(object, paste(object$drivers, collapse = "+")),
    .id = "model")
}
#' @param x object of class \code{mediation_triad}
#' @param tname target name (default \code{"target"})
#' @param mname mediator name (default \code{"mediator"})
#' @param dname driver name (default \code{"driver"})
#' @param centerline horizontal line at value (default = \code{NULL}); set to \code{NA} for no line or \code{NULL} for mean
#' @param fitline include fit line from coefficients in \code{x} if \code{TRUE}
#' @param main main title (defautl \code{tname})
#' @param colors named colors to use if \code{fitline} is \code{TRUE}
#' @param size size of text (default \code{2})
#' 
#' @rdname mediation_triad
#' @export
ggplot_mediation_triad <- function(x, 
                             tname = "target", mname = "mediator", dname = "driver",
                             centerline = NULL, fitline = FALSE,
                             main = paste(tname, "by", mname, "and", dname),
                             colors = seq_len(nrow(dat)),
                             size = 2,
                             ...) {
  
  p <- ggplot2::ggplot(x$data) +
    ggplot2::ggtitle(main)
  
  if("label" %in% names(x$data)) {
    p <- p + 
      ggplot2::aes(label = .data$label) +
      ggplot2::geom_text(size = size)
  } else {
    p <- p +
      ggplot2::geom_point(alpha = 0.2)
  }
  
  if("Sex" %in% names(x$data)) {
    p <- p +
      ggplot2::facet_wrap(~ .data$Sex)
  }
  
  # set up mediator and target.
  p <- p + 
    ggplot2::aes(.data$mediator, .data$target) +
    ggplot2::xlab(mname) +
    ggplot2::ylab(tname)

  if(is.null(centerline)) {
    centerline <- mean(x$data$target, na.rm = TRUE)
  }
  if(!is.na(centerline)) {
    p <- p +
      ggplot2::geom_hline(yintercept = centerline, linetype = "dashed", col = "grey60")
  }

  if(fitline) {
    # Add fitted model line.
    dat <- data.frame(slope = x$coef_med[x$med_name],
                      intercept = x$coef_med[x$drivers],
                      col = x$drivers,
                      row.names = x$drivers)
    if(nrow(dat) == length(colors)) {
      if(!is.null(names(colors)))
        dat$col <- names(colors)
      else
        names(colors) <- dat$col
      dat$col <- factor(dat$col, names(colors))
    }
    p <- p +
      ggplot2::geom_abline(
        ggplot2::aes(slope = .data$slope,
                     intercept = .data$intercept,
                     col = .data$col),
        data = dat)
    if(nrow(dat) == length(colors)) {
      p <- p +
        ggplot2::scale_color_manual(name = dname,
                                    values = colors)
    }
    
  } else {
    p <- p + 
      ggplot2::aes(col = .data$group) +
      ggplot2::scale_color_discrete(name = dname) +
      ggplot2::geom_smooth(method = "lm", se=FALSE, formula = "y ~ x")
  }
  p
}
#' @param object object of class \code{mediation_triad}
#' 
#' @rdname mediation_triad
#' @export
autoplot.mediation_triad <- function(object, ...) {
  ggplot_mediation_triad(object, ...)
}
#' @rdname mediation_triad
#' @export
plot.mediation_triad <- function(x, ...) {
  ggplot_mediation_triad(x, ...)
}

