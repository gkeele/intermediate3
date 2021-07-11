#' @export
#' @rdname mediation_test
plot.mediation_test <- function(x, ...)
  ggplot_mediation_test(x, ...)
#' @export
#' @rdname mediation_test
autoplot.mediation_test <- function(x, ...)
  ggplot_mediation_test(x, ...)
#' @param x object of class \code{mediation_test}
#' @param type type of plot from \code{c("pos_LR","pos_pvalue","pvalue_LR","alleles","mediator")}
#' @param main title for plot
#' @param maxPvalue maximum p-value with default \code{0.1}
#' @param local_only local only if \code{TRUE} (default \code{FALSE})
#' @param significant whow signficant if \code{TRUE} (default)
#' @param lod show as LOD plot if \code{TRUE} (default)
#' @param target_index include vertical line at target if not \code{NULL}
#' @param colors colors to use with targets and alleles
#' @param size size of points (default 2)
#' @param alpha symbol transparency (default 0.5)
#' @param show show facets by triad model or difference (default \code{"facets"})
#' @param ... additional parameters
#' 
#' @export
#' @rdname mediation_test
ggplot_mediation_test <- function(x, type = c("pos_LR","pos_pvalue","pvalue_LR","alleles","mediator"),
                               main = params$target,
                               maxPvalue = 0.1, 
                               local_only = FALSE, 
                               significant = TRUE,
                               lod = FALSE,
                               target_index = NULL,
                               colors = RColorBrewer::brewer.pal(8, "Dark2"),
                               size = 2,
                               alpha = 0.5,
                               show = c("facets","difference"),
                               ...) {
  
  show <- match.arg(show)
  if(show == "difference") {
    return(ggplot_mediation_test_diff(x, ...))
  }
  type <- match.arg(type)
  if(is.null(local_only) | !("local" %in% names(x$best)))
    local_only <- FALSE
  if(is.null(significant))
    significant <- TRUE
  
  params <- x$params
  unmediated <- params$LR_target
  index_name <- params$index_name
  
  targetFit <- x$targetFit
  x <- x$best
  
  if(lod) {
    unmediated <- unmediated / log(10)
    x$LR_mediation <- x$LR_mediation / log(10)
  }
  
  if(!("symbol" %in% names(x)))
    x <- dplyr::rename(x, symbol = "id")
  
  if(tmp <- any(is.na(x$triad))) {
    warning(paste("some triad are NA:", sum(tmp)))
  }
  
  relabel <- levels(x$triad)

  if(type != "pos_LR")
    significant <- TRUE
  if(significant) {
    x <- dplyr::filter(x, 
                       .data$pvalue <= maxPvalue)
  } else {
    relabel <- c(relabel, 
                 paste0("n.s. (p>", round(maxPvalue, 2), ")"))
    tmp <- as.character(x$triad)
    tmp[x$pvalue > maxPvalue] <- relabel[5]
    x$triad <- factor(tmp, levels = relabel)
  }
  if(!nrow(x))
    stop(paste("no pvalues below", maxPvalue))
  
  x <- dplyr::arrange(x, dplyr::desc(.data$triad))
  
  if(index_name != "index" & "index" %in% names(x)) {
    # Make sure we don't clash with column named index.
    x$index <- NULL
  }
  x <- dplyr::rename(x, pos = index_name)
  
  # For expression, use qtl_pos if not missing.
  if(!type %in% c("alleles","mediator")) {
    if(local_only)
      x <- dplyr::filter(x, local)
    else {
      if("qtl_pos" %in% names(x))
        x <- dplyr::mutate(x,
                           pos = ifelse(.data$local,
                                        .data$pos, .data$qtl_pos))
    }
    
    if(all(c("local","qtl_ct") %in% names(x))) {
      # Set up plot symbol.
      shapes <- c(17,16,2,1)
      names(shapes) <- c("distal", "local", "distal_info", "local_info")
      x <- dplyr::mutate(x,
                         shape = names(shapes)[1 + .data$local + 2 * (.data$qtl_ct > 1)])
    }
  }

  # Colors
  cols <- c(RColorBrewer::brewer.pal(4, "Dark2"), "#CCCCCC")
  names(cols) <- relabel
  
  if(!(type %in% c("alleles","mediator"))) {
    switch(type,
           pos_pvalue = {
             p <- ggplot2::ggplot(x) +
               ggplot2::aes(x = .data$pos,
                            y = -log10(.data$pvalue)) +
               ggplot2::aes(symbol = .data$symbol,
                            LR_mediation = .data$LR_mediation) +
               ggplot2::facet_grid(~ .data$triad, scales = "free_x") +
               ggplot2::xlab("Position (Mbp)") +
               ggplot2::ylab("-log10 of p-value")
             if(!is.null(target_index))
               p <- p +
                 ggplot2::geom_vline(xintercept = .data$target_index,
                                     col = "darkgrey")
           },
           pvalue_LR = {
             p <- ggplot2::ggplot(x) +
               ggplot2::aes(y = .data$LR_mediation,
                            x = -log10(.data$pvalue)) +
               ggplot2::aes(symbol = .data$symbol,
                            position = .data$pos) +
               ggplot2::facet_grid(~ .data$triad, scales = "free_x") +
               ggplot2::geom_hline(yintercept = unmediated,
                                   col = "darkgrey") +
               ggplot2::xlab("-log10 of p-value") +
               ggplot2::ylab("LR for Mediation")
           },
           pos_LR = {
             p <- ggplot2::ggplot(x) + 
               ggplot2::aes(y = .data$LR_mediation,
                            x = .data$pos) +
               ggplot2::aes(symbol = .data$symbol,
                            pvalue = .data$pvalue) +
               ggplot2::geom_hline(yintercept = unmediated,
                                   col = "darkgrey") +
               ggplot2::facet_grid(~triad, scales = "free_x") +
               ggplot2::xlab("Position (Mbp)") +
               ggplot2::ylab("LR for Mediation")
               ggplot2::scale_color_manual(values = cols)
             if(!is.null(target_index))
               p <- p +
                 ggplot2::geom_vline(xintercept = .data$target_index,
                                     col = "darkgrey")
           })
    if("biotype" %in% names(x)) {
      p <- p + ggplot2::aes(col = .data$biotype)
    }
    if("info" %in% names(x)) {
      p <- p + ggplot2::aes(info = .data$info)
    }
    if(exists("shapes")) {
      p <- p + ggplot2::geom_point(aes(shape = .data$shape),
                                   size = size, alpha = alpha) +
        ggplot2::scale_shape_manual(values = shapes)
      if("qtl_pos" %in% names(x)) {
        p <- p + ggplot2::aes(chr = .data$chr, qtl_pos = .data$qtl_pos)
      }
    } else {
      p <- p + 
        ggplot2::geom_point(size = size, alpha = alpha)
    }
    
    p + 
      ggplot2::ggtitle(main)
  } else {
    # Plot target above mediators
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow = 2)))
    
    plotfn <- function(x, type, targetCoef, layout.pos.row, ylabel = type) {
      p <- ggplot2::ggplot(allele_prep(x, type, col = colors)) +
        ggplot2::aes(x = -log10(.data$pvalue),
                     y = .data$value,
                     col = .data$geno,
                     symbol = .data$symbol) +
        ggplot2::facet_wrap(~ .data$triad, scales = "free_x") +
        ggplot2::geom_hline(data = targetCoef, 
                            ggplot2::aes(yintercept = .data$value,
                                         col = .data$geno),
                            linetype = "dashed") +
        ggplot2::geom_point(size = size) +
        ggplot2::ylab(ylabel)
      if(type == "alleles") {
        p <- p + ggplot2::ggtitle(main) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank())
      } else {
        p <- p + 
          ggplot2::theme(strip.text = ggplot2::element_blank())
      }

      print(p,
            vp = grid::viewport(layout.pos.row = layout.pos.row,
                                layout.pos.col = 1))
    }
    
    targetCoef <- target_prep(targetFit, "alleles", col = colors)
    switch(type,
           alleles = {
             plotfn(x, "mediator",
                    dplyr::mutate(targetCoef,
                                  value = (.data$value - min(.data$value)) /
                                    diff(range(.data$value))),
                    2, "mediator effects")
             plotfn(x, "alleles", targetCoef, 1, "allele effects")
           },
           mediator = {
             plotfn(x, "alleles", targetCoef, 1, "allele effects")
             plotfn(x, "mediator",
                    dplyr::mutate(targetCoef,
                                  value = (.data$value - min(.data$value)) /
                                    diff(range(.data$value))),
                    2, "mediator effects")
           })
  }
}
target_prep <- function(targetFit, type, col = RColorBrewer::brewer.pal(8, "Dark2")) {
  targetFit <- as.data.frame(t(targetFit$coef))
  codes <- LETTERS[seq_along(col)]
  m <- switch(type,
              alleles  = match(codes, names(targetFit)),
              mediator = match(paste(codes, "m", sep = "_"), names(targetFit)))
  
  col_names <- names(col)
  targetFit <- targetFit[m]
  names(targetFit) <- col_names
  out <- tidyr::pivot_longer(targetFit, names_to = "geno", values_to = "value")
  out <- dplyr::mutate(out, geno = factor(.data$geno, col_names))
  switch(type,
         alleles  = dplyr::mutate(out, value = .data$value - mean(.data$value)),
         mediator = dplyr::mutate(out, value = (.data$value - min(.data$value)) /
                                    diff(range(.data$value))))
}
allele_prep <- function(x, type, col = RColorBrewer::brewer.pal(8, "Dark2")) {
  codes <- LETTERS[seq_along(col)]
  m <- switch(type,
              alleles  = match(codes, names(x)),
              mediator = match(paste(codes, "m", sep = "_"), names(x)))
  
  col_names <- c("qtl_pos", names(col))
  names(x)[m] <- col_names
  out <- dplyr::group_by(
    tidyr::pivot_longer(
      dplyr::select(
        x,
        .data$symbol, .data$triad, .data$pvalue, .data$chr, 
        dplyr::one_of(col_names)),
      -dplyr::any_of(c("symbol", "triad", "pvalue", "chr", "qtl_pos")),
      names_to = "geno", values_to = "value"),
    .data$symbol
  )
  out <- dplyr::mutate(out, geno = factor(.data$geno, col_names))
  out <- switch(type,
                alleles  = dplyr::mutate(out, value = .data$value - mean(.data$value)),
                mediator = dplyr::mutate(out, value = (.data$value - min(.data$value))
                                         / diff(range(.data$value))))
  dplyr::arrange(
    dplyr::ungroup(out),
    .data$pvalue)
}
ggplot_mediation_test_diff <- function(x, lod = FALSE, ...) {
  dat <- x$best
  dat$LR_target <- x$params$LR_target
  
  if(lod) {
    dat <- dplyr::mutate(
      dat,
      LR_target = .data$LR_target / log(10),
      LR_mediation = .data$LR_mediation / log(10))
    lr_type <- "LOD"
  } else {
    lr_type <- "LR"
  }
  
  ggplot2::ggplot(dplyr::filter(dat, .data$pvalue <= 0.05)) +
    ggplot2::aes(.data$LR_target - .data$LR_mediation, -log10(.data$pvalue),
                 col = .data$triad, symbol = .data$symbol) +
    ggplot2::geom_point(alpha = 1, size = 3) +
    ggplot2::xlab(paste(lr_type, "target vs mediation"))
}