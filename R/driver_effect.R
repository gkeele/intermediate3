#' Driver effect on target and mediator
#' 
#' Compare driver effect on mediator, target and target adjusted by mediator.
#' 
#' @param object data frame with effects
#' @param driver_levels levels of driver
#' 
#' @export
#' 
#' @examples
#' data(Tmem68)
#' # Find mediators with significant effect
#' # Find and remove Tmem68 from mediators because it is target.
#' m <- match("Tmem68", Tmem68$annotation$symbol)
#' Tmem68$annotation[m,]
#' med_lod <- mediator_lod(mediator = Tmem68$mediator[,-m],
#'                         driver = Tmem68$qtl.geno,
#'                         annotation = Tmem68$annotation[-m,],
#'                         covar_med = NULL)
#' med_signif <- med_lod$lod >= 5
#' # Add info column.
#' med_lod$info <- paste("chr =", med_lod$chr)
#' 
#' med_test <- mediation_test(target = Tmem68$target,
#'                       mediator = Tmem68$mediator[, med_signif, drop = FALSE],
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = med_lod,
#'                       covar_tar = Tmem68$covar,
#'                       method = "double-lod-diff")
#' (sum_med <- summary(med_test))
#' 
driver_effect <- function(out, driver_levels = LETTERS[1:8]) {
  out1 <- out[, c("target","group","mediator","pvalue",
                  paste0(driver_levels, "_m"),
                  driver_levels,
                  paste0(driver_levels, "_p"))]
  out1 <- 
    dplyr::arrange(
      dplyr::ungroup(
        dplyr::mutate(
          dplyr::group_by(
            dplyr::mutate(
              tidyr::gather(
                out1,
                driver_level, effect, -target, -group, -mediator, -pvalue),
              fitType = stringr::str_remove(driver_level, "[A-H]_*"),
              driver_level = stringr::str_remove(driver_level, "_[a-z]"),
              fitType = ifelse(fitType == "", "a", fitType),
              fitType = c(m="Mediator", p="Target", a="Adjusted")[fitType],
              fitType = factor(fitType, c("Mediator", "Target", "Adjusted"))),
            target, group, mediator, fitType),
          effect = effect - mean(effect))),
      pvalue)
  class(out1) <- c("driver_effect", class(out1))
  out1
}
#' @export
ggplot_driver_effect <- function(out1,
                              colors = qtl2::CCcolors,
                              max_facet = 12) {
  udriver <- unique(out1$driver_level)
  if(length(colors) != length(udriver)) {
    colors <- seq_along(udriver)
    names(colors) <- udriver
  }
  driver_names <- names(colors)
  names(driver_names) <- udriver
  out1 <- 
    dplyr::mutate(
      dplyr::filter(
        dplyr::mutate(
          out1,
          mediator = paste0(mediator, " (", signif(pvalue, 2), ")")),
        mediator %in% unique(mediator)[seq_len(max_facet)]),
      driver_level = factor(driver_names[driver_level], names(colors)),
      mediator = factor(mediator, unique(mediator)))
  ggplot2::ggplot(out1) +
    ggplot2::aes(fitType, effect, color = driver_level, group = driver_level) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size=2) +
    ggplot2::facet_wrap(~ mediator) +
    ggplot2::ggtitle(paste(out1$target[1], "for group", out1$group[1])) +
    ggplot2::scale_color_manual(name = "driver",
                                values = colors)
  
}
#' @export
autoplot.driver_effect <- function(...) ggplot_driver_effect(...)
