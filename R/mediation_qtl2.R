# Mediation tests for qtl2 data across index
#
#' Test mediation across set of indexed drivers for `qtl2` data
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param annotation optional annotation data frame for mediators
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param genoprob genoprob object of class [qtl2::calc_genoprob()]
#' @param map list of map positions
#' @param drop_lod drop in `LOD` (= `LR/log(10)`) to define index set
#' @param query_variant function to query variant database
#' @param cores number of cores to use
#' @param target_scan optional object from [qtl2::scan1snps()] for target (created if missing)
#' @param ... additional parameters for [mediation_index()]
#'
#' @importFrom qtl2 genoprob_to_snpprob get_common_ids scan1snps top_snps
#' @importFrom dplyr arrange desc distinct filter mutate rename
#' 
#' @return Object of class `mediation_qtl2`, which inherits from class [mediation_index()]
#' 
#' @export
#'
mediation_qtl2 <- function(target, mediator,
                           annotation, covar_tar, covar_med, kinship,
                           genoprobs, map,
                           drop_lod = 1.5, min_lod = 3,
                           query_variant,
                           cores = 1, target_scan, ...) {

  chr_id <- names(genoprobs)
  if(length(chr_id) > 1) {
    warning("only first chromosome used")
    chr_id <- chr_id[1]
    genoprobs <- subset(genoprobs, chr = chr_id)
  }
  start <- min(map[[chr_id]])
  end <- max(map[[chr_id]])
  
  # Find peak for target to get SNP distribution pattern (sdp)
  if(missing(target_scan)) {
    target_scan <-
      qtl2::scan1snps(
        genoprobs = genoprobs[,chr_id],
        map = map, 
        pheno = target,
        kinship = kinship,
        addcovar = addcovar,
        chr = chr_id, start = start, end = end,
        query_func = query_variant,
        cores = cores,
        keep_all_snps = FALSE)
  }

  # Get alleles in genoprobs  
  prob_alleles <- attr(genoprobs, "alleles")
  # Get top SNPs (and their reflection)
  ts <- 
    dplyr::filter(
      qtl2::top_snps(
        target_scan$lod,
        target_scan$snpinfo),
      lod >= max(lod) - drop_lod,
      lod >= min_lod)
  ts_sdp <- unique(ts$sdp)
  ts_sdp <- unique(c(ts_sdp, 2 ^ length(prob_alleles) - 1 - ts_sdp))
  
  # Get SNP probabilities for SNPs with same sdp as target peak.
  m <- which(target_scan$snpinfo$sdp %in% ts_sdp)
  snpinfo <- target_scan$snpinfo[m,, drop = FALSE]
  snpinfo$index <- seq_len(nrow(snpinfo))

  driver_med <-
    qtl2::genoprob_to_snpprob(
      genoprobs,
      snpinfo)[[1]]

  # Create annotation from snpinfo with driver names.
  annot_driver <- dplyr::rename(snpinfo, id = "snp_id")
  annot_driver$driver_names <- dimnames(driver_med)[[3]]

  # Reduce to common IDs
  m <- qtl2::get_common_ids(
    target,
    covar_tar,
    covar_med,
    mediator,
    kinship,
    complete.cases = TRUE)
  target <- target[m,, drop = FALSE]
  mediator <- mediator[m,, drop = FALSE]
  kinship <- kinship[m,m]
  covar_tar <- covar_tar[m,, drop=FALSE]
  covar_med <- covar_med[m,, drop=FALSE]
  
  # Find SNPs with joint LOD within 1.5 of peak.
  med_joint <- 
    dplyr::filter(
      mj <- mediation_joint(
        target,
        mediator,
        NULL,
        annot_driver,
        covar_tar,
        covar_med,
        kinship,
        driver_med),
      LR >= max(LR) - drop_lod * log(10))
  
  # Mediation test indexed by SNPs identified with med_joint.
  med_index <- 
    mediation_index(
      target,
      mediator,
      NULL,
      annotation,
      covar_tar,
      covar_tar,
      kinship,
      driver_med[,, med_joint$id, drop = FALSE],
      driver_index = med_joint$pos, ...)
  
  # Add SNP distribution pattern (sdp).
  m <- match(med_joint$id, med_index$best$id)
  med_index$best$sdp <- med_joint$sdp[m]
  ts_sdp <- unique(med_index$best$sdp)
  ts_sdp <- c(ts_sdp, 2 ^ length(prob_alleles) - 1 - ts_sdp)

  ## Add in all SNPS to best
  # Get all SNPs in interval with same sdp.
  topsnps <- 
    qtl2::index_snps(
      map,
      dplyr::filter(
        query_variant(chr_id, start, end),
        pos >= min(med_index$best$pos),
        pos <= max(med_index$best$pos),
        sdp %in% ts_sdp))
  
  # Match up index for joining.
  m <- match(med_index$best$id, topsnps$snp_id)
  med_index$best$index <- topsnps$index[m]
  
  # Group snp_ids that have same index
  topsnps <-
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::filter(
            topsnps,
            index %in% med_index$best$index),
          sdp, index),
        id = ifelse(n() == 1, # id is snp_id or range of snp_id's
                        snp_id, 
                        paste(snp_id[1], n(), snp_id[n()], sep = "_")),
        pos_start = pos[1],
        pos_end = pos[n()],
        ensembl_gene = paste(unique(ensembl_gene), collapse = ","),
        consequence = paste(unique(consequence), collapse = ",")))

  # Join topsnps to best to get all SNPs.
  med_index$best <- 
    dplyr::mutate(
      dplyr::right_join(
        dplyr::select(
          med_index$best,
          -chr, -pos, -sdp, -id),
        topsnps,
        by = "index"),
      sdp = pmin(sdp, 2 ^ length(prob_alleles) - 1 - sdp),
      pattern = sdp_to_pattern(sdp, prob_alleles))
  
  med_index$joint <- mj
  med_index$map <- map[[chr_id]]
  
  med_index$params$target_LR <- ts$lod[1] * log(10)
  med_index$params$target_index <- ts$pos[1]
  
  class(med_index) <- c("mediation_qtl2", class(med_index))

  med_index
}
#' @export
#' @rdname mediation_qtl2
plot.mediation_qtl2 <- function(x, ...)
  ggplot_mediation_qtl2(x, ...)
#' @export
#' @rdname mediation_qtl2
autoplot.mediation_qtl2 <- function(x, ...)
  ggplot_mediation_qtl2(x, ...)
#' @export
#' @rdname mediation_qtl2
ggplot_mediation_qtl2 <- function(x, response = c("pvalue","IC"), 
                                  pattern_name = "pattern", ...) {
  best <- 
    dplyr::rename(
      dplyr::filter(
        x$best,
        pos_end > pos_start),
      index1 = "pos_start",
      index2 = "pos_end")
  if(!(pattern_name %in% names(best)))
    best$pattern <- "black"
  
  # Fix up x$best for plot as mediation_index object
  x$best <-
    dplyr::select(
      tidyr::gather(
        x$best,
        place, pos, dplyr::contains("pos_")),
      -place)
  

  response <- match.arg(response)
  p <- ggplot_mediation_index(x, response, pattern_name = pattern_name, ...)
  switch(
    response,
    pvalue = { 
      p <- p + 
        ggplot2::geom_segment(
          aes(x = index1, 
              xend = index2, 
              y = -log10(pvalue),
              yend = -log10(pvalue),
              col = get(pattern_name)),
          data = best,
          inherit.aes = FALSE)
      },
    IC = {
      p <- p + 
        ggplot2::geom_segment(
          aes(x = index1, 
              xend = index2, 
              y = IC,
              yend = IC,
              col = get(pattern_name)),
          data = best,
          inherit.aes = FALSE)
      })
  p
}
