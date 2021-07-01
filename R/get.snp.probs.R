#' Get the SNP probs nearest to the given chr and position.
#'
#' @param chr character indicating the chromosome to use. This must match
#'      one of the names of genoprobs and map.
#' @param pos numeric value in Mb indicating the position on the chromosome.
#' @param genoprobs list containing 3-dimensional arrays of allele probs, one 
#'            for each chromosome. In qtl2 format.
#'            Samples in rows, 8 founders in columns, markers in slices.
#' @param map list containing numeric vectors of marker positions in Mb.
#' @param query_fxn function for querying the CC_SNP data base. Obtain this
#'            by calling qtl2:::create_variant_query_func().
#'            
#' @return an n x 1 numeric matrix.
#'                          
#' @export
get.snp.probs = function(chr, pos, genoprobs, map, query_fxn) {
  
  # Get genoprobs closest to the requested position.
  gp = get_allele_probs(chr, pos, genoprobs, map)
  
  # Get SNPs in a window around the requested position.
  variants = query_fxn(chr, pos - 0.25, pos + 0.25)
  
  # Get the SNP closest to the requested position.
  wh = which.min(abs(variants$pos - pos))
  variants = variants[wh,]
  
  strains = c("A_J",       "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ",
              "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ",     "WSB_EiJ")
  alleles = as.numeric(variants[,strains])
  sp = t(tcrossprod(alleles, gp)) - 1.0
  
  return(sp)
  
} # get.snp.probs()

