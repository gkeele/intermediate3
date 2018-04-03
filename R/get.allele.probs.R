#' Get the haplotype probs nearest to the given chr and position.
#'
#' @param chr character indicating the chromosome to use. This must match one of the names of genoprobs and map.
#' @param pos numeric value in Mb indicating the position on the chromosome.
#' @param genoprobs qtl2-style list of genoprobs containing 3-dimensional arrays of allele probs, one for each chromosome. Samples in rows, 8 founders in columns, markers in slices.
#' @param map qtl2-style list containing numeric vectors of marker positions in Mb.
#' @return n x 8 numeric matrix.
#' @examples
#' ap = get.allele.probs(chr = "1", pos = 5.0, genoprobs = genoprobs, map = map)
#' @export
get.allele.probs = function(chr, pos, genoprobs, map) {
  
  if(is.null(pos) || chr == "") {
    stop("Chr argument must contain a valid chromosome ID.")
  }
  
  if(is.null(pos) || pos == "") {
    stop("Pos argument must contain a valid chromosome position.")
  }
  
  if(!is.numeric(pos)) {
    stop("Pos argument must be a numeric value in Mb.")
  }
  
  if(is.null(genoprobs)) {
    stop("Genoprobs argument must contain a qtl2-style genoprobs object.")
  }  
  
  if(is.null(map)) {
    stop("Map argument must contain a qtl2-style map object.")
  }
  
  # Get the map data for the current chromosome.
  chr.map = map[[chr]]
  if(is.null(chr.map)) {
    stop(paste("Chr", chr, "not found in map."))
  }
  
  # Get the index of the marker closest to the position.
  wh = which.min(abs(chr.map - pos))
  
  # Get the genoprobs on the chromosome at the closest marker.
  gp = genoprobs[[chr]][,,wh]
  
  return(gp)
  
} # get.allele.probs()
