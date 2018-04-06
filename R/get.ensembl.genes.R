#' Get the ensembl genes on a given set of chromosomes.
#'
#' @param build integer that is the Ensembl build to use.
#' @param chr character vector of chromosomes from which to get genes.
#' @return data.frame with the genes on the given chromosomes.
#' @examples
#' ens = get_ensembl_genes(chr = "19", build = 90)
#' @export
get_ensembl_genes = function(chr = c(1:19, "X", "Y", "MT"), build = 90) {
  
  if(is.null(build) || build == "") {
    stop("Build argument must contain a valid ensembl build ID.")
  }
  
  if(is.null(chr) || build == "") {
    stop("Chr argument must contain a valid chromsome ID.")
  }
  
  hub = AnnotationHub::query(AnnotationHub::AnnotationHub(), c("ensembl", "mus musculus", "gtf", build))
  wh = grep(paste0("Mus_musculus.GRCm38.", build,".gtf"), hub$title)
  
  if(length(wh) == 0) {
    stop(paste("No Ensembl build", build, "is available from the",
               "Bioconductor Annotation Hub."))
  }
  
  ensembl = hub[[names(hub)[wh]]]
  ensembl = ensembl[ensembl$type == "gene"]
  ensembl = ensembl[which(!is.na(match(as.vector(seqnames(ensembl)), chr)))]
  
  retval = data.frame(ensembl = ensembl$gene_id,
                      symbol  = ensembl$gene_name,
                      chr     = seqnames(ensembl), 
                      start   = start(ensembl),
                      end     = end(ensembl),
                      strand  = strand(ensembl),
                      pos     = 0.5 * (start(ensembl) + end(ensembl)),
                      biotype = ensembl$gene_biotype,
                      stringsAsFactors = FALSE)
  
  return(retval)
  
} # get_ensembl_genes() 
