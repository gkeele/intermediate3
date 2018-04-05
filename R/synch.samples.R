#' Synch up sample IDs between objects.
#' 
#' If each object exists, we get the intersection of the sample IDs,
#' sort them and subset each object.
#' 
#' @param pheno matrix or data.frame containing the phenotypes. Sample IDs must be in rownames. Required.
#' @param probs numeric matrix containing the founder allele probabilities. Sample IDs must be in rownames. Required.
#' @param expr  numeric matrix containing the gene expression data. Sample IDs must be in rownames. Required.
#' @param covar numeric matrix containing the mapping covariates. Sample IDs must be in rownames. Optional.
#' 
#' @return list with three elements: pheno, probs, expr and covar. Sample IDs will all be in the same order.
#' @export
synch.samples = function(pheno, probs, expr, covar) {
  
  samples = intersect(rownames(pheno), rownames(probs))
  samples = intersect(samples, rownames(expr))
  if(!missing(covar)) {
    samples = intersect(samples, rownames(covar))
  }
  
  if(length(samples) == 0) {
    stop("There are no samples in common. Please check the rownames on all objects.")
  }

  samples = sort(samples)
  pheno = pheno[samples,,drop = FALSE]
  probs = probs[samples,,drop = FALSE]
  expr  = expr[samples,,drop = FALSE]

  if(!missing(covar)) {
    covar = covar[samples,,drop = FALSE]
  }

  message(paste("Scanning with", nrow(pheno), "samples."))
  
  return(list(pheno = pheno, probs = probs, expr = expr, covar = covar))

} # synch.samples()
