#' Helper function to extract subset chromosomes from Tmem68
#' 
#' @param object list object with mediator, annotation, and other elements
#' @param chr chromosome(s) to select (all if \code{NULL})
#' 
#' @export
Tmem_helper <- function(object, chr = NULL) {
  if(is.null(chr))
    return(object)
  
  id <- which(object$annotation$chr %in% chr)
  if(!length(id))
    return(object)

  object$mediator <- object$mediator[,id]
  object$annotation <- object$annotation[id,]
  object
}