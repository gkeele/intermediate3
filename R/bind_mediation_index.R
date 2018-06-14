#' @export
bind_mediation_index <- function(object, ...) {
  # Remake one element as mediation_index object.
  out <- object[[1]]
  
  object <- purrr::transpose(object)
  for(i in c("test","fit","fitsLR","best","joint")) {
    out[[i]] <-
      dplyr::bind_rows(
        object[[i]],
        .id = "mediator_id")
  }

  isnt <- !sapply(object$normF, is.null)
  out$normF <- {
    if(any(isnt))
      as.data.frame(t(as.data.frame(object$normF[isnt])))
    else
      NULL
  }
  object <- purrr::transpose(object$driver)
  for(i in names(object)) {
    isnt <- !sapply(object[[i]], is.null)
    object[[i]] <- {
      if(any(isnt))
        as.data.frame(t(as.data.frame(object[[i]][isnt])))
      else
        NULL
    }
  }
  out$driver <- object

  out
}