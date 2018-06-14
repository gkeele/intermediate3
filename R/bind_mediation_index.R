#' @export
bind_mediation_index <- function(object, ...) {
  # Remake one element as mediation_index object.
  out <- object[[1]]
  
  object <- purrr::transpose(object)
  # Bind elements that are data frames.
  for(i in c("test","fit","fitsLR","best","joint")) {
    out[[i]] <-
      dplyr::bind_rows(
        object[[i]],
        .id = "mediator_id")
  }
  # Bind elements that are vectors; set to NULL if empty.
  for(i in c("normF","driver")) {
    isnt <- !sapply(object[[i]], is.null)
    out[[i]] <- {
      if(any(isnt))
        as.data.frame(t(as.data.frame(object[[i]][isnt])))
      else
        NULL
    }
  }

  out
}