
#'
#' @name list_to_matrix
#' @title Convert a list to matrix form
#'
#' @param obj A 'list' object including the ids of individuals
#' @param tp A set of observation time
#' @param group_id Groups individuals belong to; this parameter can be skipped and takes the default `NULL` value
#'
#' @export
#' @rdname list_to_matrix
#'
#'
list_to_matrix <- function(obj, tp, group_id = NULL){

  # Check arguments for errors
  if(!is.list(obj)) {
    stop("Argument 'obj' must be a 'list' object.")
  }

  if(is.list(obj) && is.null(tp)) {
    stop("please input 'tp'")
  }

  observed_time <- tp
  individual_ids <- unique(unlist(obj))
  individual_number <- length(individual_ids)
  time_point_number <- length(observed_time)
  data <- matrix(0, individual_number, time_point_number)
  rownames(data) <- individual_ids
  colnames(data) <- observed_time
  if(is.null(group_id)){
    for (i in seq_len(time_point_number)){
      idx <- rownames(data) %in% obj[[observed_time[i]]]
      data[idx, i] <- 1
    }
  }else{
    for (i in seq_len(time_point_number)){
      idx <- match(obj[[observed_time[i]]], individual_ids)
      data[idx, i] <- group_id[[observed_time[i]]][idx]
    }
  }
  data <- as.matrix(data)
  colnames(data) <- observed_time
  rownames(data) <- paste('ID', round(individual_ids, 0), sep='')
  return(data)
}

