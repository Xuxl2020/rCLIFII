
#'
#' @name data.frame_to_matrix
#' @title Convert data frame form to matrix form
#'
#' @param obj A 'data.frame' object including ids, observation time, and groups (not required) of individuals
#'
#' @export
#' @rdname data.frame_to_matrix
#'
#'
#'
data.frame_to_matrix <- function(obj){

  # Check arguments for errors
  if(class(obj) != 'data.frame') {
    stop("Argument 'obj' must be a 'data.frame' object.")
  }

  input_data <- unique(obj)
  individual_ids <- unique(input_data$ID)
  input_data$Date <- input_data$Date - min(input_data$Date) + 1
  observed_time <- unique(input_data$Date)
  individual_number <- length(individual_ids)
  time_point_number <- length(observed_time)

  output_data <- matrix(0, individual_number, time_point_number)
  rownames(data) <- individual_ids
  colnames(data) <- observed_time

  if(ncol(input_data)==2){
    for (i in seq_len(time_point_number)){
      match_ID <- input_data$ID[which(input_data$Date %in% input_data$Date[i])]
      idx <- which(individual_ids %in% match_ID)
      output_data[idx, i] <- 1
    }
  }else{
    match_ID <- input_data$ID[which(input_data$Date %in% input_data$Date[i])]
    idx <- which(individual_ids %in% match_ID)
    output_data[idx, i] <- input_data$Group[idx]
  }
  return(output_data)
}

