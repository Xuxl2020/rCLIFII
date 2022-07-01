#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Bootstrap sampling
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'
#' Generate bootstrap samples from original observations by sampling individuals with replacement.
#'
#' @name lir_bootstrap
#' @title Nonparametric bootstrap sampling for animal movement data
#'
#' @param X A list or matrix containing identities of individuals identified per sampling period
#' @param tp A set of observation time
#' @param seed Random seed; the default is NULL
#' @details
#' For more details on this function, please see Efron and Tibshirani (1993).
#'
#' @return The bootstrap samples of animal movement data.
#'
#' @export
#' @rdname lir_bootstrap
#'
#' @references Efron, B., & Tibshirani, R. J. (1986). Bootstrap methods for standard errors,
#' confidence intervals, and other measures of statistical accuracy. Statistical science, 1(1), 54-75.
#'
#' @references Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. Chapman and Hall Press.
#'
#'
#' @examples
#' # Example
#' # load data
#' data(simulationA)
#' tp <- simulationA@tp
#' # if X is a list
#' list_simulation_A <- simulationA@list_simulation_A
#' bootstrap_sample <- lir_bootstrap(list_simulation_A, tp)
#' # if X is a matrix
#' matrix_simulation_A <- simulationA@matrix_simulation_A
#' bootstrap_sample <- lir_bootstrap(matrix_simulation_A, tp)
#'

lir_bootstrap <- function(X, tp, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Check arguments for errors
  if(!(class(X)[1] %in% c("list", "matrix"))) {
    stop("Argument 'X' must be a 'list' or 'matrix' object.")
  }

  if((class(X)[1] == "matrix") && (ncol(X) != length(tp))) {
    stop("Number of columns in 'X' is not equal to length from 'tp'")
  }

  if(is.list(X)){
    observed_individual <- unlist(X)
    matrix_data <- list_to_matrix(X, tp)
  }

  if(class(X)[1] == "matrix") {
    if(is.null(rownames(X))){
      # Set default individual names in a matrix to 1:p.
      observed_individual <- 1:nrow(X)
      rownames(X) <- paste('ID', observed_individual, sep='')
      colnames(X) <- tp
      matrix_data <- X
    } else {
    observed_individual <- 1:nrow(X)
    colnames(X) <- tp
    matrix_data <- X
    }
  }

  unique_observed_individual <- observed_individual[!duplicated(observed_individual)]
  len <- length(unique_observed_individual)
  bootstrap_sample <- sample(unique_observed_individual, len, replace = TRUE)
  # update sample
  matrix_bootstrap_sample <- matrix_data[unique(match(bootstrap_sample, unique_observed_individual)),]
  return(matrix_bootstrap_sample)
}

#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Block bootstrap sampling
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'
#' Observation time is divided into non-overlapping blocks; obtain one bootstrap subsample from
#' each block by sampling the periods with replacement; a complete bootstrap sample is obtained
#' by merging the bootstrap subsamples according to their order in the original time series.
#'
#'
#' @name lar_bootstrap
#' @title Nonparametric block bootstrap sampling for animal social structure data
#'
#' @param X A list or matrix containing identities of individuals and groups affiliated with them per sampling period
#' @param block_list A block list for a series of observation time
#' @param group_id Groups to which individuals belong; if X is a list, please input group_id;
#' if X is a matrix, this parameter can be skipped and and takes the default `NULL` value
#'
#' @param seed Random seed; the default is NULL
#'
#' @details
#' For more details on this function, please see the block bootstrap sampling for time series analysis.
#'
#' @return The bootstrap samples of animal social structure data.
#'
#' @export
#' @rdname lar_bootstrap
#'
#' @references Politis, D. N. (2003). The impact of bootstrap methods on time series analysis.
#' Statistical science, 18, 219-230.
#'
#'
#' @examples
#'
#' # Example
#' # load data
#' data(simulationD)
#' block_list <- simulationD@block_list
#' group_id <- simulationD@group.id
#' # if X is a list
#' list_simulation_D <- simulationD@list_simulation_D
#' bootstrap_sample <- lar_bootstrap(list_simulation_D, block_list, group_id)
#' # if X is a matrix
#' matrix_simulation_D <- simulationD@matrix_simulation_D
#' bootstrap_sample <- lar_bootstrap(matrix_simulation_D, block_list)
#'
#'
#'
lar_bootstrap <- function(X, block_list, group_id = NULL, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Check arguments for errors
  if(!(class(X)[1] %in% c("list", "matrix"))) {
    stop("Argument 'X' must be a 'list' or 'matrix' object.")
  }

  if(is.list(X) && is.null(group_id)) {
      stop('please input group_id')
  }

  if((class(X)[1] == "list") && (length(X) != length(group_id))){
    stop("length from 'X' is not equal to length from 'group_id'")
  }

  if((class(X)[1] == "matrix") && (ncol(X) != length(unlist(block_list)))) {
    stop("Number of columns in 'X' is not equal to observation time point")
  }

  if(is.list(X)){
    matrix_data <- list_to_matrix(X, unlist(block_list), group_id)
  }

  if(class(X)[1] == "matrix") {
    if(is.null(rownames(X))){
      # Set default individual names in a matrix to 1:p.
      observed_individual <- 1:nrow(X)
      rownames(X) <- paste('ID', observed_individual, sep='')
      colnames(X) <- unlist(block_list)
      matrix_data <- X
    } else {
      observed_individual <- 1:nrow(X)
      colnames(X) <- unlist(block_list)
      matrix_data <- X
    }
  }

  bootstrap_sample <- unlist(lapply(block_list, function(x)sort(sample(x, replace = TRUE))))
  idx <- match(bootstrap_sample, unlist(block_list))
  # update sample
  matrix_bootstrap_sample <- matrix_data[, idx]
  return(matrix_bootstrap_sample)
}




