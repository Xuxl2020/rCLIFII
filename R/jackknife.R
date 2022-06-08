#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# jackknife resampling
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'
#' jackknife with sampling periods within k-time-unit intervals being omitted
#' in turn, giving T/k jackknife replicates, with T is the number of sampling periods;
#'
#' @name jackknife
#' @title Nonparametric jackknife resampling for animal movement (social structure) data
#'
#' @param X A list or matrix containing animal movement or social structure data
#' @param tp A set of observation time
#' @param k An integer represents k-time-unit intervals
#'
#' @return The jackknife samples of animal movement or social structure data.
#'
#' @export
#' @rdname jackknife
#'
#' @references Quenouille, M. (1949). Approximate tests of correlation in time series. Journal of the
#' Royal Statistical Society, Soc. Ser. B, 11, 18-84.
#'
#' @references Tukey, J. W. (1958). Bias and confidence in not quite large samples (abstract), Annals of
#' Mathematical Statistics, 29, 614.
#'
#' @references Efron, B. (1979). Bootstrap methods: Another look at the jackknife. Annals of Statistics, 7, 1-26.
#'
#'
#' @examples
#' # Example
#' # load data
#' data(simulationA)
#' tp <- simulationA@tp
#' # if X is a list
#' list_simulation_A <- simulationA@list_simulation_A
#' jackknife_sample <- jackknife(list_simulation_A, tp, k=50)
#' # if X is a matrix
#' matrix_simulation_A <- simulationA@matrix_simulation_A
#' jackknife_sample <- jackknife(matrix_simulation_A, tp, k=50)
#'

jackknife <- function(X, tp, k) {

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

  # update sample
  bin <- bin_make(tp, k)
  jackknife_sample <- list()

  for (i in 1:length(bin)){
    idx <- match(bin[[i]], tp)
    jackknife_sample[[i]] <- matrix_data[,-idx]
  }

  return(jackknife_sample)

}



