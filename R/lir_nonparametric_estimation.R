
#' @title A nonparametric estimator for lagged identification rate
#'
#' @param X A list or matrix containing the identities of individuals identified during each sampling period
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period.
#' It indicates the same number of individuals identified in all sampling periods if a positive integer.
#' @param tp A set of observation time
#' @param mtau The maximum allowable lag time
#'
#' @return A list with the following elements:
#'
#' \item{tau}{a set of time lags}
#' \item{R_tau}{A nonparametric estimator of \eqn{\hat{R}(\tau)} was given by Whitehead (2007)}
#' \item{R_m}{\eqn{\sum_{i,j} \{m_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|} }
#' \item{R_n}{\eqn{\sum_{i,j} \{n_i*n_j|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}}
#' \item{mij}{number of individuals identified at both time \eqn{t_i} and \eqn{t_j}}
#' \item{nij}{\eqn{n_i}*\eqn{n_j}}
#' \item{tauij}{\eqn{\tau_{ij}=|t_i-t_j|}}
#'
#'
#' @details
#' The lagged identification rate \eqn{R(\tau)} is the probability that an individual
#' in the study area is identified now and is reidentified after a time lag of \eqn{\tau}.
#' A nonparametric estimator of \eqn{R(\tau)} was given by Whitehead (2007):
#' \deqn{\hat{R}(\tau)=\frac{\sum_{i,j} \{m_{t_i,t_j}|\tau_{ij}=\tau\}}{\sum_{i,j}\{n_i*n_j|\tau_{ij}=\tau\}},}
#' where where \eqn{m_{t_i,t_j}} is the number of individuals identified at both
#' time \eqn{t_i} and \eqn{t_j}, and \eqn{n_{t_i}} represents the number of individuals
#' identified at time \eqn{t_i}.
#'
#'
#' @export
#' @rdname lir_nonparametric_estimation
#'
#' @examples
#' # Example
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' N <- 100; n <- 10
#' X <- list()
#' for (i in tp){
#' X[[i]] <- sample(1:N, n)
#' }
#'
#' res <- lir_nonparametric_estimation(X, n, tp)
#' res$R_tau; res$tau
#'
#'
lir_nonparametric_estimation <- function(X, n, tp, mtau=1000) {

  tp <- tp-min(tp)+1

  # Check arguments for errors
  if(!(class(X)[1] %in% c("list", "matrix"))) {
    stop("Argument 'X' must be a 'list' or 'matrix' object.")
  }

  if((class(X)[1] == "matrix") && (ncol(X) != length(tp))) {
    stop("Number of columns in 'X' is not equal to length from 'tp'")
  }

  if (length(n)>1 && length(n)!=length(tp)) {
    stop("'n' or 'tp' is valid.")
  }

  if(is.list(X)){
    data <- list_to_matrix(X, tp)
  }

  if(class(X)[1] == "matrix") {
    data <- X
  }

  if (length(n) == 1) {
    n <- rep(n, length(tp))
  }

  tT <- max(tp)
  R_m <- rep(0, tT)
  R_n <- rep(0, tT)
  mij <- c()
  nij <- c()
  tauij <- c()

  k <- 1
  for (i in 1:(length(tp)-1)){
    for (j in (i+1):length(tp)){
      if(tp[j] - tp[i]<=mtau){
      mij[k] <- sum(data[,j]*data[,i]==1)
      nij[k] <- n[i]*n[j]
      tauij[k] <- tp[j] - tp[i]
      R_m[tauij[k]] <- R_m[tauij[k]] + mij[k]
      R_n[tauij[k]] <- n[i]*n[j] + R_n[tauij[k]]
      k <- k + 1
    }
  }
}
  R_tauij <- R_m[tauij]/R_n[tauij]
  tau <- tauij
  R_tau <- R_m[tau]/R_n[tau]
  R_data <- list(R_tau=R_tau, tau=tau, R_m=R_m, R_n=R_n,
                 mij=mij, nij=nij, tauij=tauij)
  return(R_data)
}

