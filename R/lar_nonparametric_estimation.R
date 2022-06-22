

#' @title A nonparametric estimator for lagged association rate
#'
#' @param X A list or matrix containing the identities of individuals within
#' study area, and the states or status of individuals during each sampling period.
#' @param tp A set of observation time.
#' @param group_id Groups of individuals. If X is a list, please input group_id. If X is a matrix,
#' this parameter can be skipped and takes the default `NULL` value.
#'
#' @return A list with the following elements:
#' \item{tau}{a set of time lags}
#' \item{gtau}{A nonparametric estimator of \eqn{\hat{g}(\tau)} was given by Whitehead (2007)}
#' \item{g_m}{\eqn{\sum_{i,j} \{A_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|} }
#' \item{g_n}{\eqn{\sum_{i,j} \{A_i|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}}
#' \item{Aij}{number of associations of pair of individuals at both time \eqn{t_i} and \eqn{t_j}}
#' \item{Ai}{number of associations of pair of individuals at time \eqn{t_i}}
#' \item{tauij}{\eqn{\tau_{ij}=|t_i-t_j|}}
#'
#'
#' @details
#' The lagged association rate \eqn{g(\tau)} is the probability that if two individuals
#' are associated now, they will still be associated after a time lag of  \eqn{\tau}.
#' A nonparametric estimator of \eqn{g(\tau)} was given by Whitehead (2007):
#' \deqn{\hat{g}(\tau)=\frac{\sum_{i,j} \{A_{t_i,t_j}|\tau_{ij}=\tau\}}{\sum_{i,j}\{A_i|\tau_{ij}=\tau\}},}
#' where where \eqn{A_{t_i,t_j}} is the number of associations of pair of individuals at both
#' time \eqn{t_i} and \eqn{t_j}, and \eqn{A_{t_i}} represents the number of associations
#' of pair of individuals at time \eqn{t_i}.
#'
#'
#' @export
#' @rdname lar_nonparametric_estimation
#'
lar_nonparametric_estimation <- function(X, tp, group_id = NULL) {

  tp <- tp-min(tp)+1

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

  if((class(X)[1] == "matrix") && (ncol(X) != length(tp))) {
    stop("Number of columns in 'X' is not equal to observation time point")
  }

  if(is.list(X)){
    data <- list_to_matrix(X, tp, group_id)
  }

  if(class(X)[1] == "matrix") {
    data <- X
  }

  # tT <- max(tp)
  # g_m <- rep(0, tT)
  # g_n <- rep(0, tT)
  # Aij <- c()
  # Ai <- c()
  # tauij <- c()
  # p <- 1
  # for (i in 1:(length(tp)-1)) {
  #   for (j in (i+1):length(tp)){
  #     tauij[p] <- tp[j] - tp[i]
  #     num <- data[, i]*data[, j]*(data[, i]+2*data[, j])*(data[, i]+3*data[, j])
  #     Aij[p] <- sum(choose(table(num[which(num!=0)]), 2))
  #     if(sum(num!=0)!=0){
  #       obj <- data[which(num!=0), i]
  #       s <- 0
  #       for(g in 1:length(obj)){
  #         s <- s + sum(data[which(num!=0)[g]:dim(data)[1],i]==obj[g])
  #       }
  #       Ai[p] <- s - length(obj)
  #     }else{
  #       Ai[p] <- 0
  #     }
  #     g_m[tp[j]-tp[i]] <- g_m[tp[j]-tp[i]] + Aij[p]
  #     g_n[tp[j]-tp[i]] <- g_n[tp[j]-tp[i]] + Ai[p]
  #     p <- p + 1
  #   }
  # }

  tT <- max(tp)
  g_m <- rep(0, tT)
  g_n <- rep(0, tT)
  Aij <- c()
  Ai <- c()
  tauij <- c()

  p <- 1
  for (i in 1:(length(tp)-1)) {
    for (j in (i+1):length(tp)){
      tauij[p] <- tp[j] - tp[i]
      num <- data[, i]*data[, j]*(data[, i]+2*data[, j])*(data[, i]+3*data[, j])
      num <- data[, i]*data[, j]
      Aij[p] <- 2*sum(choose(table(num[which(num!=0)]), 2))
      Ai[p] <- sum(choose(table(data[which(data[, i]!=0), i]), 2))
      g_m[tp[j]-tp[i]] <- g_m[tp[j]-tp[i]] + Aij[p]
      g_n[tp[j]-tp[i]] <- g_n[tp[j]-tp[i]] + Ai[p]
      p <- p + 1
    }
  }

  gtauij <- g_m[tauij]/g_n[tauij]
  tau <- unique(tauij)
  g_tau <- g_m[tau]/g_n[tau]

  g_data <- list(tau=tau, g_tau=g_tau, g_m=g_m, g_n=g_n,
                Aij=Aij, Ai=Ai, tauij=tauij)

  return(g_data)
}


