
#'
#' @name lir_estimation_c
#' @title Estimation of variance inflation factor in QAIC under LIR models
#'
#' @param R_m \eqn{\sum_{i,j} \{m_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}
#' @param R_n \eqn{\sum_{i,j} \{n_i*n_j|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}
#' @param mij Number of individuals identified at both time t_i and t_j
#' @param nij ni*nj where ni is the number of individuals identified at time t_i
#' @param tauij Time lag between time t_i and t_j
#' @param mtau Maximum allowable lag time
#'
#' @export
#' @rdname lir_estimation_c
#'
#'
#'
lir_estimation_c <- function(R_m, R_n, mij, nij, tauij, mtau = 1000){
  # Estimation of c is c as calculated for the most general model (Model3 for LIR models).
  theta <- Model3(mij, nij, tauij, mtau)$par
  Rtau_C <- function(tauij){theta[1]*exp(-theta[2]*tauij) + theta[3]}
  unique_tau <- unique(tauij)
  real_Rm <- R_m[tauij][match(unique(tauij), tauij)]
  est_Rm <- R_n[tauij][match(unique(tauij), tauij)]*Rtau_C(unique_tau)

  # Categories were lumped so that all time lag categories contained an expected value of at least six.
  comb <- group_comb(real_Rm, est_Rm, unique_tau)
  comb_real_Rm <- comb$comb_real_val
  comb_est_Rm <- comb$comb_est_val

  df <- length(comb_real_Rm)
  c <- (comb_real_Rm-comb_est_Rm)^2 %*% matrix(1/comb_est_Rm)/(df-3-1)
  c <- ifelse(c>1 & c<3, c, 1)
}

#'
#'
#'
#' @name lar_estimation_c
#' @title Estimation of variance inflation factor in QAIC under LAR models
#'
#' @param g_m \eqn{\sum_{i,j} \{A_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}
#' @param g_n \eqn{\sum_{i,j} \{A_i|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}
#' @param Aij Number of observed associated individuals (in pairs) observed at both time t_i and t_j
#' @param Ai  Number of observed associated individuals (in pairs) at time t_i
#' @param tauij Time lag between time t_i and t_j
#' @param mtau Maximum allowable lag time
#'
#' @export
#' @rdname lar_estimation_c
#'
#'
lar_estimation_c <- function(g_m, g_n, Aij, Ai, tauij, mtau = 1000){
  # Estimation of c is c as calculated for the most general model (Model6 for LAR models).
  theta <- Model6(Aij, Ai, tauij, mtau)$par
  g_tau <- function(tauij){((1-theta[1])*exp(-theta[2]*tauij)+theta[1])*exp(-theta[3]*tauij)}
  unique_tau <- unique(tauij)
  real_gm <- g_m[tauij][match(unique(tauij), tauij)]
  est_gm <- g_n[tauij][match(unique(tauij), tauij)]*g_tau(unique_tau)
  # Categories were lumped so that all time lag categories contained an expected value of at least six.
  comb <- group_comb(real_gm, est_gm, unique_tau)
  comb_real_gm <- comb$comb_real_val
  comb_est_gm <- comb$comb_est_val

  df <- length(comb_real_gm)
  c <- (comb_real_gm-comb_est_gm)^2 %*% matrix(1/comb_est_gm)/(df-3-1)
  c <- ifelse(c>1 & c<3, c, 1)
}


#'
#'
#' Categories were lumped so that all time lag categories contained an expected value of at least six.
#'
#' @name group_comb
#' @title Summarize categories
#'
#' @param real_val A series of real values
#' @param est_val A series of estimates
#' @param tauij Time lag between time t_i and t_j
#'
#' @export
#' @rdname group_comb
#'
#'
group_comb <- function(real_val, est_val, tauij){
  unique_tau <- unique(tauij)
  num_est_val <- length(est_val)
  g <- rep(NA, num_est_val)
  idx <- seq_len(num_est_val)
  a <- est_val
  for(i in idx){
    b <- idx[which(cumsum(a) > 6)[1]]
    g[idx[which(idx <= b)]] <- b
    a <- a[-c(which(idx <= b))]
    idx <- idx[-c(which(idx <= b))]
  }

  g[(1:num_est_val)[-which(g>0)]] <- range(g[which(g>0)])[2]
  group_data <- data.frame(cbind(real_val, est_val, unique_tau, g))

  new_data <- group_data %>% dplyr::group_by(g)
  new <- new_data %>% dplyr::summarise(
    comb_real_val = sum(real_val),
    comb_est_val = sum(est_val)
  )
  return(new)
}
