
#' @name lir_simulation_A/B/C
#' @title Simulate three types of animal movement patterns
#'
#' @param Z A positive integer, representing the population size within a closed area, members of which can be either inside
#' or outside the study area.
#' @param N A positive integer, representing the population size within the study area.
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period.
#' It indicates the same number of individuals identified in all sampling periods if a positive integer.
#' @param tp A set of observed time.
#' @param lambda The rate of emigration from a study area.
#' @param mu The rate of reimmigration from a study area.
#' @import ggplot2
#'
#' @details
#'
#' Type A (Closed): a population of N individuals present in the study area throughout,
#' with no immigration, emigration, birth, death.
#'
#' Type B (Permanent emigration): a population of N individuals in the study area with
#' permanent emigration at a rate of lambda per individual per time unit,
#' with departed individuals being replaced 1:1 by new individuals.
#'
#' Type C (Emigration plus reimmigration): a closed population of Z individuals,
#' members of which can be either inside or outside the study area. Individuals
#' in the study area leave the study area at a rate of lambda per individual per
#' time unit and individuals outside the study area reenter it with a probability
#' of mu per individual per time unit.
#'
#' @return Simulation data for Type A-C.
#'
#'
#' @export
#' @rdname lir_simulations
#'
#' @examples
#'
#'
#' # Example
#' # set population size and number of subsampling
#' Z <- 300; N <- 100; n <- 40
#' # set the observation time
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # set parameters
#' lambda_B <- 0.008; lambda_C <- 0.08; mu_C <- 0.04
#'
#' # simulation of Type A
#' simulation_A <- lir_simulation_A(N, n, tp)
#' lir_data <- lir_nonparametric_estimation(simulation_A, n, tp)
#' mij <- lir_data$mij
#' nij <- lir_data$nij
#' tauij <- lir_data$tauij
#' theta <- round(Model1(mij, nij, tauij)$par, 8)
#' # a plot object
#' graph_data <- data.frame(tau = lir_data$tau, R_tau = lir_data$R_tau)
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, R_tau)) +
#' geom_point() +
#'   scale_y_continuous(limits = c(0, 0.02)) +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x = 'Time Lag', y = 'Lagged Identification Rates', title = 'Simulation_A') +
#'   geom_hline(yintercept = theta, linetype = 2, color = 'blue')
#'   geom_hline(yintercept = 0.01, linetype = 3, color = 'red')
#'
#' # simulation of Type B
#' simulation_B <- lir_simulation_B(N, n, lambda_B, tp)
#' lir_data <- lir_nonparametric_estimation(simulation_B, n, tp)
#' mij <- lir_data$mij
#' nij <- lir_data$nij
#' tauij <- lir_data$tauij
#' theta <- round(Model2(mij, nij, tauij)$par, 8)
#' alpha <- theta[1]; beta <- theta[2]
#' graph_data <- data.frame(tau = lir_data$tau, R_tau = lir_data$R_tau)
#' tT <- max(tp-min(tp))
#' line_data <- data.frame(tau = seq_len(tT), R_tau = alpha*exp(-beta*seq_len(tT)))
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, R_tau)) +
#'   geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x = 'Time Lag', y = 'Lagged Identification Rates', title = 'Simulation_B') +
#'   geom_line(aes(tau, R_tau), line_data, linetype = 2, color='blue')
#'
#' # simulation of Type C
#' simulation_C <- lir_simulation_C(Z, N, n, lambda_C, mu_C, tp)
#' lir_data <- lir_nonparametric_estimation(simulation_C, n, tp)
#' mij <- lir_data$mij
#' nij <- lir_data$nij
#' tauij <- lir_data$tauij
#' theta <- round(Model3(mij, nij, tauij)$par, 8)
#' alpha <- theta[3]; beta <- theta[2]; gamma <- theta[1]
#' graph_data <- data.frame(tau = lir_data$tau, R_tau = lir_data$R_tau)
#' tT <- max(tp-min(tp))
#' line_data <- data.frame(tau = seq_len(tT), R_tau = gamma*exp(-beta*seq_len(tT))+alpha)
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, R_tau)) +
#'   geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x = 'Time Lag', y = 'Lagged Identification Rates', title = 'Simulation_C') +
#'   geom_line(aes(tau, R_tau), line_data, linetype = 2, color='blue')
#'
#'
#'

lir_simulation_A <- function(N, n, tp) {

  tp <- tp-min(tp)+1
  tT <- max(tp)

  if (length(n)>1 && length(n)!=length(tp)) {
    stop("'n' or 'tp' is valid.")
  }

  if (length(n) == 1) {
    n0 <- rep(n, length(tp))
  }

  list_data <- list()
  for(i in tp){
    list_data[[i]] <- sample(1:N, n, replace = F)
  }

  matrix_data <- list_to_matrix(list_data, tp)
  simulation_A <- matrix_data
  return(simulation_A)
}


#' @rdname lir_simulations
#' @export
#'
lir_simulation_B <- function(N, n, lambda, tp) {

  tp <- tp-min(tp)+1

  if (length(n)>1 && length(n)!=length(tp)) {
    stop("'n' or 'tp' is valid.")
  }

  if (length(n) == 1) {
    n <- rep(n, length(tp))
  }

  tT <- max(tp)
  pop <- seq(N)
  Npop <- list()
  Npop[[1]] <- pop
  for(i in 2:tT){
    unif <- stats::runif(N, 0, 1)
    Npop[[i]] <- (unif<=lambda)*unif + Npop[[i-1]]
  }

  data <- list()
  for (i in tp) {
    data[[i]] <- sample(Npop[[i]], n[which(tp==i)])
  }

  matrix_data <- list_to_matrix(data, tp)
  simulation_B <- matrix_data
  return(simulation_B)
}

#' @rdname lir_simulations
#' @export
#'
#'
lir_simulation_C <- function(Z, N, n, lambda, mu, tp) {

  tp <- tp-min(tp)+1

  if (length(n)>1 && length(n)!=length(tp)) {
    stop("'n' or 'tp' is valid.")
  }

  if (length(n) == 1) {
    n <- rep(n, length(tp))
  }

  tT <- max(tp)
  pop <- seq(N)
  pop <- seq(N)
  Z1 <- 1:N
  Z2 <- (N+1):(Z)

  ### updata
  popin_new <- list()
  popin_new[[1]] <- Z1
  popout_new <- list()
  popout_new[[1]] <- Z2
  for(i in 2:tT) {
    unif1 <- stats::runif(length(popin_new[[i-1]]))
    unif2 <- stats::runif(length(popout_new[[i-1]]))
    popin_new[[i]] <- c(popin_new[[i-1]][which(unif1>lambda)],popout_new[[i-1]][which(unif2<=mu)])
    popout_new[[i]] <- c(popin_new[[i-1]][which(unif1<=lambda)],popout_new[[i-1]][which(unif2>mu)])
  }

  ### sample
  data <- list()
  for (i in tp) {
    data[[i]] <- sample(popin_new[[i]], n[which(tp==i)])
  }

  matrix_data <- list_to_matrix(data, tp)
  simulation_C <- matrix_data
  return(simulation_C)
}
