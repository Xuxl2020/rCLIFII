#' @name models_internal
#' @title Build the models
#'
#' @param mij The number of individuals identified at both time t_i and t_j
#' @param nij The number of individuals identified at time t_i
#' @param Aij The number of observed associated individuals (in pairs) observed at both time t_i and t_j
#' @param Ai The number of observed associated individuals (in pairs) at time t_i
#' @param tauij Time lag between time t_i and t_j
#' @param mtau The maximum allowable lag time.
#' @param model Models of LIR and LAR, or model is formulated by yourself
#' @param model_cl_fun If you formulate your model, please input function to calculate the composite likelihood about your model
#' @param cl.H If you formulate your model, please input the sensitivity matrix with respect to parameters in your model
#' @param model.K If you formulate your model, please input the number of parameters in your model
#'
#' @export
#' @rdname models_internal
#'
Model1 <- function(mij, nij, tauij, mtau=1000) {
  ind = which(tauij<=mtau)
  mij <- mij[ind]
  nij <- nij[ind]
  LL <- function(theta){
    Rtau <- theta
    if (min(Rtau)<=0) return(9000000000000);
    if (max(Rtau)>=1) return(9000000000000);
    object <- sum(mij*log(Rtau)+(nij-mij)*log(1-Rtau))
    return(-object)
  }
  opt <- stats::optim(0.001, LL, method ="Brent", lower = 0.001, upper = 0.09)
  return(opt)
}

#' @export
#' @rdname models_internal
#'
Model2 <- function(mij, nij, tauij, mtau=1000) {
  ind = which(tauij<=mtau)
  mij <- mij[ind]
  nij <- nij[ind]
  tauij <- tauij[ind]
  LL <- function(theta){
    Rtau <- theta[1]*exp(-theta[2]*tauij)
    if (min(Rtau)<=0) return(9000000000000);
    if (max(Rtau)>=1) return(9000000000000);
    object <- sum(mij*log(Rtau)+(nij-mij)*log(1-Rtau))
    return(-object)
  }
  opt <- stats::optim(c(0.001, 0.001), LL)
  return(opt)
}

#' @export
#' @rdname models_internal
#'
Model3 <- function(mij, nij, tauij, mtau=1000) {
  ind = which(tauij<=mtau)
  mij <- mij[ind]
  nij <- nij[ind]
  tauij <- tauij[ind]
  LL <- function(theta){
    Rtau <- theta[1]*exp(-theta[2]*tauij) + theta[3]
    if (min(Rtau)<=0) return(9000000000000);
    if (max(Rtau)>=1) return(9000000000000);
    object <- sum(mij*log(Rtau)+(nij-mij)*log(1-Rtau))
    return(-object)
  }
  opt <- stats::optim(c(0.001, 0.001, 0.001), LL)
  return(opt)
}

#' @export
#' @rdname models_internal
#'
#'
Model4 <- function(Aij, Ai, tauij, mtau=1000){
  ind = which(tauij<=mtau)
  Aij <- Aij[ind]
  Ai <- Ai[ind]
  LL <- function(theta){
    gtau <- theta
    if (min(gtau)<=0) return(9000000000000);
    if (max(gtau)>=1) return(9000000000000);
    return(-sum(Aij*log(gtau) + (Ai-Aij)*log(1-gtau)))
  }
  opt <- stats::optim(0.01, LL, method ="Brent", lower = 0.05, upper =0.9)
  return(opt)
}

#' @export
#' @rdname models_internal
#'
Model5 <- function(Aij, Ai, tauij, mtau=1000) {
  ind = which(tauij<=mtau)
  Aij <- Aij[ind]
  Ai <- Ai[ind]
  tauij <- tauij[ind]
  LL <- function(para){
      gtau <- (1 - para[1])*exp(-para[2]*tauij) + para[1]
      if (min(gtau)<=0) return(9000000000000);
      if (max(gtau)>=1) return(9000000000000);
      return(-sum(Aij*log(gtau) + (Ai-Aij)*log(1-gtau)))
    }
  opt <- stats::optim(c(0.01, 0.01), LL)
  return(opt)
}

#' @export
#' @rdname models_internal
#'
Model6 <- function(Aij, Ai, tauij, mtau=1000) {
  ind = which(tauij<=mtau)
  Aij <- Aij[ind]
  Ai <- Ai[ind]
  tauij <- tauij[ind]
  LL <- function(para){
    gtau <- ((1-para[1])*exp(-para[2]*tauij) + para[1])*exp(-para[3]*tauij)
    if (min(gtau)<=0) return(9000000000000);
    if (max(gtau)>=1) return(9000000000000);
    return(-sum(Aij*log(gtau) + (Ai-Aij)*log(1-gtau)))
  }
  opt <- stats::optim(c(0.01, 0.01, 0.01), LL)
  return(opt)
}

#' @export
#' @rdname models_internal
#'
lir.model.res <- function(model, mij, nij, tauij, mtau, model_cl_fun = NULL, cl.H = NULL, model.K = NULL){
  if(model=='Model1') {
    est <- Model1(mij, nij, tauij, mtau)
    sder <- sum(-mij/est$par^2-(nij-mij)/((1-est$par)^2))
    H <- matrix(sder)
    K = 1
  }
  if(model=='Model2') {
    est <- Model2(mij, nij, tauij, mtau)
    alpha <- est$par[1]
    beta <- est$par[2]
    Rtau <- alpha*exp(-beta*tauij)
    Hb11 <- sum(exp(-beta*tauij)^2*(mij/Rtau^2+(nij-mij)/(1-Rtau)^2))
    Hb22 <- sum(alpha*exp(-beta*tauij)*tauij^2*(-mij/Rtau+(nij-mij)/(1-Rtau)) +
                  (alpha*exp(-beta*tauij)*tauij)^2*(mij/Rtau^2+(nij-mij)/(1-Rtau)^2))
    Hb12 <- sum(-exp(-beta*tauij)*tauij*(mij/Rtau-(nij-mij)/(1-Rtau)) -
                  tauij*alpha*exp(-beta*tauij)^2*(mij/Rtau^2+(nij-mij)/(1-Rtau)^2))
    H <- matrix(c(Hb11, Hb12, Hb12, Hb22), 2, 2)
    K = 2
  }
  if (model=='Model3') {
    est <- Model3(mij, nij, tauij, mtau)
    alpha <- est$par[1]
    beta <- est$par[2]
    gamma <- est$par[3]
    Rtau <- alpha*exp(-beta*tauij) + gamma
    Hc11 <- sum(exp(-beta*tauij)^2*( mij/Rtau^2+(nij-mij)/(1-Rtau)^2 ))
    Hc22 <- sum(alpha*exp(-beta*tauij)*tauij^2*(-mij/Rtau+(nij-mij)/(1-Rtau)) +
                  (alpha*exp(-beta*tauij)*tauij)^2*(mij/Rtau^2+(nij-mij)/(1-Rtau)^2) )
    Hc33 <- sum(mij/Rtau^2+(nij-mij)/(1-Rtau)^2)
    Hc12 <- sum(-exp(-beta*tauij)*tauij*(mij/Rtau-(nij-mij)/(1-Rtau)) -
                  tauij*alpha*exp(-beta*tauij)^2*(mij/Rtau^2+(nij-mij)/(1-Rtau)^2))
    Hc13 <- sum(exp(-beta*tauij)*( mij/Rtau^2+(nij-mij)/(1-Rtau)^2 ))
    Hc23 <- sum(-alpha*exp(-beta*tauij)*tauij*( mij/Rtau^2+(nij-mij)/(1-Rtau)^2))
    H <- matrix(c(Hc11,Hc12,Hc13,Hc12,Hc22,Hc23,Hc13,Hc23,Hc33), 3, 3)
    K = 3
  }
  if (model == 'model_cl_fun'){
    est <- model_cl_fun(mij, nij, tauij, mtau)
    H <- cl.H
    K <- model.K
  }
  par <- est$par
  val <- est$value
  res <- list(H=H, K=K, par=par, val=val)
  return(res)
}

#' @export
#' @rdname models_internal
#'
lar.model.res <- function(model, Aij, Ai, tauij, mtau, model_cl_fun = NULL, cl.H = NULL, model.K = NULL){
  if(model=='Model4') {
    est <- Model4(Aij, Ai, tauij, mtau)
    sder <- sum(-Aij/est$par^2-(Ai-Aij)/((1-est$par)^2))
    H <- matrix(sder)
    K = 1
  }
  if(model=='Model5') {
    est <- Model5(Aij, Ai, tauij, mtau)
    alpha <- est$par[1]
    beta <- est$par[2]
    Rtau <- (1-alpha)*exp(-beta*tauij)+alpha
    He11 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(1-exp(-beta*tauij))^2)
    He12 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(1-exp(-beta*tauij))*(-tauij*(1-alpha)*exp(-beta*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*tauij*exp(-beta*tauij))
    He22 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(tauij^2*(1-alpha)^2*exp(-2*beta*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*tauij^2*(1-alpha)*exp(-beta*tauij))
    H <- matrix(c(He11,He12,He12,He22),2,2)
    K = 2
  }
  if (model=='Model6') {
    est <- Model6(Aij, Ai, tauij, mtau)
    alpha <- est$par[1]
    beta <- est$par[2]
    gamma <- est$par[3]
    Rtau <- ((1-alpha)*exp(-beta*tauij)+alpha)*exp(-gamma*tauij)
    Hf11 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(1-exp(-beta*tauij))^2*exp(-2*gamma*tauij))
    Hf12 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(-tauij*(1-alpha)*(1-exp(-beta*tauij))*exp(-2*gamma*tauij-beta*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*(tauij*exp(-gamma*tauij-beta*tauij)))
    Hf13 <- -sum( (-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(-tauij*(1-exp(-beta*tauij))*((1-alpha)*exp(-beta*tauij)+alpha)*exp(-2*gamma*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*(-tauij*exp(-gamma*tauij)*(1-exp(-beta*tauij))))
    Hf22 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(tauij^2*(1-alpha)^2*exp(-2*gamma*tauij-2*beta*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*(tauij^2*(1-alpha)*exp(-gamma*tauij-beta*tauij)))
    Hf23 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(tauij^2*(1-alpha)*((1-alpha)*exp(-beta*tauij)+alpha)*exp(-2*gamma*tauij-beta*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*(tauij^2*(1-alpha)*exp(-beta*tauij-gamma*tauij)))
    Hf33 <- -sum((-Aij/Rtau^2-(Ai-Aij)/(1-Rtau)^2)*(tauij^2*((1-alpha)*exp(-beta*tauij)+alpha)^2*exp(-2*gamma*tauij))+(Aij/Rtau-(Ai-Aij)/(1-Rtau))*(tauij^2*((1-alpha)*exp(-beta*tauij)+alpha)*exp(-gamma*tauij)))
    H <- matrix(c(Hf11,Hf12,Hf13,Hf12,Hf22,Hf23,Hf13,Hf23,Hf33),3,3)
    K = 3
  }

  if (model == 'model_cl_fun'){
    est <- model_cl_fun(Aij, Ai, tauij, mtau)
    H <- cl.H
    K <- model.K
  }
  par <- est$par
  val <- est$value
  res <- list(H=H, K=K, par=par, val=val)
  return(res)
}




