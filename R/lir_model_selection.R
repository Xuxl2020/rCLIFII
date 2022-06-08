
#' @name lir_model_selection
#' @title Model selection for models of lagged identification rate
#'
#' @param X A list or matrix containing the identities of individuals identified in each sampling period
#' @param model Models of lagged identification rate, model = 'lir_1', 'lir_2', or 'lir_3'.
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period.
#' It indicates the same number of individuals identified in all sampling periods if a positive integer.
#' @param tp A set of observed time.
#' @param nboot The number of bootstrap samples desired
#' @param mtau The maximum allowable lag time.
#' @param ncores doParallel.
#' @param seed Random seed.
#' @details
#' See Akaike (1973) for Akaike information criterion (AIC);
#' See Burnhan and Anderson (2002) for Quasi-Akaike information criterion (QAIC);
#' See this paper for composite likelihood information criterion (CLIC).
#'
#' @return The values of model selection criteria.
#'
#'
#' @export
#' @rdname lir_model_selection


lir_model_selection <- function(X, n, tp, model, nboot, mtau = 1000, ncores = 4, seed = NULL){

  tp <- tp-min(tp)+1

  if (length(n)>1 && length(n)!=length(tp)) {
    stop("'n' or 'tp' is valid.")
  }

  if (length(n) == 1) {
    n <- rep(n, length(tp))
  }

  B <- as.integer(nboot)
  if (B <= 1){
    stop("nboot must be a positive integer bigger than 1")
  }

  tT <- max(tp-min(tp))
  len <- length(tp)

  lir_data <- lir_nonparametric_estimation(X, n, tp)

  R_m <- lir_data$R_m
  R_n <- lir_data$R_n
  mij <- lir_data$mij
  nij <- lir_data$nij
  tauij <- lir_data$tauij

  if(model == 'lir_1'){
    mod0 <- 'Model1'
  }
  if(model == 'lir_2'){
    mod0 <- 'Model2'
  }else{
    mod0 <- 'Model3'
  }

  model.H <- lir.model.res(mod0, mij, nij, tauij, mtau)$H
  model.est <- lir.model.res(mod0, mij, nij, tauij, mtau)$par
  model.val <- lir.model.res(mod0, mij, nij, tauij, mtau)$val
  model.K <- lir.model.res(mod0, mij, nij, tauij, mtau)$K

  if(ncores>1){
    cl <- parallel::makeCluster(ncores) # not to overload your computer
    doParallel::registerDoParallel(cl)
    RESULTS = foreach::`%dopar%`(foreach::foreach(i = seq_len(B), .combine = rbind), {
      sampboot <- lir_bootstrap(X, tp, seed)
      dat <- lir_nonparametric_estimation(sampboot, n, tp)
      mij <- dat$mij
      nij <- dat$nij
      tauij <- dat$tauij
      res <- lir.model.res(model=mod0, mij=mij, nij=nij, tauij=tauij, mtau)
      out <- res$par
      return(out)
    })
    parallel::stopCluster(cl)
  }else{
    RESULTS <- matrix(0, B, model.K)
    for(i in seq_len(B)){
      sampboot <- lir_bootstrap(X, tp, seed)
      dat <- lir_nonparametric_estimation(sampboot, n, tp)
      mij <- dat$mij
      nij <- dat$nij
      tauij <- dat$tauij
      res <- lir.model.res(model=mod0, mij=mij, nij=nij, tauij=tauij, mtau)
      out <- res$par
      RESULTS[i,] <- out
      }
  }

    dimpara <- sum(diag(model.H%*%(stats::var(RESULTS))))
    AIC <- 2*model.val + 2*model.K
    estimation_c <- lir_estimation_c(R_m, R_n, mij, nij, tauij, mtau)
    QAIC <- 2*model.val/estimation_c + 2*model.K
    CLIC <- 2*model.val + dimpara
    res <- data.frame(AIC=AIC, QAIC=QAIC, CLIC=CLIC)
    return(res)
}


