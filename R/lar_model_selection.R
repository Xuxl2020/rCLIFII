


#' @name lar_model_selection
#' @title Model selection for models of lagged association rate
#' @param X A list or matrix containing the identities of individuals within
#' study area, and the states or status of individuals during each sampling period.
#' @param model Models of lagged identification rate, model = 'lar_1', 'lar_2', or 'lar_3'.
#' @param block_list A block list for a series of observation time. For example,
#' block_list = list(c(1:5), c(51:55), c(101:105), c(501:505), c(601:605)).
#' @param group_id Groups of individuals. If X is a list, please input group_id. If X is a matrix,
#' this parameter can be skipped and takes the default `NULL` value.
#' @param nboot The number of bootstrap samples desired
#' @param mtau The maximum allowable lag time
#' @param ncores doParallel
#' @param seed Random seed
#' @details
#' See Akaike (1973) for Akaike information criterion (AIC);
#' See Burnhan and Anderson (2002) for Quasi-Akaike information criterion (QAIC);
#' See this paper for composite likelihood information criterion (CLIC).
#'
#' @return The values of model selection criteria.
#'
#'
#' @export
#' @rdname lar_model_selection


lar_model_selection <- function(X, model, block_list, nboot, group_id = NULL, mtau = 1000, ncores = 4, seed = NULL){
  B <- as.integer(nboot)
  if (B <= 1){
    stop("nboot must be a positive integer bigger than 1")
  }
  tp <- unlist(block_list)
  tT <- max(tp-min(tp))
  len <- length(tp)
  lar_data <- lar_nonparametric_estimation(X, tp, group_id)
  g_m <- lar_data$g_m
  g_n <- lar_data$g_n
  Aij <- lar_data$Aij
  Ai <- lar_data$Ai
  tauij <- lar_data$tauij
  if(model == 'lar_1'){
    mod0 <- 'Model4'
  }
  if(model == 'lar_2'){
    mod0 <- 'Model5'
  }
  if(model == 'lar_3'){
    mod0 <- 'Model6'
  }
  model.H <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$H
  model.est <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$par
  model.val <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$val
  model.K <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$K
  if(ncores>1){
    cl <- parallel::makeCluster(ncores) # not to overload your computer
    doParallel::registerDoParallel(cl)
    RESULTS = foreach::`%dopar%`(foreach::foreach(i = seq_len(B), .combine = rbind), {
      sampboot <- lar_bootstrap(X, block_list, group_id, seed)
      dat <- lar_nonparametric_estimation(sampboot, tp, group_id)
      Aij <- dat$Aij
      Ai <- dat$Ai
      tauij <- dat$tauij
      res <- lar.model.res(model=mod0, Aij=Aij, Ai=Ai, tauij=tauij, mtau)
      out <- res$par
      return(out)
    })
    parallel::stopCluster(cl)
  }else{
    RESULTS <- matrix(0, B, model.K)
    for(i in seq_len(B)){
      sampboot <- lar_bootstrap(X, block_list, group_id, seed)
      dat <- lar_nonparametric_estimation(sampboot, tp, group_id)
      Aij <- dat$Aij
      Ai <- dat$Ai
      tauij <- dat$tauij
      res <- lar.model.res(model=mod0, Aij=Aij, Ai=Ai, tauij=tauij, mtau)
      out <- res$par
      RESULTS[i,] <- out
    }
  }
  dimpara <- sum(diag(model.H%*%(stats::var(RESULTS))))
  AIC <- 2*model.val + 2*model.K
  estimation_c <- lar_estimation_c(g_m, g_n, Aij, Ai, tauij, mtau)
  QAIC <- 2*model.val/estimation_c + 2*model.K
  CLIC <- 2*model.val + dimpara
  res <- data.frame(AIC=AIC, QAIC=QAIC, CLIC=CLIC)
  return(res)
}

















