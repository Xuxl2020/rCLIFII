


#' @name lar_model_selection
#' @title Model selection for models of lagged association rate
#' @param X A list or matrix containing the identities of individuals within
#' study area, and the states or status of individuals during each sampling period
#' @param model Models of lagged identification rate, model = 'lar_1', 'lar_2', 'lar_3', or your model 'model_cl_fun'
#' @param tp A set of observed time
#' @param group_id Groups of individuals. If X is a list, please input group_id. If X is a matrix,
#' this parameter can be skipped and takes the default `NULL` value
#' @param model_cl_fun If you formulate your model, please input function to calculate the composite likelihood about your model
#' @param cl.H If you formulate your model, please input the sensitivity matrix with respect to parameters in your model
#' @param model.K If you formulate your model, please input the number of parameters in your model
#' @param method The method = 'Bootstrap', 'BBootstrap', or 'Jackknife'
#' @param bin_len An integer represents len-time-unit intervals
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


lar_model_selection <- function(X, model, method, tp,
                                mtau = 1000,
                                ncores = 4,
                                nboot = -1,
                                bin_len = -1,
                                group_id = NULL,
                                model_cl_fun = NULL,
                                cl.H = NULL,
                                model.K = NULL,
                                seed = NULL){
  tp <- unlist(tp)

  if (method == 'BBootstrap'){
    B <- as.integer(nboot)
    if (B <= 1){
      stop("nboot must be a positive integer bigger than 1")
    }
  }
  if (method == 'Jackknife'){
    num_bin <- as.integer(nboot)
    if (num_bin <= 1){
      stop("num_bin must be a positive integer bigger than 1")
    }
  }

  tT <- max(tp-min(tp))
  len <- length(tp)
  lar_data <- lar_nonparametric_estimation(X, tp, group_id)
  g_m <- lar_data$g_m
  g_n <- lar_data$g_n
  Aij <- lar_data$Aij
  Ai <- lar_data$Ai
  tauij <- lar_data$tauij

  if(model == 'lar_1'){
    model = 'Model4'
    model.H <- lar.model.res('Model4', Aij, Ai, tauij, mtau)$H
    model.est <- lar.model.res('Model4', Aij, Ai, tauij, mtau)$par
    model.val <- lar.model.res('Model4', Aij, Ai, tauij, mtau)$val
    model.K <- lar.model.res('Model4', Aij, Ai, tauij, mtau)$K
  }
  if(model == 'lar_2'){
    model = 'Model5'
    model.H <- lar.model.res('Model5', Aij, Ai, tauij, mtau)$H
    model.est <- lar.model.res('Model5', Aij, Ai, tauij, mtau)$par
    model.val <- lar.model.res('Model5', Aij, Ai, tauij, mtau)$val
    model.K <- lar.model.res('Model5', Aij, Ai, tauij, mtau)$K
  }
  if(model == 'lar_3'){
    model = 'Model6'
    model.H <- lar.model.res('Model6', Aij, Ai, tauij, mtau)$H
    model.est <- lar.model.res('Model6', Aij, Ai, tauij, mtau)$par
    model.val <- lar.model.res('Model6', Aij, Ai, tauij, mtau)$val
    model.K <- lar.model.res('Model6', Aij, Ai, tauij, mtau)$K
  }
  if (model == 'model_cl_fun') {
    model.H <- cl.H
    model.est <- model_cl_fun(Aij, Ai, tauij, mtau)$par
    model.val <- model_cl_fun(Aij, Ai, tauij, mtau)$val
    model.K <- model.K
  }

  if (method == 'Jackknife'){
    if (bin_len > 0){
      jsamples <- jackknife(X, tp, bin_len)
    }

    RESULTS <- matrix(0, length(jsamples), length(model.est))
    for(j in 1:length(jsamples)){
      tp0 = as.numeric(colnames(jsamples[[j]]))
      lar_data <- lar_nonparametric_estimation(as.matrix(jsamples[[j]]), tp0)
      Aij <- lar_data$Aij
      Ai <- lar_data$Ai
      tauij <- lar_data$tauij
      if (model == 'model_cl_fun'){
        RESULTS[j,] <- model_cl_fun(Aij, Ai, tauij, mtau)$par
      } else {
        RESULTS[j,] <- lar.model.res(model, Aij, Ai, tauij, mtau)$par
      }
    }
  }
  if (method == 'BBootstrap'){
  if(ncores>1){
    cl <- parallel::makeCluster(ncores) # not to overload your computer
    doParallel::registerDoParallel(cl)
    RESULTS = foreach::`%dopar%`(foreach::foreach(i = seq_len(B), .combine = rbind), {
      tp_list <- bin_make(tp, bin_len)
      sampboot <- lar_bootstrap(X, tp_list, group_id, seed)
      dat <- lar_nonparametric_estimation(sampboot, tp, group_id)
      Aij <- dat$Aij
      Ai <- dat$Ai
      tauij <- dat$tauij
      if (model == 'model_cl_fun'){
        res <- model_cl_fun(Aij, Ai, tauij, mtau)
      } else {
        res <- lar.model.res(model, Aij, Ai, tauij, mtau)
      }
      out <- res$par
      return(out)
    })
    parallel::stopCluster(cl)
  }else{
    RESULTS <- matrix(0, B, model.K)
    for(i in seq_len(B)){
      tp_list <- bin_make(tp, bin_len)
      sampboot <- lar_bootstrap(X, tp_list, group_id, seed)
      dat <- lar_nonparametric_estimation(sampboot, tp, group_id)
      Aij <- dat$Aij
      Ai <- dat$Ai
      tauij <- dat$tauij
      if (model == 'model_cl_fun'){
        RESULTS[j,] <- model_cl_fun(Aij, Ai, tauij, mtau)$par
      } else {
        RESULTS[j,] <- lar.model.res(model, Aij, Ai, tauij, mtau)$par
      }
      out <- RESULTS$par
      RESULTS[i,] <- out
    }
  }
  }

  dimpara <- sum(diag(model.H%*%(stats::var(RESULTS))))
  AIC <- 2*model.val + 2*model.K
  estimation_c <- lar_estimation_c(g_m, g_n, Aij, Ai, tauij, mtau)
  QAIC <- 2*model.val/estimation_c + 2*model.K
  CLIC <- model.val + dimpara
  res <- data.frame(AIC=AIC, QAIC=QAIC, CLIC=CLIC)
  return(res)
}

















