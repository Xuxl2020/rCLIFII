
#' @name lir_model_selection
#' @title Model selection for models of lagged identification rate
#'
#' @param X A list or matrix containing the identities of individuals identified in each sampling period
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period
#' It indicates the same number of individuals identified in all sampling periods if a positive integer
#' @param tp A set of observed time
#' @param model Models of lagged identification rate, model = 'lir_1', 'lir_2', 'lir_3', or formulate model by yourself 'model_cl_fun'
#' @param model_cl_fun If you formulate your model, please input function to calculate the composite likelihood about your model
#' @param cl.H If you formulate your model, please input the sensitivity matrix with respect to parameters in your model
#' @param model.K If you formulate your model, please input the number of parameters in your model
#' @param method The method = 'Bootstrap', 'BBootstrap', or 'Jackknife'
#' @param nboot The number of bootstrap samples desired
#' @param bin_len An integer represents len-time-unit intervals
#' @param mtau The maximum allowable lag time
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


lir_model_selection <- function(X, n, tp,
                                model,
                                method,
                                ncores = 4,
                                mtau = 1000,
                                nboot = -1,
                                bin_len = -1,
                                model_cl_fun = NULL,
                                cl.H = NULL,
                                model.K = NULL,
                                seed = NULL){

  tp = tp - min(tp) + 1

  if (length(n)>1 && length(n)!=length(tp)) {
    stop("'n' or 'tp' is valid.")
  }

  if (length(n) == 1) {
    n <- rep(n, length(tp))
  }

  if (method == 'Bootstrap'|method == 'BBootstrap'){
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

  lir_data <- lir_nonparametric_estimation(X, n, tp, mtau)

  R_m <- lir_data$R_m
  R_n <- lir_data$R_n
  mij <- lir_data$mij
  nij <- lir_data$nij
  tauij <- lir_data$tauij

  if(model == 'lir_1'){
    model = 'Model1'
    model.H <- lir.model.res(model, mij, nij, tauij, mtau)$H
    model.est <- lir.model.res(model, mij, nij, tauij, mtau)$par
    model.val <- lir.model.res(model, mij, nij, tauij, mtau)$val
    model.K <- lir.model.res(model, mij, nij, tauij, mtau)$K
  }
  if(model == 'lir_2'){
    model = 'Model2'
    model.H <- lir.model.res('Model2', mij, nij, tauij, mtau)$H
    model.est <- lir.model.res('Model2', mij, nij, tauij, mtau)$par
    model.val <- lir.model.res('Model2', mij, nij, tauij, mtau)$val
    model.K <- lir.model.res('Model2', mij, nij, tauij, mtau)$K
  }
  if(model == 'lir_3'){
    model = 'Model3'
    model.H <- lir.model.res('Model3', mij, nij, tauij, mtau)$H
    model.est <- lir.model.res('Model3', mij, nij, tauij, mtau)$par
    model.val <- lir.model.res('Model3', mij, nij, tauij, mtau)$val
    model.K <- lir.model.res('Model3', mij, nij, tauij, mtau)$K
  }

  if(model == 'model_cl_fun'){
    model.H <- cl.H
    model.est <- model_cl_fun(mij, nij, tauij, mtau)$par
    model.val <- model_cl_fun(mij, nij, tauij, mtau)$val
    model.K <- model.K
  }

  if(method == 'Jackknife'){
    if (bin_len > 0){
      jsamples <- jackknife(X, tp, bin_len)
    }

    RESULTS <- matrix(0, length(jsamples), length(model.est))
    for(j in 1:length(jsamples)){
      tp0 = as.numeric(colnames(jsamples[[j]]))
      lir_data <- lir_nonparametric_estimation(as.matrix(jsamples[[j]]), n[tp %in% tp0], tp0, mtau)
      mij <- lir_data$mij
      nij <- lir_data$nij
      tauij <- lir_data$tauij
      if (model == 'model_cl_fun'){
        RESULTS[j,] <- model_cl_fun(mij, nij, tauij, mtau)$par
      } else {
        RESULTS[j,] <- lir.model.res(model, mij, nij, tauij, mtau)$par
      }
    }
  }
  if (method == 'Bootstrap'| method == 'BBootstrap'){
    if(ncores>1){
      cl <- parallel::makeCluster(ncores) # not to overload your computer
      doParallel::registerDoParallel(cl)
      RESULTS = foreach::`%dopar%`(foreach::foreach(i = seq_len(B), .combine = rbind), {
        if (method == 'Bootstrap'){
          sampboot <- lir_bootstrap(X, tp, seed)
        }
        if (method == 'BBootstrap'){
          tp_list <- bin_make(tp, bin_len)
          sampboot <- lar_bootstrap(X, tp_list)
        }
        dat <- lir_nonparametric_estimation(sampboot, n, tp, mtau)
        mij <- dat$mij
        nij <- dat$nij
        tauij <- dat$tauij
        if (model == 'model_cl_fun'){
          res <- model_cl_fun(mij, nij, tauij, mtau)
        } else {
          res <- lir.model.res(model, mij, nij, tauij, mtau)
        }
        out <- res$par
        return(out)
      })
      parallel::stopCluster(cl)
    }else{
      RESULTS <- matrix(0, B, model.K)
      for(i in seq_len(B)){
        if (method == 'Bootstrap'){
          sampboot <- lir_bootstrap(X, tp)
        }
        if (method == 'BBootstrap'){
          tp_list <- bin_make(tp, bin_len)
          sampboot <- lar_bootstrap(X, tp_list)
        }
        dat <- lir_nonparametric_estimation(sampboot, n, tp, mtau)
        mij <- dat$mij
        nij <- dat$nij
        tauij <- dat$tauij
        if (model == 'model_cl_fun'){
          RESULTS[i,] <- model_cl_fun(mij, nij, tauij, mtau)$par
        } else {
          RESULTS[i,] <- lir.model.res(model, mij, nij, tauij, mtau)$par
        }
      }
    }
  }
  dimpara <- sum(diag(model.H%*%(stats::var(RESULTS))))
  AIC <- 2*model.val + 2*model.K
  estimation_c <- lir_estimation_c(R_m, R_n, mij, nij, tauij, mtau)
  QAIC <- 2*model.val/estimation_c + 2*model.K
  CLIC <- model.val + abs(dimpara)
  res <- data.frame(AIC=AIC, QAIC=QAIC, CLIC=CLIC)
  return(res)
}

