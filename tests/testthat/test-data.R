
test_that("Results of model selection for LIR", {
  Z <- 300; N <- 100; n <- 40
  # set the observation time
  tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
  # set parameters
  lambda_C <- 0.08; mu_C <- 0.04
  # ncores, nboot
  ncores <- 2; nboot <- 2
  # simulation of Type C
  simulation_C <- lir_simulation_C(Z, N, n, lambda_C, mu_C, tp)
  # model selection
  lir_model_selection(X=simulation_C, n=40, tp, model='lir_1', method='Bootstrap', nboot=nboot, mtau = 1000, ncores = ncores)
  lir_model_selection(X=simulation_C, n=40, tp, model='lir_2', method='Bootstrap', nboot=nboot, mtau = 1000, ncores = ncores)
  lir_model_selection(X=simulation_C, n=40, tp, model='lir_3', method='Bootstrap', nboot=nboot, mtau = 1000, ncores = ncores)
  # set population size and number of subsampling
  N <- 100; n <- 4; W <- 10; U <- 20
  # set the observation time
  block_list <- list(c(1:5), c(51:55), c(101:105), c(501:505), c(601:605))
  tp <- unlist(block_list)
  # set parameters
  lambda_E <- 0.008
  # ncores, nboot
  ncores <- 2; nboot <- 2
  # simulation of Type E
  simulation_E <- lar_simulation_E(N, n, W, lambda_E, tp)
  # model selection
  lar_model_selection(X=simulation_E, model='lar_1', block_list, method='Jackknife', nboot=nboot, k=50, mtau = 1000, ncores = ncores)
  lar_model_selection(X=simulation_E, model='lar_2', block_list, method='Jackknife', nboot=nboot, k=50, mtau = 1000, ncores = ncores)
  lar_model_selection(X=simulation_E, model='lar_3', block_list, method='Jackknife', nboot=nboot, k=50, mtau = 1000, ncores = ncores)
})
