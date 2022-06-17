
#'
#' The width is the range between the 2.5 and 97.5% percentiles of the estimated parameters for many (such as 100) random runs
#'
#' @name width_ci
#' @title Compute the width of confidence intervals
#'
#' @param x A set of estimated parameters
#' @param alpha Confidence level is 1-alpha
#'
#' @return The width of confidence intervals.
#'
#' @export
#' @rdname width_ci
#'
#' @examples
#' x <- rnorm(100)
#' width_ci(x)
#'
width_ci <- function(x, alpha = 0.05){
  ci <- as.numeric(stats::quantile(x, probs = c(alpha/2, 1-alpha/2)))
  width <- as.numeric(diff(ci))
  res <- data.frame(low=ci[1], high=ci[2], width=width)
  return(res)
}
