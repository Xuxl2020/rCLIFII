
#'
#' Break the observation time into different bins, with the width of the bin as bin_len
#'
#' @name bin_make
#' @title Generate the bin of observation time
#'
#' @param tp A set of observation time
#' @param bin_len The width of bin
#'
#' @return A list of different bins.
#'
#' @export
#' @rdname bin_make
#'
#' @examples
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' bin_len <- 50
#' bin_make(tp, bin_len)
#'

bin_make <- function(tp, bin_len){

  len <- max(tp)-min(tp)+1
  n_cut <- ceiling(len/bin_len)
  dat <- data.frame(tp=1:len, bin=rep(1:n_cut, each=bin_len)[1:len])
  idx <- match(tp, dat$tp)
  dat <- dat[idx,]
  tp_list <- split(dat$tp, dat$bin)
  return(tp_list)

}

