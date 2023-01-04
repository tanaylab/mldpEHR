#' Vectorized mean, similiar to \code{pmin} and \code{pmax}
#'
#' @param ... numeric vectors to average
#' @param na.rm a logical indicating whether missing values should be removed
#'
#' @return a vector with mean of \code{...} arguments
#'
#' @examples
#' psum(c(1, 2, 3), c(4, 5, 6))
#' psum(c(1, 2, 3), pi)
#' psum(c(1, NA, 3), c(4, 5, 6), na.rm = TRUE)
#'
#' @noRd
pmean <- function(..., na.rm = FALSE) {
    d <- do.call(cbind, list(...))
    res <- rowMeans(d, na.rm = na.rm)
    idx_na <- !rowMeans(!is.na(d))
    res[idx_na] <- NA
    return(res)
}
