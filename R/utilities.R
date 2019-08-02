#' Min-max scaling
#'
#' \code{scaling.min.max} rescales the range of a data set to [min, max].
#' @param data input data matrix, of size \eqn{n} by \eqn{d}.
#' @param min desired minimum of the transformed data.
#' @param max desired maximum of the transformed data.
#' @param na.rm logical indicating whether \code{NA} values should be stripped before the computation proceeds.
#' @return a matrix containing the scaled data.
#' @keywords utilities
#' @export
scaling.min.max = function(data, min = 0, max = 1) {
  data = scale(data,
               center = apply(as.matrix(data), 2, min),
               scale = diff(apply(as.matrix(data), 2, range)))
  data = data * (max - min) + min
  
  return(data)
}