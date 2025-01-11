# Deal with the univariate case separately
#' Linear Model Adequacy Test
#' @description Some description
#'
#' @param x parameter
#' @examples
#' Include example
#' @export
#' 
#' 
kernel_density <- function(x, t, data, upper = 25){
  f_dens <- kernel_density_1d_cpp(as.matrix(x), t, as.matrix(data), upper)
  return(f_dens)
}

kernel_cdf <- function(x, t, data, upper = 25){
  f_dens <- kernel_cdf_1d_cpp(as.matrix(x), t, as.matrix(data), upper)
  return(f_dens)
}

