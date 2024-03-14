#' Coefficient of Variation (CV) calculation
#' Calculates by dividing the mean with the standard deviation
#' 
#' @param x data vector from which the coefficient of variation (CV) will be calculated
#'
#' @return the coefficients of variation
#' @export 
#'
#' @examples
CV = function(x){
  coeff = sd(x, na.rm = T)/mean(x, na.rm = T)
  return(coeff)
}