#' Do oneway test setting the variation according to the Bartlett results
#'
#' @param x response variable 
#' @param y grouping variable
#' @return a table of statistical output
#' @export
#'
#' @examples
bartlett_oneway = function(x, y){
  bartlett_result <- stats::bartlett.test(x ~ y)
  
  oneway_result <- if(bartlett_result$p.value > 0.05) {
    # var.equal is TRUE if bartlett result p > 0.05
    stats::oneway.test(x ~ y, var.equal = T)
  } else { 
    # var.equal is FALSE if bartlett result p < 0.05
    stats::oneway.test(x ~ y, var.equal = F)
  } 
  return(list(bartlett_result, oneway_result))
  }
