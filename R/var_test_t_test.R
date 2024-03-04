#' Computes Two Sample t-test followed by F test 
#' F test compares variances which determines the var.equal statement in the following t-test
#' 
#' @param x response variable
#' @param y grouping variable
#'
#' @return A list of test statistics and mean values
#' @export
#'
#' @examples
var_test_t_test = function(x, y){
  var_test_results <- stats::var.test(x ~ y)
  
  t_test_results <- if(var_test_results$p.value > 0.05) {
    # var.equal is TRUE if bartlett result p > 0.05
    stats::t.test(x ~ y, var.equal = T)
  } else { 
    # var.equal is FALSE if bartlett result p < 0.05
    stats::t.test(x ~ y, var.equal = F)
  } 
  return(list(var_test_results, t_test_results))
}

