#' F test followed by t.test according to the calculated equal/unequal variations
#'
#' @param data A data.table where the var.test and t.test will be computed 
#' @param response_variable the variable to be tested  in the format "response_variable"
#' @param grouping_factor the grouping variable  in the format "grouping_factor"
#' @param log_normalise logical, wether the response variable will be log scaled or not depending on the normal distribution of the data
#'
#' @return A table of F values, degrees of freedom and P values of the computes F test and t test results
#' @export
#'
#' @examples
t_test_results <- function(data, response_variable, grouping_factor, log_normalise){
  # Check if the specified columns exist in the dataset
  if (sum(as.character(c(response_variable, grouping_factor)) %in% names(data)) < 2) {
    stop("Specified column(s) not found in the dataset.")
  }

  if (log_normalise == F) {
    results = var_test_t_test(data[[response_variable]], data[[grouping_factor]])
  } else if (log_normalise == T) {
    results = var_test_t_test(log(data[[response_variable]]), data[[grouping_factor]])
  } else {"Define log normalisation"}
  
  results_table = data.table::data.table(
    pair = c("F_test_results(F)", "t_test(t)"),
    p_value = signif(c(results[[1]]$p.value, results[[2]]$p.val), 2),
    F_value = signif(c(results[[1]]$statistic, results[[2]]$statistic), 2),
    df_value = c(paste0(results[[1]]$parameter[2], ":", results[[1]]$parameter[1]), results[[2]]$parameter)
  )
  return(results_table)
}
