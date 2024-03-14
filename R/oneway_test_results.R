#' Checks the equal variation and computes the oneway test accordingly
#' 
#'
#' @param data the dataframe where the test will be conducted
#' @param response_variable the variable to be tested  in the format "response_variable"
#' @param grouping_factor the grouping variable  in the format "grouping_factor"
#' @param log_normalise logical, wether the response variable will be log scaled or not depending on the normal distribution of the data

#' 
#' @return Gives the F values, p values and derees of freedom of the computed oneway test as well as the Barttlett test used to compute the variations
#' @export
#'
#' @examples 
oneway_test_results <- function(data, response_variable, grouping_factor, log_normalise){
  # Check if the specified columns exist in the dataset
  if (length(setdiff(c(response_variable, grouping_factor), names(data))) > 0) {
    stop("Specified column(s) not found in the dataset.")
  }

if (log_normalise == FALSE) {
  results = bartlett_oneway(data[[response_variable]], data[[grouping_factor]])
} else if (log_normalise == TRUE) {
  results = bartlett_oneway(log(data[[response_variable]]), data[[grouping_factor]])
} else {
  stop("Define log normalisation")
  }


results_table = data.table::data.table(
  pair = c("bartlett_results(K_val)", "oneway_results(F_val)"),
  p_value = signif(c(results[[1]]$p.value, results[[2]]$p.val), 2),
  F_value = signif(c(results[[1]]$statistic, results[[2]]$statistic), 2),
  df_value = c(results[[1]]$parameter, paste0(signif(results[[2]]$parameter[1],1), ":", signif(results[[2]]$parameter[2], 1)))
)
return(results_table)
}

### convert all the data$variable format to data[,"variable"] format. Try things from there onwards
