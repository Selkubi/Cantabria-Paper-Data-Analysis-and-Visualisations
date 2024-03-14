
#' Calculate mean and standard deviation of values with a grouping factor
#' @import data.table
#' 
#' @param data A data table
#' @param response_variable The response variable of which the mean+-sd calculation will be computed on
#' @param grouping_factor The grouping factor over which the mean+-sd calculation will be computed on
#'
#' @return A table of the results
#' @export
#'
#' @examples


get_mean_and_sd <- function(data, response_variable, grouping_factor){
  table <- data.table::data.table(name = aggregate(data[, response_variable, with = FALSE], by = data[, .(get(grouping_factor))], FUN = sd)[,1],
                      value = paste0(signif(aggregate(data[, response_variable, with = FALSE], by = data[, .(get(grouping_factor))], FUN = mean)[, 2], 3), "+-", 
                                     signif(aggregate(data[, response_variable, with = FALSE], by = data[, .(get(grouping_factor))], FUN = sd)[, 2], 1) ))
  setnames(table, "value", response_variable)
  return(table)
}
