#' Title
#'
#' @param data 
#' @param response_variable 
#' @param grouping_factor 
#' @param log_normalise 
#'
#' @return
#' @export
#'
#' @examples
statistics_pipeline_wrapper <- function(data, response_variable, grouping_factor, log_normalise){
  
  mean_table <- oneway_test_results(data = data, response_variable, grouping_factor, log_normalise)
  
  datasets <- list(Natural = data[alteration == "Natural"], 
                   Temperate = data[Class == "Temperate"], 
                   Mediterranean = data[Class == "Mediterranean"])
  
  if(mean_table[2]$p_value < 0.05) {
    results_table <- lapply(datasets, t_test_results,
                                          response_variable = response_variable,
                                          grouping_factor = grouping_factor,
                                          log_normalise = log_normalise)
    
    adjusted_p_values <- p.adjust(c(Natural = results_table$Natural$p_value[2], 
                                    Temperate = results_table$Temperate$p_value[2], 
                                    Mediterranean = results_table$Mediterranean$p_value[2]),
                                  method = "bonferroni", n = 3)
    }
  else if (mean_table[2]$p_value > 0.05){
    adjusted_p_values = c(Natural = "-", 
                          Temperate = "-", 
                          Mediterranean = "-")
    }
  else {stop}
  
  return(adjusted_p_values)
}
