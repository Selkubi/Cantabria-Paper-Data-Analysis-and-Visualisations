#' Title
#'
#' @param data 
#' @param response_variable 
#' @param grouping_factor 
#' @param log_normalise 
#'
#' @return A table of adjusted p values
#' @export
#'
#' @examples
statistics_pipeline_wrapper <- function(data, response_variable, grouping_factor, log_normalise){
  
  mean_table <- oneway_test_results(data, response_variable, grouping_factor, log_normalise)
  
  # if (sum(mean_table$pair == c("bartlett_results(K_val)", "oneway_results(F_val)")) == 2 ) {
  #   stop("Specified column(s) not found in the dataset.")
  # }
  
  datasets <- list(Natural = data[alteration == "Natural"], 
                   Temperate = data[Class == "Temperate"], 
                   Mediterranean = data[Class == "Mediterranean"])
  
  if(mean_table[pair == "oneway_results(F_val)"]$p_value < 0.05) {
    
    results_table <- lapply(datasets, t_test_results,
                                          response_variable = response_variable,
                                          grouping_factor = grouping_factor,
                                          log_normalise = log_normalise)
    
    adjusted_p_values <- stats::p.adjust(c(Natural = results_table$Natural$p_value[2], # t_test(t)
                                                  Temperate = results_table$Temperate$p_value[2], # t_test(t)
                                                  Mediterranean = results_table$Mediterranean$p_value[2]), # t_test(t)
                                       method = "bonferroni", n = 3)
    
    adjusted_p_values <- data.table("value" = c("p_values", "t_values", "df_value"),
               "Natural" = c(adjusted_p_values[["Natural"]], results_table$Natural[2,3:4]), # t_test(t)
               "Temperate" = c(adjusted_p_values[["Temperate"]], results_table$Temperate[2,3:4]), # t_test(t)
               "Mediterranean" = c(adjusted_p_values[["Mediterranean"]], results_table$Mediterranean[2,3:4])) # t_test(t)
  } else if(mean_table["oneway_results(F_val)"]$p_value >= 0.05){
    
    adjusted_p_values <- matrix(c("-", "-","-"), nrow = 1, ncol = 3, dimnames = list(c("t_test_p_values"), c("Natural", "Temperate", "Mediterranean")))
 
     } else {stop("Error in p values")}

  
  if(mean_table[pair == "bartlett_results(K_val)"]$p_value < 0.05) {
    
    results_table <- lapply(datasets, t_test_results,
                            response_variable = response_variable,
                            grouping_factor = grouping_factor,
                            log_normalise = log_normalise)
    
    F_results_table <- data.table("value" = c("p_values", "F_values", "df_value"),
                                    "Natural" = c(results_table$Natural[1,2:4]), # F_test_results(F)
                                    "Temperate" = c(results_table$Temperate[1,2:4]), # F_test_results(F)
                                    "Mediterranean" = c(results_table$Mediterranean[1,2:4])) # F_test_results(F)
  
    }  else if(mean_table[pair == "bartlett_results(K_val)"]$p_value >= 0.05){
    
    F_results_table <- matrix(c("-", "-","-"), nrow = 1, ncol = 3, dimnames = list(c("F_test_p_values"), 
                                                                                   c("Natural", "Temperate", "Mediterranean")))
  } else {stop}
  
  
  return(list(t_test = adjusted_p_values, F_test = F_results_table))
}
