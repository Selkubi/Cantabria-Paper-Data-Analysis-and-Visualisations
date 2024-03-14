var_test_t_test = function(data, x, y){
  
  datasets <- list(Natural = data[alteration == "Natural"], 
                   Temperate = data[Class == "Temperate"], 
                   Mediterranean = data[Class == "Mediterranean"])
  
  
  var_test_results <- stats::var.test(x ~ y)
  
  t_test_results <- if(var_test_results$p.value >= 0.05) {
    # var.equal is TRUE if bartlett result p > 0.05
    stats::t.test(x ~ y, var.equal = T)
  } else { 
    # var.equal is FALSE if bartlett result p < 0.05
    stats::t.test(x ~ y, var.equal = F)
  } 
  return(list(var_test_results, t_test_results))
}