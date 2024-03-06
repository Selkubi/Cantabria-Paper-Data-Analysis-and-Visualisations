#' Converts below detection limit
#' The function 'convert_below_detection_limit' converts the data points which are "<0,1", "<0,01", "<0,02", "<0,06" within the dataframe to their half values. 
#' 
#' @param data 
#'
#' @return the dataframe with the converted datapoints
#' @export
#'
#' @examples
convert_below_detection_limit = function(data) {
  data = gsub(x = data, pattern = c("<0,1", "<0,01", "<0,02", "<0,06"), 
              replacement = c("0.05","0.005", "0.01", "0.03") ) #Or replacement is c("0.1","0.01", "0.02", "0.06") for abslute detection limits
  return(as.numeric(data))
}