#' Common theme function for the boxplots and pca plots
#'
#' @return ggplott theme
#' @export
#'
#' @examples
#' 
theme_boxplot <- function() {
  ggplot2::theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11, color = "black"), 
        strip.placement ="outside")+
  ggplot2::theme(text = element_text(size = 11),
        axis.line.x =  element_line(color = "black"),
        axis.line.y =  element_line(),
        axis.text = element_text(size = 11, color = "black"), 
        axis.text.x = element_text(hjust = 1),
        axis.title = element_blank(), 
        legend.position = "none") 
}
