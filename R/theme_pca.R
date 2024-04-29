#' Common theme function for the boxplots and pca plots
#'
#' @return ggplott theme
#' @export
#'
#' @examples
#' 
theme_pca <- function() {
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size=11),
    axis.title = ggplot2::element_text(size=11),
    axis.text = ggplot2::element_text(size=11, color="black"), 
    legend.position = "none", 
    strip.placement ="outside", 
    strip.background = ggplot2::element_blank(), 
    strip.text = ggplot2::element_blank(),
    axis.line =  ggplot2::element_line(color="black"),
    panel.border = ggplot2::element_blank(),
    axis.line.y.right = ggplot2::element_line(color="white"),
    panel.grid = ggplot2::element_blank())
}
