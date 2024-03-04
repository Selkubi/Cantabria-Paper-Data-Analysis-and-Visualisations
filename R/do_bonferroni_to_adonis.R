#' Does bonferroni corrections on 3 pairs of permanova (adonis) results
#'
#' @param x_nat  p value of the first pairwise comparison
#' @param x_te  p value of the second pairwise comparison
#' @param x_me  p value of the third pairwise comparison
#'
#' @return A table of results of the three pairwise comparisons Bonferroni corrected, with their p values, F values and degrees of Freedomi 
#' @export
#'
#' @examples 
do_bonferroni_to_adonis = function(x_nat, x_te, x_me) {
  # compute the bonferroni corrected p values
  adjusted_p_values = p.adjust(c(x_nat["Pr(>F)"][1,], x_te["Pr(>F)"][1,], x_me["Pr(>F)"][1,]), method = "bonferroni", n = 3)
  # make a table with the corrected p values and the related F and df values
  results = data.table(pair = c("nA-nM", "nA-aA", "nM-aM"),
                       p_value = signif(adjusted_p_values, 3),
                       F_value = signif(c(x_nat[["F"]][1], x_te[["F"]][1], x_me[["F"]][1]), 3),
                       df_value = c(paste0(x_nat[["Df"]][1], ":", x_nat[["Df"]][2]),
                                    paste0(x_te[["Df"]][1], ":", x_te[["Df"]][2]),
                                    paste0(x_me[["Df"]][1], ":", x_me[["Df"]][2]))
  )
  return(results)
}