#' Convenience function for calculating group-mean covariates.
#' 
#' Creates a between-study (or between-cluster) version of the covariate in
#' question.
#' 
#' 
#' @param var The covariate cotaining the values to be group averaged.
#' @param grp The group from which the average should be calculated.
#' @return A column or vector containing the group.mean covariate.
#' @keywords robumeta
#' @examples
#' 
#' 
#' # Load data
#' data(corrdat)
#' 
#' # Create a group mean covariate 
#' age_m <- group.mean(corrdat$age, corrdat$studynum)
#' 
#' 
#' @export 
group.mean <- function(var, grp) {
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(tapply(var, grp, mean, na.rm = TRUE)[grp])
}
