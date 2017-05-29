#' Convenience function for calculating group-centered covariates.
#' 
#' Creates a within-study (or within-cluster) version of the covariate in
#' question.
#' 
#' 
#' @param var The covariate to be group centered.
#' @param grp A vector corresponding to the group identification.
#' @return A column or vector containing the group.centered covariate.
#' @keywords robumeta
#' @examples
#' 
#' 
#' # Load data
#' data(corrdat) 
#' 
#' # Create a group centered covariate 
#' males_c <- group.center(corrdat$males, corrdat$studyid)
#' 
#' 
#' @export 
group.center <- function(var, grp) {
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(var - tapply(var, grp, mean, na.rm = TRUE)[grp])
}
