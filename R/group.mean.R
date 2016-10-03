group.mean <- function(var, grp) {
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(tapply(var, grp, mean, na.rm = TRUE)[grp])
}