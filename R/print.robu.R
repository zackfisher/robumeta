#' Outputs Model Information
#' 
#' Prints relevant information from robu function.
#' 
#' 
#' @param x Object from robu class.
#' @param digits Controls the number of digits to print when printing numeric
#' values.
#' @param ...  Additional arguments to be passed to the fitting function.
#' @references
#' 
#' Hedges, L.V., Tipton, E., Johnson, M.C. (2010) Robust variance estimation in
#' meta-regression with dependent effect size estimates. \emph{Research
#' Synthesis Methods}. \bold{1}(1): 39--65. Erratum in \bold{1}(2): 164--165.
#' DOI: 10.1002/jrsm.5
#' 
#' Tipton, E. (in press) Small sample adjustments for robust variance
#' estimation with meta-regression. \emph{Psychological Methods}.
#' @keywords robu
#' @examples
#' 
#' 
#' # Load data
#' data(hierdat)
#' 
#' ### Small-Sample Corrections - Hierarchical Dependence Model
#' HierMod  <-  robu(formula = effectsize ~ binge + followup + sreport
#'                    + age, data = hierdat, studynum = studyid, 
#'                    var.eff.size = var, modelweights = "HIER", small = FALSE)
#' 
#' print(HierMod) # Output results
#' 
#' @export
print.robu <- function(x, digits = 3,...){
  
  user_weighting <- x$user_weighting
  modelweights   <- x$modelweights
  mod_info       <- x$mod_info
  
  output              <- x$reg_table
  output              <- format(output, trim = TRUE, digits = digits, 
                                scientific = FALSE)
  colnames(output)    <- c("", "Estimate","StdErr", "t-value", "dfs", "P(|t|>)", 
                           "95% CI.L","95% CI.U", "Sig")
  
  if(!user_weighting){
    
    switch(modelweights,
           
           HIER = { # Begin HIER 
             
             cat(x$mod_label, "\n")
             cat("\nModel:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
             cat(paste("Number of clusters ="), x$N, "\n")
             cat(paste("Number of outcomes ="), sum(x$k), paste("(min ="), min(x$k),
                 paste(", mean ="), format(mean(x$k), digits = 3), paste(", median ="), 
                 stats::median(x$k), paste(", max ="), max(x$k),")\n")
             cat(paste("Omega.sq ="), mod_info$omega.sq, "\n")
             cat(paste("Tau.sq ="), mod_info$tau.sq, "\n\n")
             print(output, digits = 3)
             cat("---\n")
             cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
             cat("---\n")
             cat(x$mod_notice)
             
           }, # End HIER
           
           CORR = { # Begin CORR
             
             cat(x$mod_label, "\n")
             cat("\nModel:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
             cat(paste("Number of studies ="), x$N, "\n")
             cat(paste("Number of outcomes ="), sum(x$k), paste("(min ="), min(x$k),
                 paste(", mean ="), format(mean(x$k), digits = 3), paste(", median ="), 
                 stats::median(x$k), paste(", max ="), max(x$k),")\n")
             cat(paste("Rho ="), mod_info$rho, "\n")
             cat(paste("I.sq ="), mod_info$I.2, "\n")
             cat(paste("Tau.sq ="), mod_info$tau.sq, "\n\n")
             print(output)
             cat("---\n")
             cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
             cat("---\n")
             cat(x$mod_notice)      
             
           } # End CORR
           
    ) 
    
  } else { # Begin userweights
    
    cat(x$mod_label, "\n")
    cat("\nModel:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
    cat(paste("Number of studies ="), x$N, "\n")
    cat(paste("Number of outcomes ="), sum(x$k), paste("(min ="), min(x$k),
        paste(", mean ="), format(mean(x$k), digits = 3), paste(", median ="), 
        stats::median(x$k), paste(", max ="), max(x$k),")\n")
    print(output)
    cat("---\n")
    cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
    cat("---\n")
    cat(x$mod_notice) 
    
  } 
}
