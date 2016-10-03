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
                 median(x$k), paste(", max ="), max(x$k),")\n")
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
                 median(x$k), paste(", max ="), max(x$k),")\n")
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
        median(x$k), paste(", max ="), max(x$k),")\n")
    print(output)
    cat("---\n")
    cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
    cat("---\n")
    cat(x$mod_notice) 
    
  } 
}