#' Sensitivity Analysis for Correlated Effects RVE
#' 
#' \code{sensitivity} is used to assess the impact of differing rho values on
#' the correlated effects meta-regression model.
#' 
#' 
#' @aliases sensitivity
#' @param x A dataframe containing values of rho, tau squared, coefficient
#' estimates, and standard errors.
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
#' # Correlated Effects Model
#' CorrMod   <-  robu(formula = effectsize ~ followup + males + binge + college, 
#'                    data = corrdat, studynum = studyid, var.eff.size = var, 
#'                    rho = .8, modelweights = "CORR", small = FALSE)
#' 
#' sensitivity(CorrMod) # Output sensitivity
#' 
#' @export
sensitivity <- function(x){
  
  modelweights   <- x$modelweights
  user_weighting <- x$user_weighting
  
  if(modelweights == "HIER")  
    stop("Sensitivity analysis is not available for hierarchical effects.")
  
  if(user_weighting == TRUE)  
    stop("Sensitivity analysis is not available for user specified weights.")
  
  mod_info  <- x$mod_info
  p         <- x$p
  N         <- x$N
  Xreg      <- x$Xreg
  y         <- x$y
  X         <- x$X
  data.full <- x$data.full
  X.full    <- x$X.full
  k         <- data.full$k
  k_list    <- x$k_list
  ml        <- x$ml
  term1     <- mod_info$term1
  term2     <- mod_info$term2
  small     <- x$small
  labels    <- x$labels
  mod_label <- x$mod_label
  rho.test  <- seq(0, 1, .2)
  
  rho_labels   <- c(paste("Rho = ", seq(0, 1, .2), sep=""))
  var_labels   <- rep("", 2 * (p + 1))
  var_labels[seq(1, length(var_labels), by = 2)] <- labels
  var_labels   <- c(var_labels, "Tau.sq")
  
  col2_labels  <- rep("", 2 * (p + 1))
  col2_labels[seq(1, length(col2_labels), by = 2)] <- "Coefficient"
  col2_labels[seq(2, length(col2_labels), by = 2)] <- "Std. Error"
  col2_labels  <- c(col2_labels, "Estimate")
  sen          <- data.frame(cbind(var_labels, col2_labels))
  
  for (i in (1: length(rho.test))){
    
    tau.sq1             <- term1 + rho.test[i] * term2 
    tau.sq              <- ifelse(tau.sq1 < 0, 0, tau.sq1)
    data.full$r.weights <- 1 / (as.vector(data.full$k) * 
                                  (as.vector(data.full$avg.var.eff.size)
                                   + as.vector(tau.sq)))
    W.r.big             <- diag(data.full$r.weights)  # W
    W.r                 <- by(data.full$r.weights, data.full$study, 
                              function(x) diag(x, nrow = length(x)))
    sumXWX.r            <- Reduce("+", Map(function(X, W) 
      t(X) %*% W %*% X, 
      X, W.r))
    sumXWy.r            <- Reduce("+", Map(function(X, W, y) 
      t(X) %*% W %*% y, 
      X, W.r, y))
    b.r                 <- solve(sumXWX.r) %*% sumXWy.r 
    data.full$pred.r    <- Xreg %*% b.r
    data.full$e.r       <- cbind(data.full$effect.size) - 
      data.full$pred.r
    data.full$e.r       <- as.numeric(data.full$e.r)
    sigma.hat.r         <- by(data.full$e.r, data.full$study, 
                              function(x) tcrossprod(x))
    
    if (!small) { # Begin small = FALSE 
      
      sumXWeeWX.r <- Reduce("+", Map(function(X, W, V) 
        t(X) %*% W %*% V %*% W %*% X, 
        X, W.r, sigma.hat.r))
      VR.r        <- solve(sumXWX.r) %*% sumXWeeWX.r %*% 
        solve(sumXWX.r)  
      SE          <- sqrt(diag(VR.r)) * sqrt(N / (N - (p + 1)))
      
    } else { 
      
      Q             <- solve(sumXWX.r) #
      Q.list        <- rep(list(Q), N)
      H             <- Xreg %*% Q %*% t(Xreg) %*% W.r.big 
      ImH           <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) - H
      data.full$ImH <- cbind(ImH)
      ImHj          <- by(data.full$ImH, data.full$study, 
                          function(x) as.matrix(x))
      dfS           <- c(rep(0, p + 1))
      diag_one      <- by(rep(1, nrow(X.full)), X.full$study, 
                          function(x) diag(x, nrow = length(x)))
      ImHii         <- Map(function(X, Q, W, D) 
        D - X %*% Q %*% t(X) %*% W,
        X, Q.list, W.r, diag_one)
      eigenvec <- lapply(ImHii, function(x) eigen(x)$vectors) 
      eigenval <- lapply(ImHii, function(x) eigen(x)$values)
      I        <- ImHii
      A.MBB    <- Map(function (eigenvec, eigenval, k_list) 
        eigenvec %*% diag(1/sqrt(eigenval), 
                          k_list, k_list) 
        %*% t(eigenvec),
        eigenvec, eigenval, k_list)
      A.MBB1    <- Map(function(K, A, I) 
        if (K > 1) A else matrix(sqrt(solve(I))), 
        k_list, A.MBB, I)
      A.MBB2    <- A.MBB 
      
      sumXWA.MBBeeA.MBBWX.r <- Reduce("+", Map(function(X,W,A,S) 
        t(X) %*% W %*% A %*% 
          S %*% A %*% W %*% X, 
        X, W.r, A.MBB2, 
        sigma.hat.r))
      data.full$ImH         <- ImH
      ImH                   <- lapply(split(data.full$ImH, 
                                            data.full$study), 
                                      matrix, ncol=nrow(data.full))
      giTemp                <- Map(function(I, A, W, X, Q)
        t(I) %*% A %*% W %*% X %*% Q, 
        ImHj, A.MBB2, W.r, X, Q.list) 
      dfs                   <- c(rep(0, p + 1))
      
      for (i in 1:(p + 1)) { 
        L      <- c(rep(0,p + 1))
        L[i]   <- 1
        Ll     <- rep(list(L), N)
        gi     <- Map(function(G, L) G %*% cbind(L), giTemp, Ll)
        G      <- Reduce("+", lapply(gi, function(x) tcrossprod(x)))
        B      <- solve(sqrt(W.r.big) )%*% G %*% solve(sqrt(W.r.big))
        e.val2 <- eigen(B)
        dfs[i] <- sum(e.val2$values)^2 / sum(e.val2$values^2)
      }
      
      VR.MBB1 <- solve(sumXWX.r) %*% sumXWA.MBBeeA.MBBWX.r %*% 
        solve(sumXWX.r)
      VR.r    <- VR.MBB1
      SE      <- sqrt(diag(VR.r))         
    } 
    
    vals <- c()
    temp_vals <- c()
    for (i in 1:(p + 1)) { 
      temp_vals <- c(b.r[i], SE[i])
      vals <- c(vals, temp_vals)
    }
    vals <- c(vals, tau.sq)
    vals <- format(vals, digits=3, justify="centre")
    sen <- cbind(sen, vals)
  }
  colnames(sen)     <- c(" ", " ", rho_labels)
  format(sen[,1], justify = "left")
  cat(mod_label, "\n")
  cat("Model:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
  cat(paste("Sensitivity Analysis"), "\n\n")
  print.data.frame(sen, quote = FALSE, row.names = FALSE, right = FALSE)
}
