robu     <- function(formula, data, studynum,var.eff.size, userweights,
                     modelweights = c("CORR", "HIER"), rho = 0.8, 
                     small = TRUE, ...) {
  
  # Evaluate model weighting scheme.
  modelweights <- match.arg(modelweights)

  if (modelweights == "CORR" && rho > 1 | rho < 0)  
      stop ("Rho must be a value between 0 and 1.")
  
  if (missing(userweights)){
     user_weighting = FALSE 
  } else { 
     user_weighting = TRUE
  }

  cl                       <- match.call() # Full model call
  mf                       <- match.call(expand.dots = FALSE)
  ml                       <- mf[[2]] # Model formula 
  m                        <- match(c("formula", "data", "studynum", 
                                      "var.eff.size", "userweights"), names(mf))
  mf                       <- mf[c(1L, m)] 
  mf$drop.unused.levels    <- TRUE
  mf[[1L]]                 <- as.name("model.frame")
  mf                       <- eval(mf, parent.frame()) 
  
  if(!user_weighting){ 
    
    dframe                 <- data.frame(effect.size = mf[,1],
                                        model.matrix(formula, mf),
                                        studynum = mf[["(studynum)"]],
                                        var.eff.size = mf[["(var.eff.size)"]])
    
    X.full.names           <- names(dframe)[-match(c("effect.size", 
                                                     "studynum", 
                                                     "var.eff.size"), 
                                                   names(dframe))] 
    
  } else { # Begin userweights
    
    dframe                 <- data.frame(effect.size = mf[,1],
                                        model.matrix(formula, mf), 
                                        studynum = mf[["(studynum)"]], 
                                        var.eff.size = mf[["(var.eff.size)"]],
                                        userweights = mf[["(userweights)"]])
    
    X.full.names           <- names(dframe)[-match(c("effect.size", 
                                                     "studynum", 
                                                     "userweights", 
                                                     "var.eff.size"), 
                                                   names(dframe))] 
  } # End userweights
  study_orig_id            <- dframe$studynum
  dframe$study             <- as.factor(dframe$studynum)
  dframe$study             <- as.numeric(dframe$study)
  dframe                   <- dframe[order(dframe$study),]
  k_temp                   <- as.data.frame(unclass(rle(sort(dframe$study))))
  dframe$k                 <- k_temp[[1]][ match(dframe$study, k_temp[[2]])]
  dframe$avg.var.eff.size  <- ave(dframe$var.eff.size, dframe$study)
  dframe$sd.eff.size       <- sqrt(dframe$var.eff.size)
  
  switch(modelweights, 
         
    HIER = { # Begin HIER
      
        dframe$weights <- 1 / dframe$var.eff.size
      
    }, # End HIER
    
    CORR = { # Begin CORR
      
        dframe$weights <- 1 / (dframe$k * dframe$avg.var.eff.size)
      
    } # End CORR
    
  ) 
  
  X.full           <- dframe[c("study", X.full.names)]
  
  data.full.names  <- names(dframe)[-match(c("studynum",X.full.names), 
                                          names(dframe))] 
  data.full        <- dframe[c(data.full.names)]
  k                <- data.full[ !duplicated(data.full$study), ]$k
  k_list           <- as.list(k) 
  M                <- nrow(data.full) # Number of units in analysis
  p                <- ncol(X.full) - 2 # Number of (non-intercept) covariates 
  N                <- max(data.full$study) # Number of studies
  W                <- as.matrix(by(data.full$weights, data.full$study, 
                                   function(x) diag(x, nrow = length(x)),
                                   simplify = FALSE))
  X                <- data.matrix(X.full)
  X                <- lapply(split(X[,2:(p + 2)], X[,1]), matrix, ncol = p + 1)
  y                <- by(data.full$effect.size, data.full$study, 
                         function(x) matrix(x))
  J                <- by(rep(1, nrow(X.full)), X.full$study, 
                         function(x) matrix(x, nrow = length(x), 
                                               ncol = length(x)))
  sigma            <- by(data.full$sd.eff.size, data.full$study, 
                         function(x) tcrossprod(x))
  vee              <- by(data.full$var.eff.size, data.full$study, 
                         function(x) diag(x, nrow = length(x)))
  SigmV            <- Map(function(sigma, V) 
                          sigma - V, sigma, vee)
  sumXWX           <- Reduce("+", Map(function(X, W) 
                                      t(X) %*% W %*% X, 
                                      X, W))
  sumXWy           <- Reduce("+", Map(function(X, W, y) 
                                      t(X) %*% W %*% y, 
                                      X, W, y))
  sumXWJWX         <- Reduce("+", Map(function(X, W, J) 
                                      t(X) %*% W %*% J %*% W %*% X, 
                                      X, W, J))

  Matrx_WKXX       <- Reduce("+",  
                             Map(function(X, W, k) { t(X) %*% (W / k) %*% X},  
                             X, W, k_list))  
    
  Matrx_wk_XJX_XX <- Reduce("+", 
                            Map(function(X, W, J, k) {(W / k)[1,1] * ( t(X) %*% J %*% X - t(X) %*% X) }, 
                           X, W, J, k_list))  
  
  switch(modelweights, 
    
    HIER = { # Begin HIER
        
        tr.sumJJ <- Reduce("+", Map(function(J) 
                                    sum(diag(J %*% J)), 
                                    J)) 
        sumXJX   <- Reduce("+", Map(function(X, J) 
                                    t(X) %*% J %*% X, 
                                    X, J))
        sumXWJJX <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% W %*% J %*% J %*% X, 
                                    X, W, J))
        sumXJJWX <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% J %*% J %*% W %*% X, 
                                    X, W, J))
        sumXWWX  <- Reduce("+", Map(function(X, W) 
                                    t(X) %*% W %*% W %*% X, 
                                    X, W))
        sumXJWX  <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% J %*% W %*% X, 
                                    X , W, J))
        sumXWJX  <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% W %*% J %*% X, 
                                    X, W, J))
    } # End HIER
    
  ) 
  
  b              <- solve(sumXWX) %*% sumXWy 
  Xreg           <- as.matrix(X.full[-c(1)], dimnames = NULL)
  data.full$pred <- Xreg %*% b 
  data.full$e    <- data.full$effect.size - data.full$pred 
  
  if (!user_weighting) { 
    
  switch(modelweights, 
    
    HIER = { # Begin HIER
      
      # Sigma_aj = tau.sq * J_j + omega.sq * I_j + V_j 
      # Qe is sum of squares 1
      # Qe = Sigma(T'WT)-(Sigma(T'WX)(Sigma(X'WX))^-1(Sigma(X'WT)
      # where W = V^(-1) and V = data.full$var.eff.size
      # Also, Qe = (y-xb)' W (y-xb)
      sumV <- sum(data.full$var.eff.size)
      W    <- diag(1 / data.full$var.eff.size) 
      sumW <- sum(W)
      Qe   <- t(data.full$e) %*% W %*% data.full$e
     
      # Qa is sum of squares 2
      # Qa = sum(T-XB.hat)'J(T-XB.hat)
      # where B.hat = (X'WX)^-1(X'WT)
      # Also, Qa = (y-xb)'A (y-xb), A=diag(J)
      e      <- by(data.full$e, data.full$study, function(x) matrix(x))
      sumEJE <- Reduce("+", Map(function(e, J) t(e) %*% J %*% e, e, J))
      Qa     <- sumEJE 
      
      # MoM estimators for tau.sq and omega.sq can be written as
      # omega.sq.h = A2(Qa-C1)-A1(Qe-C2) / B1A2-B2A1
      # tau.sq.h = Qe-C2/A2 - omega.sq.h(B2/A2) where
      
      # Vi = (t(X)WX)^-1
      V.i    <- solve(sumXWX)
      
      # A1 = Sigma(kj^2) - tr(V*Sigma(kj*t(Xj)*Jj*Wj*Xj)) - 
      #                    tr(V*Sigma(kj*t(Xj)*Jj*Wj*Xj)) +
      #                    tr(V*[Sigma(t(Xj)*Jj*Xj)]*V*Sigma(t(Xj)*Wj*Jj*Wj*Xj))
      # B1 = Sigma(kj)   - tr(V Sigma(t(Xj)*Jj*Wj*Xj)) -
      #                    tr(V Sigma(t(Xj)*Wj*Jj*Xj)) +
      #                    tr(V*[Sigma(t(Xj)*Jj*Xj)]*V*Sigma(t(Xj)*Wj^2*Xj)) 
      # C1 = tr(W^-1)    - tr(V*Sigma(t(X)*Jj*Xj))
      
      A1    <- tr.sumJJ  - sum(diag(V.i %*% sumXJJWX)) - 
                           sum(diag(V.i %*% sumXWJJX)) + 
                           sum(diag(V.i %*% sumXJX %*% V.i %*% sumXWJWX))
      
      B1    <- length(data.full$study) -
                           sum(diag(V.i %*% sumXWJX)) -
                           sum(diag(V.i %*% sumXJWX)) +
                           sum(diag(V.i %*% sumXJX%*%V.i %*% sumXWWX))
      C1    <- sumV - sum(diag(V.i %*% sumXJX))
     
      # A2 = tr(W) - tr(V*Sigma(t(X)*Wj*Jj*Wj*Xj))
      # B2 = tr(W) - tr(V*Sigma(t(X)*Wj^2*Xj))
      # C2 = Sigma(kj-p)
      
      A2   <- sumW - sum(diag(V.i %*% sumXWJWX)) 
      B2   <- sumW - sum(diag(V.i %*% sumXWWX)) 
      C2   <- length(data.full$study) - (p + 1) 
      
      # MoM estimator for omega.sq.h = A2(Qa-C1)-A1(Qe-C2) / B1A2-B2A1
      # Estimate of between-studies-wthin-cluster variance component
      omega.sq1  <- ((Qa - C1) * A2 - (Qe - C2) * A1) / (B1 * A2 - B2 * A1)
      omega.sq   <- ifelse(omega.sq1 < 0, 0, omega.sq1)
      
      # MoM estimators for tau.sq: Qe-C2/A2 - omega.sq.h(B2/A2)
      # Estimate of between-clusters variance component 
      tau.sq1  <- ((Qe - C2) / A2) - omega.sq  * (B2 / A2)
      tau.sq   <- ifelse(tau.sq1 < 0, 0, tau.sq1) 

      # Approximate inverse variance weights
      data.full$r.weights <- (1 / (data.full$var.eff.size + tau.sq + omega.sq))
  
      # Model info list for hierarchical effects
      mod_info            <- list(omega.sq = omega.sq, tau.sq = tau.sq)

    }, # End HIER
         
    CORR = { # Begin CORR
      
      W       <- diag (data.full$weights) 
      sumW    <- sum(data.full$weights) # Sum (k.j*w.j)
      Qe      <- t(data.full$e) %*% W %*% data.full$e 
      
      # The following components (denom, termA, termB, term1, term2)
      # are used in the calculation of the estimate of the residual 
      # variance component tau.sq.hat. 
      # Note: The effect of correlation on the estimates occurs entirely 
      # through the rho*term2 component.
      
      denom   <- sumW - sum(diag(solve(sumXWX) %*% sumXWJWX)) 
      termA   <- sum(diag(solve(sumXWX) %*% Matrx_WKXX)) #ZH_edit 
      termB   <- sum(diag(solve(sumXWX) %*% Matrx_wk_XJX_XX ))#ZH_edit 
      term1   <- (Qe - N + termA) / denom 
      term2   <- termB / denom 
      tau.sq1 <- term1 + rho * term2 
      tau.sq  <- ifelse(tau.sq1 < 0, 0, tau.sq1)
      df      <- N - termA - rho * (termB) 
      I.2.1   <- ((Qe - df) / Qe) * 100
      I.2     <- ifelse(I.2.1 < 0, 0, I.2.1)
      
      # Approximate inverse variance weights
      data.full$r.weights <- 1 / (data.full$k * 
                             (data.full$avg.var.eff.size + tau.sq))
      
      # Model info list for correlated effects
      mod_info            <- list(rho = rho, I.2 = I.2, tau.sq = tau.sq,
                                  term1 = term1, term2 = term2)
      
    } # End CORR 
    
  ) 
  
  } else { # Begin userweights
  
    data.full$r.weights <- data.full$userweights
    
    # Model info list for userweights
    mod_info            <- list(k = k, N = N, p = p, M = M)
    
  } # End userweights
  
  W.r.big          <- diag(data.full$r.weights)  # W
  W.r              <- by(data.full$r.weights, data.full$study, # Wj
                         function(x) diag(x, nrow = length(x)))
  sumXWX.r         <- Reduce("+", Map(function(X, W) 
                                      t(X) %*% W %*% X, 
                                      X, W.r))
  sumXWy.r         <- Reduce("+", Map(function(X, W, y) 
                                      t(X) %*% W %*% y, 
                                      X, W.r, y))
  b.r              <- solve(sumXWX.r) %*% sumXWy.r 
  data.full$pred.r <- Xreg %*% b.r
  data.full$e.r    <- cbind(data.full$effect.size) - data.full$pred.r
  data.full$e.r    <- as.numeric(data.full$e.r)
  sigma.hat.r      <- by(data.full$e.r, data.full$study, 
                         function(x) tcrossprod(x))
  
  if (!small) { # Begin small = FALSE
    
    sumXWeeWX.r  <- Reduce("+", Map(function(X, W, V) 
                                    t(X) %*% W %*% V %*% W %*% X, 
                                    X, W.r, sigma.hat.r))
    
    VR.r         <- solve(sumXWX.r) %*% sumXWeeWX.r %*% solve(sumXWX.r)  
    SE           <- sqrt(diag(VR.r)) * sqrt(N / (N - (p + 1)))
    t            <- b.r / SE
    dfs          <- N - (p + 1)
    prob         <- 2 * (1 - pt(abs(t), dfs))
    CI.L         <- b.r - qt(.975, dfs) * SE
    CI.U         <- b.r + qt(.975, dfs) * SE
    
  } else { # Begin small = TRUE

    Q             <- solve(sumXWX.r) # Q = (X'WX)^(-1)
    Q.list        <- rep(list(Q), N)
    H             <- Xreg %*% Q %*% t(Xreg) %*% W.r.big # H = X * Q * X' * W
    ImH           <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) - H
    data.full$ImH <- cbind(ImH)
    
    ImHj <- lapply(split(x = ImH,f =  as.factor(data.full$study)), 
                   function(x){matrix(x, ncol =M)})
    #ImHj          <- by(data.full$ImH, data.full$study, 
    #                    function(x) as.matrix(x))
        
    if (!user_weighting){
      
       Working_Matrx_E_j <- W.r
       Working_Matrx_E <- W.r.big

    } else { # Begin userweights
      
      V.big        <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) %*% 
        diag(data.full$avg.var.eff.size)
      v.j          <- by(data.full$avg.var.eff.size, data.full$study, 
                         function(x) diag(x, nrow = length(x)))
      Working_Matrx_E_j <- v.j
      Working_Matrx_E <- V.big
    } # End userweights
    
    
    A.MBB_inv_square <- Map(
      function (W_E, ImH_j) {
        tcrossprod(sqrt(W_E) %*% ImH_j %*%sqrt(Working_Matrx_E))
      },
      Working_Matrx_E_j, ImHj)
    
    eigenvec <- lapply(A.MBB_inv_square, function(x) eigen(x)$vectors) 
    eigenval <- lapply(A.MBB_inv_square, function(x) eigen(x)$values)
    
    A.MBB  <- Map(function (eigenvec, eigenval, k_list, W_E) 
      sqrt(W_E) %*% eigenvec %*% diag(1/sqrt(eigenval), k_list, k_list) %*% t(eigenvec) %*%sqrt(W_E),
      eigenvec, eigenval, k_list,Working_Matrx_E_j)
   
    sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
      t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
      X, W.r, A.MBB, sigma.hat.r)
    
    
    sumXWA.MBBeeA.MBBWX.r <- Reduce("+", sumXWA.MBBeeA.MBBWX.r) 
    giTemp                <- Map(function(I, A, W, X, Q)
                                 t(I) %*% A %*% W %*% X %*% Q, 
                                 ImHj, A.MBB2, W.r, X, Q.list)

   


     giTemp <- do.call(rbind,giTemp)
     gi_matrix <- lapply(X = 1:(p+1), FUN = function(i){ matrix(giTemp[,i], nrow = M)  })


     if (!user_weighting) {
          W.mat <- matrix(rep(1/sqrt(data.full$r.weights),times = N),nrow = M)
          B_matrix_half <- lapply(X = gi_matrix, FUN = function(gi_mat){ W.mat * gi_mat})
     }else{

          B_matrix_half <- gi_matrix
         # B_matrix_half <- lapply(X = gi_matrix, FUN = function(gi_mat){ solve(sqrt(V.big)) %*% gi_mat})
     }

     B_mat <- lapply(X = B_matrix_half, FUN = tcrossprod)

     B_trace_square <- sapply(X = B_mat, FUN = function(B){ (sum(diag(B)))^2})

     B_square_trace <- sapply(X = B_mat, FUN = function(B){sum(B * B)})

     dfs <- B_trace_square/B_square_trace

    
    
    VR.MBB1 <- solve(sumXWX.r) %*% sumXWA.MBBeeA.MBBWX.r %*% solve(sumXWX.r)
    VR.r    <- VR.MBB1
    SE      <- sqrt(diag(VR.r))
    t       <- b.r / SE
    prob    <- 2 * (1 - pt(abs(t), df = dfs)) 
    CI.L    <- b.r - qt(.975, dfs) * SE
    CI.U    <- b.r + qt(.975, dfs) * SE
    
  } # End small = TRUE
        
    reg_table           <- data.frame(cbind(b.r, SE, t, dfs, prob, CI.L, CI.U))
    #names(X.full)[2]    <- "intercept"
    labels              <- c(colnames(X.full[2:length(X.full)]))
    sig                 <- ifelse(prob < .01, "***", 
                           ifelse(prob > .01 & prob < .05, "**",
                           ifelse(prob > .05 & prob < .10, "*", "")))
    reg_table           <- cbind(labels, reg_table, sig)
    colnames(reg_table) <- c("labels", "b.r", "SE", "t", "dfs", "prob", "CI.L", 
                             "CI.U", "sig")
  
   if (!small) { # Begin small = FALSE
 
      mod_label_sm   <- ""
      mod_notice     <- ""
      
   } else { # Begin small = TRUE
     
      mod_label_sm  <- "with Small-Sample Corrections"            
      mod_notice    <- "Note: If df < 4, do not trust the results"
    
   } # End small = TRUE
       
  
   if (!user_weighting) {
  
     switch(modelweights,
           
       HIER = { # Begin HIER
             
          mod_label <- c("RVE: Hierarchical Effects Model", mod_label_sm)

        }, # End HIER
         
        CORR = { # Begin CORR
             
          mod_label <- c("RVE: Correlated Effects Model", mod_label_sm)
             
        } # End CORR
           
     ) 
         
   } else { # Begin userweights

     mod_label <- c("RVE: User Specified Weights", mod_label_sm)
         
   } # End userweights

  res <- list(data.full = data.full, X.full = X.full, reg_table = reg_table, 
              mod_label = mod_label, mod_notice = mod_notice, modelweights = 
              modelweights, mod_info = mod_info, user_weighting = 
              user_weighting, ml = ml, cl = cl, N = N, M = M, k = k, 
              k_list = k_list, p = p, X = X, y = y, Xreg = Xreg, b.r = b.r, 
              VR.r = VR.r, dfs = dfs, small = small, data = data, labels = 
              labels, study_orig_id = study_orig_id)
               
  class(res) <- "robu"
  res
}
