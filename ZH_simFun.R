library(robumeta)
library(parallel)
r_corrSMD <- function(k, n, r, tau_sq) {
  alpha <- rnorm(n = 1, sd = sqrt(tau_sq + 2 * r / n))
  
  e <- rnorm(n = k, sd = sqrt(2 * (1 - r) / n))
  cov_mat <- as.matrix(rWishart(1, df = 2 * (n - 1), Sigma = r + diag(1 - r, nrow = k))[,,1])
  sigma_sq <- diag(cov_mat) / (2 * n - 2)
  (alpha + e) / sqrt(sigma_sq)
}

# generates m clusters, each with k correlated standardized mean difference estimates

r_cluster_SMDs <- function(iterations, m, kVals, nVals, r, Isq) {
  kVals <- rep(kVals, length.out = m)
  nVals <- rep(nVals, length.out = m)
  tau_sq <- 2 * (sum(kVals / nVals) / sum(kVals)) * Isq / (1 - Isq)
  replicate(iterations, 
            unlist(mapply(r_corrSMD, k = kVals, n = nVals, 
                          MoreArgs = list(r = r, tau_sq = tau_sq), 
                          SIMPLIFY = FALSE)))
}

data_gen <- function(m, kVals, nVals, r, Isq, X, seed = NULL) {

  # expand kVals and nVals if necessary
  kVals <- rep(kVals, length.out = m)
  nVals <- rep(nVals, length.out = m)
  
  # subset the design matrix
  choose_studies <- (0:(10 * m - 1) %% NROW(X)) + 1
  X <- X[choose_studies,]
  choose_effects <- unlist(lapply(kVals, function(k) 1:10 <= k))
  X <- X[choose_effects,]
  
  # simulate data
  if(!is.null(seed)) set.seed(seed)
  y_mat <- as.vector(r_cluster_SMDs(iterations =1, m, kVals, nVals, r, Isq))
  n <- rep(nVals, kVals)
  cluster <- factor(rep(1:m, kVals))

  Vars <- 2 / n + y_mat / (4 * n - 4)
  weights <- as.vector(1 / tapply(Vars, cluster, sum)[cluster])
  
  as.data.frame(cbind(y_mat,X, Vars, weights,cluster))
}

#Added 6 points distribution

robu_Wild    <-function(formula, data, studynum,var.eff.size, userweights,
                        modelweights = c("CORR", "HIER"), rho = 0.8, 
                        small = TRUE, Replication_num =1000, coef_index=1, six_pts = F, ...) {
  
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
    
    diag_one      <- by(rep(1, M), X.full$study, 
                        function(x) diag(x, nrow = length(x)))
    ImHii         <- Map(function(X, Q, W, D) 
      D - X %*% Q %*% t(X) %*% W,
      X, Q.list, W.r, diag_one)
    
    if (!user_weighting){
      
      switch(modelweights, 
             
             HIER = { # Begin HIER
               
               # inside = Wj^(-1/2) * (I-Hjj) * Wj^(-3/2)
               inside   <- Map(function(W, I) 
                 solve(sqrt(W)) %*% I %*% solve(sqrt(W)^3),
                 W.r, ImHii)
               I        <- inside
               eigenvec <- lapply(inside, function(x) eigen(x)$vectors) 
               eigenval <- lapply(inside, function(x) eigen(x)$values)
               
             }, # End HIER
             
             CORR = { # Begin CORR
               
               eigenvec <- lapply(ImHii, function(x) eigen(x)$vectors) 
               eigenval <- lapply(ImHii, function(x) eigen(x)$values)
               I        <- ImHii
               
             } # End CORR
             
      ) 
      
    } else { # Begin userweights
      
      V.big        <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) %*% 
        diag(data.full$avg.var.eff.size)
      V.big.list   <- rep(list(V.big), N)
      v.j          <- by(data.full$avg.var.eff.size, data.full$study, 
                         function(x) diag(x, nrow = length(x)))
      v.j.sqrt     <- lapply(v.j, function (x) sqrt(x))
      inside       <- Map(function(V, I) 
        I  %*% V %*% t(I),
        V.big.list, ImHj)
      eigenvec     <- lapply(inside, function(x) eigen(x)$vectors)
      eigenval     <- lapply(inside, function(x) eigen(x)$values)
      I            <- inside
      
    } # End userweights
    
    eigenval_inv_sqrt <- lapply(X = eigenval, FUN = function(x){ifelse(x<10^-10, 0, 1/sqrt(x))})
    
    A.MBB  <- Map(function (eigenvec, eigenval_inv_sqrt, k_list) 
      eigenvec %*% 
        diag(eigenval_inv_sqrt, k_list, k_list) %*% t(eigenvec),
      eigenvec, eigenval_inv_sqrt, k_list)
    A.MBB1 <- Map(function(K, A, I) 
      if (K > 1) A else matrix(sqrt(solve(I))), 
      k_list, A.MBB, I)
    
    if (!user_weighting){
      
      switch(modelweights, 
             
             HIER = { # Begin HIER
               
               A.MBB2                <- Map(function(W, A) 
                 solve(sqrt(W)) %*% A %*% solve(sqrt(W)),
                 W.r, A.MBB1) 
               sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
                 t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
                 X, W.r, A.MBB2, sigma.hat.r)
             }, # End HIER
             
             CORR = { # Begin CORR
               
               A.MBB2                <- A.MBB1
               sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
                 t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
                 X, W.r, A.MBB2, sigma.hat.r)
             } # End CORR
             
      ) 
      
    } else { # Begin userweights
      
      A.MBB2                <- Map(function(V, A) 
        V %*% A,
        v.j.sqrt, A.MBB1) 
      sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
        t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
        X, W.r, A.MBB2, sigma.hat.r)
    } # End userweights
    
    
    
    
    sumXWA.MBBeeA.MBBWX.r <- Reduce("+", sumXWA.MBBeeA.MBBWX.r) 
    giTemp                <- Map(function(I, A, W, X, Q)
      t(I) %*% A %*% W %*% X %*% Q, 
      ImHj, A.MBB, W.r, X, Q.list)
    
    
    
    
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
  
  
  
  ##########################################################################################
  
  ###Data needed to generate bootstrap samples
  Study_id <- X.full[,1]
  X_matrx <- Xreg
  Y_matrx <- data.full$effect.size
  Y_hat_matrx <- data.full$pred.r
  coef_est <- b.r
  residual_matrx <- data.full$e.r
  residual_by_group <- by(data = residual_matrx, INDICES = Study_id,FUN = function(x){x})
  group_num <- length(residual_by_group)
  outcome_name <- as.character(ml[[2]])
  
  
  ###Function to generate bootstrap samples
  
  Get_Bootdata <- function(){
    if(six_pts == F){
      
      RandomVector_list <-as.list(rbinom(n = group_num,size = 1,0.5)*2 -1)
      
    }else{
      
      Sixth_points <- c(- sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2))
      RandomVector_list <- Sixth_points[ceiling(runif(n =group_num) *6)]
    }
    
    
    Boot_residual_by_group <- Map(f = function(residual, randomVec){residual*randomVec}, 
                                  residual_by_group, RandomVector_list)
    
    Boot_residual <- unlist(Boot_residual_by_group)
    
    Boot_b <- b.r
    Boot_b[coef_index] <- 0
    Boot_Y <- Xreg %*% Boot_b + Boot_residual
    Boot_data <- data
    Boot_data[outcome_name] <-Boot_Y
    return(Boot_data)
  }
  #Get Bootstrap data list
  Boot_data_list <- replicate(n = Replication_num, expr = Get_Bootdata(),simplify = F)
  
  #Function to calculate bootstrap result
  Get_BootResult <- function(Boot_data){
    cl[[1]] <-  quote(robu)
    cl$data <-quote(Boot_data)
    new_result <- eval(cl)
    return(new_result)
  }
  #Get bootstrap result list.
  BootResult_list <- lapply(X = Boot_data_list,FUN = Get_BootResult)
  
  Boot_t_val_list <- lapply(X = BootResult_list, FUN = function(x){x$reg_table["t"][coef_index,]})
  
  Boot_t_val_vec <-unlist(Boot_t_val_list)
  
  p_val_Wild <- sum(abs(Boot_t_val_vec) >abs(t[coef_index,]))/Replication_num
  
  res <- list(data.full = data.full, X.full = X.full, reg_table = reg_table, 
              mod_label = mod_label, mod_notice = mod_notice, modelweights = 
                modelweights, mod_info = mod_info, user_weighting = 
                user_weighting, ml = ml, cl = cl, N = N, M = M, k = k, 
              k_list = k_list, p = p, X = X, y = y, Xreg = Xreg, b.r = b.r, 
              VR.r = VR.r, dfs = dfs, small = small, data = data, labels = 
                labels, study_orig_id = study_orig_id)
  
  class(res) <- "robu"
  res
  return(list(CRVE_result = res, 
              Wild_boot_t_val_vec = Boot_t_val_vec,
              Wild_boot_p_val = p_val_Wild))
}



# design matrix
X1a <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
X2a <- c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0)
X3a <- c(-1.00270563046404, -1.00270563046404, -1.0027056304 6404, -1.00270563046404, -1.00270563046404, -1.00270563046404, -1.00270563046404, -1.00270563046404, -1.00270563046404, -1.00270563046404, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.726181921265885, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, -0.318437830881828, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.392635406184748, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.817420360941905, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.568773376668874, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, 0.11943092553188, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.495922972116935, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.420088914385043, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, -0.127873519955423, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, 1.48853945847427, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.110662703146518, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, -0.998292141102654, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.878111535736483, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.86555970026061, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 0.0675938578282997, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, 1.0784881856417, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -1.08667333136972, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, -0.886201982800008, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885, 0.509457246647885)
X4a <- c(-0.294746498126055, 0.061176475934361, 1.4878607081176, 0.0943764314860022, 2.53605642667157, 0.788042145116689, 0.429165296817196, 2.25382626613944, 0.458196996590399, 1.68670754205783, 1.4719089241765, 1.12505222839121, -0.00971640776314253, -0.0110725162292895, 0.349798471404283, 2.73141224814032, 0.0547264680075376, 1.53496305022157, 0.313043866127804, 0.335604172060825, 1.59199379738799, 0.984214065955563, 0.0867736362637365, 0.378356921226552, 1.07320973481322, 0.557733365051443, 1.19042610779747, 1.22834962431115, 0.515395722330571, 0.431435141611246, 0.156248314818041, 1.15524529166603, 0.861921460763251, 1.73161732961643, 0.960351621298908, 1.73914743336814, 1.95928922531143, 0.932181838572626, 0.077484215175638, 1.25457443203552, 1.23798951196419, 2.49291479732462, 0.616961826790627, 0.957789283694637, -0.151601159039986, 1.04307651995346, 0.703813784300677, 1.76864878245583, 0.832235630358274, 0.680355616043406, 0.153440997762995, 0.319863181727738, 1.24794695313152, 2.1408648656992, 0.237579384114682, 1.17145563851091, 1.6243942111374, 1.57281579466851, 1.4420228413386, 1.04918031557555, 1.12180196669322, 0.640214593376593, 0.168340772133958, 2.0334736556839, 1.08265707390848, -0.106770695307287, 0.423182543441357, 1.24010008534131, 1.72387375857721, 1.57294369985468, 0.471584119924622, 0.307182583534007, -0.530155069482054, 0.0406314350954493, 1.55641658316176, 0.749768294018018, 0.438987803408893, 1.89875600995333, 1.14384706937116, 0.538063492671569, -0.0675479630064653, 0.723026683470194, 0.251538555166737, 0.961884404406095, 1.02065275320635, 1.71024843996069, 0.982120102514053, 0.0698976121433963, -0.465323204793089, 1.29531563806573, 0.106808176086012, 2.75662644278713, 0.140181280572896, 1.57401361905011, 1.99679075172621, 2.14763240038146, 0.420542510755185, 1.70487889977921, 0.665415437082732, -0.0368977193623724, 0.209934427773169, 2.66767085034218, 0.826386967608597, 1.80450550896599, 0.942684448823455, 0.123515536632673, 0.926274586939407, -0.112105382000053, 2.22278366565292, 0.851647218387409, 0.241209618824241, 1.49914087707787, 1.75382623887113, 0.0422482169931735, 1.37191829287146, 0.687126081945874, 2.28891718761916, 1.33002191732935, 0.228269982635948, 0.886001816394112, 1.34886445477324, 1.19715517027587, 1.49737117826425, 0.27421129794853, 1.15397589852662, 0.485795220656077, -0.100740861035327, 0.38264822146087, 1.25562999297498, 0.841031157080616, 1.72253611062841, 0.942114308281235, 1.46360593749851, 1.96290317258425, 1.05780338815015, 0.153118026188228, 0.407299766060794, 1.41732321217371, 0.918504045677449, 0.778639147908423, 1.6434044922428, 1.05381514948917, 0.37251438830644, 0.496001837337652, 1.50827452839003, 0.705146498336872, -0.289159686407903, 0.158578189186895, 0.855755001406733, 1.00698888767942, 0.148171875072054, -0.356965624444167, 1.31587909249981, 1.28418619783321, 1.48885123553784, 1.53099373325098, 0.687417543016071, -0.0263748905841946, 0.640814979369598, 0.911041982910652, 1.69266472745245, 1.02406889748813, 0.895270095628613, 0.460879309953268, 1.39336796123866, 1.29982146857978, 1.18983761538061, 2.72823462053557, 1.01858398875457, 1.38510897022177, 0.496061069861388, 1.15274101492898, 0.427656125283214, 1.62015402915792, 0.0234868640236752, 0.342868664237198, 2.80444499279049, 0.762453910596334, 0.803700772026355, 1.03348767458784, 1.27432560491757, 2.09737220486058, 0.363962165530275, -0.224908254027832, 1.72375012622472, 0.829744442556771, 1.42544149931656, 1.46412764849452, 1.09717762836835, 1.29932720295109, 1.20354170667481, 1.42240970089826, 1.5086985967493, 1.75179914337685, 1.91579458543882, 0.127984913087748, 0.407648716334308, 1.20161847912127, 0.715919793964415, 0.409403555686791)
X5.10 <- c(0.71, 7.81, 0.33, 1.13, 1.34, 0.5, 0.59, 2.72, 2.87, 3.24, 0.68, 1.92, 3.77, 6.32, 27.83, 2.07, 0.49, 0.51, 1.02, 5.83, 3.94, 1.38, 0.73, 1.79, 0.13, 0.62, 0.4, 0.6, 1.11, 0.47, 0.54, 1.66, 1.47, 0.11, 1.47, 0.49, 0.09, 2.44, 0.35, 0.31, 0.75, 1.39, 2.72, 0.23, 0.66, 1.94, 0.21, 1.18, 6.27, 0.3, 0.08, 0.64, 4.7, 0.76, 1.16, 0.22, 5.22, 0.73, 2.15, 5.34, 1.31, 0.8, 0.36, 0.41, 0.15, 1.34, 2.48, 0.13, 1.76, 2.13, 0.41, 0.13, 6.07, 0.66, 0.96, 3.13, 1.98, 0.76, 1.77, 1.07, 3.27, 0.58, 0.51, 4.13, 0.76, 0.42, 0.61, 0.95, 2.87, 5.57, 5.56, 0.07, 0.22, 0.49, 1.98, 6.11, 0.37, 0.26, 1.35, 1.4, 1.41
           ,1.77, 1.38, 1.03, 0.57, 0.76, 2.11, 1.9, 1.05, 0.67, 0.59, 0.68, 0.41, 2.59, 1.99, 0.76, 0.33, 3.35, 0.32, 1.08, 0.52, 0.36, 0.51, 4.26, 1.2, 0.7, 0.11, 0.12, 0.13, 1.23, 2.52, 1.96, 0.28, 0.43, 1.06, 1.57, 2.54, 2.93, 5.27, 0.53, 3.66, 0.38, 0.54, 0.78, 1.17, 2.14, 1.21, 0.68, 0.51, 0.48, 1.34, 10.91, 0.97, 0.5, 0.91, 0.17, 0.39, 0.49, 5.48, 0.16, 0.77, 2.75, 3.45, 1.59, 0.37, 3.65, 1.69, 0.48, 0.48, 1.85, 1.16, 1.91, 0.23, 0.23, 0.92, 2.18, 1.24, 1.12, 0.59, 0.52, 0.44, 0.65, 1.78, 0.37, 1.45, 2.11, 0.63, 0.39, 1.72, 0.19, 0.36, 2.1, 0.58, 0.95, 1.19, 0.8, 0.53, 2.12, 1.09, 4.36)

# design matrix. 

X1b<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
X1c<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
X1d<- c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
X1e<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)
X1f<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

#X_full <- cbind(intercept=1,X1a,X2a,X3a,X4a, X5.10) #200 total effect sizes in 20 studies each with 10
X_full <- cbind(X1a,X2a,X3a,X4a, X5.10, X1b, X1c, X1d, X1e, X1f) 


Simdata_list<- lapply(X = 1:200, FUN = function(x){   data_gen(m = 20,
                                                     kVals = 10,
                                                     nVals = 30,
                                                     r = .5,
                                                     Isq = 1/3,
                                                     X = X_full,
                                                     seed = NULL)})

