#' Fitting Robust Variance Meta-Regression Models
#' 
#' \code{robu} is used to meta-regression models using robust variance
#' estimation (RVE) methods. \code{robu} can be used to estimate correlated and
#' hierarchical effects models using the original (Hedges, Tipton and Johnson,
#' 2010) and small-sample corrected (Tipton, 2013) RVE methods. In addition,
#' \code{robu} contains options for fitting these models using user-specified
#' weighting schemes (see the Appendix of Tipton (2013) for a discussion of
#' non- efficient weights in RVE).
#' 
#' 
#' @aliases robu CORR HIER USER
#' @param formula An object of class \code{"formula"}. A typical
#' meta-regression formula will look similar to \code{y ~ x1 + x2...}, where
#' \code{y} is a vector of effect sizes and \code{x1 + x2...} are (optional)
#' user-specified covariates. An intercept only model can be specified with
#' \code{y ~ 1} and the intercept can be ommitted as follows \code{y ~ -1
#' +...}.
#' @param data A data frame, list or environment or an object coercible by
#' as.data.frame to a data frame.
#' @param studynum A vector of study numbers to be used in model fitting.
#' \code{studynum} must be a numeric or factor variable that uniquely
#' identifies each study.
#' @param var.eff.size A vector of user-calculated effect-size variances.
#' @param rho User-specified within-study effect-size correlation used to fit
#' correlated (\code{modelweights = "CORR"}) effects meta-regression models.
#' The value of \code{rho} must be between 0 and 1. The default value for
#' \code{rho} is 0.8.  \code{rho} is not specified for hierarchical
#' (\code{modelweights = "HIER"}) effects models.
#' @param modelweights User-specified model weighting scheme.  The two two
#' avialable options are \code{modelweights = "CORR"} and \code{modelweights =
#' "HIER"}.  The default is \code{"CORR"}. See Hedges, Tipton and Johnson
#' (2010) and Tipton (2013) for extended explanations of each weighting scheme.
#' @param userweights A vector of user-specified weights if non-efficient
#' weights are of interest. Users interested in non-efficient weights should
#' see the Appendix of Tipton (2013) for a discussion of the role of
#' non-efficient weights in RVE).
#' @param small \code{small = TRUE} is used to fit the meta-regression models
#' with the small- sample corrections for both the residuals and degrees of
#' freedom, as detailed in Tipton (2013). Users wishing to use the original RVE
#' estimator must specify \code{small = FALSE} as the corrected estimator is
#' the default option.
#' @param ...  Additional arguments to be passed to the fitting function.
#' @return
#' 
#' \item{output}{ A data frame containing some combination of the robust
#' coefficient names and values, standard errors, t-test value, confidence
#' intervals, degrees of freedom and statistical significance.  }
#' 
#' \item{n}{The number of studies in the sample \code{n}}.
#' 
#' \item{k}{The number of effect sizes in the sample \code{k}}.
#' 
#' \item{k descriptives}{the minimum \code{min.k}, mean \code{mean.k}, median
#' \code{median .k}, and maximum \code{max.k} number of effect sizes per study.
#' }
#' 
#' \item{tau.sq.}{ \code{tau.sq} is the between study variance component in the
#' correlated effects meta-regression model and the between-cluster variance
#' component in the hierarchical effects model. \code{tau.sq} is calculated
#' using the method-of-moments estimator provided in Hedges, Tipton, and
#' Johnson (2010).  For the correlated effects model the method-of-moments
#' estimar depends on the user-specified value of rho.  }
#' 
#' \item{omega.sq.}{ \code{omega.sq} is the between-studies-within-cluster
#' variance component for the hierarchical effects meta-regression model.
#' \code{omega.sq} is calculated using the method-of-moments estimator provided
#' in Hedges, Tipton, and Johnson (2010) erratum.  }
#' 
#' \item{I.2}{ \code{I.2} is a test statistics used to quantify the amount of
#' variability in effect size estimates due to effect size heterogeneity as
#' opposed to random variation.  }
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
#' # Small-Sample Corrections - Hierarchical Dependence Model
#' HierModSm <-  robu(formula = effectsize ~ binge + followup + sreport
#'                    + age, data = hierdat, studynum = studyid, 
#'                    var.eff.size = var, modelweights = "HIER", small = TRUE)
#' 
#' print(HierModSm) # Output results
#' 
#' @export
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
                                        stats::model.matrix(formula, mf),
                                        studynum = mf[["(studynum)"]],
                                        var.eff.size = mf[["(var.eff.size)"]])
    
    X.full.names           <- names(dframe)[-match(c("effect.size", 
                                                     "studynum", 
                                                     "var.eff.size"), 
                                                   names(dframe))] 
    
  } else { # Begin userweights
    
    dframe                 <- data.frame(effect.size = mf[,1],
                                        stats::model.matrix(formula, mf), 
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
  dframe$avg.var.eff.size  <- stats::ave(dframe$var.eff.size, dframe$study)
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
      data.full$r.weights <- (1 / (as.vector(data.full$var.eff.size) + 
                                     as.vector(tau.sq) + 
                                     as.vector(omega.sq)))
  
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
      data.full$r.weights <- 1 / (as.vector(data.full$k) * 
                             (as.vector(data.full$avg.var.eff.size) + 
                                as.vector(tau.sq)))
      
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
    prob         <- 2 * (1 - stats::pt(abs(t), dfs))
    CI.L         <- b.r - stats::qt(.975, dfs) * SE
    CI.U         <- b.r + stats::qt(.975, dfs) * SE
    
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
      
      Working_Matrx_E <- diag(1/data.full$r.weights)  #1/W
      Working_Matrx_E_j <- by(data.full$r.weights, data.full$study, # Wj
                              function(x) diag(1/x, nrow = length(x))) #1/W_j

      
      switch(modelweights, 
             HIER = {
               # Inside Matrix = E_j^0.5 * ImH_j *E * t(ImH_j) * E_j^0.5
               # In this case, the formula can be simplified to 
               # Inside Matrix = E_j^0.5 * ImH_jj * E_j^1.5
               InsideMatrx_list <-   Map(
                 function (W_E_j, ImH_jj) {
                   sqrt(W_E_j) %*% ImH_jj %*% (W_E_j^1.5)
                 },
                 Working_Matrx_E_j, ImHii)
               
               eigenres_list <- lapply(InsideMatrx_list, function(x) eigen(x))
               eigenval_list <- lapply(eigenres_list, function(x) x$values)
               eigenvec_list <- lapply(eigenres_list, function(x) x$vectors)
               
               A.MBB  <- Map(function (eigenvec, eigenval, k, W_E_j) {
                                eigenval_InvSqrt <- ifelse(eigenval< 10^-10, 0, 1/sqrt(eigenval)) # Pseudo_Inverse
                                sqrt(W_E_j) %*% eigenvec %*% diag(eigenval_InvSqrt, k, k) %*% t(eigenvec) %*%sqrt(W_E_j) # Pseudo_Inverse
                             },
                             eigenvec_list, 
                             eigenval_list, 
                             k_list,
                             Working_Matrx_E_j)
               
             },
             
             CORR = {
               # In this case, the formula can be simplified to 
               # A_MBB =  ImH_jj ^ (-0.5)
               
               eigenres_list <- lapply(ImHii, function(x) eigen(x))
               eigenval_list <- lapply(eigenres_list, function(x) x$values)
               eigenvec_list <- lapply(eigenres_list, function(x) x$vectors)
               
               A.MBB  <- Map(function (eigenvec, eigenval, k, W_E_j) {
                               eigenval_InvSqrt <- ifelse(eigenval< 10^-10, 0, 1/sqrt(eigenval)) # Pseudo_Inverse
                               eigenvec %*% diag(eigenval_InvSqrt, k, k) %*% t(eigenvec) # Pseudo_Inverse
                            },
                            eigenvec_list, 
                            eigenval_list, 
                            k_list,
                            Working_Matrx_E_j)
             })
      
      
    } else { # Begin userweights
      
      V.big        <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) %*% 
        diag(data.full$avg.var.eff.size)
      v.j          <- by(data.full$avg.var.eff.size, data.full$study, 
                         function(x) diag(x, nrow = length(x)))
      v.j.sqrt_list     <- lapply(v.j, function (x) sqrt(x))
      
      Working_Matrx_E_j <- v.j
      Working_Matrx_E <- V.big
      
      InsideMatrx_list <-   Map(
        function (ImH_j) {
           ImH_j %*% Working_Matrx_E %*% t(ImH_j)
        },
        ImHj)
      
      eigenres_list <- lapply(InsideMatrx_list, function(x) eigen(x))
      eigenval_list <- lapply(eigenres_list, function(x) x$values)
      eigenvec_list <- lapply(eigenres_list, function(x) x$vectors)
      
      A.MBB  <- Map(function (eigenvec, eigenval, k, v.j.sqrt) {
                      eigenval_InvSqrt <- ifelse(eigenval< 10^-10, 0, 1/sqrt(eigenval)) # Pseudo_Inverse
                      v.j.sqrt %*% eigenvec %*% diag(eigenval_InvSqrt, k, k) %*% t(eigenvec)  # Pseudo_Inverse
                    },
                    eigenvec_list, 
                    eigenval_list, 
                    k_list,
                    v.j.sqrt_list)
      
      
    } # End userweights
    
    sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
      t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
      X, W.r, A.MBB, sigma.hat.r)
    
    
    
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
    prob    <- 2 * (1 - stats::pt(abs(t), df = dfs)) 
    CI.L    <- b.r - stats::qt(.975, dfs) * SE
    CI.U    <- b.r + stats::qt(.975, dfs) * SE
    
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
