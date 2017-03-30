#' Predict method for robumeta model fits.
#' 
#' \code{predict.robu} obtains predictions and prediction confidence intervals of a fitted robumeta model object for a constraints vector.
#' 
#' @param  object A fitted robumeta model object.
#' @param  constraints_vector A contrast vector to test. It the vector contains NAs, \code{predict.robu} will remove the corresponding covariates from the original data, and refit a new robumeta model.The prediction and prediction interval will be estimated based on that new model. 
#' @param  level Tolerance/confidence level.
#' @return \code{predict.robu }   produces a vector of prediction and confidence interval with element names \code{prediction, se, df, t_score, lowerBound,} and \code{upperBound}.
#' @param \code{prediction}   the predicted value
#' @param \code{se}  standard error of predicted mean
#' @param \code{t_score} The t-statistic calculated based on the predicted mean.
#' @param \code{df}  The small sample corrected degrees of freedom of the distribution of the t-statistic. 
#' @param \code{lowerBound} lower bound of the confidence interval for the predicted mean.
#' @param \code{upperBound} upper bound of the confidence interval for the predicted mean.
#' @examples
#' library(robumeta)
#' library(clubSandwich)
#' 
#' robu_mod <- robu(LOR1 ~ study_design + duration + service_hrs, 
#'                  data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, 
#'                  modelweights = "HIER",small = T)
#' predict(object = robu_mod,contrast_vector = c(1,1,0,0,0),level = 0.95)




predict.robu <- function(object, contrast_vector, level = 0.95){
  
  if (!inherits(object, "robu")) 
    warning("calling predict.robu(<fake-lm-object>) ...")
  
  
  original_data <- object$Xreg
  combined_data <- rbind(object$Xreg, contrast_vector)
  
  #remove columns with NAs.
  combined_data_noNA <- combined_data[,!colSums(is.na(combined_data))]
  ####
  #Check if combined_data_noNA is NULL
  ####
  
  #construct new formula and new data frame.
  original_data_noNA <- combined_data_noNA[- nrow(combined_data_noNA),]
  newformula <- as.formula(paste0("effect.size ~",paste0(colnames(original_data_noNA),collapse = " + "), "-1"))
  
  data.full <- as.data.frame(cbind(original_data_noNA,
                                   effect.size = object$data.full$effect.size, 
                                   var.eff.size=object$data.full$var.eff.size,
                                   study = object$data.full$study))
  
  #Build new robu model.
  new_modelcall <- object$cl
  new_modelcall$formula <- newformula
  new_modelcall$data <- quote(data.full)
  new_modelcall$studynum <- quote(study)
  new_modelcall$var.eff.size <- quote(var.eff.size)
  
  new_robu_model <- eval(new_modelcall)
  
  #Calculate outcome prediction
  test_vec <- t(combined_data_noNA[nrow(combined_data_noNA),])
  coef_est <- new_robu_model$b.r #estimated coefficients
  g_hat <- test_vec %*% coef_est
  
  #calculate F-test with small sample correction
  Test_result <- Wald_test(new_robu_model, constraints = test_vec, vcov="CR2")
  
  #Retrive t statistics and confidence interval. 
  df <- Test_result$df 
  Fvalue <- Test_result$F
  
  tcrit <- qt(1-(1-level)/2,df)
  SE_Yhat <- g_hat/sqrt(Fvalue)
  LB <- g_hat - tcrit*SE_Yhat
  UB <- g_hat + tcrit*SE_Yhat
  
  res <- c(prediction = g_hat,
           se = SE_Yhat,
           df = df,
           t_score = sqrt(Fvalue),
           lowerBound = LB,
           upperBound = UB
  )
  return(res)
}
