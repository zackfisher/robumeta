#' Prediction method for a robumeta object.
#' 
#' \code{predict.robu} produces the predicted mean and confidence interval of a fitted robumeta model object given a prediction vector.
#' 
#' @param  object A fitted robumeta model object.
#' @param  pred.vector A prediction vector containing the new covariate values.
#' @param  level Confidence level.
#' @param  ... Additional arguments to predict.
#'
#' @return \code{prediction} the predicted value based on the prediction vector.
#' @return \code{se} The standard error for the  predicted mean.
#' @return \code{t} The t-statistic calculated based on the predicted mean.
#' @return \code{df} The small sample corrected degrees of freedom of the distribution of the t-statistic. 
#' @return \code{lower} The lower bound of the confidence interval for the predicted mean.
#' @return \code{upper} The upper bound of the confidence interval for the predicted mean.
#'
#' @details 
#' \itemize{
#'    \item{\code{intercept}} {
#'      If an intercept is included in the robumeta model, 
#'      the first element should always be 1, representing the intercept, 
#'      followed by the covariate values in appropriate order. 
#'      If the robumeta model does not have an intercept, the prediction 
#'      vector should begin with the first covariate value.
#'    }
#'    \item{\code{variable}} {
#'      For continuous variables, use the variable value as the 
#'      corresponding element value in \code{ pred.vector}. 
#'      For a categorical variable the original variable value should 
#'      be transformed to match the coding system used in the robumeta 
#'      model (e.g. dummy coding, deviation coding, etc.).
#'    }
#'    \item{\code{NA}} {
#'      If the vector contains NAs, \code{predict.robu} will remove the 
#'      corresponding covariates from the original data, and refit a new 
#'      robumeta model. The prediction and confidence interval will be 
#'      estimated based on the new model. 
#'    }
#' }
#'
#' \preformatted{
#'   robu_mod <- robu(LOR1 ~ study_design + duration + service_hrs, 
#'                    data = dropoutPrevention, 
#'                    studynum = studyID, 
#'                    var.eff.size = varLOR, 
#'                    modelweights = "HIER",
#'                    small = TRUE)
#' }
#' 
#' In this robumeta model, the first covariate is a categorical variable 
#' that contains three levels: "Matched" (33 percent, dummy code: 00), 
#' "Randomized"(24 percent, 01) and "non-match non-randomized"(43 percent, 
#' 10). The corresponding prediction vector begins with 1 (intercept), 
#' and followed by 0, 0, the dummy code for "Matched". The last two elements
#' are 38 and 5, the values for duration and sevice_hrs.
#'  
#' \preformatted{                                          
#'   predict(object = robu_mod, pred.vector = c(1,0,0,38,5),level = 0.95)
#' }
#' 
#' If we do not know the value of duration, the prediction vector should 
#' be c(1,0,0,NA,5). predict.robu() will refit a new model without the 
#' covariate duration, and the prediction will be based on it.
#' 
#' \preformatted{
#'   predict(object = robu_mod, pred.vector = c(1,0,0,NA,5),level = 0.95)
#' }
#' 
#' @export 
predict.robu <- function(object, pred.vector, level = 0.95, ...){
  
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    stop("clubSandwich needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!inherits(object, "robu")) 
    warning("calling predict.robu(<fake-lm-object>) ...")

  original_data <- object$Xreg
  combined_data <- rbind(object$Xreg, pred.vector)
  
  #remove columns with NAs.
  combined_data_noNA <- combined_data[,!colSums(is.na(combined_data))]
  ####
  #Check if combined_data_noNA is NULL
  ####
  
  #construct new formula and new data frame.
  original_data_noNA <- combined_data_noNA[- nrow(combined_data_noNA),]
  newformula <- stats::as.formula(paste0("effect.size ~",paste0(colnames(original_data_noNA),collapse = " + "), "-1"))
  
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
  Test_result <- clubSandwich::Wald_test(new_robu_model, constraints = test_vec, vcov="CR2")
  
  #Retrieve t statistics and confidence interval. 
  df <- Test_result$df 
  Fvalue <- Test_result$F
  
  tcrit <- stats::qt(1-(1-level)/2,df)
  SE_Yhat <- g_hat/sqrt(Fvalue)
  LB <- g_hat - tcrit*SE_Yhat
  UB <- g_hat + tcrit*SE_Yhat
  
  res <- c(prediction = g_hat,
           se = SE_Yhat,
           df = df,
           t = sqrt(Fvalue),
           lower = LB,
           upper = UB
  )
  return(res)
}
