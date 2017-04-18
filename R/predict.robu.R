#' Predict method for robumeta model fits.
#' 
#' \code{predict.robu} obtains the predicted outcome mean and confidence interval of a fitted robumeta model object given a prediction vector.
#' 
#' @param  object A fitted robumeta model object.
#' @param  pred_vector A prediction vector for hypothesis test. This contains the covariate values we have. The conditional mean value, as well as its confidence interval, of the outcome (dependent) variable will be predicted given the values.
#' 
#'  INTERCEPT: If intercept is included in the robumeta model, the first element should always be 1, which represents the intercept, and followed by the covariate values in order. If the robumeta model does non have an intercept, the prediction vector should begin with the first covariate value.
#' 
#'  VARIABLE TYPE: For a continuous variable, we just use the variable value as the corresponding element value in the vector. For a categorical variable, the original variable value should be transformed by the coding system identical to the one used in the robumeta model.(Dummy coding, deviation coding, etc.)
#'  
#'  NA VALUES: If the vector contains NAs, \code{predict.robu} will remove the corresponding covariates from the original data, and refit a new robumeta model.The prediction and confidence interval will be estimated based on that new model. 
#'  
#' @param  level Tolerance/confidence level.
#'
#' @return \code{predict.robu}   produces a vector of prediction and confidence interval with element names \code{prediction, se, df, t_score, lowerBound,} and \code{upperBound}.
#' @return \code{prediction:}   the predicted value based on the prediction vector.
#' @return \code{se:}  standard error of predicted mean
#' @return \code{t_score:} The t-statistic calculated based on the predicted mean.
#' @return \code{df:}  The small sample corrected degrees of freedom of the distribution of the t-statistic. 
#' @return \code{lowerBound:} lower bound of the confidence interval for the predicted mean.
#' @return \code{upperBound:} upper bound of the confidence interval for the predicted mean.
#' @examples
#' require(robumeta)
#' require(clubSandwich)
#'
#' ## Fit a robumeta model by robu() function. 
#' 
#' 
#' 
#' robu_mod <- robu(LOR1 ~ study_design + duration + service_hrs, 
#'                  data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, 
#'                  modelweights = "HIER",small = T)
#' 
#' ## Construct the prediction vector for the predicted mean given 
#' ## study_design = "Matched", duration = 38 and service_rs = 5.
#' ##
#' ## In this robumeta model, the first covariate is a categorical variable that contains three levels:
#' ## "Matched" (33%, dummy code: 00), "Randomized"(24%, 01) and "non-match non-randomized"(43%, 10).
#' ##
#' ## The corresponding prediction vector begins with 1 (intercept), and followed by 0 0, the dummy code for "Mateched".
#' ## The last two elements are 38 and 5, which are the values for duration and sevice_hrs.
#'                                            
#' predict(object = robu_mod,pred_vector = c(1,0,0,38,5),level = 0.95)
#' 
#' ## If we do not know the value of duration, the prediction vector should be c(1,0,0,NA,5). 
#' ## predict.robu() will refit a new model without the covariate duration, and the prediction will be based on it.
#' 
#' predict(object = robu_mod,pred_vector = c(1,0,0,NA,5),level = 0.95)
#' 




predict.robu <- function(object, pred_vector, level = 0.95){
  
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    stop("clubSandwich needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!inherits(object, "robu")) 
    warning("calling predict.robu(<fake-lm-object>) ...")

  original_data <- object$Xreg
  combined_data <- rbind(object$Xreg, pred_vector)
  
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
  Test_result <- clubSandwich::Wald_test(new_robu_model, constraints = test_vec, vcov="CR2")
  
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
