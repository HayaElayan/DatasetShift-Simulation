require(shiny)
require(ggplot2)
require(scales)
require(boot)
require(magrittr)
require(dplyr)
require(MASS)
require(Matrix)
require(ROCR)
require(MLmetrics)
require(Metrics)
require(matrixcalc)
require(DescTools)
require(predtools)
require(tidyr)
require(waiter)
require(shinyjs)
require(bslib)
require(future)
require(promises)
require(shinyBS)


plan(sequential) 
invisible(future(NULL))
# Create the matrix

CORR_MATRIX <- matrix(c(
  1,   0.4, 0.4, 0,   0.4,
  0.4, 1,   0.5, 0,   0.3,
  0.4, 0.5, 1,   0,   0.6,
  0,   0,   0,   1,   0.4,
  0.4, 0.3, 0.6, 0.4, 1
), nrow = 5, byrow = TRUE)

predictors_names <- c("X1", "X2", "X3", "X4", "X5")
colnames(CORR_MATRIX) <-predictors_names
binary_predictors <- c('X3','X4', 'X5')
lambda_hyperparam_values <- c(seq(0, 1, by = 0.1), 0.05, 0.08, 0.15, 0.25, 0.03)
folds_num=4

means_target <- c(0,0, 0,0,0)
binary_prev_target <- c(0.3, 0.4, 0.1)
B_coeffcients_target <-  c(-3.1, 0.5, 0.3, 0.7, 0.8, 0.2)
#iteration<-1


generate_development <- function( N_samples, target_split, binary_predictors,
                                  means_target, means_delta, corr_matrix, variances_delta,
                                  binary_prev_target, binary_prev_delta, predictors_names, 
                                  B_coeffcients_target, B_coeffcients_delta, iteration){
  
  source_n <- round(100000*(1-target_split)) 
  target_n <- round(100000*target_split)
  
  source <- generate_covariates_source(source_n, binary_predictors,
                                       means_target, means_delta, corr_matrix, variances_delta,
                                       binary_prev_target, binary_prev_delta, predictors_names, iteration)
  
  target <- generate_covariates_target(target_n, binary_predictors, means_target, corr_matrix,
                                       binary_prev_target, predictors_names, iteration, "target")
  
  
  source_with_outcome <- generate_outcome_source(source, B_coeffcients_target,  B_coeffcients_delta,iteration)
  target_with_outcome <- generate_outcome_target(target, B_coeffcients_target, iteration, "target")
  source_with_outcome$isSource <- 1
  target_with_outcome$isSource <- 0
  
  source_n <- round(N_samples*(1-target_split)) 
  target_n <- round(N_samples*target_split)
  
  source_with_outcome <- source_with_outcome[sample(nrow(source_with_outcome), size = source_n, replace = FALSE), ]
  target_with_outcome<- target_with_outcome[sample(nrow(target_with_outcome), size = target_n, replace = FALSE), ]
  
  full_data <- list(source=source_with_outcome, target=target_with_outcome)
  return(full_data)
}


generate_covariates_source<-function(N, binary_predictors, means, means_delta, corr_mat,
                                     variances_delta , binary_prev, binary_prev_delta, colnames, iteration){
  set.seed(123454*iteration)
  
  corr_mat <- covarinace_matrix(corr_mat, variances_delta)
  means <- (means+means_delta)
  binary_prev <- binary_prev+binary_prev_delta
  
  
  data <- as.data.frame(mvrnorm(n=N,
                                mu=means,
                                Sigma=corr_mat ))
  colnames(data) <- colnames
  
  to_binary_variables <- mapply(to_binray, data[, binary_predictors], binary_prev, iteration ) #apply binary conversion function
  data[, binary_predictors ] <- as.data.frame(to_binary_variables) #pass result to dataframe
  
  return(data)
}


generate_covariates_target<-function(N, binary_predictors, means, corr_mat,
                                     binary_prev, colnames, iteration, dataset_type){
  if(dataset_type=="target"){
    set.seed(123455*iteration) #target
  }else{
    set.seed(123456*iteration) #validation
  }
  
  data <- as.data.frame(mvrnorm(n=N,
                                mu=means,
                                Sigma=corr_mat ))
  colnames(data) <- colnames
  
  to_binary_variables <- mapply(to_binray, data[, binary_predictors], binary_prev, iteration ) #apply binary conversion function
  data[, binary_predictors ] <- as.data.frame(to_binary_variables) #pass result to dataframe
  
  return(data)
}


generate_covariates_target<-function(N, binary_predictors, means, corr_mat,
                                     binary_prev, colnames, iteration, dataset_type){
  if(dataset_type=="target"){
    set.seed(12345*iteration) #target
  }else{
    set.seed(12346*iteration) #validation
  }
  
  data <- as.data.frame(mvrnorm(n=N,
                                mu=means,
                                Sigma=corr_mat ))
  colnames(data) <- colnames
  
  to_binary_variables <- mapply(to_binray, data[, binary_predictors], binary_prev, iteration ) #apply binary conversion function
  data[, binary_predictors ] <- as.data.frame(to_binary_variables) #pass result to dataframe
  
  return(data)
}

generate_outcome_source <- function(covariates, coeffs_target, coeffs_delta, iteration ){
  set.seed(12344*iteration)
  
  mat <- covariates  %>% data.matrix()
  
  coeffs <- coeffs_target + coeffs_delta 
  
  alpha <- coeffs[1]
  
  xb <- alpha + (mat%*%coeffs[-1])
  p <-  1/(1 + exp(-xb))
  y <- rbinom(n = nrow(covariates), size = 1, prob = p)
  #print(table(y)/nrow(covariates))
  
  drop_cols <- which(coeffs[-1] == 0)
  #print(drop_cols)
  if(length(drop_cols) > 0){
    covariates<-covariates %>% dplyr::select(-all_of(drop_cols))
  }
  covariates$Y <- y
  return(covariates)
}


generate_outcome_target <- function(covariates, coeffs_target, iteration, dataset_type){
  if(dataset_type=="target"){
    set.seed(12345*iteration) #target
  }else{
    set.seed(12346*iteration) #validation
  }
  
  mat <- covariates  %>% data.matrix()
  
  coeffs <- coeffs_target
  
  alpha <- coeffs[1]
  
  xb <- alpha + (mat%*%coeffs[-1])
  p <-  1/(1 + exp(-xb))
  y <- rbinom(n = nrow(covariates), size = 1, prob = p)
  
  drop_cols <- which(coeffs[-1] == 0)
  
  if(length(drop_cols) > 0){
    covariates<-covariates %>% dplyr::select(-all_of(drop_cols))
  }
  covariates$Y <- y
  return(covariates)
}


covarinace_matrix <- function(corr_mat, delta_v){
  diag(corr_mat) <- diag(corr_mat)+delta_v
  std <- sqrt(diag(corr_mat))
  cm <- corr_mat * outer(std, std)
  return(cm)
}


brier_score <- function(pred, obs) mean((pred - obs)^2, na.rm = TRUE)



to_binray <- function(variable, prevalence, iteration){
  set.seed(13567*iteration)
  # Determine threshold value
  threshold <- quantile(variable, 1 - prevalence)
  # Convert continuous variable to binary
  binary_variable <- ifelse(variable >= threshold, 1, 0)
  return(binary_variable)
}


generate_validation <- function( N_samples=10000, binary_predictors,
                                 means_target, corr_matrix,
                                 binary_prev_target, predictors_names, 
                                 B_coeffcients_target){
  iteration <- 1
  target <- generate_covariates_target(N_samples , binary_predictors, means_target, corr_matrix,
                                       binary_prev_target, predictors_names, iteration, "validation")
  
  target_with_outcome <- generate_outcome_target(target, B_coeffcients_target, iteration, "validation")
  target_with_outcome$isSource <- 0
  
  return(target_with_outcome)
}

develop_cpm <- function(dev){
  model <- glm(Y ~ X1 + X2 + X3 + X4 + X5, data = dev, family = binomial)
  return(model)
}


intercept_calibration <- function(source, target) {
  target <- target %>% dplyr::select(-c(isSource))
  source <- source %>% dplyr::select(-c(isSource))
  train_data <- rbind(source, target)
  
  source_model <- glm(Y ~ X1 + X2 + X3 + X4 + X5, data=train_data, family = 'binomial', x=TRUE, y=TRUE)
  
  target$lp <- predict(source_model, newdata=target) #get Lp of source model on target data
  calibrated_model <- glm(Y~offset(lp), data=target, family='binomial',x=T, y=T)#update intercept only
  
  
  return(list(calibrated_model=calibrated_model, source_model=source_model))
  
}

logistic_calibration <-function(source, target) {
  
  target <- target %>% dplyr::select(-c(isSource))
  source <- source %>% dplyr::select(-c(isSource))
  train_data <- rbind(source, target)
  
  source_model <- glm(Y ~ X1 + X2 + X3 + X4 + X5, data=train_data, family = 'binomial', x=TRUE, y=TRUE)
  
  target$lp <- predict(source_model, newdata=target) #get Lp of source model on target data
  calibrated_model <- glm(Y~lp, data=target, family='binomial',x=T, y=T)#update all model's coeff
  
  return(list(calibrated_model=calibrated_model, source_model=source_model))
} 



propensity_weighting_with_lambda<- function(source, target){
  temp <- rbind(source, target)
  
  membership_model <- glm(isSource~X1+X2+X3+X4+X5, data=temp, family = 'binomial')
  
  
  score <- predict(membership_model,type="response", newdata = source) #predict p(isSource=1|X)
  
  propensity_weight <-  (1 - score)/ score #p(p(isSource=0|X)/ p(isSource=1|X))
  #print(propensity_weight)
  
  propensity_weight <- propensity_weight*(nrow(source)/ nrow(target))
  
  propensity_weight[propensity_weight >= 1]= 1
  source$propensity_weight <- propensity_weight
  target$propensity_weight <- 1
  
  lambda <- find_lambda(source, target)
  print(lambda)
  weights_adjusted <- c(source$propensity_weight*lambda,target$propensity_weight)
  return(weights_adjusted)
}


find_lambda <- function(source, target){
  hyperparam_values <- c(seq(0, 1, by = 0.1), 0.05, 0.08, 0.15, 0.25, 0.03)
  mse_results <- numeric(length(hyperparam_values))
  
  for (i in seq_along(hyperparam_values)) {
    lambda <- hyperparam_values[i]
    
    
    folds <- custom_cv(source, target) 
    mse_values <- numeric(length(folds))
    
    for (j in seq_along(folds)) {
      train_source_indices <- folds[[j]]$train_source
      train_target_indices <- folds[[j]]$train_target
      test_indices <- folds[[j]]$test
      
      train_source_data <- source[train_source_indices, ] 
      train_target_data <- target[train_target_indices, ] 
      test_data <- target[test_indices, ] 
      
      source_target_weights <- propensity_weighting_x(train_source_data, train_target_data, lambda)
      
      result <- weighted_logistic_regression_x(train_source_data, train_target_data, test_data, weight = source_target_weights)
      mse_values[j] <- calculate_mse(true_labels = test_data$Y, pred_probs = result)
      
    }
    
    mse_results[i] <- mean(mse_values)
  }
  
  
  optimal_param <- hyperparam_values[which.min(mse_results)]
  return(optimal_param)
}



propensity_weighting_x<- function(source, target, lambda){
  temp <- rbind(source, target)
  
  membership_model <- glm(isSource~X1+X2+X3+X4+X5, data=temp, family = 'binomial')
  
  score <- predict(membership_model,type="response", newdata = source) #predict p(isSource=1|X)
  
  
  propensity_weight <-  (1 - score)/ score #p(p(isSource=0|X)/ p(isSource=1|X))
  
  propensity_weight <- propensity_weight*(nrow(source)/ nrow(target))
  
  propensity_weight[propensity_weight >= 1]= 1
  source$propensity_weight <- propensity_weight
  target$propensity_weight <- 1
  
  weights <- c(source$propensity_weight*lambda,target$propensity_weight )
  
  return(weights)
  
}


# Define custom cross-validation function
custom_cv <- function(source, target) {
  indices_source <- sample(rep(1:4, length.out = nrow(source)))
  indices_target <- sample(rep(1:4, length.out = nrow(target)))
  
  folds_list <- lapply(1:4, function(i) {
    train_source_indices <- which(indices_source != i)
    train_target_indices <- which(indices_target != i)
    test_indices <- which(indices_target == i)
    list(train_source = train_source_indices, train_target= train_target_indices, test = test_indices)
  })
  return(folds_list)
}


# Define function to calculate mean squared error
calculate_mse <- function(true_labels, pred_probs) {
  mse <- MLmetrics::MSE(y_pred = pred_probs, y_true = true_labels)
  return(mse)
}

weighted_logistic_regression_x <- function(source, target, test_data, weight) {
  train_data <- rbind(source, target)
  
  model <- glm(Y~X1+X2+X3+X4+X5+isSource, data = train_data, weights = weight, family = binomial)
  
  pred_probs <- predict(model, newdata = test_data, type = "response")
  
  return( pred_probs)
  
}

weighted_LR <- function(source, target, ps_weights) {
  train_data <- rbind(source, target)
  
  model <- glm(Y~X1+X2+X3+X4+X5+isSource, data = train_data, weights = ps_weights, family = binomial)
  #print(summary(model))
  return( model)
  
}




# Function to simulate and model data
simulate_data <-
  function(N_SAMPLES = 1000,  target_split = 0.5, X0_b_delta = 0, 
           X1_m_delta = 0, X1_v_delta = 0,  X1_b_delta = 0,  
           X2_m_delta = 0, X2_v_delta = 0, X2_b_delta = 0,
           X3_p_delta = 0, X3_b_delta = 0, 
           X4_p_delta = 0, X4_b_delta = 0,
           X5_p_delta = 0,  X5_b_delta = 0, cpm_method) {
    
    variances_delta <- c(X1_v_delta, X2_v_delta, 0, 0, 0)#input
    means_delta <- c(X1_m_delta, X2_m_delta, 0, 0, 0) #input
    binary_prev_delta <- c(X3_p_delta, X4_p_delta, X5_p_delta)#input
    B_coeffcients_delta <-
      c(X0_b_delta, X1_b_delta,   X2_b_delta,  X3_b_delta,  X4_b_delta,  X5_b_delta)#input
    val_results <- list()
    cal_plots <- list()
    test <- c()
    
    for(iteration in 1:50){
      source_target <-
        generate_development(
          N_SAMPLES,
          target_split,
          binary_predictors,
          means_target,
          means_delta,
          CORR_MATRIX,
          variances_delta,
          binary_prev_target,
          binary_prev_delta,
          predictors_names,
          B_coeffcients_target,
          B_coeffcients_delta,
          iteration
        )
      dev <- rbind(source_target$source, source_target$target)
      if(iteration == 1){
        test <-
          generate_validation(
            10000,
            binary_predictors,
            means_target,
            CORR_MATRIX,
            binary_prev_target,
            predictors_names,
            B_coeffcients_target
          )
      }
      print(cpm_method)
      if(cpm_method == "Full-data Logistic Regression"){
        # Fit model
        model <- develop_cpm(dev)
        all_results <- CPM_validation_results(model, test)
        val_results[[iteration]] <- all_results$val_results
        cal_plots[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Target-dataset only logistic regression"){
        # Fit model
        model <- develop_cpm(source_target$target)
        all_results <- CPM_validation_results(model, test)
        val_results[[iteration]] <- all_results$val_results
        cal_plots[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Ancillary-dataset only logistic regression"){
        # Fit model
        model <- develop_cpm(source_target$source)
        all_results <- CPM_validation_results(model, test)
        val_results[[iteration]] <- all_results$val_results
        cal_plots[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Memebership-based Weighting"){
        # Fit model
        weights <- propensity_weighting_with_lambda(source_target$source, source_target$target)
        model <- weighted_LR(source_target$source, source_target$target, weights)
        all_results <- proposed_model_validation_results(model, test)
        val_results[[iteration]] <- all_results$val_results
        cal_plots[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Intercept Recalibration on all data"){
        # Fit model
        model <- intercept_calibration(source_target$source, source_target$target)
        all_results <- calibrated_model_validation_results(model, "intercept_only", test)
        val_results[[iteration]] <- all_results$val_results
        cal_plots[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Logistic Recalibration on all data"){
        # Fit model
        model <- logistic_calibration(source_target$source, source_target$target)
        all_results <- calibrated_model_validation_results(model, "logistic", test)
        val_results[[iteration]] <- all_results$val_results
        cal_plots[[iteration]] <- all_results$cal_plot
      }
    }
    # Shifted test data
    return(list(
      dev = dev,
      train = source_target$source,
      target = source_target$target,
      val_results = val_results,
      cal_plots = cal_plots)
    )
  }


simulate_noshift <-
  function(N_SAMPLES = 1000,  target_split = 0.5, cpm_method) {
    variances_delta <- c(0, 0, 0, 0, 0)#input
    means_delta <- c(0, 0, 0, 0, 0) #input
    binary_prev_delta <- c(0, 0, 0)
    B_coeffcients_delta <-
      c(0, 0, 0, 0, 0, 0)#input
    val_results_noshift <- list()
    cal_plots_noshift <- list()
    
    for(iteration in 1:50){
      source_target <-
        generate_development(
          N_SAMPLES,
          target_split,
          binary_predictors,
          means_target,
          means_delta,
          CORR_MATRIX,
          variances_delta,
          binary_prev_target,
          binary_prev_delta,
          predictors_names,
          B_coeffcients_target,
          B_coeffcients_delta,
          iteration
        )
      dev <- rbind(source_target$source, source_target$target)
      if(iteration == 1){
        test <-
          generate_validation(
            10000,
            binary_predictors,
            means_target,
            CORR_MATRIX,
            binary_prev_target,
            predictors_names,
            B_coeffcients_target
          )
      }
      
      
      print(cpm_method)
      if(cpm_method == "Full-data Logistic Regression"){
        # Fit model
        model <- develop_cpm(dev)
        all_results <- CPM_validation_results(model, test)
        val_results_noshift[[iteration]] <- all_results$val_results
        cal_plots_noshift[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Target-dataset only logistic regression"){
        # Fit model
        model <- develop_cpm(source_target$target)
        all_results <- CPM_validation_results(model, test)
        val_results_noshift[[iteration]] <- all_results$val_results
        cal_plots_noshift[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Ancillary-dataset only logistic regression"){
        # Fit model
        model <- develop_cpm(source_target$source)
        all_results <- CPM_validation_results(model, test)
        val_results_noshift[[iteration]] <- all_results$val_results
        cal_plots_noshift[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Memebership-based Weighting"){
        # Fit model
        weights <- propensity_weighting_with_lambda(source_target$source, source_target$target)
        model <- weighted_LR(source_target$source, source_target$target, weights)
        all_results <- proposed_model_validation_results(model, test)
        val_results_noshift[[iteration]] <- all_results$val_results
        cal_plots_noshift[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Intercept Recalibration on all data"){
        # Fit model
        model <- intercept_calibration(source_target$source, source_target$target)
        all_results <- calibrated_model_validation_results(model, "intercept_only", test)
        val_results_noshift[[iteration]] <- all_results$val_results
        cal_plots_noshift[[iteration]] <- all_results$cal_plot
      }
      
      if(cpm_method == "Logistic Recalibration on all data"){
        # Fit model
        model <- logistic_calibration(source_target$source, source_target$target)
        all_results <- calibrated_model_validation_results(model, "logistic", test)
        val_results_noshift[[iteration]] <- all_results$val_results
        cal_plots_noshift[[iteration]] <- all_results$cal_plot
      }
    }
    # Shifted test data
    return(list(
      dev = dev,
      train = source_target$source,
      target = source_target$target,
      val_results_noshift = val_results_noshift,
      cal_plots_noshift = cal_plots_noshift)
    )
  }

#bootstrap validation results
CPM_validation_results <- function(model, validation){
  if(all(is.na(model))){
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    val_results[1,1] <- NA
    val_results[1,2] <- NA
    val_results[1,3] <- NA
    val_results[1,4] <- NA
    
    Sm.full <- rep(NA, 100)
    return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
  }else{
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    
    pr_val <- predict(model, type="response", newdata = validation) # predict probabilities 
    lp_val <- predict(model, type="link", newdata = validation ) # predict lp type=link
    

    #val_cstat_model <- #roc(Y ~ pr_val,data=validation)
    #val_results[1,1] <- val_cstat_model$auc
    
    pred <- prediction(pr_val, validation$Y)
    perf <- performance(pred, "tpr", "fpr")    
    auc <- performance(pred, "auc")@y.values[[1]]
    val_results[1,1] <- auc
    
    val_results[1,4] <- brier_score(validation$Y, pred=pr_val) 
    tryCatch(
      { 
        val_citl_model <- glm(Y ~ offset(lp_val),family=binomial, data=validation)
        val_results[1,2] <- summary(val_citl_model)$coefficients[1,1]
        
        val_cslope_model <- glm(Y ~ lp_val,family=binomial(link='logit'), data=validation)
        val_results[1,3] <- summary(val_cslope_model)$coefficients[2,1]
        
        
        Sm       <- lowess(pr_val, validation$Y, iter = 0)
        #pp.full  <- seq(min(pr_val), max(pr_val), length = 100) #xxxx
        pp.full<- seq(0.001, 0.99, length=100)
        Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y #yyyyy
        
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }, error=function(e){
        val_results[1,2] <- NA
        val_results[1,3] <- NA
        
        Sm.full <- rep(NA, 100)
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }
      
    )
  }
}

# Define your own replacement function
alpha <- function(col, alpha = 1) {
  # Convert color to RGB and apply transparency
  rgb.matrix <- grDevices::col2rgb(col, alpha = TRUE) / 255
  grDevices::rgb(rgb.matrix[1,], rgb.matrix[2,], rgb.matrix[3,], alpha = alpha)
}


calibrated_model_validation_results <- function(model, type, validation){
  
  val_results <- matrix(nrow = 1,ncol = 4)
  colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
  pr_val <- nrow(validation)
  lp_val <- nrow(validation)
  
  validation$lp <- predict(model$source_model, newdata = validation, type = "link")
  #print(table(validation$Y, validation$lp))
  
  if(type=='intercept_only'){
    pr_val <- predict(model$calibrated_model, type="response", newdata = validation, offset = validation$lp) # predict probabilities 
    lp_val <- predict(model$calibrated_model, newdata = validation, offset = validation$lp ) # predict lp type=link
    #print(lp_val)
  }else{
    pr_val <- predict(model$calibrated_model, type="response", newdata = validation, lp = validation$lp) # predict probabilities 
    lp_val <- predict(model$calibrated_model, newdata = validation, lp = validation$lp ) # predict lp type=link
  }
  
  
  # calculate performance of the  model in the validation sample
  #val_cstat_model <- roc(Y ~ pr_val,data=validation)
  #val_results[1,1] <- val_cstat_model$auc
  
  pred <- prediction(pr_val, validation$Y)
  perf <- performance(pred, "tpr", "fpr")    
  auc <- performance(pred, "auc")@y.values[[1]]
  val_results[1,1] <- auc
  
  val_citl_model <- glm(Y ~ offset(lp_val),family=binomial, data=validation)
  
  
  val_cslope_model <- glm(Y ~ lp_val,family=binomial(link='logit'), data=validation)
  
  
  val_results[1,4] <- brier_score(validation$Y, pred=pr_val) #brier_score(model$calibrated_model)
  
  tryCatch(
    {
      val_results[1,2] <- summary(val_citl_model)$coefficients[1,1]
      val_results[1,3] <- summary(val_cslope_model)$coefficients[2,1]
      
      Sm <- lowess(pr_val, validation$Y, iter = 0)
      pp.full<- seq(0.001, 0.99, length=100)
      Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y #yyyyy
      return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      
    }, error=function(e){
      val_results[1,2] <- NA
      val_results[1,3] <- NA
      Sm.full <- rep(NA, 100)
      return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
    }
  )
}


proposed_model_validation_results <- function(model, validation){
  if(all(is.na(model))){
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    val_results[1,1] <- NA
    val_results[1,2] <- NA
    val_results[1,3] <- NA
    val_results[1,4] <- NA
    
    Sm.full <- rep(NA, 100)
    return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
  }else{
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    
    pr_val <- predict(model, type="response", newdata = validation) # predict probabilities 
    lp_val <- predict(model, newdata = validation ) # predict lp type=link
    
    # print(pr_val[1:10])
    # print(lp_val[1:10])
    # calculate performance of the  model in the validation sample
    #val_cstat_model <- roc(Y ~ pr_val,data=validation)
    #val_results[1,1] <- val_cstat_model$auc
    
    pred <- prediction(pr_val, validation$Y)
    perf <- performance(pred, "tpr", "fpr")    
    auc <- performance(pred, "auc")@y.values[[1]]
    val_results[1,1] <- auc
    
    val_results[1,4] <- brier_score(validation$Y, pred=pr_val)#brier_score(model)
    
    tryCatch(
      { 
        val_citl_model <- glm(Y ~ offset(lp_val),family=binomial, data=validation)
        val_results[1,2] <- summary(val_citl_model)$coefficients[1,1]
        
        val_cslope_model <- glm(Y ~ lp_val,family=binomial(link='logit'), data=validation)
        val_results[1,3] <- summary(val_cslope_model)$coefficients[2,1]
        
        
        Sm <- lowess(pr_val, validation$Y, iter = 0)
        #pp.full  <- seq(min(pr_val), max(pr_val), length = 100) #xxxx
        pp.full<- seq(0.001, 0.99, length=100)
        Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y #yyyyy
        
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }, error=function(e){
        val_results[1,2] <- NA
        val_results[1,3] <- NA
        
        Sm.full <- rep(NA, 100)
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }
      
    )
  }
}



get_smoothed_mean <- function(pp_sm_list, sc_name){
  smoothed_lines <- list()
  
  model_name <- pp_sm_list[[1]]['Model']
  
  for(i in 1:length(pp_sm_list)){
    
    Sm.full  <-  pp_sm_list[[i]]$sm[[1]]
    
    smoothed_lines[[i]] <- Sm.full
  }
  
  
  smoothed_matrix <- do.call(cbind, smoothed_lines)
  smoothed_mean <- apply(smoothed_matrix, 1, median, na.rm = TRUE) #rowMeans(smoothed_matrix, na.rm = TRUE)
  
  return(list(smoothed_mean=smoothed_mean, model_name=model_name, sc_name=sc_name))
}


get_val_data <- function(model_df){
  temp_df <- data.frame(AUC= get_ci(model_df$AUC), CITL=get_ci(model_df$CITL),
                        CSLOPE=get_ci(model_df$CSLOPE), BrierScore=get_ci(model_df$BrierScore))
  return(temp_df)
}

get_ci <- function(data){
  #data <- sort(data)
  #print(data)
  data <- data[!is.na(data)]
  
  lower_band <- quantile(data, probs = 0.025, na.rm = TRUE)
  upper_band <- quantile(data, probs = 0.975, na.rm = TRUE)
  return(list(median = median(data), lower_band= lower_band, upper_band = upper_band, divergence = length(data)))
}



val_metric_plot <- function(df) { 
  df <- do.call(rbind, df)
  df <- get_val_data(df)
  
  # Build a data frame with Metric, Estimate, and Divergence columns
  summary_table <- tibble(
    Metric = c("AUC", "Calibration-in-the-large", "Calibration slope", "Brier Score"),
    Estimate = c(
      sprintf("%.3f (%.3f, %.3f)", df$AUC.median, df$AUC.lower_band, df$AUC.upper_band),
      sprintf("%.3f (%.3f, %.3f)", df$CITL.median, df$CITL.lower_band, df$CITL.upper_band),
      sprintf("%.3f (%.3f, %.3f)", df$CSLOPE.median, df$CSLOPE.lower_band, df$CSLOPE.upper_band),
      sprintf("%.3f (%.3f, %.3f)", df$BrierScore.median, df$BrierScore.lower_band, df$BrierScore.upper_band)
    ),
    Divergence = c(
      sprintf("%d/50 iterations", df$AUC.divergence),
      sprintf("%d/50 iterations", df$CITL.divergence),
      sprintf("%d/50 iterations", df$CSLOPE.divergence),
      sprintf("%d/50 iterations", df$BrierScore.divergence)
    )
  ) %>%
    mutate(Metric = paste0("<b>", Metric, "</b>")) %>%
    rename(`Estimate (95% CI)` = Estimate)
  
  return(summary_table)
}


calibration_instability_plot<- function(pp_sm_list) {
  smoothed_lines <- list()
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "Predicted", ylab = "Observed")
  
  abline(a = 0, b = 1, col = 'black')  # Add ideal line
  
  pp.full <- seq(0.001, 0.99, length = 100)
  for (i in 1:length(pp_sm_list)) {
    Sm.full <- pp_sm_list[[i]]$sm
    lines(pp.full, Sm.full, col = alpha('grey', 0.4))  # Add lines for single iteration
    rug(Sm.full, side = 2, col = alpha('dark green', 0.4))  # Add rug plot on y-axis
    smoothed_lines[[i]] <- Sm.full
  }
  
  smoothed_matrix <- do.call(cbind, smoothed_lines)
  smoothed_mean <- apply(smoothed_matrix, 1, median, na.rm = TRUE)
  smoothed_quantiles <- apply(smoothed_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  lines(pp.full, smoothed_mean, col = "dark orange", lwd = 1)  # Add median line
  lines(pp.full, smoothed_quantiles[1, ], col = alpha("blue", 0.8), lwd = 1, lty = "dashed")  # Add lower band
  lines(pp.full, smoothed_quantiles[2, ], col = alpha("dark blue", 0.8), lwd = 1, lty = "dashed")  # Add upper band
  legend("bottomright", legend = c("Ideal", "Single iteration", "Median", "Lower band", "Upper band"),
         col = c('black', alpha('grey', 0.4), "dark orange", alpha("blue", 0.8), alpha("dark blue", 0.8)),
         lty = c("solid", "solid", "solid", "dashed", "dashed"),  bty = "n")
  
}



#bootstrap validation results
CPM_validation_results <- function(model, validation){
  if(all(is.na(model))){
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    val_results[1,1] <- NA
    val_results[1,2] <- NA
    val_results[1,3] <- NA
    val_results[1,4] <- NA
    
    Sm.full <- rep(NA, 100)
    return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
  }else{
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    
    pr_val <- predict(model, type="response", newdata = validation) # predict probabilities 
    lp_val <- predict(model, type="link", newdata = validation ) # predict lp type=link
    
    # print(pr_val[1:10])
    # print(lp_val[1:10])
    # calculate performance of the  model in the validation sample
    #val_cstat_model <- roc(Y ~ pr_val,data=validation)
    #val_results[1,1] <- val_cstat_model$auc
    pred <- prediction(pr_val, validation$Y)
    perf <- performance(pred, "tpr", "fpr")    
    auc <- performance(pred, "auc")@y.values[[1]]
    val_results[1,1] <- auc
    
    val_results[1,4] <- brier_score(validation$Y, pred=pr_val) 
    tryCatch(
      { 
        val_citl_model <- glm(Y ~ offset(lp_val),family=binomial, data=validation)
        val_results[1,2] <- summary(val_citl_model)$coefficients[1,1]
        
        val_cslope_model <- glm(Y ~ lp_val,family=binomial(link='logit'), data=validation)
        val_results[1,3] <- summary(val_cslope_model)$coefficients[2,1]
        
        
        Sm       <- lowess(pr_val, validation$Y, iter = 0)
        #pp.full  <- seq(min(pr_val), max(pr_val), length = 100) #xxxx
        pp.full<- seq(0.001, 0.99, length=100)
        Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y #yyyyy
        
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }, error=function(e){
        val_results[1,2] <- NA
        val_results[1,3] <- NA
        
        Sm.full <- rep(NA, 100)
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }
      
    )
  }
}



calibrated_model_validation_results <- function(model, type, validation){
  
  val_results <- matrix(nrow = 1,ncol = 4)
  colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
  pr_val <- nrow(validation)
  lp_val <- nrow(validation)
  
  validation$lp <- predict(model$source_model, newdata = validation, type = "link")
  #print(table(validation$Y, validation$lp))
  
  if(type=='intercept_only'){
    pr_val <- predict(model$calibrated_model, type="response", newdata = validation, offset = validation$lp) # predict probabilities 
    lp_val <- predict(model$calibrated_model, newdata = validation, offset = validation$lp ) # predict lp type=link
    #print(lp_val)
  }else{
    pr_val <- predict(model$calibrated_model, type="response", newdata = validation, lp = validation$lp) # predict probabilities 
    lp_val <- predict(model$calibrated_model, newdata = validation, lp = validation$lp ) # predict lp type=link
  }
  
  
  # calculate performance of the  model in the validation sample
  #val_cstat_model <- roc(Y ~ pr_val,data=validation)
  #val_results[1,1] <- val_cstat_model$auc
  pred <- prediction(pr_val, validation$Y)
  perf <- performance(pred, "tpr", "fpr")    
  auc <- performance(pred, "auc")@y.values[[1]]
  val_results[1,1] <- auc
  
  val_citl_model <- glm(Y ~ offset(lp_val),family=binomial, data=validation)
  
  
  val_cslope_model <- glm(Y ~ lp_val,family=binomial(link='logit'), data=validation)
  
  
  val_results[1,4] <- brier_score(validation$Y, pred=pr_val) #brier_score(model$calibrated_model)
  
  tryCatch(
    {
      val_results[1,2] <- summary(val_citl_model)$coefficients[1,1]
      val_results[1,3] <- summary(val_cslope_model)$coefficients[2,1]
      
      Sm <- lowess(pr_val, validation$Y, iter = 0)
      pp.full<- seq(0.001, 0.99, length=100)
      Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y #yyyyy
      return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      
    }, error=function(e){
      val_results[1,2] <- NA
      val_results[1,3] <- NA
      Sm.full <- rep(NA, 100)
      return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
    }
  )
}


proposed_model_validation_results <- function(model, validation){
  if(all(is.na(model))){
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    val_results[1,1] <- NA
    val_results[1,2] <- NA
    val_results[1,3] <- NA
    val_results[1,4] <- NA
    
    Sm.full <- rep(NA, 100)
    return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
  }else{
    val_results <- matrix(nrow = 1,ncol = 4)
    colnames(val_results) <- c('AUC', 'CITL', 'CSLOPE', 'BrierScore') #'MER' = mean estimated risk
    
    
    pr_val <- predict(model, type="response", newdata = validation) # predict probabilities 
    lp_val <- predict(model, newdata = validation ) # predict lp type=link
    
    # print(pr_val[1:10])
    # print(lp_val[1:10])
    # calculate performance of the  model in the validation sample
    #val_cstat_model <- roc(Y ~ pr_val,data=validation)
    #val_results[1,1] <- val_cstat_model$auc
    pred <- prediction(pr_val, validation$Y)
    perf <- performance(pred, "tpr", "fpr")    
    auc <- performance(pred, "auc")@y.values[[1]]
    val_results[1,1] <- auc
    
    val_results[1,4] <- brier_score(validation$Y, pred=pr_val)#brier_score(model)
    
    tryCatch(
      { 
        val_citl_model <- glm(Y ~ offset(lp_val),family=binomial, data=validation)
        val_results[1,2] <- summary(val_citl_model)$coefficients[1,1]
        
        val_cslope_model <- glm(Y ~ lp_val,family=binomial(link='logit'), data=validation)
        val_results[1,3] <- summary(val_cslope_model)$coefficients[2,1]
        
        
        Sm <- lowess(pr_val, validation$Y, iter = 0)
        #pp.full  <- seq(min(pr_val), max(pr_val), length = 100) #xxxx
        pp.full<- seq(0.001, 0.99, length=100)
        Sm.full  <- approx(Sm, xout = pp.full, ties = mean)$y #yyyyy
        
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }, error=function(e){
        val_results[1,2] <- NA
        val_results[1,3] <- NA
        
        Sm.full <- rep(NA, 100)
        
        return(list(val_results = as.data.frame(val_results), cal_plot = list(sm = Sm.full)) )
      }
      
    )
  }
}



get_smoothed_mean <- function(pp_sm_list, sc_name){
  smoothed_lines <- list()
  
  model_name <- pp_sm_list[[1]]['Model']
  
  for(i in 1:length(pp_sm_list)){
    
    Sm.full  <-  pp_sm_list[[i]]$sm[[1]]
    
    smoothed_lines[[i]] <- Sm.full
  }
  
  
  smoothed_matrix <- do.call(cbind, smoothed_lines)
  smoothed_mean <- apply(smoothed_matrix, 1, median, na.rm = TRUE) #rowMeans(smoothed_matrix, na.rm = TRUE)
  
  return(list(smoothed_mean=smoothed_mean, model_name=model_name, sc_name=sc_name))
}


get_val_data <- function(model_df){
  temp_df <- data.frame(AUC= get_ci(model_df$AUC), CITL=get_ci(model_df$CITL),
                        CSLOPE=get_ci(model_df$CSLOPE), BrierScore=get_ci(model_df$BrierScore))
  return(temp_df)
}

get_ci <- function(data){
  #data <- sort(data)
  #print(data)
  data <- data[!is.na(data)]
  
  lower_band <- quantile(data, probs = 0.025, na.rm = TRUE)
  upper_band <- quantile(data, probs = 0.975, na.rm = TRUE)
  return(list(median = median(data), lower_band= lower_band, upper_band = upper_band, divergence = length(data)))
}



val_metric_plot <- function(df) { 
  df <- do.call(rbind, df)
  df <- get_val_data(df)
  
  # Build a data frame with Metric, Estimate, and Divergence columns
  summary_table <- tibble(
    Metric = c("AUC", "Calibration-in-the-large", "Calibration slope", "Brier Score"),
    Estimate = c(
      sprintf("%.3f (%.3f, %.3f)", df$AUC.median, df$AUC.lower_band, df$AUC.upper_band),
      sprintf("%.3f (%.3f, %.3f)", df$CITL.median, df$CITL.lower_band, df$CITL.upper_band),
      sprintf("%.3f (%.3f, %.3f)", df$CSLOPE.median, df$CSLOPE.lower_band, df$CSLOPE.upper_band),
      sprintf("%.3f (%.3f, %.3f)", df$BrierScore.median, df$BrierScore.lower_band, df$BrierScore.upper_band)
    ),
    Divergence = c(
      sprintf("%d/50 iterations", df$AUC.divergence),
      sprintf("%d/50 iterations", df$CITL.divergence),
      sprintf("%d/50 iterations", df$CSLOPE.divergence),
      sprintf("%d/50 iterations", df$BrierScore.divergence)
    )
  ) %>%
    mutate(Metric = paste0("<b>", Metric, "</b>")) %>%
    rename(`Estimate (95% CI)` = Estimate)
  
  return(summary_table)
}


calibration_instability_plot<- function(pp_sm_list) {
  smoothed_lines <- list()
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "Predicted", ylab = "Observed")
  
  abline(a = 0, b = 1, col = 'black')  # Add ideal line
  
  pp.full <- seq(0.001, 0.99, length = 100)
  for (i in 1:length(pp_sm_list)) {
    Sm.full <- pp_sm_list[[i]]$sm
    lines(pp.full, Sm.full, col = alpha('grey', 0.4))  # Add lines for single iteration
    rug(Sm.full, side = 2, col = alpha('dark green', 0.4))  # Add rug plot on y-axis
    smoothed_lines[[i]] <- Sm.full
  }
  
  smoothed_matrix <- do.call(cbind, smoothed_lines)
  smoothed_mean <- apply(smoothed_matrix, 1, median, na.rm = TRUE)
  smoothed_quantiles <- apply(smoothed_matrix, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  lines(pp.full, smoothed_mean, col = "dark orange", lwd = 1)  # Add median line
  lines(pp.full, smoothed_quantiles[1, ], col = alpha("blue", 0.8), lwd = 1, lty = "dashed")  # Add lower band
  lines(pp.full, smoothed_quantiles[2, ], col = alpha("dark blue", 0.8), lwd = 1, lty = "dashed")  # Add upper band
  legend("bottomright", legend = c("Ideal", "Single iteration", "Median", "Lower band", "Upper band"),
         col = c('black', alpha('grey', 0.4), "dark orange", alpha("blue", 0.8), alpha("dark blue", 0.8)),
         lty = c("solid", "solid", "solid", "dashed", "dashed"),  bty = "n")
  
}


# Calibration plot function
calibration_plot <- function(cal_plots) {
  
  calibration_instability_plot(cal_plots)
}

validation_results_table <- function(val_results){
  val_metric_plot(val_results)
  
}


calibration_plot_noshift <- function(cal_plots_noshift) {
  
  calibration_instability_plot(cal_plots_noshift)
}

validation_results_table_noshift <- function(val_results_noshift){
  val_metric_plot(val_results_noshift)
  
}

ui <- navbarPage(title = "🔬Heterogeneity Simulator",
                 
                 # --- Info Tab ---
                 tabPanel(
                   "ℹ️ About",
                   fluidPage(
                     withMathJax(),
                     tags$h2("Understanding Data Heterogeneity Impact on Model Performance", style = "margin-top: 20px; color: #2C3E50; font-weight: bold; text-align: center;"),
                     tags$h3("Simulating Distribution Shifts and Their Effects on Predictive Models", style = "margin-top: 20px; color: #2C3E50; font-weight: bold; text-align: center;"),
                     
                     tags$hr(),
                     tags$div(
                       tags$ul(
                         tags$li(HTML(
                           "This app provides an interactive simulation designed to explore the impact of heterogeneity within the development 
                           dataset on predictive model performance due to data distribution shifts, which occur when differences arise 
                            in the underlying distribution of the data. 
                                  <br> <br>
                                  The development dataset in this simulation consists of two datasets:
                                  <br>
                           <br> 1) Target dataset is data sampled from a specific target population where the model is intended to be deployed. 
                           This dataset is generated using a set of fixed parameters. <br>
                            <br> 2) Ancillary dataset is additional and potentially related data sampled from an ancillary population from a different time point or location. 
                           Users generate this dataset under different distribution shift scenarios, to simulate heterogeneity between the target and the ancillary datasets.<br>"
                         )),
                         
                         tags$h3("Distribution Shift Scenarios")
                         ,
                         tags$li(
                           HTML(
                             "Users can simulate the following main scenarios:
              <br><strong>1. No data shift:</strong> Generate the ancillary population from same distribution as the target population.
              <br><strong>2. Case-mix shift:</strong> Generate the ancillary population from a distribution that differs from the target
              population by a shift in mean, variance, and/or event rate (prevalence) of any predictor.
              <br><strong>3. Predictor-outcome association shift:</strong> Generate the ancillary population by shifting the coefficients
              (β) of any predictor.
              <br><strong>4. Event rate shift:</strong> Generate the ancillary population with a different intercept value.
              <br><strong>5. Combination of shifts:</strong> Apply multiple shifts simultaneously to assess compounded effects."
                           )
                         ),
                         tags$h3(
                           HTML(
                             "Predictor-generation models")
                         ),
                         
                         tags$li(
                           HTML(
                             "<strong>Predictor-generation model for the target population</strong> includes 5 variables:
         <br>
         \\[
           X_{\\text{target}} = (X_1, X_2, X_3, X_4, X_5)^T \\sim MVN(\\mu_{\\text{target}}, \\Sigma_{\\text{target}})
           \\]
          <br>Where the mean vector \\(\\mu_{\\text{target}}\\) is set to 0, and the covariance matrix \\(\\Sigma_{\\text{target}}\\) is:
          
          <br><br>
          <table border='1' cellpadding='5' cellspacing='0' style='border-collapse: collapse; text-align: center; margin: auto;'>
          <tr>
          <th></th><th>X₁</th><th>X₂</th><th>X₃</th><th>X₄</th><th>X₅</th>
          </tr>
          <tr><th>X₁</th><td>1</td><td>0.4</td><td>0.4</td><td>0</td><td>0.4</td></tr>
          <tr><th>X₂</th><td>0.4</td><td>1</td><td>0.5</td><td>0</td><td>0.3</td></tr>
          <tr><th>X₃</th><td>0.4</td><td>0.5</td><td>1</td><td>0</td><td>0.6</td></tr>
          <tr><th>X₄</th><td>0</td><td>0</td><td>0</td><td>1</td><td>0.4</td></tr>
          <tr><th>X₅</th><td>0.4</td><td>0.3</td><td>0.6</td><td>0.4</td><td>1</td></tr>
          </table>
          <br><br>
          Predictors <code>X₃, X₄, X₅</code> are dichotomized using thresholds based on the variable’s distribution.
          Values above the threshold are assigned 1, and those below 0, with prevalences:
          <code>p₃ = 0.3</code>, <code>p₄ = 0.4</code>, <code>p₅ = 0.1</code>.<br> <br>"
                           )
                         ),
                         
                         tags$li(
                           HTML(
                             "<strong>Predictor-generation model for the ancillary population</strong>:
                  <br>
         \\[
           X_{\\text{ancillary}} = MVN(\\mu_{\\text{target}}+ \\Delta_{\\mu}, \\Sigma_{\\text{target}} +  \\Delta_{\\Sigma})
           \\] <br>      <br> where \\( \\Delta \\) refers to the user-specified shift."
                           )
                         ),
                         
                         tags$h3(
                           HTML(
                             "Outcome-generation models")
                         ),
                         
                         tags$li(
                           HTML(
                             "The following <strong>outcome-generation model for the target population</strong> is uded to generate a binary outcome {Y}:
                   <br>
                   \\[
                   P(Y = 1 \\mid X) = \\pi_{\\text{target}} = \\left[1 + \\exp\\left( - \\left( \\beta_{0,\\text{target}} + \\sum_{k=1}^{5} \\beta_{k,\\text{target}} X_k \\right) \\right) \\right]^{-1}
                   \\]
                   <br>
                   With coefficients:
                   <code>β₀ = -3.1, β₁ = 0.5, β₂ = 0.3, β₃ = 0.7, β₄ = 0.8, β₅ = 0.2</code>  <br>  <br>"
                           )
                         ),
                         
                         tags$li(
                           HTML(
                             "<strong>The outcome-generation model for the ancillary population</strong>:
               <br>
               \\[
               P(Y = 1 \\mid X) = \\pi_{\\text{ancillary}} =
               \\left[ 1 + \\exp\\left( - \\left( (\\beta_{0,\\text{target}} + \\Delta\\beta_{0}) +
               \\sum_{k=1}^{5} (\\beta_{k,\\text{target}} + \\Delta\\beta_{k}) X_k \\right) \\right) \\right]^{-1}
               \\]
               <br>      <br> where \\( \\Delta \\) refers to the user-specified shift. <br> <br>
              "
                           )
                         ) ,
                         tags$li("Throughout all the simulation iterations across all simulation scenarios, a large population-level datasets is going to be generated—comprising a total of 100,000 observations— for both the target and ancillary populations. 
                                 The goal is to simulate overarching populations from which users can draw random samples for developing the model."
                         ),
                         tags$h3("Development dataset"),
                         tags$li(
                           "Users can draw the development dataset using combination of the target and ancillary populations subsets,
                            with varying proportions of splits."
                         ), 
                         tags$h3("Validation dataset"),
                         tags$li("The model performance is evaluated on a large validation dataset (10,000), 
                                 drawn from the same data-generating mechanism as the target population."
                         ) ,
                         tags$li("Each simulation run performs 50 iterations."),
                         tags$h3("Survey", style = "color: darkred;"),
                         tags$li("After exploring the app, please participate in this survey", style = "color: darkred;", 
                                 tags$a(href="https://tinyurl.com/simulatorappsurvey", "https://tinyurl.com/simulatorappsurvey.") ,tags$br(),
                                 "The purpose of this survey is to gather your feedback on the app. 
                           Your participation will help us understand how effective the tool is for learning, how engaging and easy to use it is,
                                and how it could be improved.", style = "color: darkred;" ),
                         tags$br(),
                         
                         p(" © 2025 Haya Elayan. All rights reserved.")
                       ),
                       style = "font-size: 15px; max-width: 800px; margin: auto;"
                     )
                     
                   )
                 ),
                 # --- Simulation Tab ---
                 tabPanel(
                   "📊 Simulation",
                   fluidPage(
                     useShinyjs(),
                     #theme = shinytheme("darkly"), 
                     tags$head(
                       # Load FontAwesome for icon
                       tags$link(
                         rel = "stylesheet",
                         href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"
                       ),
                       # Initialize Bootstrap 3 tooltips on document ready and Shiny updates
                       tags$script(HTML("
          $(document).on('shiny:connected shiny:inputchanged shiny:value', function() {
            $('[data-toggle=\"tooltip\"]').tooltip();
          });
        "))
                     ),
                     tags$div(
                       tags$h2("Simulate and Analyze Data",
                               style = "color: #2C3E50; font-weight: bold; text-align: left; margin-top: 20px;")
                     ),
                     tags$div(
                       tags$h4(
                         "Adjust simulation parameters and visualize model behavior.",
                         style = "color: #7F8C8D; text-align: left; margin-bottom: 20px;"
                       )
                     ),
                     titlePanel("Simulation Parameters"),
                     sidebarLayout(
                       sidebarPanel(
                         tags$label(
                           HTML("Development dataset sample size"),
                           tags$i(
                             class = "fas fa-info-circle",
                             style = "margin-left: 5px; color: #007BFF; cursor: pointer;",
                             `data-toggle` = "tooltip",
                             title = "Development dataset =  Ancillary dataset + Target dataset"
                           )
                         ),
                         div( style = "margin-bottom: 20px",
                              sliderInput(
                                "N_SAMPLES", NULL, min = 500, max = 9500, value = 1000,  step = 500
                              )),
                         
                         tags$label(
                           HTML("
                            Proportion split of the development dataset:<br>
                            <b>Target dataset proportion</b>"),
                           tags$i(
                             class = "fas fa-info-circle",
                             style = "margin-left: 5px; color: #007BFF; cursor: pointer;",
                             `data-toggle` = "tooltip",
                             title = "This is the proportion of the data allocated to the target dataset. The rest goes to the ancillary dataset = 1 - (Target dataset proportion)."
                           )
                         ),
                         div(style = "margin-bottom: 20px", sliderInput("target_split", NULL, min = 0.1, max = 0.9, value = 0.5, step = 0.01)),
                         tags$label(
                           HTML("Model development method:"),
                           tags$i(
                             class = "fas fa-info-circle",
                             style = "margin-left: 5px; color: #007BFF; cursor: pointer;",
                             `data-toggle` = "tooltip",
                             title = "See 'Methods' tab above for more information about each method."
                           )
                         ),
                         selectInput(
                           "cpm_method",
                           NULL,
                           choices = c("Full-data Logistic Regression","Target-dataset only logistic regression",
                                       "Ancillary-dataset only logistic regression", "Memebership-based Weighting",
                                       "Intercept Recalibration on all data", "Logistic Recalibration on all data")
                         )
                         ,
                         tags$label(
                           HTML("Predictors:")),
                         
                         tags$head(
                           tags$style(HTML("
      .panel-title > a:after {
        content: '\\25B6'; /* ▶ right arrow */
        float: right;
        margin-right: 10px;
        transition: transform 0.3s ease;
      }

      .panel-title > a[aria-expanded='true']:after {
        content: '\\25BC'; /* ▼ down arrow */
      }
    "))
                         ),
                         bsCollapsePanel(HTML("<strong>Intercept</strong>"),
                                         sliderInput( "X0_b_delta", "Intercept Shift",  min = -3, max = 3, value = 0,  step = 0.1 )),
                         
                         tags$head(
                           tags$style(HTML("
      .panel-title > a:after {
        content: '\\25B6'; /* ▶ right arrow */
        float: right;
        margin-right: 10px;
        transition: transform 0.3s ease;
      }

      .panel-title > a[aria-expanded='true']:after {
        content: '\\25BC'; /* ▼ down arrow */
      }
    ")) ),
                         bsCollapsePanel(HTML("<strong>X1 (Continuous)</strong>"),
                                         sliderInput("X1_m_delta", "X1 Mean Shift", min = -10, max = 10, value = 0, step = 0.1),
                                         sliderInput("X1_v_delta", "X1 Variance Shift", min = -10, max = 10, value = 0, step = 0.1),
                                         sliderInput("X1_b_delta", "X1 Beta Shift", min = -3, max = 3, value = 0, step = 0.1)
                         )
                         
                         
                         ,
                         tags$head(
                           tags$style(HTML("
      .panel-title > a:after {
        content: '\\25B6'; /* ▶ right arrow */
        float: right;
        margin-right: 10px;
        transition: transform 0.3s ease;
      }

      .panel-title > a[aria-expanded='true']:after {
        content: '\\25BC'; /* ▼ down arrow */
      }
    "))),
                         bsCollapsePanel(HTML("<strong>X2 (Continuous)</strong>"),
                                         sliderInput(
                                           "X2_m_delta",  "X2 Mean Shift", min = -10, max = 10,   value = 0, step = 0.1  ),
                                         sliderInput(
                                           "X2_v_delta",   "X2 Variance Shift",  min = -10,    max = 10,    value = 0,  step = 0.1 ),
                                         sliderInput(
                                           "X2_b_delta", "X2 Beta Shift", min = -3, max = 3, value = 0, step = 0.1 )
                         ),
                         
                         
                         tags$head(
                           tags$style(HTML("
      .panel-title > a:after {
        content: '\\25B6'; /* ▶ right arrow */
        float: right;
        margin-right: 10px;
        transition: transform 0.3s ease;
      }

      .panel-title > a[aria-expanded='true']:after {
        content: '\\25BC'; /* ▼ down arrow */
      }
    "))),
                         bsCollapsePanel(HTML("<strong>X3 (Binary)</strong>"),
                                         sliderInput(
                                           "X3_p_delta", "X3 Prevalence Shift",  min = -0.3,  max = 0.7, value = 0, step = 0.1),
                                         sliderInput(   "X3_b_delta", "X3 Beta Shift", min = -3, max = 3, value = 0,  step = 0.1)
                         ),
                         
                         tags$head(
                           tags$style(HTML("
      .panel-title > a:after {
        content: '\\25B6'; /* ▶ right arrow */
        float: right;
        margin-right: 10px;
        transition: transform 0.3s ease;
      }

      .panel-title > a[aria-expanded='true']:after {
        content: '\\25BC'; /* ▼ down arrow */
      }
    "))),
                         bsCollapsePanel(HTML("<strong>X4 (Binary)</strong>"),
                                         sliderInput(
                                           "X4_p_delta", "X4 Prevalence Shift", min = -0.4, max = 0.6, value = 0,  step = 0.1
                                         ),
                                         sliderInput(
                                           "X4_b_delta", "X4 Beta Shift", min = -3, max = 3, value = 0, step = 0.1
                                         )
                         ),
                         
                         tags$head(
                           tags$style(HTML("
      .panel-title > a:after {
        content: '\\25B6'; /* ▶ right arrow */
        float: right;
        margin-right: 10px;
        transition: transform 0.3s ease;
      }

      .panel-title > a[aria-expanded='true']:after {
        content: '\\25BC'; /* ▼ down arrow */
      }
    "))),
                         bsCollapsePanel(HTML("<strong>X5 (Binary)</strong>"),
                                         sliderInput("X5_p_delta", "X5 Prevalence Shift", min = -0.1, max = 0.9,  value = 0, step = 0.1
                                         ),
                                         sliderInput(
                                           "X5_b_delta", "X5 Beta Shift",  min = -3, max = 3, value = 0, step = 0.1)
                         ),
                         
                         div(
                           style = "display: flex; flex-wrap: wrap",  # gap adds space between buttons
                           
                           actionButton("simulate", "Simulate Data", style = "margin-bottom: 20px; margin-top: 20px; margin-right:10px"),
                           actionButton("reset", "Reset Inputs", icon = icon("undo"), style = "margin-bottom: 20px; margin-top: 20px")
                         ),
                         
                         hidden(div(
                           id = "sim-spinner",
                           tags$h5("Running Simulation..."),
                           tags$i(class = "fas fa-hourglass-half fa-spin fa-2x")
                         ))
                       ),
                       
                       mainPanel(
                         waiter::use_waiter(),
                         conditionalPanel(
                           condition = "output.plotReady == true",
                           
                           # No Shift Section
                           tags$label(
                             h3(span("No Data Shift Scenario", style = "color: #2E8B57; font-weight: bold;"),
                                tags$i(
                                  class = "fas fa-info-circle",
                                  style = "margin-left: 5px; color: #007BFF; cursor: pointer;",
                                  `data-toggle` = "tooltip",
                                  title = "The ancillary population is generated from same distribution as the target population"
                                ))
                           ),
                           
                           div(
                             style = "display: flex; align-items: flex-start; gap: 20px; margin-bottom: 40px;",
                             div(
                               id = "table-container_noshift",
                               style = "width: 54%; min-height: 300px; margin: 0; padding: 0;",
                               tableOutput("results_summary_noshift")
                             ),
                             div(
                               id = "calPlot-container_noshift",
                               style = "width: 49%; margin-top: -80px; padding: 0;",
                               plotOutput("calPlot_noshift", height = "370px")
                             )
                           ),
                           
                           # Shift Section
                           h3(span("Selected Scenario", style = "color: #0072B2; font-weight: bold;")),
                           div(
                             style = "display: flex; align-items: flex-start; gap: 20px;",
                             div(
                               id = "table-container",
                               style = "width: 54%; min-height: 300px; margin: 0; padding: 0;",
                               tableOutput("results_summary")
                             ),
                             div(
                               id = "calPlot-container",
                               style = "width: 49%; margin-top: -80px; padding: 0;",
                               plotOutput("calPlot", height = "370px")
                             )
                           )
                         )
                       )
                     )
                   )
                 ),
                 tabPanel(
                   "🛠️ Methods",
                   fluidPage(
                     withMathJax(),
                     
                     h4("This simulation focuses on logistic regression-based models to estimate the probability of the binary outcome \\(Y\\)."),
                     tags$li(HTML("<strong>Develop model on full data</strong>: <br>This method uses the standard unweighted logistic regression to develop a model on both target and ancillary data. <br> <br><br>")),
                     tags$li(HTML("<strong>Develop model on target only</strong>: <br>This method uses the standard unweighted logistic regression to develop a model on the target data only. <br> <br><br>")),
                     tags$li(HTML("<strong>Develop model on ancillary only</strong>: <br>This method uses the standard unweighted logistic regression to develop a model on the ancillary data only. <br> <br><br>")),
                     
                     tags$li(
                       HTML(
                         "<strong>Model Re-calibration: Updating the intercept only</strong>: <br> 
                         One of discrete model updating methods that uses the new available data to update the model is Re-calibration by updating the intercept only, which intends to correct “calibration-in-the-large”, 
                         i.e make the average predicted probability equal to the observed overall event rate. This method is implemented by develop the model on all data and update on target only,
                         a logistic regression model is developed on both ancillary and target (full development dataset), referred as `full model', then the the linear predictor 
                         of the `full model' is updated  on target only to obtain \\(LP_\\text{target}\\). Next, a logistic regression 
                         is fit on the target with fixing \\(LP_\\text{target}\\) as offset variable (the regression coefficient is constrained to be 1).<br> <br><br>"
                       )
                     ),
                     tags$li(
                       HTML(
                         "<strong>Model Re-calibration: Logistic calibration</strong>: <br> 
                         Another method of discrete model updating that uses the new available data to update the model is Re-calibration by Logistic Calibration, which intends to correct both the intercept and the overall calibration slope. 
                         This method is implemented by developing the model on all data and update on target only,
                         a logistic regression model is developed on both ancillary and target (full development dataset), referred as `full model', then the the linear predictor 
                         of the `full model' is updated  on target only to obtain \\(LP_\\text{target}\\). Next, a logistic regression 
                         is fit on the target with setting \\(LP_\\text{target}\\) as the only variable.<br> <br><br>"
                       )
                       
                     ),
                     tags$li(
                       HTML(
                         "<strong>Membership-based recalibration method: Develop on ancillary and target</strong><br>
                      This method aims to develop a weighted model by giving individuals from the ancillary data high
                      weights when they are relevant to target and low weights when they are less relevant to target. 
                      The weights are defined based on the relatedness of the individual samples from the ancillary
                      data compared with the target data using a probabilistic similarity metric called the 
                      'membership model' that is based on the propensity score, thus adjusting the distribution of 
                      the ancillary data to correct for distributional changes.<br>
                      <br>
                      First, the membership propensity score \\(PS\\) is utilized, which is defined as the conditional probability 
                      of a randomly chosen individual in the study population belonging to the ancillary group, 
                      given their observed values from a set of predictor variables 
                      \\(\\boldsymbol{X}_i = (X_{i1}, \\ldots, X_{iK})\\), where \\(K\\) denotes 
                      the number of measured predictor variables, such that:
                      <br><br>
                      \\[
                      PS_i = P(R = 1 \\mid \\mathbf{X}_i), \\quad \\text{for } i \\in \\{1, \\ldots, N_{\\text{ancillary}}\\}
                      \\]
                      <br>
                      
                      <br>
                      A binary logistic regression model (membership model) is used to estimate \\(PS\\), where the outcome of the membership model is 0 for 
                      individuals of target and 1 for individuals of ancillary. 
                      Then, a weight \\(w_i\\) is produced for an individual \\(i\\) from ancillary dataset by dividing the Membership
                      propensity score of target by the Membership propensity score of ancillary. 
                      To account for differences in sample sizes between the target and ancillary datasets, 
                      the weight is adjusted by multiplying it by the ratio of sample sizes (i.e., the sample size of the ancillary dataset divided by the sample size of the 
                      target dataset). This adjustment ensures that each dataset contributes proportionally to the final weights, 
                      reflecting the relative sizes of the datasets. The weight is calculated as follows: 
                      \\[
                        w_i = \\min\\left( 
                        \\left( \\frac{P(R = 0 \\mid \\mathbf{X}_i)}{P(R = 1 \\mid \\mathbf{X}_i)} 
                        \\times \\frac{N_{\\text{ancillary}}}{N_{\\text{target}}} \\right),\\ 1 
                        \\right) \\times \\lambda
                        \\]
                                          
                      where \\(\\lambda\\) is an overall `forgetting factor' to further down-weight the ancillary, estimated as a hyperparameter using 
                      4-fold cross-validation. Adding \\(\\lambda\\) is crucial to correct for predictor-outcome association shift, if any, as
                      the Membership propensity score only correct for case-mix shift.  
                      After that, a weighted regression model is used with \\(K\\) set of predictors to weight each individual from ancillary in the 
                      log likelihood based on their relatedness to target, while weighting individuals from target with weight equal to 1. 
                      Let individual \\(i\\) and \\(w_i\\) is the corresponding individual weight from weights set that was calculated based on the
                      importance of each individual. 
                      The weights are then included in the log-likelihood \\(LL\\) of the model as follows:
                      \\[
                        \\text{LL} = \\sum_{i=1}^{N} \\left(Y_i \\times \\log(P_i) + (1 - Y_i) \\times \\log(1 - P_i)\\right) \\times w_i
                        \\]
                        where \\(P\\) is the predicted probability for a binary outcome. 
                      <br>
                      <br>This method  is implemented by Recalibration, where \\(PS\\) is estimateed in the weighting scheme given the observed values from 
                      a set of predictor variables \\(\\boldsymbol{X}_i = (X_{i1}, \\ldots, X_{iK})\\) only. Then after developing the weighted model, 
                      the calibration in the large (CITL) of the weighted model is adjustedthrough adding a dummy variable 
                      for the target membership \\(R\\). The dummy variable \\(R\\) indicates the absence or presence 
                      of the effect of target membership that causes a shift for the outcome, where it takes a value of 0 if the samples
                      belong to target, and a value of 1 if the samples belong to ancillary. The dummy variable is added as an additional
                      predictor in the developed model. 
                      <br>
                      Note: This method reuqire longer time to run compared to other methods due to the cross-validation process.
                      "
                       )
                     )
                   )
                 )
                 
)


server <- function(input, output, session) {
  sim_data <- reactiveVal()
  noshift_data <- reactiveVal()
  hide("sim-spinner")
  
  output$plotReady <- reactive({
    !is.null(sim_data()) && !is.null(noshift_data())
  })
  outputOptions(output, "plotReady", suspendWhenHidden = FALSE)
  error_notification_id <- NULL
  
  observeEvent(input$simulate, {
    sim_data(NULL)
    noshift_data(NULL)
    
    if (!is.null(error_notification_id)) {
      removeNotification(error_notification_id)
      error_notification_id <<- NULL
    }
    
    # Capture inputs outside the future (good!)
    N_SAMPLES <- input$N_SAMPLES
    target_split <- input$target_split
    X0_b_delta <- input$X0_b_delta
    X1_m_delta <- input$X1_m_delta
    X1_v_delta <- input$X1_v_delta
    X1_b_delta <- input$X1_b_delta
    X2_m_delta <- input$X2_m_delta
    X2_v_delta <- input$X2_v_delta
    X2_b_delta <- input$X2_b_delta
    X3_p_delta <- input$X3_p_delta
    X3_b_delta <- input$X3_b_delta
    X4_p_delta <- input$X4_p_delta
    X4_b_delta <- input$X4_b_delta
    X5_p_delta <- input$X5_p_delta
    X5_b_delta <- input$X5_b_delta
    cpm_method <- input$cpm_method
    
    show("sim-spinner")
    
    # --- FUTURE 1: simulate_data ---
    fut1 <- future({
      
      message(Sys.time(), ": Starting simulate_data")
      
      result <- simulate_data(
        N_SAMPLES = N_SAMPLES,
        target_split = target_split,
        X0_b_delta = X0_b_delta,
        X1_m_delta = X1_m_delta,
        X1_v_delta = X1_v_delta,
        X1_b_delta = X1_b_delta,
        X2_m_delta = X2_m_delta,
        X2_v_delta = X2_v_delta,
        X2_b_delta = X2_b_delta,
        X3_p_delta = X3_p_delta,
        X3_b_delta = X3_b_delta,
        X4_p_delta = X4_p_delta,
        X4_b_delta = X4_b_delta,
        X5_p_delta = X5_p_delta,
        X5_b_delta = X5_b_delta,
        cpm_method = cpm_method
      )
      
      message(Sys.time(), ": Finished simulate_data")
      result
    }, seed = TRUE , packages = c("tidyr", "promises", "Matrix"))
    
    # --- FUTURE 2: simulate_noshift ---
    fut2 <- future({
      
      message(Sys.time(), ": Starting simulate_noshift")
      
      result <- simulate_noshift(
        N_SAMPLES = N_SAMPLES,
        target_split = target_split,
        cpm_method = cpm_method
      )
      
      message(Sys.time(), ": Finished simulate_noshift")
      result
    }, seed = TRUE, packages = c("tidyr", "promises", "Matrix"))
    
    # --- Handle both futures ---
    promise_all(fut1 = fut1, fut2 = fut2) %...>% (function(results) {
      sim_data(results$fut1)
      noshift_data(results$fut2)
      
      hide("sim-spinner")
      enable("simulate")
      
      # Remove any old error message if present
      if (!is.null(error_notification_id)) {
        removeNotification(error_notification_id)
        error_notification_id <<- NULL
      }
      
    }) %...!% (function(e) {
      print(e)

      hide("sim-spinner")
      enable("simulate")
      
      # Remove any old notification first
      if (!is.null(error_notification_id)) {
        removeNotification(error_notification_id)
      }
      
      # Show a persistent error message
      id <- showNotification(
        "⚠️ Failed to compute covariance matrix. Please choose different values for shifts.",
        #e$message,
        type = "error",
        duration = NULL  # persistent until removed
      )
      
      error_notification_id <<- id
    })
  })
  
  
  observeEvent(input$reset, {
    # (unchanged)
    updateNumericInput(session, "N_SAMPLES", value = 1000)
    updateSliderInput(session, "target_split", value = 0.5)
    updateSliderInput(session, "X0_b_delta", value = 0)
    updateSliderInput(session, "X1_m_delta", value = 0)
    updateSliderInput(session, "X1_v_delta", value = 0)
    updateSliderInput(session, "X1_b_delta", value = 0)
    updateSliderInput(session, "X2_m_delta", value = 0)
    updateSliderInput(session, "X2_v_delta", value = 0)
    updateSliderInput(session, "X2_b_delta", value = 0)
    updateSliderInput(session, "X3_p_delta", value = 0)
    updateSliderInput(session, "X3_b_delta", value = 0)
    updateSliderInput(session, "X4_p_delta", value = 0)
    updateSliderInput(session, "X4_b_delta", value = 0)
    updateSliderInput(session, "X5_p_delta", value = 0)
    updateSliderInput(session, "X5_b_delta", value = 0)
    updateSelectInput(session, "cpm_method", selected = "Full-data Logistic Regression")
    
    # Also clear any existing error notification
    if (!is.null(error_notification_id)) {
      removeNotification(error_notification_id)
      error_notification_id <<- NULL
    }
  })
  
  # --- Outputs ---
  output$calPlot <- renderPlot({
    req(sim_data())
    calibration_plot(sim_data()$cal_plots)
  })
  
  output$results_summary <- renderTable({
    req(sim_data())
    validation_results_table(sim_data()$val_results)
  }, sanitize.text.function = function(x) x)
  
  output$calPlot_noshift <- renderPlot({
    req(noshift_data())
    calibration_plot_noshift(noshift_data()$cal_plots_noshift)
  })
  
  output$results_summary_noshift <- renderTable({
    req(noshift_data())
    validation_results_table_noshift(noshift_data()$val_results_noshift)
  }, sanitize.text.function = function(x) x)
}



# Run the app
shinyApp(ui = ui, server = server)
