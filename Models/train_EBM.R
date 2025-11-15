# models/train_EBM.R

train_EBM <- function(dt, y) {
  
  ebm <- interpret$ExplainableBoostingRegressor(n_jobs = 1L, 
                                                max_rounds = 1000L,
                                                objective = "poisson_deviance",
                                                max_bins = 100L,
                                                learning_rate = 0.06,
                                                early_stopping_rounds = 5L,
                                                early_stopping_tolerance = 1e-4)
  
  X_train <- dt
  y_train <- as.numeric(y)

  X_py <- r_to_py(X_train)
  y_py <- r_to_py(y_train)

  ebm$fit(X_py, y_py)
  return(list(model = ebm))
}

predict_EBM <- function(model, vdt, target = "ClaimNb") {
  X_test <- data.matrix(vdt)
  y_test <- as.numeric(vdt[[target]])
  
  
  y_pred <- model$predict(r_to_py(X_test))
  y_pred_r <- py_to_r(y_pred)
  
  return(list(
    predictions = y_pred_r,
    actuals = y_test
  ))
}
