library(reticulate)
library(xgboost)
optuna <- import("optuna")

train_XGB_BO <- function(df_train, df_validation, n_trials){

    hy_params <- optimise_hyperparams(df_train, df_validation, n_trials)

    model <- xgboost::xgb.train(
        data = xgboost::xgb.DMatrix(data.matrix(df_train[,-1]),label = df_train$ClaimNb), 
        params = hy_params,
        nrounds = 1000
    )

    return(model)

    
}

optimise_hyperparams <- function(df_train, df_validation, n_trials){
    # Function for optimisation
    objective_function <- function(trial){
        # Set params
        params <- list(
            objective = "count:poisson",
            eval_metric = "poisson-nloglik",
            learning_rate = trial$suggest_float("learning_rate", 0.01, 0.1),
            max_depth = trial$suggest_int("max_depth", 2, 10),
            subsample = trial$suggest_float("subsample", 0.6, 1),
            colsample_bytree = trial$suggest_float("colsample_bytree", 0.6, 1.0)
        )

        # Train model
        model <- xgboost::xgb.train(
            data = xgboost::xgb.DMatrix(data.matrix(df_train[,-1]),label = df_train$ClaimNb), 
            params = params,
            nrounds = 1000
        )

        # Validate model
        preds <- as.vector(
            predict(
                model,
                xgboost::xgb.DMatrix(data.matrix(df_validation[,-1])),
                type="response"
            )
        )
        pd = poisson_deviance( # from init.R
            y_true = df_validation$ClaimNb,
            y_pred = preds
        )
        return(pd)

    }

    # Run optimisation
    study <- optuna$create_study(direction = "minimize") # minimise poisson deviance. Optuna uses TPE by default
    study$optimize(objective_function, n_trials = n_trials)

    # Get best params
    best_params <- study$best_params
    result <- list(
        objective = "count:poisson",  # Fixed parameter
        eval_metric = "poisson-nloglik",  # Fixed parameter
        learning_rate = best_params$learning_rate,
        max_depth = best_params$max_depth,
        subsample = best_params$subsample,
        colsample_bytree = best_params$colsample_bytree
    )

    return(result)
}