# Distilltrees
# Distill knowledge of addative tree models (GBM) into generalised linear models (GLM)
# Ref: https://detralytics.com/wp-content/uploads/2023/10/Detra-Note_Additive-tree-ensembles.pdf

library(dplyr)
library(xgboost)


extract_bins_from_xgb <- function(model, feature_names){
  
  # Get tree dump
  tree_df <- xgb.model.dt.tree(model = model)
  
  bins <- list()
  for (feat in feature_names){
    # Extract split points for the feature
    splits <- tree_df %>%
      filter(Feature == feat, !is.na(Split)) %>%
      pull(Split) %>%
      unique() %>%
      sort()
    
    bins[[feat]] <- splits
  }
  return(bins)
}

bin_features <- function(data, bins, feature_names, feature_indices){
  data_binned <- data
  
  for (i in seq_along(feature_names)){
    feat <- feature_names[i]
    feat_idx <- as.character(feature_indices[i])
    
    if (length(bins[[feat_idx]]) > 0){
      data_binned[[paste0(feat, "_binned")]] <- as.numeric(cut(
        data[[feat]],
        breaks = c(-Inf, bins[[feat_idx]], Inf),
        labels = FALSE
      ))
    } else {
      # No splits, keep original
      data_binned[[paste0(feat, "_binned")]] <- data[[feat]]
    }
  }
  return(data_binned)
}



fit_distilltree <- function(data,
                            vdt,
                            num_features = c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density"),
                            cat_features = c("VehBrand", "VehGas", "Region")) {
  # Prepare data
  y <- data$ClaimNb
  
  formula_vars <- c(num_features, cat_features)
  
  # Convert to factors
  data_prep <- data
  for (feat in cat_features) {
    data_prep[[feat]] <- as.factor(data_prep[[feat]])
  }
  
  # Convert to factors for validation data
  vdt_prep <- vdt
  vdt_x_val <- vdt$x_val
  for (feat in cat_features) {
    vdt_x_val[[feat]] <- as.factor(vdt_x_val[[feat]])
  }
  vdt_prep$x_val <- vdt_x_val[, c(num_features, cat_features)]
  
  
  # 1. Fit XGB using train_XGBoost
  xgb_model <- train_XGBoost(
    dt = data_prep[, c(num_features, cat_features)],
    y = y,
    vdt = vdt_prep,
    nrounds = 2, # Number of boosting rounds, limited to 3 for testing and memory usage
    max_depth = 1 # Max tree depth, limited for testing and memory usage
  )
  
  # 2. Extract bins from trees
  num_feature_indices <- 0:(length(num_features) - 1)
  bins <- extract_bins_from_xgb(xgb_model, as.character(num_feature_indices))
  
  data_binned <- bin_features(data_prep, bins, num_features, num_feature_indices)
  
  binned_features <- paste0(num_features, "_binned")
  
  # 3. Fit GLM per tree on binned features
  n_trees <- xgb_model$best_iteration
  if (is.null(n_trees)) {
    n_trees <- xgb_model$niter
  }
  
  glms <- list()
  glm_formula_binned <- as.formula(paste("y ~", paste(c(binned_features, cat_features), collapse = " + ")))
  
  # Fit one GLM per tree
  glm_data = cbind(data_binned, y = y)
  for (tree_idx in 1:n_trees) {
    glm_tree <- glm(
      formula = glm_formula_binned,
      data = glm_data,
      family = poisson()
    )
    glms[[tree_idx]] <- glm_tree
  }
  
  # Aggregate predictions from GLMs
  glm_preds_list <- lapply(glms, function(glm_tree) {
    predict(glm_tree, newdata = data_binned, type = "response")
  })
  glm_preds_agg <- Reduce("+", glm_preds_list) / length(glms)
  
  # 4. Autocalibration - fit final GLM
  final_glm <- glm(
    y ~ offset(log(glm_preds_agg)),
    family = poisson()
  )
  
  
  # Store model
  model <- list(
    xgb_model = xgb_model,
    bins = bins,
    glms = glms,
    final_glm = final_glm,
    num_features = num_features,
    cat_features = cat_features,
    num_feature_indices = num_feature_indices
  )
  
  class(model) <- "distilltree"
  
  return(model)
}


predict.distilltree <- function(model, data, type = "response") {
  
  # Prepare data
  data_prep <- data
  for (feat in model$cat_features) {
    data_prep[[feat]] <- as.factor(data_prep[[feat]])
  }
  
  
  # Bin the data
  data_binned <- bin_features(data_prep, model$bins, model$num_features, model$num_feature_indices)
  
  # Aggregate predictions from all tree GLMs
  glm_preds_list <- lapply(model$glms, function(glm_tree) {
    predict(glm_tree, newdata = data_binned, type = "response")
  })
  glm_preds_agg <- Reduce("+", glm_preds_list) / length(model$glms)
  
  # Apply final calibration
  final_preds <- glm_preds_agg * exp(coef(model$final_glm))[1]
  
  return(as.vector(final_preds))
}