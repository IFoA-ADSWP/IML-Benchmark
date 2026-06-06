fit_localglmnet_encoder <- function(df, feature_cols) {
  stopifnot(is.data.frame(df), length(feature_cols) > 0)

  x <- df[, feature_cols, drop = FALSE]
  is_num <- vapply(x, is.numeric, logical(1))

  num_cols <- names(is_num)[is_num]
  cat_cols <- names(is_num)[!is_num]

  numeric_stats <- lapply(num_cols, function(col) {
    v <- as.numeric(x[[col]])
    mu <- mean(v, na.rm = TRUE)
    sdv <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    list(mean = mu, sd = sdv)
  })
  names(numeric_stats) <- num_cols

  categorical_levels <- lapply(cat_cols, function(col) {
    vals <- as.character(x[[col]])
    vals[is.na(vals)] <- "__NA__"
    sort(unique(vals))
  })
  names(categorical_levels) <- cat_cols

  list(
    feature_cols = feature_cols,
    num_cols = num_cols,
    cat_cols = cat_cols,
    numeric_stats = numeric_stats,
    categorical_levels = categorical_levels
  )
}


enforce_localglmnet_seen_levels <- function(dt, feature_cols, train_idx, holdout_idx) {
  x <- dt[, feature_cols, drop = FALSE]
  is_num <- vapply(x, is.numeric, logical(1))
  cat_cols <- names(is_num)[!is_num]

  if (length(cat_cols) == 0) {
    return(list(train_idx = train_idx, holdout_idx = holdout_idx))
  }

  train_idx <- as.integer(train_idx)
  holdout_idx <- as.integer(holdout_idx)

  for (col in cat_cols) {
    vals <- as.character(x[[col]])
    vals[is.na(vals)] <- "__NA__"

    levels_all <- unique(vals)
    levels_train <- unique(vals[train_idx])
    missing_levels <- setdiff(levels_all, levels_train)

    if (length(missing_levels) == 0) next

    for (lvl in missing_levels) {
      candidate <- holdout_idx[vals[holdout_idx] == lvl][1]
      if (is.na(candidate)) {
        stop(paste0("Could not enforce seen levels for column '", col, "'."))
      }
      holdout_idx <- setdiff(holdout_idx, candidate)
      train_idx <- c(train_idx, candidate)
    }
  }

  list(train_idx = sort(unique(train_idx)), holdout_idx = sort(unique(holdout_idx)))
}


assert_localglmnet_no_unseen_levels <- function(dt, feature_cols, train_idx, check_idx) {
  if (is.null(check_idx) || length(check_idx) == 0) return(invisible(TRUE))

  x <- dt[, feature_cols, drop = FALSE]
  is_num <- vapply(x, is.numeric, logical(1))
  cat_cols <- names(is_num)[!is_num]

  for (col in cat_cols) {
    vals <- as.character(x[[col]])
    vals[is.na(vals)] <- "__NA__"

    levels_train <- unique(vals[train_idx])
    levels_check <- unique(vals[check_idx])
    unseen <- setdiff(levels_check, levels_train)
    if (length(unseen) > 0) {
      stop(paste0(
        "Unseen categorical levels in '", col, "' for provided split: ",
        paste(unseen, collapse = ", ")
      ))
    }
  }

  invisible(TRUE)
}


transform_localglmnet_input <- function(df, encoder) {
  stopifnot(is.data.frame(df), is.list(encoder))

  x <- df[, encoder$feature_cols, drop = FALSE]

  num_mat <- NULL
  if (length(encoder$num_cols) > 0) {
    num_mat <- sapply(encoder$num_cols, function(col) {
      stats <- encoder$numeric_stats[[col]]
      v <- as.numeric(x[[col]])
      v[is.na(v)] <- stats$mean
      (v - stats$mean) / stats$sd
    }, simplify = "matrix")

    if (is.null(dim(num_mat))) {
      num_mat <- matrix(num_mat, ncol = 1)
      colnames(num_mat) <- encoder$num_cols
    }
  }

  cat_mat <- NULL
  if (length(encoder$cat_cols) > 0) {
    cat_blocks <- lapply(encoder$cat_cols, function(col) {
      levels_train <- encoder$categorical_levels[[col]]
      vals <- as.character(x[[col]])
      vals[is.na(vals)] <- "__NA__"
      if (any(!(vals %in% levels_train))) {
        stop(paste0("Unseen level found in column '", col, "'. Ensure split keeps all levels in train."))
      }

      block <- outer(vals, levels_train, FUN = "==") * 1
      colnames(block) <- paste0(col, "_", make.names(levels_train))
      block
    })
    cat_mat <- do.call(cbind, cat_blocks)
  }

  x_mat <- cbind(num_mat, cat_mat)
  storage.mode(x_mat) <- "double"
  colnames(x_mat) <- make.names(colnames(x_mat), unique = TRUE)
  x_mat
}


build_localglmnet_feature_label_map <- function(encoder) {
  raw_names <- c()

  if (length(encoder$num_cols) > 0) {
    raw_names <- c(raw_names, encoder$num_cols)
  }

  if (length(encoder$cat_cols) > 0) {
    cat_raw <- unlist(lapply(encoder$cat_cols, function(col) {
      levels_train <- encoder$categorical_levels[[col]]
      paste0(col, "_", make.names(levels_train))
    }), use.names = FALSE)
    raw_names <- c(raw_names, cat_raw)
  }

  clean_names <- make.names(raw_names, unique = TRUE)
  stats::setNames(raw_names, clean_names)
}


resolve_localglmnet_plot_features <- function(feature_names, plot_vars = NULL, encoder = NULL) {
  if (is.null(plot_vars) || length(plot_vars) == 0) {
    return(feature_names)
  }

  selected <- c()
  for (v in plot_vars) {
    if (v %in% feature_names) {
      selected <- c(selected, v)
      next
    }

    if (!is.null(encoder)) {
      if (v %in% encoder$num_cols && v %in% feature_names) {
        selected <- c(selected, v)
      }

      if (v %in% encoder$cat_cols) {
        selected <- c(selected, grep(paste0("^", v, "_"), feature_names, value = TRUE))
      }
    }
  }

  selected <- unique(selected)
  if (length(selected) == 0) {
    stop("No plotting features matched plot_vars.")
  }

  selected
}


pretty_localglmnet_feature_label <- function(feature_name, feature_label_map = NULL) {
  raw <- feature_name
  if (!is.null(feature_label_map) && feature_name %in% names(feature_label_map)) {
    raw <- feature_label_map[[feature_name]]
  }

  label <- raw
  if (grepl("_", raw)) {
    parts <- strsplit(raw, "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      var_name <- parts[1]
      level <- paste(parts[-1], collapse = "_")
      level <- gsub("__NA__", "[missing]", level, fixed = TRUE)
      label <- paste0(var_name, " = ", level)
    }
  }

  label <- gsub("__NA__", "[missing]", label, fixed = TRUE)
  label
}


prepare_localglmnet_data <- function(df,
                                     target_col,
                                     feature_cols = NULL,
                                     train_idx = NULL,
                                     val_idx = NULL,
                                     test_idx = NULL,
                                     train_fraction = 0.6,
                                     val_fraction = 0.2,
                                     test_fraction = 0.2,
                                     seed = 2) {
  stopifnot(is.data.frame(df), target_col %in% colnames(df))

  if (is.null(feature_cols)) {
    feature_cols <- setdiff(colnames(df), target_col)
  }

  keep <- stats::complete.cases(df[, c(feature_cols, target_col), drop = FALSE])
  dt <- df[keep, c(feature_cols, target_col), drop = FALSE]

  y <- as.numeric(dt[[target_col]])
  if (any(!is.finite(y))) stop("Target contains non-finite values.")
  if (any(y < 0)) stop("Poisson target must be non-negative.")

  n <- nrow(dt)

  if (is.null(train_idx) || is.null(test_idx)) {
    set.seed(seed)
    idx_all <- sample.int(n, size = n, replace = FALSE)

    n_test <- floor(test_fraction * n)
    n_val <- floor(val_fraction * n)

    holdout_idx <- idx_all[seq_len(n_test + n_val)]
    train_idx <- setdiff(idx_all, holdout_idx)

    enforced <- enforce_localglmnet_seen_levels(
      dt = dt,
      feature_cols = feature_cols,
      train_idx = train_idx,
      holdout_idx = holdout_idx
    )

    train_idx <- enforced$train_idx
    holdout_idx <- enforced$holdout_idx

    test_idx <- holdout_idx[seq_len(min(n_test, length(holdout_idx)))]
    val_idx <- setdiff(holdout_idx, test_idx)
  }

  if (!is.null(train_idx) && is.null(val_idx) && val_fraction > 0) {
    set.seed(seed)
    n_val_from_train <- floor(val_fraction * length(train_idx))
    val_idx <- sample(train_idx, size = n_val_from_train, replace = FALSE)
    train_idx <- setdiff(train_idx, val_idx)
  }

  assert_localglmnet_no_unseen_levels(dt, feature_cols, train_idx, val_idx)
  assert_localglmnet_no_unseen_levels(dt, feature_cols, train_idx, test_idx)

  encoder <- fit_localglmnet_encoder(dt[train_idx, feature_cols, drop = FALSE], feature_cols)

  X_L <- transform_localglmnet_input(dt[train_idx, feature_cols, drop = FALSE], encoder)
  X_T <- transform_localglmnet_input(dt[test_idx, feature_cols, drop = FALSE], encoder)
  X_V <- NULL
  Y_V <- NULL
  D_V <- NULL
  if (!is.null(val_idx) && length(val_idx) > 0) {
    X_V <- transform_localglmnet_input(dt[val_idx, feature_cols, drop = FALSE], encoder)
    Y_V <- y[val_idx]
    D_V <- dt[val_idx, feature_cols, drop = FALSE]
  }

  list(
    X_L = X_L,
    Y_L = y[train_idx],
    D_L = dt[train_idx, feature_cols, drop = FALSE],
    X_V = X_V,
    Y_V = Y_V,
    D_V = D_V,
    X_T = X_T,
    Y_T = y[test_idx],
    D_T = dt[test_idx, feature_cols, drop = FALSE],
    feature_names = colnames(X_L),
    feature_label_map = build_localglmnet_feature_label_map(encoder),
    encoder = encoder,
    train_idx = train_idx,
    val_idx = val_idx,
    test_idx = test_idx
  )
}


build_localglmnet <- function(q,
                              hidden_units = c(20, 15, 10),
                              hidden_activation = "tanh",
                              optimizer = keras3::optimizer_nadam(),
                              loss = "poisson") {
  stopifnot(length(q) == 1, q > 0)

  design <- keras3::layer_input(shape = c(q), dtype = "float32", name = "Design")

  attention <- design
  for (u in hidden_units) {
    attention <- keras3::layer_dense(attention, units = u, activation = hidden_activation)
  }

  attention <- keras3::layer_dense(
    attention,
    units = q,
    activation = "linear",
    name = "Attention"
  )

  response <- keras3::layer_multiply(list(design, attention))
  response <- keras3::layer_dense(
    response,
    units = 1,
    activation = "exponential",
    name = "Response"
  )

  model <- keras3::keras_model(inputs = design, outputs = response)
  keras3::compile(model, loss = loss, optimizer = optimizer)
  model
}


fit_localglmnet <- function(df,
                            target_col,
                            feature_cols = NULL,
                            train_idx = NULL,
                            val_idx = NULL,
                            test_idx = NULL,
                            train_fraction = 0.6,
                            val_fraction = 0.2,
                            test_fraction = 0.2,
                            seed = 2,
                            hidden_units = c(20, 15, 10),
                            hidden_activation = "tanh",
                            optimizer = keras3::optimizer_nadam(),
                            loss = "poisson",
                            epochs = 100,
                            batch_size = 5000,
                            validation_split = 0.2,
                            patience = 10,
                            verbose = 1) {
  prep <- prepare_localglmnet_data(
    df = df,
    target_col = target_col,
    feature_cols = feature_cols,
    train_idx = train_idx,
    val_idx = val_idx,
    test_idx = test_idx,
    train_fraction = train_fraction,
    val_fraction = val_fraction,
    test_fraction = test_fraction,
    seed = seed
  )

  model <- build_localglmnet(
    q = ncol(prep$X_L),
    hidden_units = hidden_units,
    hidden_activation = hidden_activation,
    optimizer = optimizer,
    loss = loss
  )

  fit_args <- list(
    object = model,
    x = prep$X_L,
    y = prep$Y_L,
    epochs = epochs,
    batch_size = batch_size,
    callbacks = list(
      keras3::callback_early_stopping(
        monitor = "val_loss",
        patience = patience,
        restore_best_weights = TRUE
      )
    ),
    verbose = verbose
  )

  if (!is.null(prep$X_V) && !is.null(prep$Y_V)) {
    fit_args$validation_data <- list(prep$X_V, prep$Y_V)
  } else {
    fit_args$validation_split <- validation_split
  }

  history <- do.call(keras3::fit, fit_args)

  c(prep, list(model = model, history = history))
}


plot_localglmnet_attention <- function(model,
                                       X,
                                       feature_names = colnames(X),
                                       plot_vars = NULL,
                                       encoder = NULL,
                                       raw_df = NULL,
                                       sample_size = 5000,
                                       ncol = 3) {
  stopifnot(length(feature_names) == ncol(X), !is.null(encoder))

  if (is.null(raw_df)) {
    stop("Please pass raw_df (original, unencoded rows) to plot with original variable names.")
  }
  if (!is.data.frame(raw_df)) {
    stop("raw_df must be a data.frame aligned row-wise with X.")
  }
  if (nrow(raw_df) != nrow(X)) {
    stop("raw_df and X must have the same number of rows.")
  }

  candidate_vars <- c(encoder$num_cols, encoder$cat_cols)
  if (is.null(plot_vars) || length(plot_vars) == 0) {
    selected_vars <- candidate_vars
  } else {
    selected_vars <- intersect(plot_vars, candidate_vars)
  }
  if (length(selected_vars) == 0) {
    stop("No plot_vars matched original variables.")
  }

  col_groups <- lapply(selected_vars, function(v) {
    if (v %in% encoder$num_cols) {
      v
    } else {
      grep(paste0("^", v, "_"), feature_names, value = TRUE)
    }
  })
  names(col_groups) <- selected_vars
  col_groups <- col_groups[vapply(col_groups, length, integer(1)) > 0]
  selected_vars <- names(col_groups)

  attention_layer <- keras3::keras_model(
    inputs = model$input,
    outputs = model$get_layer("Attention")$output
  )
  attention_weights <- predict(attention_layer, X, verbose = 0)

  set.seed(42)
  sample_idx <- sample(seq_len(nrow(X)), size = min(sample_size, nrow(X)), replace = FALSE)

  plot_blocks <- lapply(selected_vars, function(v) {
    grp <- col_groups[[v]]
    grp_idx <- match(grp, feature_names)
    x_raw <- raw_df[[v]][sample_idx]

    if (v %in% encoder$num_cols && length(grp_idx) == 1) {
      att <- attention_weights[sample_idx, grp_idx]
      data.frame(
        feature_name = v,
        feature_value = as.numeric(x_raw),
        attention_value = as.numeric(att),
        is_numeric = TRUE
      )
    } else {
      att <- rowSums(attention_weights[sample_idx, grp_idx, drop = FALSE] * X[sample_idx, grp_idx, drop = FALSE])
      x_lbl <- as.character(x_raw)
      x_lbl[is.na(x_lbl)] <- "[missing]"
      data.frame(
        feature_name = v,
        feature_value = x_lbl,
        attention_value = as.numeric(att),
        is_numeric = FALSE
      )
    }
  })
  plot_df <- do.call(rbind, plot_blocks)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = feature_value, y = attention_value)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, position = ggplot2::position_jitter(width = 0.1, height = 0)) +
    ggplot2::facet_wrap(~feature_name, scales = "free", ncol = ncol) +
    ggplot2::labs(
      title = "Regression Attentions",
      x = "Original feature values",
      y = "Regression attention"
    ) +
    ggplot2::theme_minimal()
}


plot_localglmnet_contributions <- function(model,
                                           X,
                                           feature_names = colnames(X),
                                           plot_vars = NULL,
                                           encoder = NULL,
                                           raw_df = NULL,
                                           sample_size = 5000,
                                           ncol = 3) {
  stopifnot(length(feature_names) == ncol(X), !is.null(encoder))

  if (is.null(raw_df)) {
    stop("Please pass raw_df (original, unencoded rows) to plot with original variable names.")
  }
  if (!is.data.frame(raw_df)) {
    stop("raw_df must be a data.frame aligned row-wise with X.")
  }
  if (nrow(raw_df) != nrow(X)) {
    stop("raw_df and X must have the same number of rows.")
  }

  candidate_vars <- c(encoder$num_cols, encoder$cat_cols)
  if (is.null(plot_vars) || length(plot_vars) == 0) {
    selected_vars <- candidate_vars
  } else {
    selected_vars <- intersect(plot_vars, candidate_vars)
  }
  if (length(selected_vars) == 0) {
    stop("No plot_vars matched original variables.")
  }

  col_groups <- lapply(selected_vars, function(v) {
    if (v %in% encoder$num_cols) {
      v
    } else {
      grep(paste0("^", v, "_"), feature_names, value = TRUE)
    }
  })
  names(col_groups) <- selected_vars
  col_groups <- col_groups[vapply(col_groups, length, integer(1)) > 0]
  selected_vars <- names(col_groups)

  attention_layer <- keras3::keras_model(
    inputs = model$input,
    outputs = model$get_layer("Attention")$output
  )
  attention_weights <- predict(attention_layer, X, verbose = 0)

  set.seed(42)
  sample_idx <- sample(seq_len(nrow(X)), size = min(sample_size, nrow(X)), replace = FALSE)

  plot_blocks <- lapply(selected_vars, function(v) {
    grp <- col_groups[[v]]
    grp_idx <- match(grp, feature_names)
    x_raw <- raw_df[[v]][sample_idx]

    contribution <- rowSums(attention_weights[sample_idx, grp_idx, drop = FALSE] * X[sample_idx, grp_idx, drop = FALSE])

    if (v %in% encoder$num_cols) {
      data.frame(
        feature_name = v,
        feature_value = as.numeric(x_raw),
        contribution = as.numeric(contribution),
        is_numeric = TRUE
      )
    } else {
      x_lbl <- as.character(x_raw)
      x_lbl[is.na(x_lbl)] <- "[missing]"
      data.frame(
        feature_name = v,
        feature_value = x_lbl,
        contribution = as.numeric(contribution),
        is_numeric = FALSE
      )
    }
  })
  plot_df <- do.call(rbind, plot_blocks)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = feature_value, y = contribution)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, color = "black", position = ggplot2::position_jitter(width = 0.1, height = 0)) +
    ggplot2::geom_hline(yintercept = 0, color = "orange", linewidth = 0.8) +
    ggplot2::facet_wrap(~feature_name, scales = "free", ncol = ncol) +
    ggplot2::labs(
      title = "Feature Contributions",
      x = "Original feature values",
      y = "Contribution"
    ) +
    ggplot2::theme_minimal()
}


fr_data <- dt_list$fre_mtpl2_freq %>% select(-IDpol)
fr_features <- setdiff(colnames(fr_data), "ClaimNb")

# Example run on French MTPL data using helpers from init.R
if (exists("dt_list") && exists("poisson_deviance") && exists("multiple_lift")) {
  set.seed(2)
  n_all <- nrow(fr_data)
  idx_all <- sample.int(n_all, size = n_all, replace = FALSE)

  n_train <- floor(0.6 * n_all)
  n_val <- floor(0.2 * n_all)

  train_idx <- idx_all[seq_len(n_train)]
  val_idx <- idx_all[seq.int(n_train + 1, n_train + n_val)]
  test_idx <- idx_all[seq.int(n_train + n_val + 1, n_all)]

  fr_fit <- fit_localglmnet(
    df = fr_data,
    target_col = "ClaimNb",
    feature_cols = fr_features,
    train_idx = train_idx,
    val_idx = val_idx,
    test_idx = test_idx,
    epochs = 100,
    batch_size = 5000,
    validation_split = 0,
    patience = 10,
    verbose = 1
  )

  xgb_fit <- train_XGBoost(
    dt = fr_data[train_idx, fr_features, drop = FALSE],
    y = fr_data$ClaimNb[train_idx],
    vdt = list(
      x_val = fr_data[val_idx, fr_features, drop = FALSE],
      y_val = fr_data$ClaimNb[val_idx]
    )
  )

  pred_test <- as.numeric(predict(fr_fit$model, fr_fit$X_T, verbose = 0)[,1])
  pred_xgb_test <- as.numeric(predict(xgb_fit, data.matrix(fr_data[test_idx, fr_features, drop = FALSE])))

  pd_test <-          poisson_deviance(y_true = fr_fit$Y_T, y_pred = pred_test)
  pd_xgb_test <-      poisson_deviance(y_true = fr_fit$Y_T, y_pred = pred_xgb_test)
  homogenous_model <- poisson_deviance(y_true = fr_fit$Y_T, y_pred = mean(fr_fit$Y_L))
  pinball <- 1 - pd_test / homogenous_model
  pinball_xgb <- 1 - pd_xgb_test / homogenous_model
  
  cat(sprintf("LocalGLMnet test Poisson deviance: %.6f\n", pd_test))
  cat(sprintf("Pinball score vs homogeneous baseline: %.6f\n", pinball))
  cat(sprintf("XGBoost test Poisson deviance: %.6f\n", pd_xgb_test))
  cat(sprintf("XGBoost pinball score vs homogeneous baseline: %.6f\n", pinball_xgb))

  plot_vars <- c("VehPower", "DrivAge", "BonusMalus", "Density")

  lift_plot <- multiple_lift(
    y_true = fr_fit$Y_T,
    y_pred_df = data.frame(LocalGLMnet = pred_test, XGBoost = pred_xgb_test),
    tiles = 10
  )
  print(lift_plot)

  attention_plot <- plot_localglmnet_attention(
    model = fr_fit$model,
    X = fr_fit$X_T,
    feature_names = fr_fit$feature_names,
    plot_vars = plot_vars,
    encoder = fr_fit$encoder,
    raw_df = fr_fit$D_T
  )
  print(attention_plot)

  contribution_plot <- plot_localglmnet_contributions(
    model = fr_fit$model,
    X = fr_fit$X_T,
    feature_names = fr_fit$feature_names,
    plot_vars = plot_vars,
    encoder = fr_fit$encoder,
    raw_df = fr_fit$D_T
  )
  print(contribution_plot)
}



