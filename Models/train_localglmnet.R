fit_localglmnet_encoder <- function(df, feature_cols) {
  stopifnot(is.data.frame(df), length(feature_cols) > 0)

  x <- df[, feature_cols, drop = FALSE]
  is_num <- vapply(x, is.numeric, logical(1))

  num_cols <- names(is_num)[is_num]
  cat_cols <- names(is_num)[!is_num]

  numeric_stats <- lapply(num_cols, function(col) {
    v <- as.numeric(x[[col]])
    mu <- mean(v, na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    sdv <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    list(mean = mu, sd = sdv)
  })
  names(numeric_stats) <- num_cols

  categorical_levels <- lapply(cat_cols, function(col) {
    vals <- as.character(x[[col]])
    vals[is.na(vals)] <- "__NA__"
    lvls <- sort(unique(vals))
    unique(c(lvls, "__UNSEEN__"))
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
      vals[!(vals %in% levels_train)] <- "__UNSEEN__"

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


prepare_localglmnet_data <- function(train_df,
                                     test_df,
                                     val_df = NULL,
                                     target_col = "ClaimNb",
                                     feature_cols = NULL) {
  stopifnot(is.data.frame(train_df), is.data.frame(test_df), target_col %in% colnames(train_df), target_col %in% colnames(test_df))

  if (is.null(val_df)) {
    val_df <- test_df
  }

  if (!is.data.frame(val_df) || !(target_col %in% colnames(val_df))) {
    stop("val_df must be a data.frame containing target_col.")
  }

  if (is.null(feature_cols)) {
    feature_cols <- setdiff(colnames(train_df), target_col)
  }

  missing_train <- setdiff(feature_cols, colnames(train_df))
  missing_val <- setdiff(feature_cols, colnames(val_df))
  missing_test <- setdiff(feature_cols, colnames(test_df))
  if (length(missing_train) > 0 || length(missing_val) > 0 || length(missing_test) > 0) {
    stop("All split data frames must contain the same feature_cols.")
  }

  y_train <- as.numeric(train_df[[target_col]])
  y_val <- as.numeric(val_df[[target_col]])
  y_test <- as.numeric(test_df[[target_col]])

  if (any(!is.finite(y_train)) || any(!is.finite(y_val)) || any(!is.finite(y_test))) {
    stop("Target contains non-finite values in one or more splits.")
  }
  if (any(y_train < 0) || any(y_val < 0) || any(y_test < 0)) {
    stop("Poisson target must be non-negative in all splits.")
  }

  encoder <- fit_localglmnet_encoder(train_df[, feature_cols, drop = FALSE], feature_cols)

  X_L <- transform_localglmnet_input(train_df[, feature_cols, drop = FALSE], encoder)
  X_V <- transform_localglmnet_input(val_df[, feature_cols, drop = FALSE], encoder)
  X_T <- transform_localglmnet_input(test_df[, feature_cols, drop = FALSE], encoder)

  list(
    X_L = X_L,
    Y_L = y_train,
    D_L = train_df[, feature_cols, drop = FALSE],
    X_V = X_V,
    Y_V = y_val,
    D_V = val_df[, feature_cols, drop = FALSE],
    X_T = X_T,
    Y_T = y_test,
    D_T = test_df[, feature_cols, drop = FALSE],
    feature_names = colnames(X_L),
    feature_label_map = build_localglmnet_feature_label_map(encoder),
    encoder = encoder
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


fit_localglmnet <- function(train_df,
                            test_df,
                            val_df = NULL,
                            target_col = "ClaimNb",
                            feature_cols = NULL,
                            hidden_units = c(20, 15, 10),
                            hidden_activation = "tanh",
                            optimizer = keras3::optimizer_nadam(),
                            loss = "poisson",
                            epochs = 100,
                            batch_size = 5000,
                            patience = 10,
                            verbose = 1) {
  prep <- prepare_localglmnet_data(
    train_df = train_df,
    test_df = test_df,
    val_df = val_df,
    target_col = target_col,
    feature_cols = feature_cols
  )

  model <- build_localglmnet(
    q = ncol(prep$X_L),
    hidden_units = hidden_units,
    hidden_activation = hidden_activation,
    optimizer = optimizer,
    loss = loss
  )

  history <- keras3::fit(
    object = model,
    x = prep$X_L,
    y = prep$Y_L,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(prep$X_V, prep$Y_V),
    callbacks = list(
      keras3::callback_early_stopping(
        monitor = "val_loss",
        patience = patience,
        restore_best_weights = TRUE
      )
    ),
    verbose = verbose
  )

  c(prep, list(model = model, history = history))
}


predict_localglmnet <- function(model_obj, new_df) {
  X_new <- transform_localglmnet_input(new_df[, model_obj$encoder$feature_cols, drop = FALSE], model_obj$encoder)
  as.numeric(predict(model_obj$model, X_new, verbose = 0)[, 1])
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
  if (is.null(raw_df) || !is.data.frame(raw_df) || nrow(raw_df) != nrow(X)) {
    stop("raw_df must be a data.frame aligned row-wise with X.")
  }

  candidate_vars <- c(encoder$num_cols, encoder$cat_cols)
  selected_vars <- if (is.null(plot_vars) || length(plot_vars) == 0) candidate_vars else intersect(plot_vars, candidate_vars)
  if (length(selected_vars) == 0) stop("No plot_vars matched original variables.")

  col_groups <- lapply(selected_vars, function(v) {
    if (v %in% encoder$num_cols) v else grep(paste0("^", v, "_"), feature_names, value = TRUE)
  })
  names(col_groups) <- selected_vars
  col_groups <- col_groups[vapply(col_groups, length, integer(1)) > 0]
  selected_vars <- names(col_groups)

  attention_idx <- which(vapply(model$model$layers, function(l) identical(l$name, "Attention"), logical(1)))
  if (length(attention_idx) != 1) {
    stop("Could not uniquely locate 'Attention' layer in localglmnet model.")
  }
  attention_layer <- keras3::keras_model(inputs = model$model$input, outputs = model$model$layers[[attention_idx]]$output)
  attention_weights <- predict(attention_layer, X, verbose = 0)

  set.seed(42)
  sample_idx <- sample(seq_len(nrow(X)), size = min(sample_size, nrow(X)), replace = FALSE)

  plot_blocks <- lapply(selected_vars, function(v) {
    grp <- col_groups[[v]]
    grp_idx <- match(grp, feature_names)
    x_raw <- raw_df[[v]][sample_idx]

    if (v %in% encoder$num_cols && length(grp_idx) == 1) {
      att <- attention_weights[sample_idx, grp_idx]
      data.frame(feature_name = v, feature_value = as.numeric(x_raw), value = as.numeric(att), is_numeric = TRUE)
    } else {
      att <- rowSums(attention_weights[sample_idx, grp_idx, drop = FALSE] * X[sample_idx, grp_idx, drop = FALSE])
      x_lbl <- as.character(x_raw)
      x_lbl[is.na(x_lbl)] <- "[missing]"
      data.frame(feature_name = v, feature_value = x_lbl, value = as.numeric(att), is_numeric = FALSE)
    }
  })
  plot_df <- do.call(rbind, plot_blocks)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = feature_value, y = value)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, position = ggplot2::position_jitter(width = 0.1, height = 0)) +
    ggplot2::facet_wrap(~feature_name, scales = "free", ncol = ncol) +
    ggplot2::labs(title = "Regression Attentions", x = "Original feature values", y = "Regression attention") +
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
  if (is.null(raw_df) || !is.data.frame(raw_df) || nrow(raw_df) != nrow(X)) {
    stop("raw_df must be a data.frame aligned row-wise with X.")
  }

  candidate_vars <- c(encoder$num_cols, encoder$cat_cols)
  selected_vars <- if (is.null(plot_vars) || length(plot_vars) == 0) candidate_vars else intersect(plot_vars, candidate_vars)
  if (length(selected_vars) == 0) stop("No plot_vars matched original variables.")

  col_groups <- lapply(selected_vars, function(v) {
    if (v %in% encoder$num_cols) v else grep(paste0("^", v, "_"), feature_names, value = TRUE)
  })
  names(col_groups) <- selected_vars
  col_groups <- col_groups[vapply(col_groups, length, integer(1)) > 0]
  selected_vars <- names(col_groups)

  attention_idx <- which(vapply(model$model$layers, function(l) identical(l$name, "Attention"), logical(1)))
  if (length(attention_idx) != 1) {
    stop("Could not uniquely locate 'Attention' layer in localglmnet model.")
  }
  attention_layer <- keras3::keras_model(inputs = model$model$input, outputs = model$model$layers[[attention_idx]]$output)
  attention_weights <- predict(attention_layer, X, verbose = 0)

  set.seed(42)
  sample_idx <- sample(seq_len(nrow(X)), size = min(sample_size, nrow(X)), replace = FALSE)

  plot_blocks <- lapply(selected_vars, function(v) {
    grp <- col_groups[[v]]
    grp_idx <- match(grp, feature_names)
    x_raw <- raw_df[[v]][sample_idx]
    contribution <- rowSums(attention_weights[sample_idx, grp_idx, drop = FALSE] * X[sample_idx, grp_idx, drop = FALSE])

    if (v %in% encoder$num_cols) {
      data.frame(feature_name = v, feature_value = as.numeric(x_raw), value = as.numeric(contribution), is_numeric = TRUE)
    } else {
      x_lbl <- as.character(x_raw)
      x_lbl[is.na(x_lbl)] <- "[missing]"
      data.frame(feature_name = v, feature_value = x_lbl, value = as.numeric(contribution), is_numeric = FALSE)
    }
  })
  plot_df <- do.call(rbind, plot_blocks)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = feature_value, y = value)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, color = "black", position = ggplot2::position_jitter(width = 0.1, height = 0)) +
    ggplot2::geom_hline(yintercept = 0, color = "orange", linewidth = 0.8) +
    ggplot2::facet_wrap(~feature_name, scales = "free", ncol = ncol) +
    ggplot2::labs(title = "Feature Contributions", x = "Original feature values", y = "Contribution") +
    ggplot2::theme_minimal()
}
