# models/train_ebm2.R

train_EBM2 <- function(dt,
                       y,
                       max_rounds = 1000L,
                       early_stopping_rounds = 10L,
                       learning_rate = 0.06,
                       random_state = 2L,
                       n_jobs = 1L) {
  if (!requireNamespace("ebm", quietly = TRUE)) {
    stop("R package 'ebm' is required. Install via install.packages('ebm').")
  }

  train_df <- as.data.frame(dt)
  train_df$ClaimNb <- as.numeric(y)

  fit <- ebm::ebm(
    formula = ClaimNb ~ .,
    data = train_df,
    objective = "poisson_deviance",
    max_rounds = as.integer(max_rounds),
    early_stopping_rounds = as.integer(early_stopping_rounds),
    learning_rate = learning_rate,
    random_state = as.integer(random_state),
    n_jobs = as.integer(n_jobs)
  )

  list(model = fit)
}


predict_EBM2 <- function(model, vdt, target = "ClaimNb") {
  new_df <- as.data.frame(vdt)
  y_test <- as.numeric(new_df[[target]])

  preds <- as.numeric(
    predict(model, newdata = new_df, type = "response")
  )

  list(
    predictions = preds,
    actuals = y_test
  )
}
