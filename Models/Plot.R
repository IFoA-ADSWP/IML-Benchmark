library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang)
library(viridis)

# Plot partial effect for DrivAge from a trained mgcv::bam GAM (Poisson, log link)
plot_GAM_DrivAge <- function(model, data,
                             var = "DrivAge",
                             n = 300,
                             plot_response_scale = FALSE,
                             se_mult = 2,
                             jitter_partial_resid = FALSE) {
  # Safety checks
  if (is.null(model)) stop("model is NULL")
  if (!var %in% names(data)) stop(glue::glue("{var} not found in data"))
  
  # Ensure numeric variable in data (use your to_num if you have it)
  if (!is.numeric(data[[var]])) {
    data[[var]] <- as.numeric(as.character(data[[var]]))
  }
  
  # prediction grid
  grd <- data.frame(!!sym(var) := seq(min(data[[var]], na.rm=TRUE),
                                      max(data[[var]], na.rm=TRUE),
                                      length.out = n))
  
  # If model has factor terms used in training, set them to reference/mode present in the model
  # model$xlevels often stores levels for factors used by predict; align them
  if (!is.null(model$xlevels)) {
    for (fv in names(model$xlevels)) {
      if (!fv %in% names(grd)) {
        # if factor exists in original data, set to first level of the trained model (safe default)
        grd[[fv]] <- factor(model$xlevels[[fv]][1], levels = model$xlevels[[fv]])
      } else {
        grd[[fv]] <- factor(grd[[fv]], levels = model$xlevels[[fv]])
      }
    }
  }
  
  # Get term predictions (type = "terms") for grid and se
  term_pred <- predict(model, newdata = grd, type = "terms", se.fit = TRUE)
  term_fit_mat <- term_pred$fit
  term_se_mat  <- term_pred$se.fit
  
  # Identify the column corresponding to the DrivAge smooth
  term_names_full <- colnames(term_fit_mat)
  term_index <- grep(var, term_names_full, ignore.case = TRUE)[1]
  if (is.na(term_index)) stop("Could not find the smooth term for DrivAge in model terms.")
  term_name <- term_names_full[term_index]
  
  s_fit <- as.numeric(term_fit_mat[, term_index])
  s_se  <- as.numeric(term_se_mat[, term_index])
  
  plot_df <- tibble(
    !!sym(var) := grd[[var]],
    s_fit = s_fit,
    s_se  = s_se,
    s_lower = s_fit - se_mult * s_se,
    s_upper = s_fit + se_mult * s_se
  )
  
  # Partial residuals at observed rows
  terms_obs <- predict(model, type = "terms")            # term contributions for observed data
  if (is.null(colnames(terms_obs))) stop("predict(..., type='terms') returned unexpected result.")
  # find column index in observed matrix
  obs_index <- grep(var, colnames(terms_obs), ignore.case = TRUE)[1]
  partial_term_obs <- as.numeric(terms_obs[, obs_index])
  resp_resids <- residuals(model, type = "response")     # response residuals
  partial_resid_val <- partial_term_obs + resp_resids
  
  partial_df <- tibble(
    !!sym(var) := data[[var]],
    partial_resid = partial_resid_val
  )
  
  # Transform to response scale if requested (Poisson + log-link)
  if (plot_response_scale) {
    # On response scale the total predicted mean multiplier from smooth is exp(s_fit)
    plot_df <- plot_df %>%
      mutate(
        mu_fit = exp(s_fit),
        mu_lower = exp(s_lower),
        mu_upper = exp(s_upper)
      )
    # For partial residuals: convert approximate partial resid to multiplicative effect
    partial_df <- partial_df %>% mutate(mu_partial = exp(partial_resid))
  }
  
  # Base ggplot
  if (!plot_response_scale) {
    p <- ggplot() +
      geom_ribbon(data = plot_df, aes_string(x = var, ymin = "s_lower", ymax = "s_upper"),
                  fill = "grey80", alpha = 0.6) +
      geom_line(data = plot_df, aes_string(x = var, y = "s_fit", color = var), size = 1.1) +
      geom_point(data = partial_df, aes_string(x = var, y = "partial_resid"),
                 alpha = 0.45, size = 1.6,
                 position = if (jitter_partial_resid) position_jitter(height = 0.0, width = 0.0) else position_identity()) +
      geom_rug(data = data, aes_string(x = var), sides = "b", alpha = 0.25) +
      labs(x = var,
           y = "Partial effect s(DrivAge) (log link)",
           title = glue::glue("GAM partial effect for {var} (link scale)"),
           subtitle = "Shaded band = approx ±2 SE; points = partial residuals") +
      scale_color_viridis_c(option = "D") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
  } else {
    p <- ggplot() +
      geom_ribbon(data = plot_df, aes_string(x = var, ymin = "mu_lower", ymax = "mu_upper"),
                  fill = "grey80", alpha = 0.6) +
      geom_line(data = plot_df, aes_string(x = var, y = "mu_fit", color = var), size = 1.1) +
      geom_point(data = partial_df, aes_string(x = var, y = "mu_partial"),
                 alpha = 0.45, size = 1.6,
                 position = if (jitter_partial_resid) position_jitter(height = 0.0, width = 0.0) else position_identity()) +
      geom_rug(data = data, aes_string(x = var), sides = "b", alpha = 0.25) +
      labs(x = var,
           y = "Multiplicative effect on expected ClaimNb (exp(s))",
           title = glue::glue("GAM partial effect for {var} (response scale)"),
           subtitle = "Shaded band = approx ±2 SE transformed via exp(); points = exp(partial residual)") +
      scale_color_viridis_c(option = "D") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
  }
  
  print(p)
  invisible(list(plot_df = plot_df, partial_df = partial_df, term_name = term_name))
}
